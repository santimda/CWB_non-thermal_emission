#!/usr/bin/env python3
"""
Check the default Apep model against observational constraints.
- X-ray: F_10-30 keV = (3.9±1.2)e-13 erg/s/cm2 from del Palacio+ 2023
- Radio: S_2GHz ~ 120 ± 15 mJy (del Palacio+ 2022, 2026, and references therein)
"""

import json
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTROLS = ROOT / "controls.f90"
RESULTS = Path(__file__).resolve().parent / "results"
F10_RANGE = (2.7e-13, 5.1e-13)
S2_RANGE = (105.0, 135.0)


# Replace one logical control while preserving the surrounding source text.
def set_control(text, name, value):
    pattern = rf"^(\s*logical,\s*parameter\s*::\s*{name}\s*=\s*).*$"
    updated, count = re.subn(pattern, rf"\g<1>{value}", text, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"Could not update control {name}")
    return updated


# Build and execute the current model, returning its terminal output.
def run_model():
    subprocess.run(["make", "run.x"], cwd=ROOT, check=True)
    command = ["./run.x"]
    try:
        return subprocess.run(
            command, cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout
    except subprocess.CalledProcessError as error:
        print(
            f"Command failed with exit status {error.returncode}: {' '.join(command)}",
            file=sys.stderr,
        )
        if error.stdout:
            print("--- stdout ---", file=sys.stderr)
            print(error.stdout, file=sys.stderr, end="")
        if error.stderr:
            print("--- stderr ---", file=sys.stderr)
            print(error.stderr, file=sys.stderr, end="")
        raise


# Extract the first numeric value matched by a regular expression.
def extract(pattern, output):
    match = re.search(pattern, output, flags=re.DOTALL)
    if match is None:
        raise RuntimeError(f"Could not find model output matching {pattern!r}")
    return float(match.group(1))


# Run the baseline model and validate its X-ray and radio constraints.
def main():
    RESULTS.mkdir(exist_ok=True)
    original_controls = CONTROLS.read_text()
    baseline_output = ""

    try:
        controls = original_controls
        controls = set_control(controls, "absgg", "on")
        controls = set_control(controls, "rad_ic1", "on")
        controls = set_control(controls, "rad_ic2", "on")
        controls = set_control(controls, "rad_syn1", "on")
        controls = set_control(controls, "rad_syn2", "on")
        controls = set_control(controls, "abssyn", "on")
        controls = set_control(controls, "maps", "on")
        controls = set_control(controls, "convolve", "on")
        CONTROLS.write_text(controls)
        baseline_output = run_model()
    finally:
        CONTROLS.write_text(original_controls)
        subprocess.run(["make", "run.x"], cwd=ROOT, check=True)

    f10 = extract(
        r"Absorption-corrected fluxes.*?F_10-30keV.*?=\s+([0-9.Ee+-]+)",
        baseline_output,
    )
    s2 = extract(r"Post-conv:.*?S\(2GHz\)=\s+([0-9.Ee+-]+)", baseline_output)
    checks = {
        "F_10-30keV": {"value": f10, "range": F10_RANGE, "passed": F10_RANGE[0] <= f10 <= F10_RANGE[1]},
        "S_2GHz_mJy": {"value": s2, "range": S2_RANGE, "passed": S2_RANGE[0] <= s2 <= S2_RANGE[1]},
    }
    (RESULTS / "baseline_constraints.json").write_text(json.dumps(checks, indent=2) + "\n")
    (RESULTS / "baseline_run.log").write_text(baseline_output)

    for name, check in checks.items():
        status = "OK" if check["passed"] else "ERROR"
        print(f"{status}: {name} = {check['value']:.6g}, expected {check['range']}")
    if not all(check["passed"] for check in checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
