#!/usr/bin/env python3
"""Run the Apep Mdot exploration used by the radio SED figure."""

import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SYSTEM_PARAMETERS = ROOT / "system_parameters.f90"
CONTROLS = ROOT / "controls.f90"
RESULTS = Path(__file__).resolve().parent / "results"
MDOTS = np.arange(1.5, 4.01, 0.5)
S2GHz_RANGE = (105.0, 135.0)
OUTPUT_FILES = ("rad_syn_WCR.dat", "windemiff.dat", "windemiff1.dat", "windemiff2.dat")
GG_TABLES = (ROOT / "output" / "opac_gg1.dat", ROOT / "output" / "opac_gg2.dat")


# Replace one Fortran real parameter while preserving the source structure.
def replace_parameter(text, name, value):
    pattern = rf"^\s*real\(dp\),\s*parameter\s*::\s*{name}\b.*$"
    updated, count = re.subn(pattern, value, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"Could not update parameter {name}")
    return updated


# Replace one logical control while preserving the surrounding source text.
def set_control(text, name, value):
    pattern = rf"^(\s*logical,\s*parameter\s*::\s*{name}\s*=\s*).*$"
    updated, count = re.subn(pattern, rf"\g<1>{value}", text, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"Could not update control {name}")
    return updated


# Format a Python value using Fortran double-precision exponent notation.
def fortran_real(value):
    return f"{value:.8e}".replace("e", "d")


# Read a Fortran double-precision parameter from the source text.
def parameter_value(text, name):
    match = re.search(
        rf"parameter\s*::\s*{name}\s*=\s*([0-9.]+)d([+-]?[0-9]+)", text
    )
    if match is None:
        raise RuntimeError(f"Could not read parameter {name}")
    return float(f"{match.group(1)}e{match.group(2)}")


# Build and execute one model case, returning its terminal output.
def run_model():
    subprocess.run(["make", "run.x"], cwd=ROOT, check=True, capture_output=True, text=True)
    return subprocess.run(
        ["./run.x"], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout


# Copy the model outputs needed by the radio SED figure into the results folder.
def save_outputs(suffix):
    for filename in OUTPUT_FILES:
        source = ROOT / "output" / filename
        if not source.exists():
            raise RuntimeError(f"Expected model output is missing: {source}")
        shutil.copy2(source, RESULTS / f"{source.stem}_{suffix}{source.suffix}")


# Extract the post-convolution 2 GHz flux from one model run.
def extract_s2(output):
    match = re.search(r"Post-conv:.*?S\(2GHz\)=\s+([0-9.Ee+-]+)", output, flags=re.DOTALL)
    if match is None:
        raise RuntimeError("Could not find post-convolution S(2GHz) in model output")
    return float(match.group(1))


# Run the six-case Mdot exploration and write its compact summary.
def main():
    RESULTS.mkdir(exist_ok=True)
    original_system = SYSTEM_PARAMETERS.read_text()
    original_controls = CONTROLS.read_text()
    missing_tables = [str(path) for path in GG_TABLES if not path.exists()]
    if missing_tables:
        raise RuntimeError(
            "The Mdot scan requires gamma-gamma tables from the baseline run: "
            + ", ".join(missing_tables)
        )

    mdot0 = parameter_value(original_system, "Mdot1")
    frac_nt0 = parameter_value(original_system, "frac_NT")
    eta_b0 = parameter_value(original_system, "eta_B1")
    constant_nt = mdot0 * frac_nt0

    for pattern in ("rad_syn_WCR_M*.dat", "windemiff_M*.dat", "windemiff1_M*.dat", "windemiff2_M*.dat"):
        for path in RESULTS.glob(pattern):
            path.unlink()

    summary = []
    try:
        controls = set_control(original_controls, "rad_ic1", "off")
        controls = set_control(controls, "rad_ic2", "off")
        controls = set_control(controls, "absgg", "off")
        controls = set_control(controls, "protons", "off")
        CONTROLS.write_text(controls)

        for mdot in MDOTS:
            mdot_physical = mdot * 1.0e-5
            frac_nt = constant_nt / mdot_physical
            eta_b = eta_b0 * (mdot_physical / mdot0) ** (-0.98)
            text = replace_parameter(
                original_system,
                "Mdot1",
                f"   real(dp), parameter :: Mdot1 = {mdot:.1f}d-5*M_sun_yr",
            )
            text = replace_parameter(
                text,
                "frac_NT",
                f"   real(dp), parameter :: frac_NT = {fortran_real(frac_nt)}",
            )
            text = replace_parameter(
                text,
                "eta_B1",
                f"   real(dp), parameter :: eta_B1 = {fortran_real(eta_b)}, eta_B2 = eta_B1",
            )
            SYSTEM_PARAMETERS.write_text(text)

            output = run_model()
            suffix = f"M{mdot:.1f}"
            save_outputs(suffix)
            s2 = extract_s2(output)
            in_range = S2GHz_RANGE[0] <= s2 <= S2GHz_RANGE[1]
            if not in_range:
                print(f"WARNING: {suffix} S(2GHz) = {s2:.3f} mJy, expected {S2GHz_RANGE}")
            summary.append({"Mdot1_1e-5": mdot, "frac_NT": frac_nt, "eta_B1": eta_b, "S_2GHz_mJy": s2, "in_expected_range": in_range})
    finally:
        SYSTEM_PARAMETERS.write_text(original_system)
        CONTROLS.write_text(original_controls)
        subprocess.run(["make", "run.x"], cwd=ROOT, check=True)

    (RESULTS / "mdot_scan_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Completed {len(summary)} Mdot cases.")


if __name__ == "__main__":
    main()
