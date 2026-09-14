#!/usr/bin/env python3
"""Regression test for the numerical summary printed by the default run."""

import json
import math
import pathlib
import re
import shutil
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]
REFERENCE = pathlib.Path(__file__).with_name("reference_terminal_output.json")
CONFIGURATION = pathlib.Path(__file__).parent / "cases" / "default"


PATTERNS = {
    "synchrotron_unabsorbed": r"L_syn\(unabs\s+([0-9.Ee+-]+)",
    "synchrotron_ssa": r"L_syn_SSA=\s+([0-9.Ee+-]+)",
    "ic_unabsorbed": r"L_IC\(un\)=\s+([0-9.Ee+-]+)",
    "ic_absorbed": r"L_IC\(abs\)=\s+([0-9.Ee+-]+)",
    "bremsstrahlung": r"L_Brem=\s+([0-9.Ee+-]+)",
    "thermal_bolometric": (
        r"L_X bolometric \(S1\)\s*=\s+([0-9.Ee+-]+)"
        r".*?L_X bolometric \(S2\)\s*=\s+([0-9.Ee+-]+)"
    ),
    "unabsorbed_flux": (
        r"Unabsorbed fluxes.*?F_3-10keV.*?=\s+([0-9.Ee+-]+)"
        r".*?F_10-30keV.*?=\s+([0-9.Ee+-]+)"
        r".*?F_0\.1-100GeV.*?=\s+([0-9.Ee+-]+)"
        r".*?F_0\.1-10TeV.*?=\s+([0-9.Ee+-]+)"
    ),
    "absorbed_flux": (
        r"Absorption-corrected fluxes.*?F_3-10keV.*?=\s+([0-9.Ee+-]+)"
        r".*?F_10-30keV.*?=\s+([0-9.Ee+-]+)"
        r".*?F_0\.1-100GeV.*?=\s+([0-9.Ee+-]+)"
        r".*?F_0\.1-10TeV.*?=\s+([0-9.Ee+-]+)"
    ),
    "radio_flux": r"S\(2GHz\)=\s+([0-9.Ee+-]+) mJy",
    "star1_free_free": r"S1_ff\(2GHz\)\s+([0-9.Ee+-]+) mJy",
    "star2_free_free": r"S2_ff\(2GHz\)\s+([0-9.Ee+-]+) mJy",
}


def values(pattern, output):
    matches = re.findall(pattern, output, flags=re.DOTALL)
    return [float(value) for match in matches for value in (match if isinstance(match, tuple) else (match,))]


def main():
    reference = json.loads(REFERENCE.read_text())
    with tempfile.TemporaryDirectory(prefix="cwb-regression-") as temporary_directory:
        test_root = pathlib.Path(temporary_directory)
        for source in ROOT.glob("*.f90"):
            shutil.copy2(source, test_root / source.name)
        for filename in ("Makefile", "Makefile.depend"):
            shutil.copy2(ROOT / filename, test_root / filename)
        shutil.copy2(ROOT / "gauntff.dat", test_root / "gauntff.dat")
        for filename in ("controls.f90", "system_parameters.f90", "global.f90"):
            shutil.copy2(CONFIGURATION / filename, test_root / filename)

        subprocess.run(
            ["make"],
            cwd=test_root,
            check=True,
            capture_output=True,
            text=True,
        )
        result = subprocess.run(
            ["./run.x"],
            cwd=test_root,
            check=True,
            capture_output=True,
            text=True,
        )

    failures = []
    for name, pattern in PATTERNS.items():
        actual = values(pattern, result.stdout)
        expected = reference[name]
        if len(actual) != len(expected):
            failures.append(f"{name}: expected {len(expected)} values, found {len(actual)}")
            continue
        for index, (actual_value, expected_value) in enumerate(zip(actual, expected), start=1):
            if not math.isclose(actual_value, expected_value, rel_tol=2e-3, abs_tol=1e-30):
                failures.append(
                    f"{name}[{index}]: expected {expected_value:g}, found {actual_value:g}"
                )

    if failures:
        print("Terminal-output regression failed:")
        print("\n".join(f"- {failure}" for failure in failures))
        return 1

    print("Terminal-output regression passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
