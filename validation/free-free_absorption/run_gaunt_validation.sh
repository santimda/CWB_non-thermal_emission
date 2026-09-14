#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
VALIDATION_DIR="$ROOT_DIR/validation/free-free_absorption"
BUILD_DIR="$ROOT_DIR/build"

cd "$ROOT_DIR"
mkdir -p "$VALIDATION_DIR/results"
make build/constants.o build/controls.o build/system_parameters.o build/global.o build/opac_ff.o

gfortran -ffree-line-length-none -fbackslash -J"$BUILD_DIR" -I"$BUILD_DIR" \
  "$VALIDATION_DIR/gaunt_driver.f90" \
  "$BUILD_DIR/constants.o" "$BUILD_DIR/controls.o" \
  "$BUILD_DIR/system_parameters.o" "$BUILD_DIR/global.o" "$BUILD_DIR/opac_ff.o" \
  -o "$VALIDATION_DIR/gaunt_driver.x"

"$VALIDATION_DIR/gaunt_driver.x"
python3 "$VALIDATION_DIR/plot_gaunt_factors.py"
