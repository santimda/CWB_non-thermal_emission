#!/usr/bin/env bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$CASE_DIR/../.." && pwd)"
PYTHON="${PYTHON:-python3}"

run_case() {
  local work="$CASE_DIR/work"

  rm -rf "$work"
  mkdir -p "$work"
  cp "$ROOT"/*.f90 "$ROOT"/Makefile "$ROOT"/Makefile.depend "$ROOT"/gauntff.dat "$work"/
  cp "$CASE_DIR/ffa_map_driver.f90" "$work"/

  sed -i \
    -e 's/integer, parameter :: ml = 200/integer, parameter :: ml = 100/' \
    -e 's/integer, parameter :: mnu = 800/integer, parameter :: mnu = 3/' \
    -e 's/integer, parameter :: mnu_data = 4/integer, parameter :: mnu_data = 3/' \
    "$work/system_parameters.f90" "$work/global.f90"

  sed -i \
     -e 's/logical, parameter :: sed_all_calc = on/logical, parameter :: sed_all_calc = off/' \
     -e 's/logical, parameter :: radio_calc = off/logical, parameter :: radio_calc = on/' \
    "$work/controls.f90"

  sed -i \
    -e 's/nu_data = (\/1.4d9, 2.2d9, 5.5d9, 9.d9\/)/nu_data = (\/0.3d9, 1.d9, 5.d9\/)/' \
    "$work/main.f90"

  sed -i \
    -e 's/call vector_log(nu_min_Hz, nu_max_Hz, nu_int, nu)/nu = (\/0.3d9, 1.d9, 5.d9\/)/' \
    -e 's/real(dp) :: nuc_max, nuc_min, nu_int, nu_min_Hz, nu_max_Hz/real(dp) :: nuc_max, nuc_min, nu_min_Hz, nu_max_Hz/' \
    "$work/CD_rad_syn.f90"

  # The map driver overwrites x(1:ml), so use the analytic CD vertex.
  sed -i 's/x_v = x(1)/x_v = D\/(1.d0 + sqrt(eta))/g' "$work/CD_abssyn.f90"

  make -C "$work"

  local objects=()
  local object
  for object in "$work"/build/*.o; do
    if [[ "$object" == "$work/build/main.o" ]]; then
      continue
    fi
    objects+=("$object")
  done
  gfortran -fopenmp -ffree-line-length-none -fbackslash -O3 -J"$work/build" \
    -I"$work/build" "$work/ffa_map_driver.f90" "${objects[@]}" \
    -llapack -lblas -o "$work/ffa_maps.x"
  (cd "$work" && ./ffa_maps.x > validation_run.log)
  "$PYTHON" "$CASE_DIR/plot_ffa.py" \
    --work-dir "$work"
  rm -rf "$work"
}

run_case
