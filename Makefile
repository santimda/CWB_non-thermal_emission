DEBUG  = N

FC     = gfortran -fopenmp
FFLAGS = -ffree-line-length-none -fbackslash
LIBS   = -llapack -lblas
BUILD_DIR = build

ifeq ($(DEBUG), Y)
  FFLAGS += -Wall -Wno-tabs -Wextra -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
else
  FFLAGS += -Wall -Wno-tabs -Wextra -fPIC -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops
#  FFLAGS += -Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-none -fbackslash
endif

# Production sources, ordered from low-level modules to the program entry point.
SOURCES = constants.f90 \
          controls.f90 \
          system_parameters.f90 \
          global.f90 \
          initialize.f90 \
          orbit.f90 \
          CD.f90 \
          CD_thermo.f90 \
          opac_ff.f90 \
          opac_gg.f90 \
          BB.f90 \
          WCR_thermal.f90 \
          CD_dist_e.f90 \
          CD_dist_p.f90 \
          CD_rot.f90 \
          CD_rad_syn.f90 \
          CD_rad_ic.f90 \
          CD_rad_br.f90 \
          CD_rad_pp.f90 \
          CD_abssyn.f90 \
          fluxes.f90 \
          WCR_conv.f90 \
          main.f90
OBJECTS = $(SOURCES:%.f90=$(BUILD_DIR)/%.o)

APEP_PYTHON ?= python3

.PHONY: clean test test-unit apep

run.x: $(OBJECTS) | output
	$(FC) $(OBJECTS) $(LIBS) -o run.x

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

output:
	mkdir -p output

.SUFFIXES:
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) -c $(FFLAGS) -J$(BUILD_DIR) -o $@ $<

test: test-unit
	python3 tests/test_terminal_output.py

test-unit: build/test_numerical_utilities
	./build/test_numerical_utilities

apep: run.x
	$(APEP_PYTHON) python_scripts/check_apep_baseline.py
	$(APEP_PYTHON) python_scripts/explore_parameter.py
	$(APEP_PYTHON) -m python_scripts.plot_radio_seds

build/test_numerical_utilities: tests/test_numerical_utilities.f90 build/constants.o build/controls.o build/system_parameters.o build/global.o | build
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) $^ -o $@

clean:
	rm -f *.x *.o *.mod *.bin *.log *~
	rm -rf ./build ./output

include Makefile.depend
