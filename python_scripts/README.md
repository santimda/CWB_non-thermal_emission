# Apep Mdot Exploration

Run the scripts from the project root. The workflow requires Python 3, NumPy,
Pandas, SciPy, and Matplotlib. LaTeX is used by Matplotlib for the publication
fonts; set `text.usetex` to `False` in `python_scripts/plot_utils.py` if LaTeX
is not installed.

1. Check the default Apep normalization with IC and gamma-gamma table
   generation enabled:

   ```text
   python python_scripts/check_apep_baseline.py
   ```

   This writes `baseline_constraints.json` and `baseline_run.log` under
   `python_scripts/results/`. The checks are
   `F_10-30 keV = (3.9 +/- 1.2)e-13 erg/s/cm2` and
   `S(2GHz) = 120 +/- 15 mJy`.

2. Run the six-case Mdot exploration:

   ```text
   python python_scripts/explore_parameter.py
   ```

   The scan covers `Mdot1 = 1.5, 2.0, ..., 4.0` in units of
   `1e-5 Msun/yr`, keeps the constant-X-ray scaling of `frac_NT` and
   `eta_B1`, and runs with `rad_ic1=off` and `rad_ic2=off`. All cases are
   retained; values outside `105-135 mJy` at 2 GHz generate warnings. The
   scan stores a compact JSON summary rather than one log file per run.

3. Run the radio SED plotting script:

   ```text
   python -m python_scripts.plot_radio_seds
   ```

   The script prints the fitted spectral indices for the ALMA B3, B6, and B7
   bands, prints the reduced chi-squared for each Mdot model using the plotted
   radio data and total uncertainties, then writes
   `python_scripts/results/SEDs_Mdot.pdf`.

4. Plot the absorbed, beam-convolved radio maps:

   ```text
   python -m python_scripts.plot_maps
   ```

   The script discovers the available convolved maps in `output/` and writes
   `python_scripts/results/maps_absorbed.pdf`. It supports one, two, or three
   available frequencies. Each frequency has its own color scale and colorbar,
   with contours at 3%, 10%, 30%, 50%, 70%, and 90% of that map's maximum. The displayed beam uses
   `sigma_x=5.6 mas` and `sigma_y=11.3 mas`, matching `WCR_conv.f90`.

The scripts restore `system_parameters.f90` and `controls.f90` after each
workflow. The baseline run generates the gamma-gamma tables and restores the
   source configuration. The Mdot scan runs with `absgg=off`, IC disabled, and
   `protons=off`. Plotting style settings are shared
with the validation figures through `python_scripts/plot_utils.py`.
