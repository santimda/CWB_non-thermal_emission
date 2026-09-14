# Free-Free Absorption Validation

This validation uses the production Fortran pipeline in an isolated work
directory. It is intended to expose geometric behavior to help in 
the interpretation of the f-f absorption dependence on the observing angle.

The validation case uses the Apep parameters from `system_parameters.f90`:

- Star 1: `Mdot1 = 2.5e-5 Msun/yr`, `vinf1 = 3350 km/s`.
- Star 2: `vinf2 = 2350 km/s`, with `Mdot2` derived from `eta`.
- `eta = 0.32`.

The maps use `0.3`, `1`, and `5 GHz`, with observer angles `psi=75` and
`psi=105` degrees. The map ranges are `x/D=[-1.5,2.5]` and `y/D=[-2,2]`.

Run from the project root:

```text
validation/free-free_absorption/run_validation.sh
```

The script copies the production sources into `work/`, changes only the
validation parameters in that copy, builds a three-frequency (`mnu=3`)
validation executable with `ml=100`, and calls the production
`total_tau_calc` routine once for each map row. The generated Fortran
transmission maps are then used to produce
`free-free_transmission.pdf`. 
The temporary `work/` directory is removed after the plot is generated.

The map deriver uses `z=0`, a `100x100` grid, and the `crossing_point` root selection.
The crossing point is derived by approximating the wind-interaction geometry analytically as the
positive-x branch of a hyperbola.
A parabolic approximation was explored during development, but the hyperbola
is a better representation of the numerical CD. This
approximation is used only to determine which stellar wind absorbs the
synchrotron emission along its path.

For the plot, points with `y < 0` between the positive-x
hyperbola branch and the line of sight through `(x_v, 0)` are masked. In the
plotted plane, `x_line = x_v + y/abs(tan(psi))` for `psi < 90` degrees and
`x_line = x_v - y/abs(tan(psi))` for `psi > 90` degrees. This is a
presentation-only mask, as the Fortran opacity calculations assume photons 
emitted from the CD.
The validation work directory is ignored by git.

## Gaunt-factor validation

Run the production Gaunt-factor diagnostic with:

```text
validation/free-free_absorption/run_gaunt_validation.sh
```

This calls `gaunt_ff_calc` for the two wind temperatures and charges from
`system_parameters.f90` on 50 logarithmically spaced frequencies from
`1e8` to `1e12 Hz`. The Fortran values are written to
`results/gaunt_fortran.dat`, and `plot_gaunt_factors.py` compares them with
the analytical approximation documented in `opac_ff.f90` in
`results/gaunt_factors.pdf`. The agreement between the tabulated gaunt 
factor and analytical approximations is excellent, although at submm frequencies
deviations start to occur, for which it is safest to use the tabulated values.

