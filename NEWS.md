# TmCalculator 1.0.9

## Bug fixes

* Corrected two entropy values in `DNA_NN_SantaLucia_2004`. The `TA/AT` step
  was `-20.4` (duplicated from the adjacent `AT/TA` row) and should be `-21.3`;
  the `GG/CC` step was `-19.0` and should be `-19.9`. Both values are shared
  with `DNA_NN_Allawi_1998`, which already carried them correctly, so the two
  tables were internally inconsistent. Tm values computed with
  `DNA_NN_SantaLucia_2004` will change slightly for sequences containing these
  steps.

## New parameter sets

* Added 19 nearest-neighbor parameter sets derived by melting-temperature
  optimization (Weber and colleagues, UFMG), bringing the total from 8 to 27.
  Unlike the existing sets, these are fitted directly at a stated sodium
  concentration and are intended to replace salt correction rather than be
  corrected.

  - DNA: `DNA_NN_Weber_2015` (1020 mM) and a salt series
    `DNA_NN_Weber_OW04_69`, `_119`, `_220`, `_621`, `_1020`.
  - RNA: `RNA_NN_Weber_VIF_*` and `RNA_NN_Weber_FIF_*` at 71, 121, 221, 621
    and 1021 mM. The VIF (variable initiation factor) sets gave better
    cross-validation in the source study.
  - RNA/DNA hybrid: `RNA_DNA_NN_Weber_2019_FT` and `_VH` (1000 mM) and
    `RNA_DNA_NN_Weber_2019_LS` (100 mM).

  Values were taken from the parameter files distributed with VarGibbs 5.0 at
  full precision rather than transcribed from the published tables. The
  transcription was validated by confirming that the reference files shipped
  alongside them reproduce the existing `DNA_NN_Allawi_1998`,
  `DNA_NN_Sugimoto_1996`, `RNA_NN_Xia_1998`, `RNA_NN_Freier_1986` and all 87
  rows of `DNA_IMM_Peyret_1999` exactly.

## Salt handling

* Parameter sets now carry a `salt_mM` attribute when they were fitted at a
  specific sodium concentration. `tm_nn()` uses it to avoid double-counting the
  ionic contribution:

  - if the requested `Na` matches the concentration the set was fitted at,
    salt correction is skipped;
  - if it does not, the correction is applied and a warning names both
    concentrations.

  Existing parameter sets carry no such attribute, so this is a no-op for all
  previous usage.

* `salt_method` gains a `"none"` option to disable salt correction explicitly.
  The documentation previously stated that `NA` would do this, but
  `match.arg()` rejected it.

* `tm_nn()` and `tm_calculate()` now report `Salt correction applied` (logical)
  and `Parameter set fitted at [Na+] (mM)` in the returned `options`, so the
  automatic skip is visible rather than silent.

## Documentation

* `tm_nn()` gained a `@return` section; the return value was previously
  undocumented.
* Removed `"SantaLucia1998-2"` from the documented `salt_method` options in
  `tm_nn()` and `tm_calculate()`. It was listed in the help pages but absent
  from the function's accepted values, so following the documentation produced
  an error.
* Added a "Choosing a parameter set" section to `?tm_nn` and a "Salt handling"
  section to `?tm_calculate`.
* `@details` for `tm_calculate()` now explains when each of the three methods
  is appropriate, including the note that nearest-neighbor parameters are
  calibrated on short duplexes, so values computed over long sequences or
  fixed-width genomic windows are best read as a relative measure of local
  thermodynamic stability rather than as an absolute experimental Tm.
* `tm_nn()` and `tm_calculate()` now cross-reference each other via `@seealso`.