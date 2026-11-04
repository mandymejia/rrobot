## Test environments

* Mac aarch64-apple-darwin23, R 4.6.1

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

None.

## Tests

Passes all the tests in `tests/testthat.R`.

## Previous submission

Flavor: r-devel-linux-x86_64-debian-gcc
Check: Rd contents, Result: NOTE
  Rd files without \usage:
    'B.Rd' 'M.Rd' 'RD_obj.Rd' 'RD_org_obj.Rd' 'alpha.Rd' 'binwidth.Rd'
    'boot_quant.Rd' 'cov_mcd.Rd' 'cutoff.Rd' 'dist.Rd' 'dot-...Rd'
    'imp_data.Rd' 'imp_datasets.Rd' 'impute_method.Rd' 'ind_incld.Rd'
    'k.Rd' 'method.Rd' 'method_univOut.Rd' 'mode.Rd' 'quantile.Rd'
    'thr.Rd' 'thresh_result.Rd' 'trans.Rd' 'verbose.Rd' 'w.Rd' 'x.Rd'
  \arguments should not be documented without \usage.

* These have been fixed.