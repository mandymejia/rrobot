## Test environments

* Mac x86_64-apple-darwin17.0, R 4.5.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

None.

## Tests

Passes all the tests in `tests/testthat.R`.

## Previous submission

  If there are references describing the methods in your package, please
  add these in the description field of your DESCRIPTION file in the form
  authors (year) <doi:...>
  authors (year, ISBN:...)
  or if those are not available: <https:...>
  with no space after 'doi:', 'https:' and angle brackets for
  auto-linking. (If you want to add a title as well please put it in
  quotes: "Title")
  For more details:
  <https://contributor.r-project.org/cran-cookbook/description_issues.html#references>

> Added a reference of our pre-print to the end of the DESCRIPTION. 

  Please add \value to .Rd files regarding exported methods and explain
  the functions results in the documentation. Please write about the
  structure of the output (class) and also what the output means. (If a
  function does not return a value, please document that too, e.g.
  \value{No return value, called for side effects} or similar)
  For more details:
  <https://contributor.r-project.org/cran-cookbook/docs_issues.html#missing-value-tags-in-.rd-files>
  Missing Rd-tags:
        summary.F_result.Rd: \value
        summary.MI_boot_result.Rd: \value
        summary.MI_result.Rd: \value
        summary.SI_boot_result.Rd: \value
        summary.SI_result.Rd: \value

> All exported methods have \value documentation now. 

  Please add small executable examples in your Rd-files to illustrate the
  use of the exported function but also enable automatic testing.

> Done! 

  Please fix and resubmit.