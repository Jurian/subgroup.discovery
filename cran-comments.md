## Test environments
* local OS X install, R 3.4.0
* ubuntu 12.04 (on travis-ci), R 3.4.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is an update for the previous CRAN release and is not backwards compatible

* The note is due to the use of GNU make, which is requried by RcppParallel
* As far as I know, this cannot be remedied and is otherwise not an issue
* Please see https://github.com/RcppCore/RcppParallel/issues/23 for why

## Reverse dependencies

None

