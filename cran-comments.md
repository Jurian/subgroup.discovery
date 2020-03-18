## Test environments
* local OS X install, R 3.6.2
* ubuntu 12.04 (on travis-ci), R 3.6.2
* fedora 28, R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes

## Version 0.3

* This is a major update for release 0.2.x and is not backwards compatible

* The note is due to the use of GNU make, which is requried by RcppParallel. As far as I know, this cannot be remedied and is otherwise not an issue. Please see https://github.com/RcppCore/RcppParallel/issues/23.

## Version 0.3.2

* The minimal R version has been set to 3.0.0, so this package can also be tested on some architectures on CRAN which use an older version

* The issue with ASAN / UBSAN warnings in CRAN submissions should be resolved thanks to update 5.0.0 from RcppParallel package

## Reverse dependencies

None

