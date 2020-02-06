## Test environments
* local OS X install, R 3.6.2
* ubuntu 12.04 (on travis-ci), R 3.6.2
* fedora 28, R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is an update for the previous CRAN release and is not backwards compatible

* The note is due to the use of GNU make, which is requried by RcppParallel. As far as I know, this cannot be remedied and is otherwise    not an issue. Please see https://github.com/RcppCore/RcppParallel/issues/23 for why

* Sadly, the package fails on Debian Linux, R-devel, GCC ASAN/UBSAN. This seems to be an issue that is caused by TBB which is used by      RcppParallel (a dependency of this package). According to this issue https://github.com/RcppCore/RcppParallel/issues/36 there is         currently no workaround, as the TBB team has to fix this. Since RcppParallel has this issue too and has been accepted to CRAN, I think   the error can be ignored for this package as well.

* Before this package had trouble building on Solaris. Therefore this package has been checked using rhub::check_on_solaris(), Oracle Solaris 10, x86, 32 bit, R-patched (experimental). This gave 0 errors, 0 warnings and 0 notes. 

## Reverse dependencies

None

