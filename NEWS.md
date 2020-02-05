subgroup.discovery v0.3.0 (Release date: 2020-01-30)
==============

Changes:

* Moved all computation to c++ code, resulting package is now much faster
* Finding the optimal peel is done in a parallel map-reduce fashion
* There is a dedicated function called prim.data.prepare() to facilitate the proper setup of input data
* The covering strategies have been removed, package is now in line with common practices, i.e. you build a model using prim() and then evaluate it using predict().
* To find observations which fall into a certain box, use the prim.box.index() function


subgroup.discovery v0.3.1 (Release date: 2020-02-06)
==============

Changes:

* Removed all namespaces except Rcpp to comply with best practices (especially not using std)
* Changed the pattern of sub-box creation so at no point there is a uninstantiated value

