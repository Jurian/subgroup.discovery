subgroup.discovery v0.3.0 (Release date: 2020-01-30)
==============

Changes:

* Moved all computation to c++ code, resulting package is now much faster
* Finding the optimal peel is done in a parallel map-reduce fashion
* There is a dedicated function called prim.data.prepare() to facilitate the proper setup of input data
* The covering strategies have been removed, package is now in line with common practices, i.e. you build a model using prim() and then evaluate it using predict().
* To find observations which fall into a certain box, use the prim.box.index() function

