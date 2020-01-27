## Subgroup Discovery

## Motivation

This package was developed to assist in discovering interesting subgroups in multi-dimensional data. 

## Description

The PRIM implementation is based on the 1998 paper "Bump hunting in high-dimensional data" by Jerome H. Friedman and Nicholas I. Fisher. PRIM involves finding a set of "rules" which combined imply unusually large (or small) values of some other target variable. Specifically one tries to find a set of sub-regions in which the target variable is substantially larger than overall mean. 

The objective of bump hunting in general is to find regions in the input (attribute/feature) space with relatively high (low) values for the target variable. The regions are described by simple rules of the type if: {condition-1 & ... & condition-n} then: estimated target value. Given the data (or a subset of the data), the goal is to produce a box B within which the target mean is as large as possible. There are many problems where finding such regions is of considerable practical interest.  

## Contributors

Developed by Jurian Baas, part of master thesis in Artificial Intelligence.

Special thanks for contributions such as suggestions and bug-fixing to Dr. A.J. Feelders, Utrecht University, Department of Information and Computing Sciences

## References

1. Friedman, Jerome H., and Nicholas I. Fisher. "Bump hunting in high-dimensional data." Statistics and Computing 9.2 (1999): 123-143.

## Citation

To cite this package, use citation("subgroup.discovery") in R

## Licence

This package is licenced under GPL-3

    Copyright (C) 2020  Jurian Baas

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/.

