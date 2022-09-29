# c-u-curve

This repository contains code for calculating the c-u-curve and some example applications.

The c-u-curve method can be used to analyse, classify and compare dynamical systems of arbitrary dimension by the two key features uncertainty and complexity.
The method is applicable to uni- and multivariate data sets, and accepts both deterministic and probabilistic value representations.

For detail information see the related publication by Ehret and Dey (2022):
Ehret, U., and P. Dey (2022), Technical note: c-u-curve: A method to analyse, classify and compare dynamical systems by uncertainty and complexity, Hydrol. Earth Syst. Sci. Discuss., 2022, 1-12.


## Requisites

MATLAB 9.9

## Usage

See test_c_u_curve.m

## Files

* test_c_u_curve.m            contains example applications for single- and multivariate data sets, and for deterministic and probabilisitc value representations
* f_c_u_curve.f               function returning all values to plot the c_u_curve (uncertainty and complexity values for user-selected time slicing schemes)
* f_entropy_anyd_fast         function returning the (joint) entropy of an 1-to-any-dimensional discrete (binned) frequency distribution
* c_u_curve_application_EhretDey2022.m contains all code to reproduce the results in Ehret and Dey (2022)
* c_u_curve_application_EhretDey2022.mat contains all data used in c_u_curve_application_EhretDey2022.m

## Contact

Uwe Ehret | uwe.ehret@kit.edu
