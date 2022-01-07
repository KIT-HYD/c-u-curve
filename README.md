# c-u-curve

This repository contains code for calculating the c-u-curve and some example applications.

The c-u-curve method can be used to analyse, classify and compare dynamical systems of arbitrary dimension by the two key features uncertainty and complexity.
The method is applicable to uni- and multivariate data sets, and accepts both deterministic and probabilistic value representations.

For detail information see the related publication by Ehret and Dej (2022).
ADD REF

## Requisites

MATLAB 9.9

## Usage

See test_c_u_curve.m

## Files

* test_c_u_curve.m            contains example applications for single- and multivariate data sets, and for deterministic and probabilisitc value representations
* f_c_u_curve.f               function returning all values to plot the c_u_curve (uncertainty and complexity values for user-selected time slicing schemes)
* f_entropy_anyd_fast         function returning the (joint) entropy of an 1-to-any-dimensional discrete (binned) frequency distribution

## Contact

Uwe Ehret | uwe.ehret@kit.edu
