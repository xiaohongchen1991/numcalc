# Release versions

version 1.0.0:

Provide a set of simple C++ APIs to do 1D or 2D numerical differentiation using finite difference method.

# Quick start

This is a simple C++ library for numerical differentiation using finite difference method. It supports taking partial derivative of 1D and 2D functions and can be extended into 3D as well.

## Basics

    #include "differentiation.hpp"
    
Given 1D function f(x) = exp(x), to compute its first derivative at x = 0, do:

    double f_x = numcalc::diff<numcal::f_x>(exp, 0.0);

For its second derivative, do

    double f_xx = numcalc::diff<numcal::f_xx>(exp, 0.0);

The computation of the partial derivatives of 2D functions is also supported. Given a 2D function f(x, y),

    auto f = [](double x, double y){return exp(x + 2.0*y);};

the first order derivatives at coordinates (1, 1) can be computed by

    double f_x = numcalc::diff<numcalc::f_x>(f, 1.0, 1.0);
    double f_y = numcalc::diff<numcalc::f_y>(f, 1.0, 1.0);

and the second order derivatives can be computed by

    double f_xx = numcalc::diff<numcalc::f_xx>(f, 1.0, 1.0);
    double f_xy = numcalc::diff<numcalc::f_xy>(f, 1.0, 1.0);
    double f_yy = numcalc::diff<numcalc::f_yy>(f, 1.0, 1.0);

## Linking

This is a header only library.

## Requirements

The library only requires a C++ compiler that supports C++14.

## TODO list

* Extend diff() to 3D cases
* Add numerical integration