/*
MIT License

Copyright (c) 2019 Xiaohong Chen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef NUMCALC_DIFFERENTIATION_HPP
#define NUMCALC_DIFFERENTIATION_HPP

/*
 * Performs second-order finite central difference
 */

#include "numcalc_config.hpp"

#include <cmath>
#include <limits>
#include <iostream>
#include <cstdlib>

namespace numcalc {

enum diff_pattern{f_x, f_y, f_z, f_xx, f_xy, f_yy};

namespace detail {
    template<class Real>
    Real adjust_step_size(Real x, Real h)
    {
        using std::numeric_limits;
        // Redefine h so that x + h is representable. Not using this trick leads to large error.
        // The compiler flag -ffast-math evaporates these operations . . .
        Real temp = x + h;
        h = temp - x;
        // Handle the case x + h == x:
        if (h == 0)
        {
            h = std::nextafter(x, (numeric_limits<Real>::max)()) - x;
        }
        return h;
    }

    template <diff_pattern>
    struct diff_tag {};

    // second-order central difference for f_x.
    // f_x \approx (f[x+h] - f[x-h])/(2*h)
    // h ~ (3*eps)^1/3
    template<class F, class Real>
    Real diff(const F f, Real x, const diff_tag<f_x>&)
    {
        using std::sqrt;
        using std::pow;
        using std::abs;
        using std::numeric_limits;

        const Real eps = (numeric_limits<Real>::epsilon)();
        Real h = pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));
        h = adjust_step_size(x, h);

        Real d = f(x + h) - f(x - h);

        return d / (2 * h);
    }

    // second-order central difference for f_xx.
    // f_xx \approx (f[x+h] - 2f[x] + f[x-h])/(h*h)
    // h ~ (48*eps)^1/4
    template<class F, class Real>
    Real diff(const F f, Real x, const diff_tag<f_xx>&)
    {
        using std::sqrt;
        using std::pow;
        using std::abs;
        using std::numeric_limits;

        const Real eps = (numeric_limits<Real>::epsilon)();
        Real h = pow(48 * eps, static_cast<Real>(1) / static_cast<Real>(4));
        h = adjust_step_size(x, h);

        Real d = f(x + h) - 2.0*f(x) + f(x - h);

        return d / (h * h);
    }

    // numerical differentiation for f_x(x, y)
    template<class F, class Real>
    Real diff(const F f, Real x, Real y, const diff_tag<f_x>&)
    {
        return diff([f, y](Real x){return f(x, y);}, x, diff_tag<f_x>());
    }

    // numerical differentiation for f_y(x, y)
    template<class F, class Real>
    Real diff(const F f, Real x, Real y, const diff_tag<f_y>&)
    {
        return diff([f, x](Real y){return f(x, y);}, y, diff_tag<f_x>());
    }

    // second-order central difference for f_xy.
    // f_xy \approx (f[x+h, y+h] - f[x+h, y-h] - f[x-h, y+h] + f[x-h, y-h])/(4*h*h)
    // h ~ (3*eps)^1/3
    template<class F, class Real>
    Real diff(const F f, Real x, Real y, const diff_tag<f_xy>&)
    {
        using std::sqrt;
        using std::pow;
        using std::abs;
        using std::numeric_limits;

        const Real eps = (numeric_limits<Real>::epsilon)();
        Real h = pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(4));
        h = adjust_step_size(x, h);

        Real d = f(x + h, y + h) - f(x - h, y + h) - f(x + h, y - h) + f(x - h, y - h);

        return d / (4 * h * h);
    }

    // numerical differentiation for f_xx(x, y)
    template<class F, class Real>
    Real diff(const F f, Real x, Real y, const diff_tag<f_xx>&)
    {
        return diff([f, y](Real x){return f(x, y);}, x, diff_tag<f_xx>());
    }

    // numerical differentiation for f_yy(x, y)
    template<class F, class Real>
    Real diff(const F f, Real x, Real y, const diff_tag<f_yy>&)
    {
        return diff([f, x](Real y){return f(x, y);}, y, diff_tag<f_xx>());
    }

    // numerical differentiation for f_x(x, y, z)
    template<class F, class Real>
    Real diff(const F f, Real x, Real y, Real z, const diff_tag<f_x>&)
    {
        return diff([f, y, z](Real x){return f(x, y, z);}, x, diff_tag<f_x>());
    }

    // overloadings for undefined tags
    template<class F, class Real, class Tag>
    Real diff(const F f, Real x, const Tag&)
    {
        std::cerr << "Unrecognized differentiation pattern!" << "\n";
        std::exit(EXIT_FAILURE);
    }

    template<class F, class Real, class Tag>
    Real diff(const F f, Real x, Real y, const Tag&)
    {
        std::cerr << "Unrecognized differentiation pattern!" << "\n";
        std::exit(EXIT_FAILURE);
    }

    template<class F, class Real, class Tag>
    Real diff(const F f, Real x, Real y, Real z, const Tag&)
    {
        std::cerr << "Unrecognized differentiation pattern!" << "\n";
        std::exit(EXIT_FAILURE);
    }
} // detail

template<diff_pattern P, class F, class Real>
inline Real diff(const F f, Real x)
{
    return detail::diff(f, x, detail::diff_tag<P>());
}

template<diff_pattern P, class F, class Real>
inline Real diff(const F f, Real x, Real y)
{
    return detail::diff(f, x, y, detail::diff_tag<P>());
}

template<diff_pattern P, class F, class Real>
inline Real diff(const F f, Real x, Real y, Real z)
{
    return detail::diff(f, x, y, z, detail::diff_tag<P>());
}

} // end of namecalc

#endif
