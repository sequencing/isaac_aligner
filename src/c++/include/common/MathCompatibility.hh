/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file MathCompatibility.hh
 **
 ** \brief  Compatibility layer for math-related constructs.
 **
 ** \author Come Raczy
 **/

#ifndef iSAACCOMMON_MATH_COMPATIBILITY_HH
#define iSAACCOMMON_MATH_COMPATIBILITY_HH

#include <cmath>

#include "config.h"

#ifndef HAVE_FLOORF
inline float floorf(float x)
{
    return static_cast<float> (floor(x));
}
#define HAVE_FLOORF
#endif

#ifndef HAVE_ROUND
inline double round(double x)
{
    return (x - floor(x) < 0.5) ? floor(x) : ceil(x);
}
#define HAVE_ROUND
#endif

#ifndef HAVE_ROUNDF
inline float roundf(float x)
{
    return (x - floorf(x) < 0.5f ? floorf(x) : ceil(x));
}
#define HAVE_ROUNDF
#endif

#ifndef HAVE_POWF
inline float powf(float x, float y)
{
    return static_cast<float> (pow(x, y));
}
#define HAVE_POWF
#endif

#ifndef HAVE_ERF
#include <boost/math/special_functions/erf.hpp>
inline double erf(double x)
{
    return boost::math::erf(x);
}
#define HAVE_ERF
#endif

#ifndef HAVE_ERFF
#include <boost/math/special_functions/erf.hpp>
inline float erff(float x)
{
    return boost::math::erf(x);
}
#define HAVE_ERFF
#endif

#ifndef HAVE_ERFC
#include <boost/math/special_functions/erf.hpp>
inline double erfc(double x)
{
    return boost::math::erfc(x);
}
#define HAVE_ERFC
#endif

#ifndef HAVE_ERFCF
#include <boost/math/special_functions/erf.hpp>
inline float erfc(float x)
{
    return boost::math::erfc(x);
}
#define HAVE_ERFCF
#endif

#endif // #ifndef iSAACCOMMON_MATH_COMPATIBILITY_HH
