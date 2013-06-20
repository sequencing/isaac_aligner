/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file FastIo.hh
 **
 ** \brief Fast IO routines for integers and fixed width floating points.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_FAST_IO_HH
#define iSAAC_COMMON_FAST_IO_HH

#include <limits>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/bool.hpp>

#include "common/MathCompatibility.hh"

namespace isaac
{
namespace common
{

/**
 ** @brief Fast output of a float value with a fixed number of decimal positions.
 **
 ** The number is written as an optional '-' followed by the decimal value. If
 ** the template parameter 'decimals == 0', then the output is an integer (no
 ** decimal point).
 **
 ** When the number does not fit within the given 'maxWidth', the characters are
 ** truncated on the left, including the '-' if any. If 'decimals+2 > maxWidth'
 ** the behavior is undefined.
 **
 ** When the width of the output value is less than 'minWidth', the output is
 ** padded on the left with ' '.
 **
 ** @param buffer: the output buffer. Must have enough space available to store
 ** max(decimals+2, maxWidth) characters.
 **
 ** @param floatValue: the value to output into the buffer.
 **/
template<int decimals, int minWidth, int maxWidth>
int sprintFloat(char *buffer, const float originalValue)
{
    // Compile-time sanity checks over the template parameters
    BOOST_MPL_ASSERT_MSG( (0 <= decimals), THE_NUMBER_OF_DECIMALS_MUST_NOT_BE_NEGATIVE, () );
    BOOST_MPL_ASSERT_MSG( (minWidth <= maxWidth), THE_MAX_WIDTH_MUST_NOT_BE_LESSER_THAN_MIN_WIDTH, () );
    BOOST_MPL_ASSERT_MSG( (decimals+2 <= maxWidth), THE_MAX_WIDTH_MUST_BE_GREATER_THAN_DECIMALS_PLUS_ONE, () );
    // Synopsis: the digits are first written in the reverse order,
    // then the padding is added, and finally the string is reversed.
    //
    // multiply to expose the necessary decimal positions after rounding
    char * const begin = buffer;
    float floatValue = originalValue;
    for (unsigned int i = 0; decimals != i; ++i)
    {
        floatValue *= 10.0;
    }
    assert(std::numeric_limits<long>::min() < floatValue
            && std::numeric_limits<long>::max() > floatValue);
    const float rounding = (originalValue < -0.0 ? -0.5f : 0.5f);
    long value = static_cast<long> (floatValue + rounding);
    const bool negative = (0 > value);
    // work with positive values to avoid the modulo problems
    if (negative)
        value = -value;
    // write all the digits after the decimal point
    for (unsigned int i = 0; decimals != i; ++i)
    {
        *buffer++ = '0' + (value % 10);
        value /= 10;
    }
    // write the decimal point iif there are decimals
    if (0 != decimals)
    {
        *buffer++ = '.';
    }
    // write the digit for the units
    *buffer++ = '0' + (value % 10);
    value /= 10;
    // write the remaining digits, if there is enough space
    char * const end = begin + maxWidth;
    while (end > buffer && 0 != value)
    {
        *buffer++ = '0' + (value % 10);
        value /= 10;
    }
    // write the '-' sign, if any, and if there is enough space
    if (end > buffer && negative)
        *buffer++ = '-';
    // add the padding if necessary
    char * const minEnd = begin + minWidth;
    if (minEnd > buffer)
    {
        memset(buffer, ' ', minEnd - buffer);
        buffer = minEnd;
    }
    std::reverse(begin, buffer); // faster than hand-crafted reverse
    return buffer - begin;
}

/**
 ** \brief append an unsigned integer to a string.
 **
 ** Performance: in the order of 100 million digits / second if the string has
 ** enough pre-allocated capacity.
 **
 ** Note: recursive implementation is about 5% slower
 **/
template <typename ContainerT, typename NumberT>
inline  void appendUnsignedNumber(ContainerT &s, NumberT value)
{
    const size_t initialLength = s.size();
    while (10 <= value)
    {
        s.push_back('0' + (value % 10));
        value /= 10;
    }
    s.push_back('0' + (value % 10));
    std::reverse(s.begin() + initialLength, s.end());
}

/**
 ** \brief append an unsigned integer to a string.
 **
 ** Performance: in the order of 100 million digits / second if the string has
 ** enough pre-allocated capacity.
 **
 ** Note: recursive implementation is about 5% slower
 **/
template <typename ContainerT>
inline  void appendUnsignedInteger(ContainerT &s, unsigned value)
{
    appendUnsignedNumber(s, value);
}

/**
 ** @brief Fast and portable output of an unsigned integer into a stream.
 **
 ** Drops all the usual formatting options to the benefit of speed.
 **
 ** Can be instantiated on unsigned versions of char, short, int and
 ** long, and generally on any type supporting '<=', '%', '/=', '+'
 ** and defining 'digits10 in std::numeric_limits.
 **/
template<class T>
std::ostream &putUnsignedInteger(std::ostream &os, T value)
{
    // generate a compilation error when instantiated on signed types
    BOOST_MPL_ASSERT_MSG(boost::is_unsigned<T>::value,
            SIGNED_TYPES_ARE_NOT_ALLOWED_FOR_putUnsignedInteger, (T));

    char begin[1 + std::numeric_limits<T>::digits10];
    char *buffer = begin;
    while (10 <= value)
    {
        *buffer++ = '0' + (value % 10);
        value /= 10;
    }
    *buffer++ = '0' + value;
    std::reverse(begin, buffer);
    return os.write(begin, buffer - begin);
}

/**
 ** @brief Fast and portable output of an unsigned integer into a stream.
 **
 ** Can be instantiated on (signed versions of) char, short, int and long.
 **
 ** @see putUnsignedInteger
 **/
template<class T>
std::ostream &putInteger(std::ostream &os, T value)
{
typedef    typename boost::make_unsigned<T>::type Unsigned;
    if (0 > value)
    {
        os.put('-');
        value = -value;
    }
    return putUnsignedInteger<Unsigned>(os, static_cast<Unsigned>(value));
}

/**
 ** @brief Read an unsigned integer.
 **
 ** Assumes that the first available character is the first digit of the
 ** integer. Reads until a non digit character is found. That non digit
 ** character is discarded if the "discardSeparator" parameter is set.
 **/
template <class T>
std::istream &getUnsignedInteger(std::istream &is, T &value,
        const bool discardSeparator = false)
{
    // generate a compilation error when instantiated on signed types
    BOOST_MPL_ASSERT_MSG( boost::is_unsigned<T>::value, SIGNED_TYPES_ARE_NOT_ALLOWED_FOR_putUnsignedInteger , (T) );
    value = 0;
    char c;
    while(is.get(c) && '0' <= c && '9' >= c)
    {
        value *= 10;
        value += (c-'0');
    }
    if (is && ! discardSeparator) is.putback(c);
    return is;
}

/**
 ** Read an integer (possibly signed).
 **
 ** Similar to getUnsignedInteger, except that the first character is either
 ** the '-' sign or the first digit.
 **/
template <class T>
std::istream &getInteger(std::istream &is, T &value,
        const bool discardSeparator = false)
{
    value = 0;
    char c;
    const bool isNegative = ('-' == is.peek());
    if (isNegative)
    is.get(c);
    typedef typename boost::make_unsigned<T>::type Unsigned;
    Unsigned unsignedValue;
    getUnsignedInteger<Unsigned>(is, unsignedValue, discardSeparator);
    if (isNegative)
    {
        value = -unsignedValue;
    }
    else
    {
        value = unsignedValue;
    }
    return is;
}

/**
 ** @brief Write a boolean as a N (false) or Y (true)
 **
 ** example:
 ** \code     putBool<'1','0'>(os, true); \endcode
 **/
template<char Y, char N>
std::ostream &putBool(std::ostream &os, bool b)
{
    return os.put(b? Y : N);
}

/**
 ** @brief Read a 0/1 value and convert it into a boolean
 **
 ** If the stream fails to produce 1 character (eof, fail or bad), no
 ** further action is taken. Else, the b is set to true/false if the
 ** character read is Y/N. The stream is set to fail if the character
 ** read is neither Y nor N.
 **/
template<char Y, char N>
std::istream &getBool(std::istream &is, bool &b)
{
    char c;
    if (is.get(c))
    {
        if ((Y == c || N == c))
        {
            b = (Y == c);
        }
        else
        {
            is.clear(std::ios_base::failbit);
        }
    }
    return is;
}

/**
 ** \brief Support for big endian systems
 **/
template <int DataSize>
inline char *reorderBytes(char *input)
{
    BOOST_MPL_ASSERT_RELATION(DataSize, >, 0);
#ifdef ISAAC_IS_BIG_ENDIAN
    std::reverse(input, input + DataSize);
#endif
    return input;
}

/**
 ** \brief Identification of the type of an integer based on its size.
 **
 ** Sizes that are not registered will default to a void type and
 ** should result into a compilation error.
 **
 ** There will be a compilation error if any of these sizes are equal
 ** to each other.
 **/
template <int DataSize>
struct SignedIntegerTypeTraits
{
    typedef void type;
};

template<> struct SignedIntegerTypeTraits<sizeof(char)>
{
    typedef char type;
};
template<> struct SignedIntegerTypeTraits<sizeof(short)>
{
    typedef short type;
};
template<> struct SignedIntegerTypeTraits<sizeof(int)>
{
    typedef int type;
};

/**
 ** \brief Identification of the type of an unsigned integer based on its size.
 **/
template <int DataSize>
struct UnsignedIntegerTypeTraits
{
    typedef typename SignedIntegerTypeTraits<DataSize>::type SignedType;
    typedef typename boost::mpl::eval_if<boost::is_void<SignedType>,boost::mpl::identity<void>, typename boost::make_unsigned<SignedType> >::type type;
};

/**
 ** \brief Identification of the type of a decimal number based on its size.
 **
 ** Sizes that are not registered will default to a void type and
 ** should result into a compilation error.
 **
 ** There will be a compilation error if any of these sizes are equal
 ** to each other.
 **/
template <int DataSize>
struct DecimalNumberTypeTraits
{
    typedef void type;
};

template<> struct DecimalNumberTypeTraits<sizeof(float)>
{
    typedef float type;
};

template<> struct DecimalNumberTypeTraits<sizeof(double)>
{
    typedef double type;
};

/**
 ** \brief Read an integer of a given length in bytes.
 **
 ** Assumes that the stream provides the data little end first. The
 ** bytes are reordered if the system is big endian (as defined in
 ** config.h at configuration time).
 **
 ** \tparam the number of bytes used to store the data
 **
 ** \return the value read (undefined when the read fails)
 **/
template <int DataSize>
typename SignedIntegerTypeTraits<DataSize>::type readInteger(std::istream &is)
{
    typedef typename SignedIntegerTypeTraits<DataSize>::type DataType;
    BOOST_MPL_ASSERT((boost::is_integral<DataType>));
    BOOST_MPL_ASSERT_RELATION(sizeof(DataType), ==, DataSize);
    DataType result = 0;
    char * const buffer = reinterpret_cast<char *>(&result);
    is.read(buffer, DataSize);
    reorderBytes<sizeof(DataType)>(buffer);
    return result;
}

/**
 ** \brief Read a signed decimal number (floating point).
 **
 ** \return the value read (0. when the read fails)
 **/
template <int DataSize>
typename DecimalNumberTypeTraits<DataSize>::type readDecimalNumber(std::istream &is)
{
    typedef typename DecimalNumberTypeTraits<DataSize>::type DataType;
    BOOST_MPL_ASSERT((boost::is_floating_point<DataType>));
    BOOST_MPL_ASSERT_RELATION(sizeof(DataType), ==, DataSize);
    DataType result = 0;
    char * const buffer = reinterpret_cast<char *>(&result);
    is.read(buffer, DataSize);
    reorderBytes<sizeof(DataType)>(buffer);
    return result;
}

inline float readFloat(std::istream &is)
{
    return readDecimalNumber<4>(is);
}

inline double readDouble(std::istream &is)
{
    return readDecimalNumber<8>(is);
}

/**
 ** \brief Read an unsigned integer of a given length in bytes.
 **
 ** \tparam the number of bytes used to store the data
 **
 ** \return the value read (0 when the read fails)
 **/
template <int DataSize>
typename UnsignedIntegerTypeTraits<DataSize>::type readUnsignedInteger(std::istream &is)
{
    typedef typename UnsignedIntegerTypeTraits<DataSize>::type DataType;
    return static_cast<DataType>(readInteger<DataSize>(is));
}

template <unsigned int DataSize>
std::ostream &writeUnsignedInteger(std::ostream &os, const typename SignedIntegerTypeTraits<DataSize>::type &value)
{
    typedef typename UnsignedIntegerTypeTraits<DataSize>::type DataType;
    BOOST_MPL_ASSERT_RELATION(sizeof(DataType), ==, DataSize);
    typename UnsignedIntegerTypeTraits<DataSize>::type reoderedValue(value);
    char * const buffer = reinterpret_cast<char *>(&reoderedValue);
    reorderBytes<sizeof(DataType)>(buffer);
    os.write(buffer, DataSize);
    return os;
}

template <unsigned int DataSize>
std::ostream &writeInteger(std::ostream &os, const typename SignedIntegerTypeTraits<DataSize>::type &value)
{
    typedef typename SignedIntegerTypeTraits<DataSize>::type DataType;
    BOOST_MPL_ASSERT_RELATION(sizeof(DataType), ==, DataSize);
    typename SignedIntegerTypeTraits<DataSize>::type reoderedValue(value);
    char * const buffer = reinterpret_cast<char *>(&reoderedValue);
    reorderBytes<sizeof(DataType)>(buffer);
    os.write(buffer, DataSize);
    return os;
}

template <unsigned int DataSize>
std::ostream &writeDecimalNumber(std::ostream &os, const typename DecimalNumberTypeTraits<DataSize>::type &value)
{
    typedef typename DecimalNumberTypeTraits<DataSize>::type DataType;
    BOOST_MPL_ASSERT((boost::is_floating_point<DataType>));
    BOOST_MPL_ASSERT_RELATION(sizeof(DataType), ==, DataSize);
    typename DecimalNumberTypeTraits<DataSize>::type reoderedValue(value);
    char * const buffer = reinterpret_cast<char *>(&reoderedValue);
    reorderBytes<sizeof(DataType)>(buffer);
    os.write(buffer, DataSize);
    return os;
}

inline std::ostream &writeFloat(std::ostream &os, const float &value)
{
    return writeDecimalNumber<4>(os, value);
}

inline std::ostream &writeDouble(std::ostream &os, const double &value)
{
    return writeDecimalNumber<8>(os, value);
}


} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_FAST_IO_HH
