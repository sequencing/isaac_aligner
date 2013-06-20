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

#ifndef iSAAC_COMMON_ENDIANNESS_HH
#define iSAAC_COMMON_ENDIANNESS_HH

namespace isaac
{
namespace common
{

template <typename T, typename IterT> IterT extractLittleEndian(const IterT p, T &result)
{
#ifdef ISAAC_IS_BIG_ENDIAN
    // TODO: untested
    const unsigned char *cp = reinterpret_cast<const unsigned char*>(p);
    result = T(*cp) + *(cp + 1) * 256U + *(cp + 2) * 256U * 256U + *(cp + 3) * 256U * 256U * 256U;
#endif
    result = *reinterpret_cast<const T*>(&*p);
    return p + sizeof(T);
}

template <typename T, typename IterT> T extractLittleEndian(const IterT p)
{
#ifdef ISAAC_IS_BIG_ENDIAN
    // TODO: untested
    const unsigned char *cp = reinterpret_cast<const unsigned char*>(p);
    return T(*cp) + *(cp + 1) * 256U + *(cp + 2) * 256U * 256U + *(cp + 3) * 256U * 256U * 256U;
#endif
    return *reinterpret_cast<const T*>(&*p);
}

} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_ENDIANNESS_HH
