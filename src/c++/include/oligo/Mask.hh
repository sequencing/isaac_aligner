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
 ** \file Mask.hh
 **
 ** Masking operations on oligomers.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_OLIGO_MASK_HH
#define iSAAC_OLIGO_MASK_HH

#include <string>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "oligo/Kmer.hh"

namespace isaac
{
namespace oligo
{

inline unsigned int getMaskCount(const unsigned int maskWidth) {return (1 << maskWidth);}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_MASK_HH
