/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/downloads/sequencing/licenses/>.
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

#include "oligo/Kmer.hh"

namespace isaac
{
namespace oligo
{

inline unsigned int getMaskCount(const unsigned int maskWidth) {return (1 << maskWidth);}

const unsigned int defaultMaskWidth = 6;

inline std::string getSortedReferenceFileName(const boost::filesystem::path &outputPrefix,
                                              const std::string &permutationName, const Kmer mask,
                                              const std::string &outputSuffix)
{
    return (boost::format("%s%s-%02d%s") %
            outputPrefix.string() % permutationName % mask % outputSuffix).str();
}

inline boost::filesystem::path getMaskPath(const unsigned int mask, const boost::format &maskPathFormat, unsigned int pathCount = 4)
{
    boost::format format = maskPathFormat;
    return (1 == format.expected_args())
        ? (format % (1 + (mask % pathCount))).str()
        : ".";
}

} // namespace oligo
} // namespace isaac

#endif // #ifndef iSAAC_OLIGO_MASK_HH
