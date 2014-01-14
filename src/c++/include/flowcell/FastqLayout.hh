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
 ** \file BamLayout.hh
 **
 ** Specialization of Layout for bam flowcell.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_FASTQ_LAYOUT_HH
#define iSAAC_FLOWCELL_FASTQ_LAYOUT_HH

#include "flowcell/Layout.hh"

namespace isaac
{
namespace flowcell
{

namespace fastq
{
    static const unsigned READ_NUMBER_MAX = 2;


    void getFastqFilePath(
        const boost::filesystem::path &baseCallsPath,
        const unsigned lane,
        const unsigned read,
        const bool compressed,
        boost::filesystem::path &result);
};

struct FastqFilePathAttributeTag
{
    typedef boost::filesystem::path value_type;
    friend std::ostream &operator << (std::ostream &os, const FastqFilePathAttributeTag &tag){return os << "FastqFilePathAttributeTag";}
};

template<>
const boost::filesystem::path & Layout::getLaneReadAttribute<Layout::Fastq, FastqFilePathAttributeTag>(
    const unsigned lane, const unsigned read, boost::filesystem::path &result) const;

template<>
inline boost::filesystem::path Layout::getLongestAttribute<Layout::Fastq, FastqFilePathAttributeTag>() const
{
    boost::filesystem::path filtersFilePath;
    getLaneReadAttribute<Layout::Fastq, FastqFilePathAttributeTag>(laneNumberMax_, fastq::READ_NUMBER_MAX, filtersFilePath);
    return filtersFilePath;
}



} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_FASTQ_LAYOUT_HH
