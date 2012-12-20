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
 ** \file UnsortedAlignmentBamAdapter.hh
 **
 ** Implements a translation interface required for serializing UnsortedAlignment into bam
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_UNSORTED_ALIGNMENT_BAM_ADAPTER_HH
#define iSAAC_BUILD_UNSORTED_ALIGNMENT_BAM_ADAPTER_HH

#include "UnsortedAlignment.hh"

namespace isaac
{
namespace build
{

class UnsortedAlignmentBamAdapter
{
    const UnsortedAlignment& alignment_;
public:
    UnsortedAlignmentBamAdapter(const UnsortedAlignment& alignment) : alignment_(alignment){}

    std::string readName() const { return alignment_.read_name; }

    std::vector<unsigned> cigar() const {
        return std::vector<unsigned>(alignment_.cigar,
                                     alignment_.cigar + sizeof(alignment_.cigar) / sizeof(alignment_.cigar[0]));
    }

    int seqLen() const { return sizeof(alignment_.seq) * 2; }

    std::vector<unsigned char> seq() const {
        return std::vector<unsigned char>(alignment_.seq,
                                          alignment_.seq + sizeof(alignment_.seq) / sizeof(alignment_.seq[0]));
    }

    std::vector<unsigned char> qual() const {
        return std::vector<unsigned char>(alignment_.qual,
                                          alignment_.qual + sizeof(alignment_.qual) / sizeof(alignment_.qual[0]));
    }

    int         refId()     const { return 0; }
    int         pos()       const { return alignment_.alignmentPos_; }
    unsigned    mapq()      const { return alignment_.sras_;}
    unsigned    flag()      const { return 0;}
    int         nextRefId() const { return 1; }
    int         nextPos()   const { return alignment_.alignmentPos_;}
    int         tlen()      const { return alignment_.templateLength_; }

};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_UNSORTED_ALIGNMENT_BAM_ADAPTER_HH
