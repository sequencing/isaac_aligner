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
 ** \file Read.hh
 **
 ** \brief Data associated to a cluster: sequence and quality strings for all
 ** the reads in the cluster.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_READ_HH
#define iSAAC_ALIGNMENT_READ_HH

#include <string>
#include <iostream>

#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
/**
 ** \brief Encapsulates the sequence and quality string for a read.
 **/
class Read
{
public:
    /// Default constructor to enable use in containers
    //explicit Read(unsigned index = 0) : index_(index) {}
    Read(const unsigned maxReadLength, const unsigned index) : index_(index), /*beginCyclesMasked_(0), */endCyclesMasked_(0)
    {
        forwardSequence_.reserve(maxReadLength);
        reverseSequence_.reserve(maxReadLength);
        forwardQuality_.reserve(maxReadLength);
        reverseQuality_.reserve(maxReadLength);
    }

    /// Copy constructor to preserve the reserverd capacity during container push_back
    Read(const Read &read)
        : index_(read.index_)
        , forwardSequence_(read.forwardSequence_)
        , reverseSequence_(read.reverseSequence_)
        , forwardQuality_(read.forwardQuality_)
        , reverseQuality_(read.reverseQuality_)
        /*, beginCyclesMasked_(read.beginCyclesMasked_)*/
        , endCyclesMasked_(read.endCyclesMasked_)
    {
        // keep the size of pre-allocated buffers
        forwardSequence_.reserve(read.forwardSequence_.capacity());
        reverseSequence_.reserve(read.reverseSequence_.capacity());
        forwardQuality_.reserve(read.forwardQuality_.capacity());
        reverseQuality_.reserve(read.reverseQuality_.capacity());
    }
    
    const std::vector<char> &getStrandSequence(bool reverse) const {return reverse ? reverseSequence_ : forwardSequence_;}
    const std::vector<char> &getStrandQuality(bool reverse) const {return reverse ? reverseQuality_ : forwardQuality_;}

    const std::vector<char> &getForwardSequence() const {return forwardSequence_;}
    const std::vector<char> &getReverseSequence() const {return reverseSequence_;}
    const std::vector<char> &getForwardQuality() const {return forwardQuality_;}
    const std::vector<char> &getReverseQuality() const {return reverseQuality_;}
    unsigned getLength() const {return forwardSequence_.size();}
    unsigned getBeginCyclesMasked() const { return 0; /*return beginCyclesMasked_;*/}
    unsigned getEndCyclesMasked() const {return endCyclesMasked_;}
    unsigned getIndex() const {return index_;}
    void maskCyclesFromEnd(unsigned cycles) {endCyclesMasked_ = cycles;}

    /// storing Read objects in the vector requires this operator although it is not expected to be executed at runtime
    Read &operator=(const Read &read)
    {
        ISAAC_ASSERT_MSG(false, "Read objects are not supposed to be reassigned");
        return *this;
    }

    /// Decode BCL data and initialize the sequencd and quality string appropriately
    void decodeBcl(std::vector<char>::const_iterator bclBegin, std::vector<char>::const_iterator bclEnd, unsigned index);
    template<class InpuT> friend InpuT& operator >>(InpuT &input, Read &read);
private:
    const unsigned index_;

    std::vector<char> forwardSequence_;
    std::vector<char> reverseSequence_;
    std::vector<char> forwardQuality_;
    std::vector<char> reverseQuality_;
    /// number of cycles masked at the start of the read.
    //unsigned beginCyclesMasked_;
    /// number of cycles masked at the end of the read.
    unsigned endCyclesMasked_;
};

std::ostream &operator<<(std::ostream &os, const Read &read);

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_READ_HH
