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
 ** \file PackedFragmentBuffer.hh
 **
 ** Helper.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_PACKED_FRAGMENT_BUFFER_HH
#define iSAAC_BUILD_PACKED_FRAGMENT_BUFFER_HH

#include "alignment/BinMetadata.hh"
#include "build/FragmentIndex.hh"

namespace isaac
{
namespace build
{

/**
 * \brief Helper to accees fragments stored in a contiguous byte vector
 */
class PackedFragmentBuffer : std::vector<char>
{
public:
    struct Index
    {
        Index(
            const isaac::reference::ReferencePosition pos,
            unsigned long dataOffset,
            unsigned long mateDataOffset_,
            const unsigned *cigarBegin,
            const unsigned *cigarEnd):
                pos_(pos), dataOffset_(dataOffset), mateDataOffset_(dataOffset),
                cigarBegin_(cigarBegin), cigarEnd_(cigarEnd)
        {}

        Index(const FStrandFragmentIndex &idx, const io::FragmentAccessor &fragment) :
            pos_(idx.fStrandPos_), dataOffset_(idx.dataOffset_), mateDataOffset_(idx.mateDataOffset_),
            cigarBegin_(fragment.cigarBegin()), cigarEnd_(fragment.cigarEnd()){}
        Index(const RStrandOrShadowFragmentIndex &idx, const io::FragmentAccessor &fragment) :
            pos_(idx.fStrandPos_), dataOffset_(idx.dataOffset_), mateDataOffset_(idx.mateDataOffset_),
            cigarBegin_(fragment.cigarBegin()), cigarEnd_(fragment.cigarEnd()){}
        Index(const SeFragmentIndex &idx, const io::FragmentAccessor &fragment) :
            pos_(idx.fStrandPos_), dataOffset_(idx.dataOffset_), mateDataOffset_(idx.dataOffset_),
            cigarBegin_(fragment.cigarBegin()), cigarEnd_(fragment.cigarEnd()){}

        bool hasMate() const
        {
            return mateDataOffset_ != dataOffset_;
        }

        unsigned getBeginClippedLength() const
        {
            ISAAC_ASSERT_MSG(cigarBegin_ != cigarEnd_, "Unexpected empty CIGAR");
            std::pair<unsigned, alignment::Cigar::OpCode> operation = alignment::Cigar::decode(*cigarBegin_);
            if (alignment::Cigar::SOFT_CLIP == operation.second)
            {
                return operation.first;
            }
            return 0;
        }

        /**
         * \brief Returns unadjusted position if it is adjusted due to a soft clipping
         */
        isaac::reference::ReferencePosition getUnclippedPosition() const
        {
            return pos_ - getBeginClippedLength();
        }

        // temporary storage for fragment.fStrandPosition_. Not guaranteed to be up to date.
        // ensure it is synchronized with fragment.fStrandPosition_ before using it
        isaac::reference::ReferencePosition pos_;
        unsigned long dataOffset_;
        // same as adataOffset_ for single-ended
        unsigned long mateDataOffset_;

        const uint32_t *cigarBegin_;
        const uint32_t *cigarEnd_;
    };


    using std::vector<char>::front;
    using std::vector<char>::size;
    using std::vector<char>::begin;
    using std::vector<char>::end;

    void resize(const alignment::BinMetadata& bin)
    {
        std::vector<char>::resize(bin.getDataSize());
    }

    void unreserve()
    {
        std::vector<char>().swap(*this);
    }

    io::FragmentAccessor &getFragment(const Index& fragmentIndex)
        {return getFragment(fragmentIndex.dataOffset_);}

    io::FragmentAccessor &getMate(const Index& fragmentIndex)
        {return getFragment(fragmentIndex.mateDataOffset_);}

    const io::FragmentAccessor &getFragment(const Index& fragmentIndex) const
        {return getFragment(fragmentIndex.dataOffset_);}

    io::FragmentAccessor &getFragment(const FragmentIndex &fragmentIndex)
        {return getFragment(fragmentIndex.dataOffset_);}

    const io::FragmentAccessor &getFragment(const FragmentIndex &fragmentIndex) const
        {return getFragment(fragmentIndex.dataOffset_);}

    const io::FragmentAccessor &getMate(const FragmentIndex& fragmentIndex) const
        {return getFragment(fragmentIndex.mateDataOffset_);}

    io::FragmentAccessor &getFragment(
        const unsigned long offset)
    {
        io::FragmentAccessor &fragment =
            *reinterpret_cast<io::FragmentAccessor*>(&at(offset));
        return fragment;
    }

    const io::FragmentAccessor &getFragment(
        const unsigned long offset) const
    {
        const io::FragmentAccessor &fragment =
            *reinterpret_cast<const io::FragmentAccessor*>(&at(offset));
        return fragment;
    }

    static unsigned long getMemoryRequirements(const alignment::BinMetadata& bin)
    {
        return bin.getDataSize();
    }

    bool orderForBam(const Index &left, const Index &right) const
    {
        if (left.pos_ < right.pos_)
        {
            return true;
        }
        else if (left.pos_ == right.pos_)
        {
            const io::FragmentAccessor &leftFragment = getFragment(left);
            const io::FragmentAccessor &rightFragment = getFragment(right);
            const unsigned long leftGlobalClusterId = leftFragment.tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + leftFragment.clusterId_;
            const unsigned long rightGlobalClusterId = rightFragment.tile_ * INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE + rightFragment.clusterId_;

            return
                leftGlobalClusterId < rightGlobalClusterId ||
                (leftGlobalClusterId == rightGlobalClusterId &&
                    // ensure singleton and shadow stay together and shadow follows the singleton
                    (leftFragment.flags_.unmapped_ < rightFragment.flags_.unmapped_ ||
                        (leftFragment.flags_.unmapped_ == rightFragment.flags_.unmapped_ &&
                            // SAAC-378 ensure second read comes after first to generate consistent bam files between different runs
                            leftFragment.flags_.secondRead_ < rightFragment.flags_.secondRead_
                        )
                    )
                );
        }

        return false;
    }
};

inline std::ostream & operator <<(std::ostream &os, const PackedFragmentBuffer::Index &index)
{
    return os << "PackedFragmentBuffer::Index(" <<
        index.pos_ << "," << index.dataOffset_ << "do " << index.mateDataOffset_ << "mdo, " <<
        alignment::Cigar::toString(index.cigarBegin_, index.cigarEnd_) << ")";
}

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_PACKED_FRAGMENT_BUFFER_HH
