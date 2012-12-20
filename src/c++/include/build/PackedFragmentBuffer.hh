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
 ** \file PackedFragmentBuffer.hh
 **
 ** Helper.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_PACKED_FRAGMENT_BUFFER_HH
#define iSAAC_BUILD_PACKED_FRAGMENT_BUFFER_HH

#include "alignment/BinMetadata.hh"
#include "io/Fragment.hh"
#include "io/FragmentIndex.hh"

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
            const reference::ReferencePosition pos,
            unsigned long dataOffset,
            unsigned long mateDataOffset_,
            const unsigned *cigarBegin,
            const unsigned *cigarEnd):
                pos_(pos), dataOffset_(dataOffset), mateDataOffset_(dataOffset),
                cigarBegin_(cigarBegin), cigarEnd_(cigarEnd)
        {}

        Index(const io::FStrandFragmentIndex &idx, const io::FragmentAccessor &fragment) :
            pos_(idx.fStrandPos_), dataOffset_(idx.dataOffset_), mateDataOffset_(idx.mateDataOffset_),
            cigarBegin_(fragment.cigarBegin()), cigarEnd_(fragment.cigarEnd()){}
        Index(const io::RStrandOrShadowFragmentIndex &idx, const io::FragmentAccessor &fragment) :
            pos_(idx.fStrandPos_), dataOffset_(idx.dataOffset_), mateDataOffset_(idx.mateDataOffset_),
            cigarBegin_(fragment.cigarBegin()), cigarEnd_(fragment.cigarEnd()){}
        Index(const io::SeFragmentIndex &idx, const io::FragmentAccessor &fragment) :
            pos_(idx.fStrandPos_), dataOffset_(idx.dataOffset_), mateDataOffset_(idx.dataOffset_),
            cigarBegin_(fragment.cigarBegin()), cigarEnd_(fragment.cigarEnd()){}

        bool hasMate() const
        {
            return mateDataOffset_ != dataOffset_;
        }

        // temporary storage for fragment.fStrandPosition_. Not guaranteed to be up to date.
        // ensure it is synchronized with fragment.fStrandPosition_ before using it
        reference::ReferencePosition pos_;
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

    io::FragmentAccessor &getFragment(const io::FragmentIndex &fragmentIndex)
        {return getFragment(fragmentIndex.dataOffset_);}

    const io::FragmentAccessor &getFragment(const io::FragmentIndex &fragmentIndex) const
        {return getFragment(fragmentIndex.dataOffset_);}

    const io::FragmentAccessor &getMate(const io::FragmentIndex& fragmentIndex) const
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
            // ensure singleton and shadow stay togeter and shadow follows the singleton
            if (leftFragment.clusterId_ < rightFragment.clusterId_)
            {
                return true;
            }
            else if (leftFragment.clusterId_ == rightFragment.clusterId_)
            {
                return leftFragment.flags_.unmapped_ < rightFragment.flags_.unmapped_;
            }
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
