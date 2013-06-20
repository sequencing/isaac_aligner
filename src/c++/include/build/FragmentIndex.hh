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
 ** \file FragmentIndex.hh
 **
 ** \brief Defines io structures for pre-bam bin fragment indexes.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_FRAGMENT_INDEX_HH
#define iSAAC_BUILD_FRAGMENT_INDEX_HH

#include "io/Fragment.hh"
#include "oligo/Nucleotides.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace build
{

// Use this as a way to convert from tile local cluster id to a cluster id that is supposedly unique within flowcell.
static const unsigned long INSANELY_HIGH_NUMBER_OF_CLUSTERS_PER_TILE = 1000000000;


struct FragmentIndex
{
    FragmentIndex() : dataOffset_(0), mateDataOffset_(0){}
/*
    FragmentIndex(const unsigned long dataOffset, const unsigned long mateDataOffset):
        dataOffset_(dataOffset), mateDataOffset_(mateDataOffset){}
*/
    FragmentIndex(reference::ReferencePosition fStrandPos):
        fStrandPos_(fStrandPos), dataOffset_(0), mateDataOffset_(0)
    {}

    reference::ReferencePosition fStrandPos_;

    unsigned long dataOffset_;
    // if matches to dataOffset_, means that mate information is not accessible. Either single-ended or mate is in a different bin.
    unsigned long mateDataOffset_;
};

inline std::ostream &operator <<(std::ostream& os, const FragmentIndex& idx)
{
    return os << "FragmentIndex(" <<
        idx.fStrandPos_ << ", " <<
        idx.dataOffset_ << "do, " <<
        idx.mateDataOffset_ << "mdo " <<
        ")" << &idx;
}

//binary layout for non-paired fragment indexes
struct SeFragmentIndex : public FragmentIndex
{
    SeFragmentIndex(){}
    SeFragmentIndex(reference::ReferencePosition fStrandPos):
        FragmentIndex(fStrandPos)
    {}
};
BOOST_STATIC_ASSERT(24 == sizeof(SeFragmentIndex));

//binary layout for unaligned (not mapped) fragment indexes. Note that shadows of the pair are stored with \see {RStrandOrShadowFragmentIndex}
struct NmFragmentIndex : public FragmentIndex
{
    NmFragmentIndex():
        FragmentIndex(reference::ReferencePosition(reference::ReferencePosition::NoMatch))
    {}
};
BOOST_STATIC_ASSERT(24 == sizeof(NmFragmentIndex));


struct FragmentIndexMate
{
    FragmentIndexMate(){}

    FragmentIndexMate(const bool shadow,
                      const bool reverse,
                      const unsigned storageBin,
                      const io::FragmentIndexAnchor anchor) :
        info_(shadow, reverse, storageBin),
        anchor_(anchor){}

    union Info
    {
        Info(): value_(0) {}
        Info(const bool shadow, const bool reverse, const unsigned storageBin) {
            fields_.shadow_ = shadow;
            fields_.reverse_ = reverse;
            fields_.storageBin_ = storageBin;
        }
        unsigned value_;
        struct Fields
        {
            // set to 1 if mate is not aligned
            unsigned shadow_:1;

            // set to 1 if mate is r-strand aligned
            unsigned reverse_:1;

            // bin where the mate is stored. If mate is r-stranded, dupe-detect only across fragments that
            // have same mate storageBin_. This will ensure choice consistency for both fragments in template
            unsigned storageBin_:30;
        } fields_;
    } info_;

    io::FragmentIndexAnchor anchor_;
};

inline std::ostream &operator <<(std::ostream& os, const FragmentIndexMate::Info::Fields& mateInfoFields)
{
    return os << "FragmentIndexMate::Info::Fields(" <<
        mateInfoFields.shadow_ << "|" <<
        mateInfoFields.reverse_ << "|" <<
        mateInfoFields.storageBin_ << //"|" <<
        ")";
}

inline std::ostream &operator <<(std::ostream& os, const FragmentIndexMate& mate)
{
    return os << "FragmentIndexMate(" <<
        mate.anchor_ << ", " <<
        mate.info_.fields_ <<
        ")";
}

// base binary layout for an end of a pair
struct PairEndIndex : public FragmentIndex
{
    PairEndIndex() {}
    PairEndIndex(reference::ReferencePosition fStrandPos,
                 const FragmentIndexMate &mate,
                 const unsigned long duplicateClusterRank):
                     FragmentIndex(fStrandPos), mate_(mate), duplicateClusterRank_(duplicateClusterRank)
    {}

    FragmentIndexMate mate_;
    unsigned long duplicateClusterRank_;
};

inline std::ostream &operator <<(std::ostream& os, const PairEndIndex& idx)
{
    return os << "PairEndIndex(" <<
        idx.fStrandPos_ << ", " <<
        idx.mate_ << ", " <<
        idx.duplicateClusterRank_ << "dcr, " <<
        idx.dataOffset_ << "do, " <<
        idx.mateDataOffset_ << "mdo " <<
        ")" << &idx;
}

//binary layout for forward-strand fragment indexes
struct FStrandFragmentIndex : public PairEndIndex
{
    FStrandFragmentIndex(){}
    FStrandFragmentIndex(reference::ReferencePosition fStrandPos,
                         const FragmentIndexMate &mate,
                         const unsigned long duplicateClusterRank):
                             PairEndIndex(fStrandPos, mate, duplicateClusterRank)
    {}
};
BOOST_STATIC_ASSERT(48 == sizeof(FStrandFragmentIndex));

inline std::ostream &operator <<(std::ostream& os, const FStrandFragmentIndex& idx)
{
    return os << "FStrandFragmentIndex(" <<
        idx.fStrandPos_ << ", " <<
        idx.mate_ << ", " <<
        idx.duplicateClusterRank_ << "dcr, " <<
        idx.dataOffset_ << "do, " <<
        idx.mateDataOffset_ << "mdo " <<
        ")" << &idx;
}

// binary layout for reverse-strand and shadow fragment indexes
struct RStrandOrShadowFragmentIndex : public PairEndIndex
{
    io::FragmentIndexAnchor anchor_;

    RStrandOrShadowFragmentIndex(){}
    RStrandOrShadowFragmentIndex(
        reference::ReferencePosition fStrandPos,
        io::FragmentIndexAnchor anchor,
        const FragmentIndexMate &mate,
        const unsigned long duplicateClusterRank):
            PairEndIndex(fStrandPos, mate, duplicateClusterRank),
            anchor_(anchor)
    {}
};

BOOST_STATIC_ASSERT(56 == sizeof(RStrandOrShadowFragmentIndex));

inline std::ostream &operator <<(std::ostream& os, const RStrandOrShadowFragmentIndex& idx)
{
    return os << "RStrandOrShadowFragmentIndex(" <<
        idx.fStrandPos_ << ", " <<
        idx.anchor_ << ", " <<
        idx.mate_ << ", " <<
        idx.duplicateClusterRank_ << "dcr, " <<
        idx.dataOffset_ << "do, " <<
        idx.mateDataOffset_ << "mdo " <<
        ")" << &idx;
}


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_FRAGMENT_INDEX_HH
