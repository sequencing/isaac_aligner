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
 ** \file FragmentIndex.hh
 **
 ** \brief Defines io structures for pre-bam bin fragment indexes.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_IO_FRAGMENT_INDEX_HH
#define iSAAC_IO_FRAGMENT_INDEX_HH

#include "alignment/FragmentMetadata.hh"
#include "alignment/BamTemplate.hh"
#include "oligo/Nucleotides.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace io
{

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


/**
 * \brief In terms of duplicate detection, anchor is same for duplicate candidates
 */
union FragmentIndexAnchor
{
    FragmentIndexAnchor() : value_(0UL){}
    explicit FragmentIndexAnchor(unsigned long value) : value_(value){}
    FragmentIndexAnchor(const alignment::FragmentMetadata & mate)
    {
        if(mate.isAligned()){
            pos_ = mate.getStrandReferencePosition().getValue();
        }else{
            shadowBases_ = oligo::pack32BclBases(mate.getBclData());
        }
    }
    // Aligned reads are anchored to their lowest cycles. This means that f-stranded reads are
    // anchored to their f-strand position and r-stranded to their r-strand alignment position
    // store ReferencePosition::value as ReferencePosition itself cannot be stored in a union
    unsigned long pos_;

    // use bases for duplicate detection of shadow (unaligned) reads.
    // the CASAVA uses 8 bases of the shadow.
    // TODO: decide whether dupe-detection will benefit from using more than
    // 32 bases of the shadow for now
    unsigned long shadowBases_;

    unsigned long value_;
};

inline std::ostream &operator <<(std::ostream& os, const FragmentIndexAnchor& anchor)
{
    return os << "FragmentIndexAnchor(" << anchor.value_ << ")";
}

struct FragmentIndexMate
{
    FragmentIndexMate(){}
    /**
     * \brief custom construction of FragmentIndexMate. Should be used for unit tests only
     */
    FragmentIndexMate(const bool shadow,
                      const bool reverse,
                      const unsigned storageBin,
                      const FragmentIndexAnchor anchor) :
        info_(shadow, reverse, storageBin),
        anchor_(anchor){}

    FragmentIndexMate(const alignment::FragmentMetadata &mate, const unsigned mateStorageBin) :
        //TODO: calculate mate storage bin!
        info_(!mate.isAligned(), mate.isReverse(), mateStorageBin),
        anchor_(mate){}

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

    FragmentIndexAnchor anchor_;
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
/**
 * \brief returns a value which can be used when ranking duplicates to chose the best one
 */
inline unsigned long getTemplateDuplicateRank(const alignment::BamTemplate &templ)
{
    return static_cast<unsigned long>(templ.getQuality()) << 32 |
        (templ.getTotalReadLength() - templ.getEditDistance()) << 16 |
        templ.getAlignmentScore();
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
    
    FStrandFragmentIndex(reference::ReferencePosition fStrandPos,
                         const FragmentIndexMate &mate,
                         const alignment::BamTemplate &bamTemplate):
                             PairEndIndex(fStrandPos, mate, getTemplateDuplicateRank(bamTemplate))
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
    FragmentIndexAnchor anchor_;

    RStrandOrShadowFragmentIndex(){}
    RStrandOrShadowFragmentIndex(reference::ReferencePosition fStrandPos,
                         FragmentIndexAnchor anchor,
                         const FragmentIndexMate &mate,
                         const unsigned long duplicateClusterRank):
                             PairEndIndex(fStrandPos, mate, duplicateClusterRank),
                             anchor_(anchor)
    {}
    
    RStrandOrShadowFragmentIndex(reference::ReferencePosition fStrandPos,
                         FragmentIndexAnchor anchor,
                         const FragmentIndexMate &mate,
                         const alignment::BamTemplate &bamTemplate):
                             PairEndIndex(fStrandPos, mate, getTemplateDuplicateRank(bamTemplate)),
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


} // namespace io
} // namespace isaac

#endif // #ifndef iSAAC_IO_FRAGMENT_INDEX_HH
