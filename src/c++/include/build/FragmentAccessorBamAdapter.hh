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
 ** \file FragmentAccessorBamAdapter.hh
 **
 ** Implements a translation interface required for serializing FragmentAccessor into bam
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_FRAGMENT_ACCESSOR_BAM_ADAPTER_HH
#define iSAAC_BUILD_FRAGMENT_ACCESSOR_BAM_ADAPTER_HH

#include <bitset>
#include <iterator>

#include "bam/Bam.hh"
#include "build/BuildContigMap.hh"
#include "build/PackedFragmentBuffer.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace build
{
struct IncludeTags
{
    IncludeTags(
        bool includeAS,
        bool includeBC,
        bool includeNM,
        bool includeOC,
        bool includeRG,
        bool includeSM,
        bool includeZX,
        bool includeZY) :
            includeAS_(includeAS),
            includeBC_(includeBC),
            includeNM_(includeNM),
            includeOC_(includeOC),
            includeRG_(includeRG),
            includeSM_(includeSM),
            includeZX_(includeZX),
            includeZY_(includeZY){}
    bool includeAS_;
    bool includeBC_;
    bool includeNM_;
    bool includeOC_;
    bool includeRG_;
    bool includeSM_;
    bool includeZX_;
    bool includeZY_;
} ;


class FragmentAccessorBamAdapter
{
    const unsigned maxReadLength_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const BuildContigMap &contigMap_;
    reference::ReferencePosition pos_;
    const io::FragmentAccessor* pFragment_;
    std::vector<char> readGroupNameBuffer_;
    std::vector<char> barcodeNameBuffer_;
    std::vector<char> readNameBuffer_;
    std::vector<char> originalCigarBuffer_;
    const unsigned *cigarBegin_;
    const unsigned *cigarEnd_;
    // if read got realigned these point at the original cigar
    const unsigned *originalCigarBegin_;
    const unsigned *originalCigarEnd_;
    std::vector<unsigned char> seqBuffer_;
    std::vector<unsigned char> qualBuffer_;
    unsigned char forcedDodgyAlignmentScore_;
    flowcell::FlowcellLayoutList const &flowCellLayoutList_;
    const IncludeTags includeTags_;
    const bool pessimisticMapQ_;

public:
    FragmentAccessorBamAdapter(
        unsigned maxReadLength,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const BuildContigMap &contigMap,
        const unsigned char forcedDodgyAlignmentScore,
        const flowcell::FlowcellLayoutList &flowCellLayoutList,
        const IncludeTags includeTags,
        const bool pessimisticMapQ
        ) :
        maxReadLength_(maxReadLength), tileMetadataList_(tileMetadataList),
        barcodeMetadataList_(barcodeMetadataList),
        contigMap_(contigMap),
        pos_(0ul), pFragment_(0), cigarBegin_(0), cigarEnd_(0), originalCigarBegin_(0), originalCigarEnd_(0),
        forcedDodgyAlignmentScore_(forcedDodgyAlignmentScore),
        flowCellLayoutList_(flowCellLayoutList),
        includeTags_(includeTags),
        pessimisticMapQ_(pessimisticMapQ)
    {
        reserve();
    }

    void reserve()
    {
        readGroupNameBuffer_.resize(100);
        barcodeNameBuffer_.resize(1024);
        readNameBuffer_.resize(getMaxReadNameLength());
        // TODO: 5 is kind of arbitrary. An adjustment of how much realignment can increase the number of operation in a cigar
        originalCigarBuffer_.reserve(5 * alignment::Cigar::getMaxOpeations(maxReadLength_) * alignment::Cigar::OPERATION_CHARS_MAX);
        seqBuffer_.reserve(maxReadLength_);
        qualBuffer_.reserve(maxReadLength_);
    }

    /// For serialization of aligned fragments and shadows
    FragmentAccessorBamAdapter &operator()(
        const PackedFragmentBuffer::Index &index,
        const io::FragmentAccessor& fragment)
    {
        pos_ = index.pos_;
        pFragment_ = &fragment;
        cigarBegin_ = index.cigarBegin_;
        cigarEnd_ = index.cigarEnd_;
        originalCigarBegin_ = fragment.cigarBegin();
        originalCigarEnd_ = fragment.cigarEnd();
        return *this;
    }

    /// For serialization of unaligned fragments
    FragmentAccessorBamAdapter &operator()(
        const io::FragmentAccessor& fragment)
    {
        pos_ = reference::ReferencePosition(reference::ReferencePosition::NoMatch);
        pFragment_ = &fragment;
        cigarBegin_ = 0;
        cigarEnd_ = 0;
        originalCigarBegin_ = 0;
        originalCigarEnd_ = 0;
        return *this;
    }

    static size_t getMaxReadNameLength() {
        // TODO: get this calculated and controlled
        return 1024;
    }

    static size_t getMaxBarcodeNameLength() {
        // TODO: get this calculated and controlled
        return 1024;
    }

    const char *readName() {
        const flowcell::TileMetadata &tileMetadata = tileMetadataList_[pFragment_->tile_];

        readNameBuffer_.clear();
        const std::string &flowcellId = tileMetadata.getFlowcellId();
        readNameBuffer_.insert(readNameBuffer_.end(), flowcellId.begin(), flowcellId.end());
        readNameBuffer_.push_back(':');
        const std::string &laneString = tileMetadata.getLaneString();
        readNameBuffer_.insert(readNameBuffer_.end(), laneString.begin(), laneString.end());
        readNameBuffer_.push_back(':');
        const std::string &tileString = tileMetadata.getTileString();
        readNameBuffer_.insert(readNameBuffer_.end(), tileString.begin(), tileString.end());
        readNameBuffer_.push_back(':');
        common::appendUnsignedNumber(readNameBuffer_, pFragment_->clusterId_);
        readNameBuffer_.push_back(':');
        readNameBuffer_.push_back('0');
        readNameBuffer_.push_back('\0');
        return &readNameBuffer_.front();
    }

    typedef std::pair<const unsigned *, const unsigned *> CigarBeginEnd;

    bool isRealigned() const {
        return cigarBegin_ != originalCigarBegin_;
    }

    bam::zTag getFragmentOC() {
        static const char OC[2] = {'O', 'C'};

        if (!includeTags_.includeOC_ || !isRealigned())
        {
            return bam::zTag();
        }
        originalCigarBuffer_.clear();
        alignment::Cigar::toString(originalCigarBegin_, originalCigarEnd_, originalCigarBuffer_);
        // push terminating zero
        originalCigarBuffer_.push_back(0);
        return bam::zTag(OC, &originalCigarBuffer_.front(), &originalCigarBuffer_.back() + 1);
    }

    CigarBeginEnd cigar() const {
        return std::make_pair(cigarBegin_, cigarEnd_);
    }

    int seqLen() const {
        return pFragment_->readLength_;
    }

    static unsigned char bamBaseFromBclByte(unsigned char bclByte){
        return !oligo::isBclN(bclByte) ? 1 << (bclByte & 0x03) : 15;
    }

    static unsigned char bamBasesFromBclShort(unsigned short bclShort){
        unsigned char *pBclShort(reinterpret_cast<unsigned char*>(&bclShort));
        return bamBaseFromBclByte(pBclShort[0]) << 4 | bamBaseFromBclByte(pBclShort[1]);
    }

    static unsigned char bamQualFromBclByte(unsigned char bclByte){
        return bclByte >> 2;
    }

    const std::vector<unsigned char> &seq() {
        seqBuffer_.resize((seqLen()+1)/2, 15);
        //todo: verify if endianness is not an issue here

        const unsigned short *twoBasesBegin(reinterpret_cast<const unsigned short*>(pFragment_->basesBegin()));
        std::transform(twoBasesBegin, twoBasesBegin + seqLen() / 2, seqBuffer_.begin(), bamBasesFromBclShort);
        if (seqLen() % 2) {
            seqBuffer_.back() = bamBaseFromBclByte(*(pFragment_->basesEnd() - 1)) << 4;
        }
        return seqBuffer_;
    }

    const std::vector<unsigned char> &qual() {
        qualBuffer_.clear();
        std::transform(pFragment_->basesBegin(), pFragment_->basesEnd(), std::back_inserter(qualBuffer_), bamQualFromBclByte);
        return qualBuffer_;
    }

    int refId() const {
        return pos_.isNoMatch() ?
            -1 :
            contigMap_.getMappedContigIndex(barcodeMetadataList_.at(pFragment_->barcode_).getReferenceIndex(),
                                            pos_.getContigId());
    }

    int pos() const {
        return pos_.isNoMatch() ? -1 : pos_.getPosition();
    }

    unsigned char mapq() const {
        if (pFragment_->flags_.properPair_)
        {
            if (pFragment_->DODGY_ALIGNMENT_SCORE == pFragment_->templateAlignmentScore_)
            {
                return  forcedDodgyAlignmentScore_;
            }
            ISAAC_ASSERT_MSG(pFragment_->DODGY_ALIGNMENT_SCORE != pFragment_->alignmentScore_,
                             "Both scores must be either present or missing." << *pFragment_);
            return  std::min<unsigned>(60U,
                                       pessimisticMapQ_ ?
                                           std::min(pFragment_->alignmentScore_, pFragment_->templateAlignmentScore_) :
                                           std::max(pFragment_->alignmentScore_, pFragment_->templateAlignmentScore_));
        }
        return pFragment_->DODGY_ALIGNMENT_SCORE == pFragment_->alignmentScore_ ? (forcedDodgyAlignmentScore_) : std::min<unsigned>(60U, pFragment_->alignmentScore_);
    }

    bam::iTag getFragmentSM() const {
        static const char SM[2] = {'S', 'M'};
        return (!includeTags_.includeSM_ || io::FragmentAccessor::DODGY_ALIGNMENT_SCORE == pFragment_->alignmentScore_) ?
            bam::iTag() : bam::iTag(SM, pFragment_->alignmentScore_);
    }

    bam::iTag getFragmentAS() const {
        static const char AS[2] = {'A', 'S'};
        return (!includeTags_.includeAS_ || !pFragment_->flags_.properPair_ || io::FragmentAccessor::DODGY_ALIGNMENT_SCORE == pFragment_->templateAlignmentScore_) ?
            bam::iTag() :  bam::iTag(AS, pFragment_->templateAlignmentScore_);
    }

    /**
     * \brief Read group IDs are
     */
    bam::zTag getFragmentRG() {
        static const char RG[2] = {'R', 'G'};

        if (!includeTags_.includeRG_)
        {
            return bam::zTag();
        }

        readGroupNameBuffer_.clear();
        const flowcell::BarcodeMetadata &barcode = barcodeMetadataList_[pFragment_->barcode_];
        // barcode index is unique within the data analysis
        common::appendUnsignedInteger(readGroupNameBuffer_, barcode.getIndex());
        readGroupNameBuffer_.push_back('\0');

        return bam::zTag(RG, &readGroupNameBuffer_.front());
    }

    bam::iTag getFragmentNM() const {
        static const char NM[2] = {'N', 'M'};
        return includeTags_.includeNM_ ? bam::iTag(NM, pFragment_->editDistance_) : bam::iTag();
    }

    bam::zTag getFragmentBC() {
        static const char BC[2] = {'B', 'C'};

        if (!includeTags_.includeBC_)
        {
            return bam::zTag();
        }

        barcodeNameBuffer_.clear();
        const std::string &barcodeName =
            barcodeMetadataList_[pFragment_->barcode_].getName();

        if(!flowCellLayoutList_[tileMetadataList_[pFragment_->tile_].getFlowcellIndex()].getBarcodeCycles().size())
        {
            // Use barcode from the sample sheet
            barcodeNameBuffer_.insert(barcodeNameBuffer_.end(), barcodeName.begin(), barcodeName.end());
        }
        else
        {
            // Use barcode from fragment
            oligo::unpackKmer(
                pFragment_->barcodeSequence_,
                flowCellLayoutList_[tileMetadataList_[pFragment_->tile_].getFlowcellIndex()].getBarcodeCycles().size(),
                back_inserter(barcodeNameBuffer_));
        }

        // Null terminate barcode
        barcodeNameBuffer_.push_back('\0');
        return bam::zTag(BC, &barcodeNameBuffer_.front());
    }

    bam::iTag getFragmentZX() const {
        static const char ZX[2] = {'Z', 'X'};
        return (includeTags_.includeZX_ && pFragment_->isClusterXySet()) ? bam::iTag(ZX, pFragment_->clusterX_) :  bam::iTag();
    }

    bam::iTag getFragmentZY() const {
        static const char ZY[2] = {'Z', 'Y'};
        return (includeTags_.includeZY_ && pFragment_->isClusterXySet()) ? bam::iTag(ZY, pFragment_->clusterY_) :  bam::iTag();
    }

    unsigned flag() const {
        std::bitset<11> bs(0);
        bs.set(0, pFragment_->flags_.paired_);
        bs.set(1, pFragment_->flags_.properPair_);
        bs.set(2, pFragment_->flags_.unmapped_);
        bs.set(3, pFragment_->flags_.paired_ && pFragment_->flags_.mateUnmapped_);
        bs.set(4, pFragment_->flags_.reverse_);
        bs.set(5, pFragment_->flags_.mateReverse_);
        bs.set(6, pFragment_->flags_.paired_ && pFragment_->flags_.firstRead_);
        bs.set(7, pFragment_->flags_.paired_ && pFragment_->flags_.secondRead_);

        bs.set(9, pFragment_->flags_.failFilter_);
        bs.set(10, pFragment_->flags_.duplicate_);
        return bs.to_ulong();
    }
    int nextRefId() const {
        return pFragment_->flags_.paired_ ?
            (pFragment_->flags_.unmapped_ && pFragment_->flags_.mateUnmapped_ ?
                -1 : contigMap_.getMappedContigIndex(barcodeMetadataList_.at(pFragment_->barcode_).getReferenceIndex(),
                                                     pFragment_->mateFStrandPosition_.getContigId())) : -1;
    }

    int nextPos() const {
        return pFragment_->flags_.paired_ ?
            (pFragment_->flags_.unmapped_ && pFragment_->flags_.mateUnmapped_ ?
                            -1 : pFragment_->mateFStrandPosition_.getPosition()) : -1;
    }

    int tlen() const { return pFragment_->bamTlen_; }

    int observedLength() const { return pFragment_->observedLength_; }

};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_FRAGMENT_ACCESSOR_BAM_ADAPTER_HH
