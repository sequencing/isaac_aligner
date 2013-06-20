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
 ** \file BarcodeMetadata.hh
 **
 ** Packaging of the metadata associated to a barcode.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_BARCODE_METADATA_HH
#define iSAAC_FLOWCELL_BARCODE_METADATA_HH

#include <iostream>
#include <iterator>
#include <vector>

#include "demultiplexing/Barcode.hh"
#include "flowcell/SequencingAdapterMetadata.hh"

namespace isaac
{
namespace flowcell
{

/**
 ** \brief metadata associated to a barcode.
 **
 ** The intended usage is for barcode management in ordered collections (the index
 ** in the collectionis is associated to each barcode metadata instance).
 ** Index 0 is reserved for mapping the barcode sequences that don't match the known ones
 **
 **/
class BarcodeMetadata
{
public:
    static const unsigned INVALID_INDEX = -1U;
    static const unsigned UNMAPPED_REFERENCE_INDEX = -1U;
    static const std::string NO_INDEX_BARCODE;
    static const std::string UNKNOWN_BARCODE;
    static const std::string UNKNOWN_SAMPLE;
    static const std::string DEFAULT_PROJECT;

    /**
     * \brief default constructor creates 'unknown index' barcode
     */
    BarcodeMetadata() : flowcellIndex_(0), lane_(0), referenceIndex_(UNMAPPED_REFERENCE_INDEX), control_(false), unknown_(true), index_(INVALID_INDEX){;}

    BarcodeMetadata(const std::string &flowcellId,
                    const unsigned flowcellIndex,
                    const unsigned lane,
                    const unsigned referenceIndex,
                    const bool unknown,
                    const flowcell::SequencingAdapterMetadataList &adapters) :
                        flowcellId_(flowcellId), flowcellIndex_(flowcellIndex), lane_(0),
                        referenceIndex_(referenceIndex), control_(false),
                        unknown_(unknown), adapters_(adapters), index_(INVALID_INDEX){setLane(lane);}

    static BarcodeMetadata constructUnknownBarcode(
        const std::string &flowcellId,
        const unsigned flowcellIndex,
        const unsigned lane,
        const unsigned referenceIndex,
        const flowcell::SequencingAdapterMetadataList &adapters)
    {
        return BarcodeMetadata(flowcellId, flowcellIndex, lane, referenceIndex, true, adapters);
    }

    static BarcodeMetadata constructNoIndexBarcode(
        const std::string &flowcellId,
        const unsigned flowcellIndex,
        const unsigned lane,
        const unsigned referenceIndex,
        const flowcell::SequencingAdapterMetadataList &adapters)
    {
        BarcodeMetadata ret(flowcellId, flowcellIndex, lane, referenceIndex, false, adapters);
        ret.setSampleName("default");
        return ret;
    }

    const std::string &getFlowcellId() const {return flowcellId_;}
    void setFlowcellId(const std::string &flowcellId) {flowcellId_ = flowcellId;}
    unsigned getFlowcellIndex() const {return flowcellIndex_;}
    void setFlowcellIndex(const unsigned flowcellIndex) {flowcellIndex_ = flowcellIndex;}

    void setLane(const unsigned lane) {lane_ = lane; laneSampleName_ = "lane" + boost::lexical_cast<std::string>(lane);}
    unsigned getLane() const {return lane_;}

    void setSampleName(const std::string &sampleName) {sampleName_ = sampleName; unknown_ = false;}
    const std::string &getSampleName() const
    {
        return isUnknown() ? UNKNOWN_SAMPLE : sampleName_.empty() ? laneSampleName_ : sampleName_;
    }

    void setAdapters(const flowcell::SequencingAdapterMetadataList &adapters) {adapters_ = adapters;}
    const flowcell::SequencingAdapterMetadataList &getAdapters() const {return adapters_;}

    const std::vector<unsigned> &getComponentMismatches() const {return componentMismatches_;}
    void setComponentMismatches(const std::vector<unsigned> &componentMismatches)
    {
        componentMismatches_ = componentMismatches;
        if(getComponentsCount() > componentMismatches_.size())
        {
            componentMismatches_.insert(componentMismatches_.end(),
                                        getComponentsCount() - componentMismatches_.size(),
                                        componentMismatches_.back());
        }
    }

    const std::string& getSequence() const {return sequence_;}
    unsigned getSequenceLength() const
    {
        return std::count_if(sequence_.begin(), sequence_.end(),
                             boost::bind(&boost::cref<char>, _1) != '-');
    }
    /**
     * \brief Sets the sequence, resets the isUnknown flag
     */
    void setSequence(const std::string& sequence) {sequence_ = sequence; unknown_ = false;}
    unsigned getComponentsCount() const {return std::count(sequence_.begin(), sequence_.end(), '-') + 1;}

    void setReference(const std::string &reference) {reference_ = reference;}
    const std::string &getReference() const {return reference_;}

    unsigned getReferenceIndex() const {return referenceIndex_;}
    bool isUnmappedReference() const {return UNMAPPED_REFERENCE_INDEX == referenceIndex_;}
    void setReferenceIndex(const unsigned referenceIndex) {referenceIndex_ = referenceIndex;}

    void setDescription(const std::string &description) {description_ = description;}
    const std::string &getDescription() const {return description_;}

    void setControl(const bool control) {control_ = control;}

    void setRecipe(const std::string &recipe) {recipe_ = recipe;}
    const std::string &getRecipe() const {return recipe_;}

    void setOperator(const std::string &oprtr) {operator_ = oprtr;}
    const std::string &getOperator() const {return operator_;}

    void setProject(const std::string &project) {project_ = project;}
    const std::string &getProject() const {return project_.empty() ? DEFAULT_PROJECT : project_;}

    unsigned getIndex() const {return index_;}
    void setIndex(unsigned index) {index_ = index;}
    bool isUnknown() const {return unknown_;}
    void setUnknown()
    {
        if (!sampleName_.empty())
        {
            BOOST_THROW_EXCEPTION(common::PreConditionException(
                (boost::format("ERROR: Sample name must be empty for 'unknown' barcodes. Actual: %s, lane: %d") %
                    sampleName_ % lane_).str()));
        }
        unknown_ = true;
    }
    bool isNoIndex() const {return !isUnknown() && sequence_.empty();}

    bool isDefault() const {return isUnknown() || isNoIndex();}

    const std::string &getName() const {
        return isUnknown() ? UNKNOWN_BARCODE : isNoIndex() ? NO_INDEX_BARCODE : sequence_;
    }

private:
    // string identifier of the flowcell
    std::string flowcellId_;
    // index of the flowcell in the global list of flowcells
    unsigned flowcellIndex_;
    // lane number within the flowcell
    unsigned lane_;
    std::string sampleName_;
    std::vector<unsigned> componentMismatches_;
    std::string sequence_;
    std::string reference_;
    unsigned referenceIndex_;
    std::string description_;
    bool control_;
    bool unknown_;
    std::string recipe_;
    std::string operator_;
    std::string project_;
    std::string laneSampleName_;
    flowcell::SequencingAdapterMetadataList adapters_;
    unsigned index_;
};

typedef std::vector<BarcodeMetadata> BarcodeMetadataList;

inline std::ostream &operator<<(std::ostream &os, const BarcodeMetadata &barcodeMetadata)
{
    os << "BarcodeMetadata(" << barcodeMetadata.getFlowcellId() << "," << barcodeMetadata.getLane() << "," <<
        barcodeMetadata.getSampleName() << "," << barcodeMetadata.getName() << "," << barcodeMetadata.getReference();
    if (!barcodeMetadata.isUnmappedReference())
    {
        os << "(" << barcodeMetadata.getReferenceIndex() << ")";
    }
    else
    {
        os << "(unmapped)";
    }

    os << ", " << barcodeMetadata.getIndex() << ")";
    return os;
}

} // namespace flowcell
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_BARCODE_METADATA_HH

