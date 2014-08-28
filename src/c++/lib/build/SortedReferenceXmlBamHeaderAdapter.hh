/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file SortedReferenceXmlBamHeaderAdapter.hh
 **
 ** Implements interface required for generating bam header out of SortedReference.xml
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_SORTED_REFERENCE_XML_BAM_ADAPTER_HH
#define iSAAC_BUILD_SORTED_REFERENCE_XML_BAM_ADAPTER_HH

#include "common/Strings.hh"
#include "reference/SortedReferenceMetadata.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace build
{

template <typename IncludeContigF> class SortedReferenceXmlBamHeaderAdapter
{
    const reference::SortedReferenceMetadata &sortedReferenceMetadata_;
    const IncludeContigF &includeContig_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const std::string &sampleName_;
public:
    SortedReferenceXmlBamHeaderAdapter(
        const reference::SortedReferenceMetadata &sortedReferenceMetadata,
        const IncludeContigF &includeContig,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::string &sampleName):
        sortedReferenceMetadata_(sortedReferenceMetadata),
        includeContig_(includeContig),
        tileMetadataList_(tileMetadataList),
        barcodeMetadataList_(barcodeMetadataList),
        sampleName_(sampleName)
        {}

    class RefSequence
    {
        const reference::SortedReferenceMetadata::Contig& contig_;
    public:
        RefSequence(const reference::SortedReferenceMetadata::Contig &contig):
            contig_(contig){}
        const std::string &name() const {return contig_.name_;}
        int length() const {return contig_.totalBases_;}
        const std::string &bamSqAs() const {return contig_.bamSqAs_;}
        const std::string &bamSqUr() const {return contig_.bamSqUr_.empty() ? contig_.filePath_.string() : contig_.bamSqUr_;}
        const std::string &bamM5() const {return contig_.bamM5_;}
    };

    typedef RefSequence RefSeqType;
    typedef std::vector<reference::SortedReferenceMetadata::Contig> RefSeqsType;
    const RefSeqsType getRefSequences() const {
        const RefSeqsType ret = sortedReferenceMetadata_.getKaryotypeOrderedContigs(includeContig_);
        return ret;
    }

    struct ReadGroupType : public std::map<std::string, std::string>::value_type
    {
        ReadGroupType(const std::map<std::string, std::string>::value_type &that):
            std::map<std::string, std::string>::value_type(that){}
        const std::string &getId() const {return first;}
        const std::string &getValue() const {return second;}
    };

    std::map<std::string, std::string> getReadGroups(const std::string &bamPuFormat) const
    {
        std::map<std::string, std::string> ret;
        std::string readGroupDefinitions;
        BOOST_FOREACH(const flowcell::TileMetadata &tile, tileMetadataList_)
        {
            BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList_)
            {
                if (barcode.getSampleName() == sampleName_ &&
                    tile.getFlowcellId() == barcode.getFlowcellId() && tile.getLane() == barcode.getLane())
                {
                    // barcode index is unique within the data anaysis
                    const std::string readGroupId = boost::lexical_cast<std::string>(barcode.getIndex());
                    std::string paltformUnit = bamPuFormat;
                    common::replaceSubstring(paltformUnit, "%F", tile.getFlowcellId());
                    common::replaceSubstring(paltformUnit, "%L", tile.getLaneString());
                    common::replaceSubstring(paltformUnit, "%B", barcode.getName());

                    std::string readGroup =
                        "@RG\tID:" + readGroupId + "\tPL:ILLUMINA\tSM:" + sampleName_ + "\tPU:"+ paltformUnit;
                    ret.insert(std::make_pair(readGroupId, readGroup));
                }
            }
        }

        return ret;
    }

};

template <typename IncludeContigF> SortedReferenceXmlBamHeaderAdapter<IncludeContigF> makeSortedReferenceXmlBamHeaderAdapter(
        const reference::SortedReferenceMetadata &sortedReferenceMetadata,
        const IncludeContigF &includeContig,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::string &sampleName)
{
    return SortedReferenceXmlBamHeaderAdapter<IncludeContigF>(sortedReferenceMetadata, includeContig, tileMetadataList, barcodeMetadataList, sampleName);
}


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_SORTED_REFERENCE_XML_BAM_ADAPTER_HH
