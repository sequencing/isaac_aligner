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
 ** \file SortedReferenceXmlBamHeaderAdapter.hh
 **
 ** Implements interface required for generating bam header out of SortedReference.xml
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_SORTED_REFERENCE_XML_BAM_ADAPTER_HH
#define iSAAC_BUILD_SORTED_REFERENCE_XML_BAM_ADAPTER_HH

#include "reference/SortedReferenceXml.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace build
{

class SortedReferenceXmlBamHeaderAdapter
{
    const reference::SortedReferenceXml &sortedReferenceXml_;
    const flowcell::TileMetadataList &tileMetadataList_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const std::string &sampleName_;
public:
    SortedReferenceXmlBamHeaderAdapter(
        const reference::SortedReferenceXml &sortedReferenceXml,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const std::string &sampleName):
        sortedReferenceXml_(sortedReferenceXml),
        tileMetadataList_(tileMetadataList),
        barcodeMetadataList_(barcodeMetadataList),
        sampleName_(sampleName)
        {}

    class RefSequence
    {
        const reference::SortedReferenceXml::Contig& contig_;
    public:
        RefSequence(const reference::SortedReferenceXml::Contig &contig):
            contig_(contig){}
        const std::string &name() const {return contig_.name_;}
        int length() const {return contig_.totalBases_;}
        const std::string &bamSqAs() const {return contig_.bamSqAs_;}
        const std::string &bamSqUr() const {return contig_.bamSqUr_.empty() ? contig_.filePath_.string() : contig_.bamSqUr_;}
        const std::string &bamM5() const {return contig_.bamM5_;}
    };

    int getRefSequenceCount() const {
        return sortedReferenceXml_.getContigs().size();
    }

    typedef RefSequence RefSeqType;
    typedef std::vector<reference::SortedReferenceXml::Contig> RefSeqsType;
    RefSeqsType getRefSequences() const {
        RefSeqsType ret = sortedReferenceXml_.getKaryotypeOrderedContigs();
        return ret;
    }

    struct ReadGroupType : public std::map<std::string, std::string>::value_type
    {
        ReadGroupType(const std::map<std::string, std::string>::value_type &that):
            std::map<std::string, std::string>::value_type(that){}
        const std::string &getId() const {return first;}
        const std::string &getValue() const {return second;}
    };

    std::map<std::string, std::string> getReadGroups() const
    {
        std::map<std::string, std::string> ret;
        std::string readGroupDefinitions;
        BOOST_FOREACH(const flowcell::TileMetadata &tile, tileMetadataList_)
        {
            BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList_)
            {
                if (tile.getFlowcellId() == barcode.getFlowcellId() && tile.getLane() == barcode.getLane())
                {
                    // barcode index is unique within the data anaysis
                    const std::string readGroupId = boost::lexical_cast<std::string>(barcode.getIndex());
                    const std::string paltformUnit = tile.getFlowcellId() + ":" + tile.getLaneString() + ":" + barcode.getName();

                    std::string readGroup =
                        "@RG\tID:" + readGroupId + "\tPL:ILLUMINA\tSM:" + sampleName_ + "\tPU:"+ paltformUnit;
                    ret.insert(std::make_pair(readGroupId, readGroup));
                }
            }
        }

        return ret;
    }

};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_SORTED_REFERENCE_XML_BAM_ADAPTER_HH
