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
 ** \file BamSerializer.hh
 **
 ** Helper class for converting BinSorter data to bam.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BAM_SERIALIZER_HH
#define iSAAC_BUILD_BAM_SERIALIZER_HH

#include <boost/iostreams/filtering_stream.hpp>

#include "bam/Bam.hh"
#include "bam/BamIndexer.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/FragmentAccessorBamAdapter.hh"
#include "build/FragmentIndex.hh"
#include "build/PackedFragmentBuffer.hh"
#include "flowcell/TileMetadata.hh"


namespace isaac
{
namespace build
{

class BamSerializer
{
public:
    BamSerializer(
        const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeOutputFileIndexMap,
        const flowcell::TileMetadataList &tileMetadataList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const BuildContigMap &contigMap,
        const unsigned maxReadLength,
        const unsigned char forcedDodgyAlignmentScore,
        const flowcell::FlowcellLayoutList &flowCellLayoutList,
        const IncludeTags includeTags):
            barcodeOutputFileIndexMap_(barcodeOutputFileIndexMap),
            bamAdapter_(maxReadLength, tileMetadataList, barcodeMetadataList, contigMap, forcedDodgyAlignmentScore, flowCellLayoutList, includeTags)
    {}

    typedef void result_type;
    void operator()(const PackedFragmentBuffer::Index& idx,
                    boost::ptr_vector<boost::iostreams::filtering_ostream> &streams,
                    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts,
                    const PackedFragmentBuffer &fragmentData)
    {
        const io::FragmentAccessor &fragment = fragmentData.getFragment(idx);
        FragmentAccessorBamAdapter& adapter = bamAdapter_(idx, fragment);
        unsigned serializedLength = bam::serializeAlignment(
            streams.at(barcodeOutputFileIndexMap_.at(fragment.barcode_)), adapter);
        bam::BamIndexPart& bamIndexPart = bamIndexParts.at(barcodeOutputFileIndexMap_.at(fragment.barcode_));
        bamIndexPart.processFragment( adapter, serializedLength );
//        ISAAC_THREAD_CERR << "Serialized to bam pos_: " << idx.pos_ << " dataOffset_: " << idx.dataOffset_ << std::endl;
    }

    void operator()(const io::FragmentAccessor &fragment,
                    boost::ptr_vector<boost::iostreams::filtering_ostream> &streams,
                    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts)
    {
        FragmentAccessorBamAdapter& adapter = bamAdapter_(fragment);
        unsigned serializedLength = bam::serializeAlignment(
            streams.at(barcodeOutputFileIndexMap_.at(fragment.barcode_)), adapter);
        bam::BamIndexPart& bamIndexPart = bamIndexParts.at(barcodeOutputFileIndexMap_.at(fragment.barcode_));
        bamIndexPart.processFragment( adapter, serializedLength );
//        ISAAC_THREAD_CERR << "Serialized unaligned pos_: " << fragment << std::endl;
    }
private:
    const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeOutputFileIndexMap_;
    FragmentAccessorBamAdapter bamAdapter_;
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BAM_SERIALIZER_HH
