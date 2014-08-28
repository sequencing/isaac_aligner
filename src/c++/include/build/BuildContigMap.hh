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
 ** \file BuildContigMap.hh
 **
 ** Helper for avoiding to deal with reference contigs that don't have any mapped reads.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BUILD_CONTIG_MAP_HH
#define iSAAC_BUILD_BUILD_CONTIG_MAP_HH

#include "alignment/BinMetadata.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "reference/SortedReferenceMetadata.hh"


namespace isaac
{
namespace build
{

class BuildContigMap : std::vector<std::vector<unsigned> >
{
    static const unsigned UNMAPPED_CONTIG = -1U;
public:

    BuildContigMap(
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const alignment::BinMetadataCRefList &bins,
        const isaac::reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const bool skipEmptyBins)
    {
        reserve(sortedReferenceMetadataList.size());
        // initialize all with 0
        BOOST_FOREACH(const isaac::reference::SortedReferenceMetadata &referenceMetadata, sortedReferenceMetadataList)
        {
            push_back(std::vector<unsigned>(referenceMetadata.getContigsCount(), UNMAPPED_CONTIG));
        }

        // set reference contig value to non-zero if there is at least one record in it
        BOOST_FOREACH(const alignment::BinMetadata &bin, bins)
        {
            ISAAC_THREAD_CERR << bin << std::endl;
            BOOST_FOREACH(const flowcell::BarcodeMetadata &barcodeMetadata, barcodeMetadataList)
            {
                if (!bin.isUnalignedBin() && !barcodeMetadata.isUnmappedReference() &&
                    (!skipEmptyBins || bin.getBarcodeElements(barcodeMetadata.getIndex())))
                {
                    at(barcodeMetadata.getReferenceIndex()).at(bin.getBinStart().getContigId()) = 0;
                }
            }
        }

        for (unsigned referenceIndex = 0; sortedReferenceMetadataList.size() > referenceIndex; ++ referenceIndex)
        {
            unsigned loadedContigIndex = 0;
            BOOST_FOREACH(unsigned &referenceContigMapping, at(referenceIndex))
            {
                if (UNMAPPED_CONTIG != referenceContigMapping)
                {
                    referenceContigMapping = loadedContigIndex;
                    ++loadedContigIndex;
                }
            }
            ISAAC_THREAD_CERR << "Will load " << loadedContigIndex << " contigs for reference id " << referenceIndex << std::endl;
        }
    }

    unsigned getMappedContigIndex(const unsigned referenceIndex, const unsigned referenceContigIndex) const
    {
        return at(referenceIndex).at(referenceContigIndex);
    }

    bool isMapped(const unsigned referenceIndex, const unsigned referenceContigIndex) const
    {
        return UNMAPPED_CONTIG != at(referenceIndex).at(referenceContigIndex);
    }
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BUILD_CONTIG_MAP_HH
