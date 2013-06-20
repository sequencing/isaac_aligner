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
 ** \file SingleEndClusterExtractor.hh
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_SINGLE_END_CLUSTER_EXTRACTOR_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_SINGLE_END_CLUSTER_EXTRACTOR_HH

#include "bam/Bam.hh"
#include "flowcell/ReadMetadata.hh"
#include "bam/BamParser.hh"
#include "reference/ReferencePosition.hh"

//#pragma GCC push_options
//#pragma GCC optimize ("0")

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{
namespace bamDataSource
{

class SingleEndClusterExtractor
{
public:

    static void nothing() {}

    template <typename ClusterInsertIt, typename PfInsertIt>
    bool extractSingleRead(
        const bam::BamBlockHeader &block,
        unsigned &clusterCount,
        const flowcell::ReadMetadataList &readMetadataList,
        ClusterInsertIt &clustersIt,
        PfInsertIt &pfIt)
    {
        ISAAC_ASSERT_MSG(1 == readMetadataList.size(), "Incorrect class used to extract paired data clusters");
        const flowcell::ReadMetadata &readMetadata = readMetadataList[0];
        if ((readMetadata.getNumber()) % 2 == block.isReadOne())
        {
            clustersIt = bam::extractBcl(block, clustersIt, readMetadata);
        }
        *pfIt++ = block.isPf();
        return --clusterCount;
    }

private:
};

} // namespace bamDataSource
} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

//#pragma GCC pop_options


#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_SINGLE_END_CLUSTER_EXTRACTOR_HH
