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
 ** \file MergeReferencesWorkflow.hh
 **
 ** \brief merge multiple reference metadata into one
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_MERGE_REFERENE_WORKFLOW_HH
#define iSAAC_WORKFLOW_MERGE_REFERENE_WORKFLOW_HH

#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class MergeReferencesWorkflow: boost::noncopyable
{
private:
    const std::vector<bfs::path> &filesToMerge_;
    const bfs::path &outputFilePath_;

public:
    MergeReferencesWorkflow(
        const std::vector<bfs::path> &filesToMerge,
        const bfs::path &outputFilePath);

    void run();
};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_MERGE_REFERENE_WORKFLOW_HH
