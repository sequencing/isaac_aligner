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
 ** \file BpbToWigWorkflow.hh
 **
 ** \brief prints wig file given the bitset and contigs map
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_BPB_TO_WIG_WORKFLOW_HH
#define iSAAC_WORKFLOW_BPB_TO_WIG_WORKFLOW_HH

#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class BpbToWigWorkflow: boost::noncopyable
{
private:
    const bfs::path sortedReferenceMetadata_;
    const bfs::path inputFilePath_;
    const std::string outputFormatString_;

    reference::SortedReferenceMetadata xml_;

public:
    BpbToWigWorkflow(
        const bfs::path &sortedReferenceMetadata,
        const bfs::path &inputFilePath,
        const std::string &outputFormatString
        );

    void run();
private:
    void printWig(
        std::ifstream &bitsetFile,
        const reference::SortedReferenceMetadata::Contigs &contigs );
    void printBed(
        std::ifstream &bitsetFile,
        const reference::SortedReferenceMetadata::Contigs &contigs );

};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_BPB_TO_WIG_WORKFLOW_HH
