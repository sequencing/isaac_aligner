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
 ** \file ExtractNeighborsWorkflow.hh
 **
 ** \brief Top level component to control the neighbor extraction process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_EXTRACT_NEIGHBORS_WORKFLOW_HH
#define iSAAC_WORKFLOW_EXTRACT_NEIGHBORS_WORKFLOW_HH

#include "common/Threads.hpp"
#include "reference/Contig.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class ExtractNeighborsWorkflow: boost::noncopyable
{
private:
    const bfs::path sortedReferenceMetadata_;
    const bfs::path neighborsFilePath_;
    const bfs::path highRepeatsFilePath_;
    common::ThreadVector threads_;

    reference::SortedReferenceMetadata xml_;

public:
    ExtractNeighborsWorkflow(
        const bfs::path &sortedReferenceMetadata,
        const bfs::path &neighborsFilePath,
        const bfs::path &highRepeatsFilePath
        );

    template <typename KmerT>
    void run();

private:
    template <typename KmerT>
    void scanMaskFile(
        const reference::SortedReferenceMetadata::MaskFile &maskFile,
        const std::vector<unsigned long> &contigOffsets,
        std::vector<bool> &neighbors,
        std::vector<bool> &highRepeats);
    void dumpResults(const std::vector<bool> &neighbors, const std::vector<bool> &highRepeats);
};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_EXTRACT_NEIGHBORS_WORKFLOW_HH
