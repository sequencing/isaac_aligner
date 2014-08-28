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
 ** \file ReorderReferenceWorkflow.hh
 **
 ** \brief Top level component to control the reference reordering process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_REORDER_REFERENCE_WORKFLOW_HH
#define iSAAC_WORKFLOW_REORDER_REFERENCE_WORKFLOW_HH

#include "common/Threads.hpp"
#include "reference/Contig.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace workflow
{

namespace bfs = boost::filesystem;

class ReorderReferenceWorkflow: boost::noncopyable
{
public:
    ReorderReferenceWorkflow(
        const bfs::path &sortedReferenceMetadata,
        const bfs::path &newXmlPath,
        const bfs::path &newFaPath,
        const std::vector<std::string> &newOrder,
        const unsigned basesPerLine
        );

    void run();

private:
    const bfs::path sortedReferenceMetadata_;
    const bfs::path newXmlPath_;
    const bfs::path newFaPath_;
    const std::vector<std::string> &newOrder_;
    const unsigned basesPerLine_;
    common::ThreadVector threads_;

    reference::SortedReferenceMetadata xml_;
    bool orderByKaryotypeIndex(const reference::Contig& left, const reference::Contig& right);
    void storeContig(std::ostream &os, const reference::Contig &contig);
    void writeBase(std::ostream &os, const char base, const bool writeNewline);

};
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_REORDER_REFERENCE_WORKFLOW_HH
