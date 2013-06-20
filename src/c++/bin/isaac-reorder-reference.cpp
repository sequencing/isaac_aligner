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
 ** \file isaac-reorder-reference.cpp
 **
 ** \brief User-facing executable for changing the order of contigs in reference
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "options/ReorderReferenceOptions.hh"
#include "workflow/ReorderReferenceWorkflow.hh"

void reoderReference(const isaac::options::ReorderReferenceOptions &options)
{
    isaac::workflow::ReorderReferenceWorkflow workflow(
        options.sortedReferenceMetadata_,
        options.newXmlPath_,
        options.newFaPath_,
        options.newOrder_,
        options.basesPerLine_
        );

    workflow.run();
}

int main(int argc, char *argv[])
{
    isaac::common::run(reoderReference, argc, argv);
}

