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
 ** \file bpbToWig.cpp
 **
 ** \brief Merges multiple reference metadata files into one
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "options/MergeReferencesOptions.hh"
#include "workflow/MergeReferencesWorkflow.hh"

void mergeReferences(const isaac::options::MergeReferencesOptions &options)
{
    isaac::workflow::MergeReferencesWorkflow workflow(
        options.filesToMerge_,
        options.outputFilePath_
        );

    workflow.run();
}

int main(int argc, char *argv[])
{
    isaac::common::run(mergeReferences, argc, argv);
}

