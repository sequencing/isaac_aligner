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
 ** \brief Extract neighbor flags from kmers and store them in a file where
 **        each bit would correspond to a position in the reference
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "options/BpbToWigOptions.hh"
#include "workflow/BpbToWigWorkflow.hh"

void bpbToWig(const isaac::options::BpbToWigOptions &options)
{
    isaac::workflow::BpbToWigWorkflow workflow(
        options.sortedReferenceMetadata_,
        options.inputFilePath_,
        options.outputFormatString_
        );

    workflow.run();
}

int main(int argc, char *argv[])
{
    isaac::common::run(bpbToWig, argc, argv);
}

