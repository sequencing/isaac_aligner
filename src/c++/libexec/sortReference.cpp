/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/downloads/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file sortReference.cpp
 **
 ** The main for the reference sorter.
 **
 ** \author Come Raczy
 **/

#include "options/SortReferenceOptions.hh"
#include "reference/ReferenceSorter.hh"

void sortReference(const isaac::options::SortReferenceOptions &options);

int main(int argc, char *argv[])
{
    isaac::common::run(sortReference, argc, argv);
}

void sortReference(const isaac::options::SortReferenceOptions &options)
{
    isaac::reference::ReferenceSorter referenceSorter(
        options.maskWidth,
        options.mask,
        options.genomeFile,
        options.genomeNeighborsFile,
        options.permutationName,
        options.outFile,
        options.repeatThreshold);
    referenceSorter.run();
}
