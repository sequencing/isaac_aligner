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
 ** \file sortReference.cpp
 **
 ** The main for the reference sorter.
 **
 ** \author Come Raczy
 **/

#include "oligo/Kmer.hh"
#include "options/SortReferenceOptions.hh"
#include "reference/ReferenceSorter.hh"

template <typename KmerT>
void sortReferenceT(const isaac::options::SortReferenceOptions &options)
{
    isaac::reference::ReferenceSorter<KmerT> referenceSorter(
        options.maskWidth,
        options.mask,
        options.genomeFile,
        options.genomeNeighborsFile,
        options.outFile,
        options.repeatThreshold);
    referenceSorter.run();
}

void sortReference(const isaac::options::SortReferenceOptions &options)
{
    if (16 == options.seedLength)
    {
        sortReferenceT<isaac::oligo::ShortKmerType>(options);
    }
    else if (32 == options.seedLength)
    {
        sortReferenceT<isaac::oligo::KmerType>(options);
    }
    else if (64 == options.seedLength)
    {
        sortReferenceT<isaac::oligo::LongKmerType>(options);
    }
    else
    {
        ISAAC_ASSERT_MSG(false, "Unexpected seedLength " << options.seedLength)
    }
}

int main(int argc, char *argv[])
{
    isaac::common::run(sortReference, argc, argv);
}
