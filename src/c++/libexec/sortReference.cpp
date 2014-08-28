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
