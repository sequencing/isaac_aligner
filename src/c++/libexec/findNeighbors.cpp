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
 ** \file findNeighbors.cpp
 **
 ** \brief The main for the identification of the neighbors (32-mers with 1 or 2 mismatches).
 **
 ** \author Come Raczy
 **/

#include "options/FindNeighborsOptions.hh"
#include "reference/NeighborsFinder.hh"

template <typename KmerT>
void findNeighborsT(const isaac::options::FindNeighborsOptions &options)
{
    isaac::reference::NeighborsFinder<KmerT> neighborsFinder(
        options.parallelSort,
        options.inputFile,
        options.outputDirectory,
        options.outputFile,
        options.tempFile,
        options.jobs);
    neighborsFinder.run();
}

void findNeighbors(const isaac::options::FindNeighborsOptions &options)
{
    if (16 == options.seedLength)
    {
        findNeighborsT<isaac::oligo::ShortKmerType>(options);
    }
    else if (32 == options.seedLength)
    {
        findNeighborsT<isaac::oligo::KmerType>(options);
    }
    else if (64 == options.seedLength)
    {
        findNeighborsT<isaac::oligo::LongKmerType>(options);
    }
    else
    {
        ISAAC_ASSERT_MSG(false, "Unexpected seedLength " << options.seedLength)
    }
}

int main(int argc, char *argv[])
{
    isaac::common::run(findNeighbors, argc, argv);
}

