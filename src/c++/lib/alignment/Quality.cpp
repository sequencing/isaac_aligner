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
 ** \file Quality.cpp
 **
 ** Various functions and tables to support alignment quality.
 **
 ** \author Come Raczy
 **/

#include <vector>

#include "alignment/Quality.hh"
#include "common/MathCompatibility.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Create the lookup table associating quality scores (index) to the
 ** correcponding probability of a matching base (value).
 **/
std::vector<double> getLogMatchLookup()
{
    std::vector<double> lookup;

    // Ns in data (which come in Q0 from Bcl) need special treatment
    // or else they wreck the logProbability of the perfect alignment. Assume they are Q1 Ns
    const double nMismatch = pow(10.0, 1.0 / -10.0);
    lookup.push_back(log(1.0 - nMismatch));

    for(int i = 1; i < 100; ++i)
    {
        const double mismatch = pow(10.0, (double)i / -10.0);
        lookup.push_back(log(1.0 - mismatch));
    }
    return lookup;
}

/**
 ** \brief Create the lookup table associating quality scores (index) to the
 ** correcponding probability of a mismatching base (value).
 **/
std::vector<double> getLogMismatchLookup()
{
    std::vector<double> lookup;
    // prevent the logarithmic singularity
    lookup.push_back(log(1.0 - pow(10.0, 1.0 / -10.0)));
    for(unsigned quality = 1; quality < 100U; ++quality)
    {
        const double logMismatch = Quality::getLogMismatch(quality);
        lookup.push_back(logMismatch);
    }
    return lookup;
}

const std::vector<double> Quality::logMatchLookup = getLogMatchLookup();
const std::vector<double> Quality::logMismatchLookup = getLogMismatchLookup();

const unsigned MASK_READ_LENGTH_MIN = 35;
void trimLowQualityEnd(Read &read, const unsigned baseQualityCutoff)
{
    if (read.getLength() < MASK_READ_LENGTH_MIN)
    {
        return;
    }

    const std::vector<char> &reverse = read.getReverseQuality();
    int qscoreSum = 0;
    int peakSum = 0;
    std::vector<char>::const_iterator trimPos = reverse.begin();

    for (std::vector<char>::const_iterator it = reverse.begin();
        reverse.end() - MASK_READ_LENGTH_MIN != it; ++it)
    {
        qscoreSum += baseQualityCutoff - *it;
        if (qscoreSum < 0)
        {
            break;
        }

        if (qscoreSum > peakSum)
        {
            peakSum = qscoreSum;
            trimPos = it;
        }
    }

    read.maskCyclesFromEnd(trimPos - reverse.begin());
    ISAAC_ASSERT_MSG(read.getEndCyclesMasked() < read.getLength(), "TAda");
}

void trimLowQualityEnds(Cluster &cluster, const unsigned baseQualityCutoff)
{
    if (!baseQualityCutoff)
    {
        return;
    }

    for (unsigned readIndex = 0; cluster.getNonEmptyReadsCount() > readIndex; ++readIndex)
    {
        trimLowQualityEnd(cluster[readIndex], baseQualityCutoff);
    }
}

} // namespace alignemnt
} // namespace isaac
