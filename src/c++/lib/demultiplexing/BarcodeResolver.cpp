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
 ** \file BarcodeResolver.cpp
 **
 ** Performs translation from barcode sequences to barcode indexes. Allows for sequence mismatches.
 **
 ** \author Roman Petrovski
 **/

#include "demultiplexing/BarcodeResolver.hh"

namespace isaac
{
namespace demultiplexing
{

bool orderBySequenceAndBarcode(const Barcode &left, const Barcode &right)
{
    return left.getSequence() < right.getSequence() ||
        (left.getSequence() == right.getSequence() && left.getBarcode() < right.getBarcode());
}

/**
 * \return Returns true if the sequence and barcode match.
 *
 * \throws common::InvalidOptionException Throws an exception if the sequence matches but the
 *                                        barcode index does not. This is an indication
 *                                        of barcode collision.
 */
bool isSequenceBarcodeEqual(const flowcell::BarcodeMetadataList &allBarcodeMetadata,
                            const Barcode &left, const Barcode &right)
{
    if (left.getSequence() != right.getSequence())
    {
        return false;
    }

    if (left.getBarcode() != right.getBarcode())
    {
        BOOST_THROW_EXCEPTION(
            common::InvalidOptionException("Barcode collision detected. Barcode "
                + boost::lexical_cast<std::string>(left) + " produced from "
                    + boost::lexical_cast<std::string>(allBarcodeMetadata.at(left.getBarcode()))
                + " collides with "
                + boost::lexical_cast<std::string>(right) + " produced from "
                    + boost::lexical_cast<std::string>(allBarcodeMetadata.at(right.getBarcode()))));
    }
    return true;
}


bool orderBySequence(const Barcode &left, const Barcode &right)
{
    return left.getSequence() < right.getSequence();
}

/**
 * \return Number of single-mismatch variations for a kmer of kmerLength. Original kmer is included in the count.
 */
unsigned get1MismatchKmersCount(const unsigned kmerLength)
{
    return kmerLength * (oligo::invalidOligo + 1);
}

/**
 * \return Number of mismatch variations for a kmer of kmerLength. Duplicates are not excluded from the count.
 */
unsigned BarcodeResolver::getMismatchKmersCount(const unsigned kmerLength, unsigned maxMismatches)
{
    unsigned ret = 1;
    const unsigned oneMismatchKmersCount = get1MismatchKmersCount(kmerLength);
    while(maxMismatches--)
    {
        ret *= oneMismatchKmersCount;
    }
    return ret;
}

/**
 * \param iteration Simulation iteration index when simulating all 1-mismatch variants
 *
 * \return Pair of kmer and the edit distance from the original (0 or 1)
 */
std::pair<Kmer, unsigned> BarcodeResolver::get1MismatchKmer(
    const Kmer original,
    const unsigned kmerLength,
    const unsigned componentOffset,
    const unsigned iteration)
{
    // how many values a single base can take (in our case 5)
    const unsigned baseVariants = (oligo::invalidOligo + 1);

    const unsigned iterationMismatchPosition = componentOffset + (iteration / baseVariants) % kmerLength;
    const Kmer mismatchBase = iteration % baseVariants;
    const Kmer shiftedMask = kmerMask_ << (iterationMismatchPosition * BITS_PER_BASE);
    const Kmer invertedShiftedMask = ~shiftedMask;
    const Kmer shiftedMismatchBase = (mismatchBase << (iterationMismatchPosition * BITS_PER_BASE));
    const Kmer retKmer = (original & invertedShiftedMask) | shiftedMismatchBase;

    const Kmer originalBase = original & shiftedMask;

    return std::make_pair(retKmer, unsigned(originalBase != shiftedMismatchBase));
}

/**
 * \param iteration Simulation iteration index when simulating all 2-mismatch variants
 *
 * \return Pair of kmer and the edit distance from the original (0, 1 or 2)
 */
std::pair<Kmer, unsigned> BarcodeResolver::get2MismatchKmer(
    const Kmer original,
    const unsigned kmerLength,
    const unsigned componentOffset,
    const unsigned iteration)
{
    const unsigned level0Iteration = iteration / get1MismatchKmersCount(kmerLength);
    const unsigned level1Iteration = iteration % get1MismatchKmersCount(kmerLength);

    const std::pair<Kmer, unsigned> level0MismatchKmer = get1MismatchKmer(original, kmerLength, componentOffset, level0Iteration);
    std::pair<Kmer, unsigned> ret = get1MismatchKmer(level0MismatchKmer.first, kmerLength, componentOffset, level1Iteration);
    ret.second += level0MismatchKmer.second;

    return ret;
}

std::pair<Kmer, unsigned> generateMismatchKmer(
    const Kmer original,
    const std::vector<unsigned> &componentLengths,
    const std::vector<unsigned> &mismatchesPerComponent,
    std::vector<unsigned> &allComponentIterations,
    unsigned iteration,
    std::vector<Barcode> &result)
{
    std::pair<Kmer, unsigned> ret = std::make_pair(original, 0U);
    unsigned componentOffset = 0;
    BOOST_REVERSE_FOREACH(const unsigned &componentIterations, allComponentIterations)
    {
        const unsigned componentIndex = &componentIterations - &allComponentIterations.front();
        const unsigned thisComponentIteration = iteration % componentIterations;
        iteration /= componentIterations;
        switch (mismatchesPerComponent.at(componentIndex))
        {
        case 0:
        {
            break;
        }
        case 1:
        {
            const std::pair<Kmer, unsigned> mismatchKmer =
                BarcodeResolver::get1MismatchKmer(ret.first, componentLengths.at(componentIndex),
                                                  componentOffset, thisComponentIteration);
            ret.first = mismatchKmer.first;
            ret.second += mismatchKmer.second;
            break;
        }
        case 2:
        {
            const std::pair<Kmer, unsigned> mismatchKmer =
                BarcodeResolver::get2MismatchKmer(ret.first, componentLengths.at(componentIndex),
                                                  componentOffset, thisComponentIteration);
            ret.first = mismatchKmer.first;
            ret.second += mismatchKmer.second;
            break;
        }
        default:
        {
            ISAAC_ASSERT_MSG(false, "Only 0,1 and 2 mismatches simulation is supported");
            break;
        }
        }
        componentOffset += componentLengths.at(componentIndex);
    }
    return ret;
}

/**
 * \param mismatchesPerComponent Vector of mismatch counts to be used for each barcode component.
 *                               Expected to be enough to cover all barcode components.
 */
void BarcodeResolver::generateBarcodeMismatches(
    const flowcell::BarcodeMetadata &barcodeMetadata,
    std::vector<Barcode> &result)
{
    const std::string &sequence = barcodeMetadata.getSequence();
    ISAAC_ASSERT_MSG(!sequence.empty(), "only default barcode can have an empty sequence and it must not be passed here");

    static const oligo::Translator translator = oligo::getTranslator(true, oligo::invalidOligo);

    // decode kmer from string and note all component lengths and how many mismatch sequences can each
    // component produce given its length and command-line mismatch count
    Kmer kmer = 0;
    std::vector<unsigned> componentLengths;
    std::vector<unsigned> allComponentIterations;
    unsigned mismatchIterations = 1;
    {
        unsigned componentLength = 0;
        BOOST_FOREACH(const char base, sequence)
        {
            if ('-' != base)
            {
                kmer = (kmer << BITS_PER_BASE) | translator[base];
                ++componentLength;
            }
            else
            {
                componentLengths.push_back(componentLength);
                allComponentIterations.push_back(
                    getMismatchKmersCount(componentLength, barcodeMetadata.getComponentMismatches().at(allComponentIterations.size())));
                mismatchIterations *= allComponentIterations.back();
                componentLength = 0;
            }
        }
        ISAAC_ASSERT_MSG(componentLength, "barcode cannot end with '-' or be empty, so, last component length cannot be 0");
        componentLengths.push_back(componentLength);
        allComponentIterations.push_back(
            getMismatchKmersCount(componentLength, barcodeMetadata.getComponentMismatches().at(allComponentIterations.size())));
        mismatchIterations *= allComponentIterations.back();
    }

    result.reserve(result.size() + mismatchIterations);
    for (unsigned iteration = 0; mismatchIterations > iteration; ++iteration)
    {
        const std::pair<Kmer, unsigned> mismatchKmer = generateMismatchKmer(
            kmer, componentLengths, barcodeMetadata.getComponentMismatches(), allComponentIterations, iteration, result);
        result.push_back(Barcode(mismatchKmer.first, BarcodeId(0, barcodeMetadata.getIndex(), 0, mismatchKmer.second)));
    }
}

std::vector<Barcode> BarcodeResolver::generateMismatches(
        const flowcell::BarcodeMetadataList &allBarcodeMetadata,
        const flowcell::BarcodeMetadataList &barcodeGroup)
{
    ISAAC_ASSERT_MSG(!barcodeGroup.empty(), "Barcode list must be not empty");
    ISAAC_ASSERT_MSG(barcodeGroup.at(0).isDefault(), "The very first barcode must be the 'unknown indexes or no index' one");
    std::vector<Barcode> ret;
    ret.reserve(barcodeGroup.size()-1);
    // Don't generate anything for the 'unknown indexes' barcode.
    std::for_each(barcodeGroup.begin() + 1,
                  barcodeGroup.end(),
                  boost::bind(&BarcodeResolver::generateBarcodeMismatches, _1, boost::ref(ret)));
    ISAAC_THREAD_CERR << "Generated " << ret.size() << " mismatch barcodes " << std::endl;
    std::sort(ret.begin(), ret.end(), orderBySequenceAndBarcode);

    // Will throw common::InvalidOptionException if there are colliding sequences produced from different
    // barcodes
    ret.erase(std::unique(ret.begin(), ret.end(),
                          boost::bind(&isSequenceBarcodeEqual, boost::ref(allBarcodeMetadata), _1, _2)), ret.end());

    return ret;
}

BarcodeResolver::BarcodeResolver(
    const flowcell::TileMetadataList &allTilesMetadata,
    const flowcell::BarcodeMetadataList &allBarcodeMetadata,
    const flowcell::BarcodeMetadataList &barcodeGroup)
    : allTilesMetadata_(allTilesMetadata)
    , allBarcodeMetadata_(allBarcodeMetadata)
    , mismatchBarcodes_(generateMismatches(allBarcodeMetadata_, barcodeGroup))
    , unknownBarcodeIndex_(barcodeGroup.at(0).getIndex())
    , barcodeHits_(allBarcodeMetadata_.size())
{
}

inline std::ostream &operator << (std::ostream &os, const std::vector<unsigned> &mismatchesPerComponent)
{
    os << mismatchesPerComponent.at(0);

    BOOST_FOREACH(const unsigned &thisCompMismatches,
                  std::make_pair(mismatchesPerComponent.begin() + 1, mismatchesPerComponent.end()))
    {
        os << ":" << thisCompMismatches;
    }
    return os;
}

/**
 * \brief Updates barcode indexes with those of teh matching mismatch barcodes.
 *        Index 0 is reserved for the undetermined barcode.
 */
void BarcodeResolver::resolve(
    std::vector<Barcode> &dataBarcodes,
    demultiplexing::DemultiplexingStats &demultiplexingStats)
{
    ISAAC_THREAD_CERR << "Resolving barcodes for " << dataBarcodes.size() << " clusters against " <<
        mismatchBarcodes_.size() << " mismatch variants" << std::endl;

    unsigned long totalBarcodeHits = 0;
    std::sort(dataBarcodes.begin(), dataBarcodes.end(), orderBySequence);
    std::vector<Barcode>::const_iterator mismatchBarcodeIterator = mismatchBarcodes_.begin();

    for(std::vector<Barcode>::iterator dataBarcodeIterator = dataBarcodes.begin();
        dataBarcodes.end() != dataBarcodeIterator; ++dataBarcodeIterator)
    {
        Barcode &dataBarcode = *dataBarcodeIterator;
        // rewind sample sheet to the first data match
        while(mismatchBarcodes_.end() != mismatchBarcodeIterator &&
            dataBarcode.getSequence() > mismatchBarcodeIterator->getSequence())
        {
            ISAAC_ASSERT_MSG(dataBarcode.getBarcode() == unknownBarcodeIndex_, "Data barcodes are expected to have the index preset to 'unknown'");
            ++mismatchBarcodeIterator;
        }
/*
        if (mismatchBarcodes_.end() == mismatchBarcodeIterator)
        {
            // there will be no matches anymore. Bail out.
            break;
        }
*/
        if (dataBarcode.getSequence() == mismatchBarcodeIterator->getSequence())
        {
            // match!, set the index in data
            BarcodeId barcodeId(dataBarcode.getTile(), mismatchBarcodeIterator->getBarcode(),
                                dataBarcode.getCluster(), mismatchBarcodeIterator->getMismatches());
            dataBarcode.setBarcodeId(barcodeId);
            ++barcodeHits_.at(mismatchBarcodeIterator->getBarcode());
            ++totalBarcodeHits;
            demultiplexingStats.recordBarcode(barcodeId);
        }
        else
        {
            std::vector<Barcode>::iterator sameBarcodeIterator = dataBarcodeIterator;
            while(dataBarcodes.end() != sameBarcodeIterator &&
                sameBarcodeIterator->getSequence() == dataBarcodeIterator->getSequence())
            {
                demultiplexingStats.recordUnknownBarcode(unknownBarcodeIndex_, dataBarcode.getTile());
                ++sameBarcodeIterator;
            }

            demultiplexingStats.recordUnknownBarcodeHits(dataBarcodeIterator->getSequence(),
                                                         std::distance(dataBarcodeIterator, sameBarcodeIterator));
            dataBarcodeIterator = sameBarcodeIterator - 1;
        }

        // else dataBarcode.getSequence() < mismatchBarcodeIterator->getSequence(), no match, move on
    }

    if (!dataBarcodes.empty())
    {
        demultiplexingStats.finalizeUnknownBarcodeHits(unknownBarcodeIndex_);
    }
    ISAAC_THREAD_CERR << "Resolving barcodes done for " << dataBarcodes.size() << " clusters against " <<
        mismatchBarcodes_.size() << " mismatch variants. Found barcode hits breakdown. Total(" << totalBarcodeHits << "):"<< std::endl;

    BOOST_FOREACH(const unsigned long &barcodeHits, barcodeHits_)
    {
        if (barcodeHits)
        {
            ISAAC_THREAD_CERR << allBarcodeMetadata_.at(&barcodeHits - &barcodeHits_.front()) << ": " << barcodeHits << std::endl;
        }
    }
}

} // namespace demultiplexing
} // namespace isaac

