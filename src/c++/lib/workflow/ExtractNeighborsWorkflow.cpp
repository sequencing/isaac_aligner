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
 ** \file ExtractNeighborsWorkflow.cpp
 **
 ** \brief see ExtractNeighborsWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/BitsetSaver.hh"
#include "reference/ContigLoader.hh"
#include "reference/ReferenceKmer.hh"
#include "workflow/ExtractNeighborsWorkflow.hh"

namespace isaac
{
namespace workflow
{

ExtractNeighborsWorkflow::ExtractNeighborsWorkflow(
    const bfs::path &sortedReferenceMetadata,
    const bfs::path &neighborsFilePath,
    const bfs::path &highRepeatsFilePath
    )
    : sortedReferenceMetadata_(sortedReferenceMetadata),
      neighborsFilePath_(neighborsFilePath),
      highRepeatsFilePath_(highRepeatsFilePath),
      threads_(boost::thread::hardware_concurrency()),
      xml_(reference::loadSortedReferenceXml(sortedReferenceMetadata_))
{
}

template <typename KmerT>
void ExtractNeighborsWorkflow::run()
{
    const reference::SortedReferenceMetadata::MaskFiles &maskFiles =
        xml_.getMaskFileList(oligo::KmerTraits<KmerT>::KMER_BASES);
    if (maskFiles.empty())
    {
        BOOST_THROW_EXCEPTION(isaac::common::PreConditionException("No mask files in " + sortedReferenceMetadata_.string()));
    }

    const reference::SortedReferenceMetadata::Contigs contigs = xml_.getKaryotypeOrderedContigs();
    const std::vector<unsigned long> contigOffsets = reference::computeContigOffsets(contigs);
    std::vector<bool> neighbors(reference::genomeLength(contigs), false);
    std::vector<bool> highRepeats(highRepeatsFilePath_.empty() ? 0 : reference::genomeLength(contigs), true);

    // there could be mutliple mask widths in the xml. Just pick one.
    const unsigned maskWidth = xml_.getDefaultMaskWidth();

    BOOST_FOREACH(const reference::SortedReferenceMetadata::MaskFile &maskFile, maskFiles)
    {
        //Don't reprocess redundant mask files of different widths
        if (maskWidth == maskFile.maskWidth)
        {
            scanMaskFile<KmerT>(maskFile, contigOffsets, neighbors, highRepeats);
        }
    }

    dumpResults(neighbors, highRepeats);
}

template <typename KmerT>
void ExtractNeighborsWorkflow::scanMaskFile(
    const reference::SortedReferenceMetadata::MaskFile &maskFile,
    const std::vector<unsigned long> &contigOffsets,
    std::vector<bool> &neighbors,
    std::vector<bool> &highRepeats)
{
    if (!exists(maskFile.path))
    {
        const boost::format message = boost::format("Mask file %s does not exist: %s") % maskFile.path;
        BOOST_THROW_EXCEPTION(common::IoException(ENOENT, message.str()));
    }

    std::ifstream maskInput(maskFile.path.c_str());
    if (!maskInput)
    {
        const boost::format message = boost::format("Failed to open mask file %s for reading: %s") % maskFile.path % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    std::size_t scannedKmers = 0, maskNeighbors = 0, maskNonHighRepeats = 0;
    while(maskInput)
    {
        reference::ReferenceKmer<KmerT> referenceKmer;
        if (maskInput.read(reinterpret_cast<char *>(&referenceKmer), sizeof(referenceKmer)))
        {
            ++scannedKmers;
            const reference::ReferencePosition pos = referenceKmer.getReferencePosition();
            if (!pos.isTooManyMatch())
            {
                if (pos.hasNeighbors())
                {
                    neighbors.at(contigOffsets.at(pos.getContigId()) + pos.getPosition()) = true;
                    ++maskNeighbors;
                }

                if (!highRepeats.empty())
                {
                    highRepeats.at(contigOffsets.at(pos.getContigId()) + pos.getPosition()) = false;
                }
                ++maskNonHighRepeats;
            }
        }
    }
    if (!maskInput.eof())
    {
        const boost::format message = boost::format("Failed to scan %s to the end") % maskFile.path % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
    else
    {
        ISAAC_THREAD_CERR << "Scanning " << maskFile.path << " found " << scannedKmers << " kmers of which " <<
            maskNeighbors << " neighbors and " << (scannedKmers - maskNonHighRepeats) << " high repeats" << std::endl;
    }
}

void ExtractNeighborsWorkflow::dumpResults(const std::vector<bool> &neighbors, const std::vector<bool> &highRepeats)
{
    io::BitsetSaver neighborsSaver(neighborsFilePath_);
    neighborsSaver.save(neighbors);
    ISAAC_THREAD_CERR << "Stored " << neighbors.size() << " neighbor locations in " << neighborsFilePath_ << std::endl;
    if (!highRepeatsFilePath_.empty())
    {
        io::BitsetSaver highRepeatsSaver(highRepeatsFilePath_);
        highRepeatsSaver.save(highRepeats);
        ISAAC_THREAD_CERR << "Stored " << highRepeats.size() << " high repeats locations in " << highRepeatsFilePath_ << std::endl;
    }
}

template void ExtractNeighborsWorkflow::run<oligo::ShortKmerType>();
template void ExtractNeighborsWorkflow::run<oligo::KmerType>();
template void ExtractNeighborsWorkflow::run<oligo::LongKmerType>();

} // namespace workflow
} // namespace isaac
