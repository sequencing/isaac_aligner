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
#include "reference/ContigLoader.hh"
#include "reference/ReferenceKmer.hh"
#include "workflow/ExtractNeighborsWorkflow.hh"

namespace isaac
{
namespace workflow
{

ExtractNeighborsWorkflow::ExtractNeighborsWorkflow(
    const bfs::path &sortedReferenceXml,
    const bfs::path &outputFilePath
    )
    : sortedReferenceXml_(sortedReferenceXml),
      outputFilePath_(outputFilePath),
      threads_(boost::thread::hardware_concurrency()),
      xml_(reference::loadSortedReferenceXml(sortedReferenceXml_))
{
}

/**
 * \brief Builds a vector of global starts of contigs (all contig bases are considered)
 *        for the contigs in a given order
 */
std::vector<unsigned long> computeContigOffsets(const reference::SortedReferenceXml::Contigs &contigs)
{
    std::vector<unsigned long> ret(contigs.size());

    unsigned long lastOffset = 0UL;

    BOOST_FOREACH(const reference::SortedReferenceXml::Contig &contig, contigs)
    {
        ret.at(contig.index_) = lastOffset;
        lastOffset += contig.totalBases_;
    }

    return ret;
}

inline bool orderByKaryotypeIndex(
    const reference::SortedReferenceXml::Contig& left,
    const reference::SortedReferenceXml::Contig& right)
{
    return left.karyotypeIndex_ < right.karyotypeIndex_;
}

void ExtractNeighborsWorkflow::run()
{
    const std::vector<reference::SortedReferenceXml::MaskFile> maskFiles = xml_.getMaskFileList("ABCD");
    if (maskFiles.empty())
    {
        BOOST_THROW_EXCEPTION(isaac::common::PreConditionException("No mask files in " + sortedReferenceXml_.string()));
    }

    reference::SortedReferenceXml::Contigs contigs = xml_.getContigs();
    std::sort(contigs.begin(), contigs.end(), &orderByKaryotypeIndex);
    const std::vector<unsigned long> contigOffsets = computeContigOffsets(contigs);
    std::vector<bool> neighbors(reference::genomeLength(contigs), false);

    // there could be mutliple mask widths in the xml. Just pick one.
    const unsigned maskWidth = maskFiles.front().maskWidth;

    BOOST_FOREACH(const reference::SortedReferenceXml::MaskFile &maskFile, maskFiles)
    {
        //Don't reprocess redundant mask files of different widths
        if (maskWidth == maskFile.maskWidth)
        {
            scanMaskFile(maskFile, contigOffsets, neighbors);
        }
    }

    dumpResults(neighbors);
}

void ExtractNeighborsWorkflow::scanMaskFile(
    const reference::SortedReferenceXml::MaskFile &maskFile,
    const std::vector<unsigned long> &contigOffsets,
    std::vector<bool> &neighbors)
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

    unsigned maskNeighbors = 0;
    while(maskInput)
    {
        reference::ReferenceKmer referenceKmer;
        if (maskInput.read(reinterpret_cast<char *>(&referenceKmer), sizeof(referenceKmer)))
        {
            const reference::ReferencePosition pos = referenceKmer.getReferencePosition();
            if (!pos.isTooManyMatch() && pos.hasNeighbors())
            {
                neighbors.at(contigOffsets.at(pos.getContigId()) + pos.getPosition()) = true;
                ++maskNeighbors;
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
        ISAAC_THREAD_CERR << "Scanning " << maskFile.path << " found " << maskNeighbors << " neighbors " << std::endl;
    }
}

void ExtractNeighborsWorkflow::dumpResults(const std::vector<bool> &neighbors)
{
    std::ofstream outputFile(outputFilePath_.c_str());
    if (!outputFile)
    {
        const boost::format message = boost::format("Failed to open file %s for writing: %s") % outputFilePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    unsigned neighborCount = 0;
    char byte = 0;
    for (size_t i = 0; i < neighbors.size(); ++ i)
    {
        const unsigned shift = i % 8;
        const char positionHasNeighbors = neighbors[i];
        byte |= positionHasNeighbors << shift;
        if (7 == shift)
        {
            if (!outputFile.write(&byte, sizeof(byte)))
            {
                const boost::format message = boost::format("Failed to write byte into %s: %s") % outputFilePath_ % strerror(errno);
                BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
            }
            byte = 0;
        }
        neighborCount += positionHasNeighbors;
    }

    if (!outputFile.write(&byte, sizeof(byte)))
    {
        const boost::format message = boost::format("Failed to write final byte into %s: %s") % outputFilePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    ISAAC_THREAD_CERR << "Stored " << neighborCount << " neighbor locations in " << outputFilePath_ << std::endl;

}

} // namespace workflow
} // namespace isaac
