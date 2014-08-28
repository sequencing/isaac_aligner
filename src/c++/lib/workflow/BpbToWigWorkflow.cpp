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
 ** \file BpbToWigWorkflow.cpp
 **
 ** \brief see BpbToWigWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/BitsetSaver.hh"
#include "workflow/BpbToWigWorkflow.hh"

namespace isaac
{
namespace workflow
{

BpbToWigWorkflow::BpbToWigWorkflow(
    const bfs::path &sortedReferenceMetadata,
    const bfs::path &inputFilePath,
    const std::string &outputFormatString
    )
    : sortedReferenceMetadata_(sortedReferenceMetadata),
      inputFilePath_(inputFilePath),
      outputFormatString_(outputFormatString),
      xml_(reference::loadSortedReferenceXml(sortedReferenceMetadata_))
{
}

inline bool orderByKaryotypeIndex(
    const reference::SortedReferenceMetadata::Contig& left,
    const reference::SortedReferenceMetadata::Contig& right)
{
    return left.karyotypeIndex_ < right.karyotypeIndex_;
}

void BpbToWigWorkflow::run()
{
    std::ifstream bitsetFile(inputFilePath_.c_str());
    if (!bitsetFile)
    {
        const boost::format message = boost::format("Failed to open bitset file %s for reading: %s") %
            inputFilePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }

    const reference::SortedReferenceMetadata::Contigs &contigs = xml_.getContigs();

    if ("wig" == outputFormatString_)
    {
        printWig(bitsetFile, contigs);
    }
    else if ("bed" == outputFormatString_)
    {
        printBed(bitsetFile, contigs);
    }
    else
    {
        ISAAC_ASSERT_MSG(false, "Unknown output format: " << outputFormatString_);
    }
}

void BpbToWigWorkflow::printWig(
    std::ifstream &bitsetFile,
    const reference::SortedReferenceMetadata::Contigs &contigs )
{
    char bits = 0;
    bitsetFile.read(&bits, 1);
    std::size_t bitPos = 1;
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, contigs)
    {
        std::cout << "variableStep chrom=" << contig.name_ << std::endl;
        std::size_t contigBases = contig.totalBases_;
        while(bitsetFile && contigBases)
        {
            const std::size_t contigPosition = contig.totalBases_ - contigBases + 1;
            const unsigned value = (bits & 0x01);
            if (value)
            {
                std::cout << contigPosition << "\t" << value << std::endl;
            }
            if (!(bitPos % 8))
            {
                bitsetFile.read(&bits, 1);
            }
            else
            {
                bits >>= 1;
            }
            --contigBases;
            ++bitPos;
        }
        if (!bitsetFile && !bitsetFile.eof())
        {
            const boost::format message = boost::format("Failed to read bitset file %s: %s") % inputFilePath_ % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
    }
}

void BpbToWigWorkflow::printBed(
    std::ifstream &bitsetFile,
    const reference::SortedReferenceMetadata::Contigs &contigs )
{
    std::cout << "track graphType=bar type=bedGraph" << std::endl;

    char bits = 0;
    bitsetFile.read(&bits, 1);
    unsigned lastValue = 0;
    std::size_t lastContigPosition = 0;
    std::size_t bitPos = 1;
    BOOST_FOREACH(const reference::SortedReferenceMetadata::Contig &contig, contigs)
    {
        std::size_t contigBases = contig.totalBases_;
        while(bitsetFile && contigBases)
        {
            const std::size_t contigPosition = contig.totalBases_ - contigBases;
            const unsigned value = (bits & 0x01);
            if (value != lastValue)
            {
                if (lastValue)
                {
                    std::cout << contig.name_ << "\t" << lastContigPosition << "\t" << contigPosition << "\t" << lastValue << std::endl;
                }
                lastContigPosition = contigPosition;
                lastValue = value;
            }
            if (!(bitPos % 8))
            {
                bitsetFile.read(&bits, 1);
            }
            else
            {
                bits >>= 1;
            }
            --contigBases;
            ++bitPos;
        }
        if (!bitsetFile && !bitsetFile.eof())
        {
            const boost::format message = boost::format("Failed to read bitset file %s: %s") % inputFilePath_ % strerror(errno);
            BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
        }
        const std::size_t contigPosition = contig.totalBases_ - contigBases;
        if (lastValue)
        {
            std::cout << contig.name_ << "\t" << lastContigPosition << "\t" << contigPosition << "\t" << lastValue << std::endl;
        }
        lastContigPosition = 0;
    }
}


} // namespace workflow
} // namespace isaac
