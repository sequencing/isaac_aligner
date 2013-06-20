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
 ** \file ReoderReferenceWorkflow.cpp
 **
 ** \brief see ReoderReferenceWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/ContigLoader.hh"
#include "workflow/ReorderReferenceWorkflow.hh"

namespace isaac
{
namespace workflow
{

ReorderReferenceWorkflow::ReorderReferenceWorkflow(
    const bfs::path &sortedReferenceMetadata,
    const bfs::path &newXmlPath,
    const bfs::path &newFaPath,
    const std::vector<std::string> &newOrder,
    const unsigned basesPerLine
    )
    : sortedReferenceMetadata_(sortedReferenceMetadata),
      newXmlPath_(newXmlPath),
      newFaPath_(newFaPath),
      newOrder_(newOrder),
      basesPerLine_(basesPerLine),
      threads_(boost::thread::hardware_concurrency()),
      xml_(reference::loadSortedReferenceXml(sortedReferenceMetadata_))
{
    if (!newOrder_.empty())
    {
        std::vector<bool> present(newOrder_.size());
        BOOST_FOREACH(reference::SortedReferenceMetadata::Contig &xmlContig, xml_.getContigs())
        {
            const std::vector<std::string>::const_iterator newOrderIt = std::find(newOrder_.begin(), newOrder_.end(), xmlContig.name_);
            if (newOrder_.end() == newOrderIt)
            {
                BOOST_THROW_EXCEPTION(isaac::common::InvalidParameterException(
                    "Contig name not listed in the new order: " + xmlContig.name_));
            }
            const unsigned newKaryotypeIndex = newOrderIt - newOrder_.begin();
            present.at(newKaryotypeIndex) = true;
            xmlContig.karyotypeIndex_ = newKaryotypeIndex;
        }

        const std::vector<bool>::const_iterator firstUnused = std::find(present.begin(), present.end(), false);
        if(present.end() != firstUnused)
        {
            BOOST_THROW_EXCEPTION(isaac::common::InvalidParameterException(
                "Contig name listed in the new order not found in the reference: " + newOrder_.at(firstUnused - present.begin())));
        }
    }
    else
    {
        ISAAC_THREAD_CERR << "Preserving the existing order of contigs" << std::endl;
    }
}

bool ReorderReferenceWorkflow::orderByKaryotypeIndex(const reference::Contig& left, const reference::Contig& right)
{
    return xml_.getContigs().at(left.index_).karyotypeIndex_ < xml_.getContigs().at(right.index_).karyotypeIndex_;
}

void ReorderReferenceWorkflow::run()
{
    std::ofstream xmlOs(newXmlPath_.c_str());
    if (!xmlOs)
    {
        BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open output file: " + newXmlPath_.string()));
    }

    if (!newOrder_.empty() || !xml_.singleFileReference())
    {
        std::ofstream fastaOs(newFaPath_.c_str());
        if (!fastaOs)
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open output file: " + newFaPath_.string()));
        }

        std::vector<reference::Contig> contigs = reference::loadContigs(xml_.getContigs(), threads_);
        std::sort(contigs.begin(), contigs.end(),
                  boost::bind(&ReorderReferenceWorkflow::orderByKaryotypeIndex, this, _1, _2));

        BOOST_FOREACH(const reference::Contig &contig, contigs)
        {
            storeContig(fastaOs, contig);
        }
    }
    else
    {
        const boost::filesystem::path &srcPath = xml_.getContigs().front().filePath_;
        boost::filesystem::copy_file(srcPath, newFaPath_);
        ISAAC_THREAD_CERR << "Copied entire reference from " << srcPath << " to " << newFaPath_ << std::endl;
        BOOST_FOREACH(reference::SortedReferenceMetadata::Contig &xmlContig, xml_.getContigs())
        {
            xmlContig.filePath_ = newFaPath_;
        }
    }

    saveSortedReferenceXml(xmlOs, xml_);
}

void ReorderReferenceWorkflow::writeBase(std::ostream &os, const char base, const bool writeNewline)
{
    if (!(os << base) || (writeNewline && !(os << "\n")))
    {
        BOOST_THROW_EXCEPTION(
            isaac::common::IoException(errno, (boost::format("Failed to write data into output file %s: ") %
                newFaPath_.string()).str()));
    }
}

void ReorderReferenceWorkflow::storeContig(std::ostream &os,const reference::Contig &contig)
{
    reference::SortedReferenceMetadata::Contig &xmlContig = xml_.getContigs().at(contig.index_);

    if (!(os << ">" << xmlContig.name_ << std::endl))
    {
        BOOST_THROW_EXCEPTION(
            isaac::common::IoException(errno, (boost::format("Failed to write contig name %s into output file: ") %
                xmlContig.name_ % newFaPath_.string()).str()));
    }

    unsigned long startPos = os.tellp();
    std::vector<char>::const_iterator current = contig.forward_.begin();
    while(contig.forward_.end() != current)
    {
        writeBase(os, *current,
                  (!((current - contig.forward_.begin() + 1) % basesPerLine_)) ||
                  (contig.forward_.end() - 1 == current));
        ++current;
    }
    unsigned long endPos = os.tellp();

    xmlContig.filePath_ = newFaPath_;
    xmlContig.offset_ = startPos;
    xmlContig.size_ = endPos - startPos;

    ISAAC_THREAD_CERR << "Stored contig: " << xmlContig.name_ << std::endl;
}

} // namespace workflow
} // namespace isaac
