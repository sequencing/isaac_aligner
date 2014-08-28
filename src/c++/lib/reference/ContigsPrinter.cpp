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
 ** \file ContigsPrinter.cpp
 **
 ** Top level component for extracting contigs metadata.
 **
 ** \author Roman Petrovski
 **/

#include <boost/assert.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/algorithm/string.hpp>

#include "common/Exceptions.hh"
#include "common/MD5Sum.hh"
#include "io/FastaReader.hh"
#include "oligo/Nucleotides.hh"
#include "reference/ContigsPrinter.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

ContigsPrinter::ContigsPrinter (
    const boost::filesystem::path &originalSortedReferenceXml,
    const boost::filesystem::path &genomeFile
)
    : originalSortedReferenceXml_(originalSortedReferenceXml), genomeFile_(genomeFile)
{
}



void ContigsPrinter::run()
{
    SortedReferenceMetadata::Contigs originalContigs;
    if (!originalSortedReferenceXml_.empty())
    {
        SortedReferenceMetadata inXml(reference::loadSortedReferenceXml(originalSortedReferenceXml_));
        originalContigs = inXml.getContigs();
    }

    SortedReferenceMetadata::Contigs outputContigs;

    std::ifstream is(genomeFile_.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open reference file " + genomeFile_.string()));
    }

    long lastAcgtCount(0), lastBasesCount(0), lastContigByteStart(0), lineNumber(0), lastContigGenomicStart(0), streamPos(0);
    unsigned index = 0;
    std::string lastHeader, line;

    SortedReferenceMetadata outXml;
    isaac::common::MD5Sum md5Sum;
    // TODO: use the MultiFastaReader component to ensure consistency (particularly for the index)
    while (std::getline(is, line))
    {
        streamPos = is.tellg();
        ++lineNumber;
        if ('>' == line[0])
        {
            std::cerr << "header at:" << streamPos << "\n";
            if (!lastHeader.empty())
            {
                const std::string contigName = lastHeader.substr(1, lastHeader.find_first_of(" \t\r") - 1);
                SortedReferenceMetadata::Contigs::const_iterator originalContigIt =
                    std::find_if(originalContigs.begin(), originalContigs.end(), boost::bind(&SortedReferenceMetadata::Contig::name_, _1) == contigName);
                outXml.putContig(lastContigGenomicStart,
                                 contigName,
                                 genomeFile_,
                                 lastContigByteStart,
                                 streamPos - lastContigByteStart - (1 + line.length()), //byte length (assuming newline is one byte long)
                                 lastBasesCount, lastAcgtCount, index, index,
                                 (originalContigs.end() == originalContigIt ? "" : originalContigIt->bamSqAs_),
                                 (originalContigs.end() == originalContigIt ? "" : originalContigIt->bamSqUr_),
                                 isaac::common::MD5Sum::toHexString( md5Sum.getDigest().data, 16 ));

                ++index;
            }
            // Reset counters
            lastAcgtCount = 0;
            lastHeader = line;
            lastContigByteStart = streamPos;
            lastContigGenomicStart += lastBasesCount;
            lastBasesCount = 0;
            md5Sum.clear();
        }
        else
        {
            lastAcgtCount += std::count_if(line.begin(), line.end(), boost::bind(&oligo::getValue, _1) != oligo::invalidOligo);
            // ignore untranslated '\r'
            lastBasesCount += std::count_if(line.begin(), line.end(), boost::bind(&boost::ref<char>, _1) != '\r');
        
            // MD5 with format characters stripped out
            std::string lineStd = line;
            boost::to_upper(lineStd);
            lineStd.erase(std::remove_if(lineStd.begin(), lineStd.end(), &isspace), lineStd.end());  // isspace should remove ' ', '\n', '\r', '\t'
    
            md5Sum.update( lineStd.c_str(), lineStd.size() );
        }
    }
    if (!is.eof()) {
        BOOST_THROW_EXCEPTION(
                isaac::common::IoException(errno, "Failed while reading sequence."
                                             " Line:" + boost::lexical_cast<std::string>(lineNumber)
                                           + " Position:" + boost::lexical_cast<std::string>(is.tellg())
                                           + " Buffer:" + line + "\n"));
    }
    std::cerr << "end of stream at:" << streamPos << "\n";
    if (!lastHeader.empty())
    {
        const std::string contigName = lastHeader.substr(1, lastHeader.find_first_of(" \t\r") - 1);
        SortedReferenceMetadata::Contigs::const_iterator originalContigIt =
            std::find_if(originalContigs.begin(), originalContigs.end(), boost::bind(&SortedReferenceMetadata::Contig::name_, _1) == contigName);
        outXml.putContig(lastContigGenomicStart,
                         contigName,
                         genomeFile_,
                         lastContigByteStart,
                         streamPos - lastContigByteStart, //byte length
                         lastBasesCount, lastAcgtCount, index, index,
                         (originalContigs.end() == originalContigIt ? "" : originalContigIt->bamSqAs_),
                         (originalContigs.end() == originalContigIt ? "" : originalContigIt->bamSqUr_),
                         isaac::common::MD5Sum::toHexString( md5Sum.getDigest().data, 16 ));
        ++index;
    }

    saveSortedReferenceXml(std::cout, outXml);
}

} // namespace reference
} // namespace isaac
