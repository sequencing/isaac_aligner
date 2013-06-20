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
 ** \file SortedReferenceXml.cpp
 **
 ** SortedReference.xml helper.
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "config.h"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/SortedReferenceXml.hh"
#include "xml/XmlReader.hh"
#include "xml/XmlWriter.hh"

namespace isaac
{
namespace reference
{

void serialize(xml::XmlReader &reader, SortedReferenceMetadata::MaskFile &mf, const unsigned int version)
{
    mf.mask_ = reader("Mask")["Mask"];
    mf.path = reader.nextChildElement("File").readElementText().string();
    mf.kmers = (reader += "Kmers").nextChildElement("Total").readElementText();

    // deal with deprecated MaxPrefixRangeCount element
    if (reader.nextElementBelowLevel(4))
    {
        reader.assertName("MaxPrefixRangeCount");
    }
    else
    {
        reader.clear();
    }
}

void serialize(xml::XmlReader &reader, SortedReferenceMetadata::MaskFiles &maskFiles, const unsigned int version)
{
    while (reader.nextElementBelowLevel(3) && reader("Mask"))
    {
        SortedReferenceMetadata::MaskFile maskFile;
        serialize(reader, maskFile, version);
        maskFiles.push_back(maskFile);
    }
    reader.clear();
}

void serialize(xml::XmlReader &reader, SortedReferenceMetadata::AllMaskFiles &maskFiles, const unsigned int version)
{
    while (reader.nextElementBelowLevel(2) && reader("Masks"))
    {
        const unsigned maskWidth = reader["Width"];
        // SeedLength is optional for older xml files. Absent SeedLength is treated as 32
        const unsigned seedLength = reader.getAttribute("SeedLength", 32);
        if (!maskFiles[seedLength].empty())
        {
            BOOST_THROW_EXCEPTION(xml::XmlReaderException(std::string("Multiple Masks elements with same SeedLength are not allowed ") + reader.getCurrentDebugContext()));
        }
        serialize(reader, maskFiles[seedLength], version);
        BOOST_FOREACH(SortedReferenceMetadata::MaskFile &maskFile, maskFiles[seedLength])
        {
            maskFile.maskWidth = maskWidth;
        }
    }

    reader.clear();
}

void serialize(xml::XmlReader &reader, SortedReferenceMetadata::Contig &c, const unsigned int version)
{
    c.genomicPosition_ = reader("Contig")["Position"];
    c.index_ = reader("Contig").nextChildElement("Index").readElementText();
    if (reader.nextElementBelowLevel(2).checkName("KaryotypeIndex"))
    {
        c.karyotypeIndex_ = reader.readElementText();
        reader += "Name";
    }
    else
    {
        c.karyotypeIndex_ = c.index_;
    }
    c.name_ = reader("Name").readElementText().string();
    c.filePath_ = (reader += "Sequence").nextChildElement("File").readElementText().string();
    c.offset_ = (reader += "Offset").readElementText();
    c.size_ = (reader += "Size").readElementText();
    c.totalBases_ = (reader += "TotalBases").readElementText();
    c.acgtBases_ = (reader += "AcgtBases").readElementText();

    // loop through the elements until the parent closes
    while (reader.nextElementBelowLevel(2))
    {
        // process optional elements
        if (reader.checkName("As"))
        {
            c.bamSqAs_ = reader.readElementText().string();
        }
        else if (reader.checkName("Ur"))
        {
            c.bamSqUr_ = reader.readElementText().string();
        }
        else if (reader.checkName("M5"))
        {
            c.bamM5_ = reader.readElementText().string();
        }
    }
    reader.clear();
}

void serialize(xml::XmlReader &reader, SortedReferenceMetadata::Contigs &contigs, const unsigned int version)
{
    while (reader.nextElementBelowLevel(1) && reader("Contig"))
    {
        SortedReferenceMetadata::Contig contig;
        serialize(reader, contig, version);
        contigs.push_back(contig);
    }
    reader.clear();
}

template <>
void serialize<xml::XmlReader>(xml::XmlReader &reader, SortedReferenceMetadata &sortedReferenceMetadata, const unsigned int version)
{
    ISAAC_ASSERT_MSG(version == SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION, "Unexpected version requested: " << version);

    sortedReferenceMetadata.formatVersion_ = (reader+="SortedReference").nextChildElement("FormatVersion").readElementText();
    if (SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION < sortedReferenceMetadata.formatVersion_ ||
        SortedReferenceMetadata::OLDEST_SUPPORTED_REFERENCE_FORMAT_VERSION > sortedReferenceMetadata.formatVersion_)
    {
        BOOST_THROW_EXCEPTION(xml::XmlReaderException(
            (boost::format("Unexpected sorted reference FormatVersion: %s. FormatVersion must be in range [%d,%d]") %
            reader.getValue().string() % SortedReferenceMetadata::OLDEST_SUPPORTED_REFERENCE_FORMAT_VERSION %
            SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION).str()));
    }

    // SoftwareVersion is optional for older xml files
    if (reader++.checkName("SoftwareVersion"))
    {
        reader++;
    }

    // Contigs may not be present
    if (reader.checkName("Contigs"))
    {
        serialize(reader, sortedReferenceMetadata.contigs_, version);
        // advance if possible
        ++reader;
    }

    // Permutations may not be present
    if (reader && reader.checkName("Permutations"))
    {
        // only ABCD permutation is supported
        reader += "Permutation";
        if (reader["Name"] != "ABCD")
        {
            BOOST_THROW_EXCEPTION(xml::XmlReaderException(std::string("Only ABCD permutation masks are supported")));
        }

        serialize(reader, sortedReferenceMetadata.maskFiles_, version);
    }

    if (!sortedReferenceMetadata.maskFiles_.empty())
    {
        sortedReferenceMetadata.defaultMaskWidth_ = sortedReferenceMetadata.maskFiles_.begin()->second.at(0).maskWidth;
    }
    else
    {
        sortedReferenceMetadata.defaultMaskWidth_ = 0;
    }

    // As we were able to successfully read the file, bump format version up to the current to avoid confusion
    // when stored or merged
    sortedReferenceMetadata.formatVersion_ = SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION;
}


SortedReferenceMetadata loadSortedReferenceXml(
    std::istream &is)
{
    xml::XmlReader reader(is);
    SortedReferenceMetadata ret;
    serialize(reader, ret, SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION);
    return ret;
}


SortedReferenceMetadata loadSortedReferenceXml(
    const boost::filesystem::path &xmlPath)
{
    SortedReferenceMetadata ret;
    std::ifstream is(xmlPath.c_str());
    if (!is)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open sorted reference file " + xmlPath.string()));
    }

    return loadSortedReferenceXml(is);
}

void serialize(xml::XmlWriter &writer, const SortedReferenceMetadata::MaskFile &mf, const unsigned int version)
{
    ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Mask")
    {
        writer.writeAttribute("Mask", mf.mask_);
        writer.writeElement("File", mf.path.string());
        ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Kmers")
        {
            writer.writeElement("Total", mf.kmers);
        }
    }
}

void serialize(xml::XmlWriter &writer, const SortedReferenceMetadata::Contig &contig, const unsigned int version)
{
    ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Contig")
    {
        writer.writeAttribute("Position", contig.genomicPosition_);

        writer.writeElement("Index", contig.index_);
        writer.writeElement("KaryotypeIndex", contig.karyotypeIndex_);
        writer.writeElement("Name", contig.name_);

        ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Sequence")
        {
            writer.writeElement("File", contig.filePath_.string());
            writer.writeElement("Offset", contig.offset_);
            writer.writeElement("Size", contig.size_);
        }

        writer.writeElement("TotalBases", contig.totalBases_);
        writer.writeElement("AcgtBases", contig.acgtBases_);

        ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "BamMetadata")
        {
            ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Sq")
            {
                writer.writeElement("As", contig.bamSqAs_);
                writer.writeElement("Ur", contig.bamSqUr_);
                writer.writeElement("M5", contig.bamM5_);
            }
        }
    }
}

template <typename T>
void serialize(xml::XmlWriter &writer, const std::vector<T> &v, const unsigned int version)
{
    BOOST_FOREACH(const T &contig, v)
    {
        serialize(writer, contig, version);
    }
}

template <>
void serialize<xml::XmlWriter>(xml::XmlWriter &writer, const SortedReferenceMetadata &sortedReferenceMetadata, const unsigned int version)
{
    ISAAC_ASSERT_MSG(version == SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION, "Unexpected version requested: " << version);

    ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "SortedReference")
    {
        writer.writeElement("FormatVersion", sortedReferenceMetadata.formatVersion_);
        writer.writeElement("SoftwareVersion", iSAAC_VERSION_FULL);

        if (!sortedReferenceMetadata.contigs_.empty())
        {
            ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Contigs")
            {
                serialize(writer, sortedReferenceMetadata.contigs_, version);
            }
        }

        if (!sortedReferenceMetadata.maskFiles_.empty())
        {
            ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Permutations")
            {
                ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Permutation")
                {
                    writer.writeAttribute("Name", "ABCD");
                    BOOST_FOREACH(const SortedReferenceMetadata::AllMaskFiles::value_type &seedMaskFiles,
                                  sortedReferenceMetadata.maskFiles_)
                    {
                        ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, "Masks")
                        {
                            writer.writeAttribute("Width", seedMaskFiles.second.front().maskWidth);
                            writer.writeAttribute("SeedLength", seedMaskFiles.first);
                            serialize(writer, seedMaskFiles.second, version);
                        }
                    }
                }
            }
        }
    }
    writer.close();
}

void saveSortedReferenceXml(
    const boost::filesystem::path &xmlPath,
    const SortedReferenceMetadata &sortedReferenceMetadata)
{
    std::ofstream ofs(xmlPath.c_str());
    if (!ofs)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open sorted reference file for write: " + xmlPath.string()));
    }
    saveSortedReferenceXml(ofs, sortedReferenceMetadata);
}

void saveSortedReferenceXml(
    std::ostream &os,
    const SortedReferenceMetadata &sortedReferenceMetadata)
{
    xml::XmlWriter writer(os);
    serialize(writer, sortedReferenceMetadata, SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION);
}

} // namespace reference
} // namespace isaac

