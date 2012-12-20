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
 ** \file SortedReferenceXml.cpp
 **
 ** SortedReference.xml helper.
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/PtreeXml.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

static std::string getMasksKey(const std::string &permutationName)
{
    const std::string ret("SortedReference.Permutations.<indexed>Permutation.<Name>"
            + permutationName
            + ".<indexed>Masks");
    return ret;
}

static const std::string widthIndexAttributePrefix("<Width>");
static std::string getMaskKey(
        const std::string &permutationName,
        const unsigned int maskWidth)
{
    const std::string ret(
            getMasksKey(permutationName) + "." + widthIndexAttributePrefix + boost::lexical_cast<std::string>(maskWidth)
            + ".<indexed>Mask");
    return ret;
}

static const std::string maskIndexAttributePrefix("<Mask>");
static std::string getMaskValuePrefix(
        const std::string &permutationName,
        const unsigned int maskWidth,
        const isaac::oligo::Kmer mask)
{
    const std::string ret(
            getMaskKey(permutationName, maskWidth)
            + "." + maskIndexAttributePrefix + boost::lexical_cast<std::string>(mask)
            );
    return ret;
}

static const std::string contigIndexAttributePrefix("<Position>");

void SortedReferenceXml::putContig(
    const unsigned long genomicOffset,
    const std::string& name,
    const boost::filesystem::path &sequencePath,
    const unsigned long byteOffset,
    const unsigned long byteSize,
    const unsigned long totalBases,
    const unsigned long acgtBases,
    const unsigned index,
    const unsigned karyotypeIndex,
    const std::string &bamSqAs,
    const std::string &bamSqUr,
    const std::string &bamM5
    )
{

    const std::string contigKey("SortedReference.Contigs.<indexed>Contig." + contigIndexAttributePrefix
            + boost::lexical_cast<std::string>(genomicOffset));

    put(contigKey + ".Index", index);
    put(contigKey + ".KaryotypeIndex", karyotypeIndex);
    put(contigKey + ".Name", name);
    put(contigKey + ".Sequence.File", sequencePath.string());
    put(contigKey + ".Sequence.Offset", byteOffset);
    put(contigKey + ".Sequence.Size", byteSize);
    put(contigKey + ".TotalBases", totalBases);
    put(contigKey + ".AcgtBases", acgtBases);
    put(contigKey + ".BamMetadata.Sq.As", bamSqAs);
    put(contigKey + ".BamMetadata.Sq.Ur", bamSqUr);
    put(contigKey + ".BamMetadata.Sq.M5", bamM5);
}
void SortedReferenceXml::addMaskFile(const std::string &permutationName, const unsigned int maskWidth,
                              const isaac::oligo::Kmer mask, const boost::filesystem::path &filePath,
                                     const size_t kmers, const unsigned maxPrefixRangeCount)
{
    std::string maskValuePrefix(getMaskValuePrefix(permutationName, maskWidth, mask));

    put(maskValuePrefix + ".File", filePath.string());
    put(maskValuePrefix + ".Kmers.Total", kmers);
    put(maskValuePrefix + ".MaxPrefixRangeCount", maxPrefixRangeCount);
}

/**
 * \brief extracts the first Masks width for ABCD permutation.
 *
 * It is assumed that if other mask widths are represented in the file, they will be
 * placed after the default one
 *
 * It is assumed that ABCD permutation is always represented
 *
 */
unsigned int SortedReferenceXml::getDefaultMaskWidth() const
{
    static const std::string defaultMasksWidthKey(getMasksKey("ABCD"));
    const std::string widthIndex(get_child(defaultMasksWidthKey).begin()->first);
    return boost::lexical_cast<unsigned int>(widthIndex.substr(widthIndexAttributePrefix.length()));
}

const boost::filesystem::path SortedReferenceXml::getMaskFile(const std::string &permutationName, const unsigned int maskWidth,
             const isaac::oligo::Kmer mask) const
{
    std::string maskValuePrefix(getMaskValuePrefix(permutationName, maskWidth, mask));

    return get<std::string>(maskValuePrefix + ".File");
}

std::ostream &operator <<(std::ostream &os, const isaac::reference::SortedReferenceXml &tree)
{
    return isaac::io::serializeAsXml(os, tree);
}

std::istream &operator >> (std::istream &is, isaac::reference::SortedReferenceXml &indexedTree)
{
    boost::property_tree::read_xml<boost::property_tree::ptree>(is, indexedTree);

    const unsigned expectedFormatVersion = isaac::reference::SortedReferenceXml::currentReferenceFormatVersion_;
    const unsigned referenceFormatVersion = indexedTree.get<unsigned>("SortedReference.FormatVersion");

    if (expectedFormatVersion != referenceFormatVersion)
    {
        BOOST_THROW_EXCEPTION(common::UnsupportedVersionException(
              (boost::format("Sorted refernce format version is wrong. This sorted reference cannot be used. "
                  "Expected format version: %u. Actual: %u.") %
                  expectedFormatVersion % referenceFormatVersion).str()));
    }

    const std::vector<std::string > indexAttrs = boost::assign::list_of
        (std::string("SortedReference.Contigs.Contig.Position"))
        (std::string("SortedReference.Permutations.Permutation.Name"))
        (std::string("SortedReference.Permutations.Permutation.Masks.Width"))
        (std::string("SortedReference.Permutations.Permutation.Masks.Mask.Mask"))
            ;
    isaac::io::index(indexAttrs, indexedTree);
    return is;
}

SortedReferenceXml::Contigs SortedReferenceXml::getKaryotypeOrderedContigs() const
{
    Contigs ret = getContigs();
    std::sort(ret.begin(), ret.end(),
              boost::bind(&reference::SortedReferenceXml::Contig::karyotypeIndex_, _1) <
              boost::bind(&reference::SortedReferenceXml::Contig::karyotypeIndex_, _2));
    return ret;
}

SortedReferenceXml::Contigs SortedReferenceXml::getContigs() const
{
    using boost::property_tree::ptree;
    std::vector<bool> knownIndexFlags;
    // Check the integrity of the contigs and find out how many there are given that they are not necessarily
    // ordered by their Index
    BOOST_FOREACH(const ptree::value_type &contigNode, get_child("SortedReference.Contigs.<indexed>Contig"))
    {
        const unsigned contigIndex = contigNode.second.get<unsigned>("Index");

        if (contigIndex + 1 > knownIndexFlags.size())
        {
            knownIndexFlags.resize(contigIndex + 1, false);
        }
        if (knownIndexFlags[contigIndex])
        {
            using boost::format;
            using common::PreConditionException;
            const format message = (format("Index %d duplicated for contigs in sorted reference") % contigIndex);
            BOOST_THROW_EXCEPTION(PreConditionException(message.str()));
        }
        knownIndexFlags[contigIndex] = true;
    }
    if (knownIndexFlags.empty())
    {
        using common::PreConditionException;
        BOOST_THROW_EXCEPTION(PreConditionException("no contigs found in sorted reference"));
    }
    Contigs ret(knownIndexFlags.size());
    BOOST_FOREACH(const ptree::value_type &contigNode, get_child("SortedReference.Contigs.<indexed>Contig"))
    {
        ISAAC_ASSERT_MSG (contigIndexAttributePrefix == contigNode.first.substr(0, contigIndexAttributePrefix.length()),
                          "SortedReference.Contigs.<indexed>Contig is expected to contain only the <Position> nodes");

        const unsigned genomicPosition = boost::lexical_cast<unsigned>(contigNode.first.substr(contigIndexAttributePrefix.length()));

        const unsigned contigIndex = contigNode.second.get<unsigned>("Index");
        const unsigned karyotypeIndex = contigNode.second.get_optional<unsigned>("KaryotypeIndex").get_value_or(contigIndex);
        ISAAC_ASSERT_MSG(knownIndexFlags.at(karyotypeIndex), "Karyotype indexes must map to contig indexes");
        const std::string sequenceFile = contigNode.second.get<std::string>("Sequence.File");
        const Contig contig(contigIndex,
                            karyotypeIndex,
                            contigNode.second.get<std::string>("Name"),
                            sequenceFile,
                            contigNode.second.get<unsigned long>("Sequence.Offset"),
                            genomicPosition,
                            contigNode.second.get<unsigned long>("TotalBases"),
                            contigNode.second.get<unsigned long>("AcgtBases"),
                            contigNode.second.get<std::string>("BamMetadata.Sq.As", ""),
                            contigNode.second.get<std::string>("BamMetadata.Sq.Ur", ""),
                            contigNode.second.get<std::string>("BamMetadata.Sq.M5", ""));
        assert(ret.size() > contig.index_ && "Invalid contig index. Expected sequentially numbered starting with 0");
        ret[contig.index_] = contig;
    }
    return ret;
}


unsigned long SortedReferenceXml::getTotalKmers() const
{
    const std::string permutationName = "ABCD";
    const unsigned maskWidth = getDefaultMaskWidth();
    const unsigned maskCount = (1 << maskWidth);
    unsigned long ret = 0;
    for (unsigned mask = 0; maskCount > mask; ++mask)
    {
        const std::string maskValuePrefix = getMaskValuePrefix(permutationName, maskWidth, mask);
        ret += get<unsigned>(maskValuePrefix + ".Kmers.Total");
    }
    return ret;
}

unsigned SortedReferenceXml::getMaxPrefixRangeCount() const
{
    const std::string permutationName = "ABCD";
    const unsigned maskWidth = getDefaultMaskWidth();
    const unsigned maskCount = (1 << maskWidth);
    unsigned ret = 0;
    for (unsigned mask = 0; maskCount > mask; ++mask)
    {
        const std::string maskValuePrefix = getMaskValuePrefix(permutationName, maskWidth, mask);
        const unsigned maxPrefixRangeCount = get<unsigned>(maskValuePrefix + ".MaxPrefixRangeCount");
        if (maxPrefixRangeCount > ret)
        {
            ret = maxPrefixRangeCount;
        }
    }
    return ret;
}

std::vector<SortedReferenceXml::MaskFile> SortedReferenceXml::getMaskFileList(const std::string &permutationName) const
{
    const unsigned maskWidth = getDefaultMaskWidth();
    const std::string masksKey = getMaskKey(permutationName, maskWidth);

    std::vector<MaskFile> ret;
    BOOST_FOREACH(const boost::property_tree::ptree::value_type &mask, get_child(masksKey))
    {
        ISAAC_ASSERT_MSG (maskIndexAttributePrefix == mask.first.substr(0, maskIndexAttributePrefix.length()),
                          "Permutations.Permutation.Masks.<indexed>Mask is expected to contain only the <Mask> nodes");

        const unsigned maskValue = boost::lexical_cast<unsigned>(mask.first.substr(maskIndexAttributePrefix.length()));
        const boost::filesystem::path path =  mask.second.get<boost::filesystem::path>("File");
        ret.resize(std::max<unsigned>(maskValue + 1, ret.size()));
        ret[maskValue].path = path;
        ret[maskValue].mask = maskValue;
        ret[maskValue].maskWidth = maskWidth;
        ret[maskValue].kmers = mask.second.get<size_t>("Kmers.Total");
        ret[maskValue].maxPrefixRangeCount = mask.second.get<unsigned>("MaxPrefixRangeCount");
    }
    return ret;
}


reference::SortedReferenceXml loadSortedReferenceXml(
    const boost::filesystem::path &xmlPath)
{
    reference::SortedReferenceXml ret;
    std::ifstream is;
    is.open(xmlPath.c_str());
    if (!is)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open sorted reference desriptor file " + xmlPath.string()));
    }

    if (!(is >> ret))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to read data from sorted reference desriptor file " + xmlPath.string()));
    }

    return ret;
}

} // namespace reference
} // namespace isaac

