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
 ** \file SortedReferenceXml.hh
 **
 ** SortedReference.xml helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_SORTED_REFERENCE_XML_HH
#define iSAAC_REFERENCE_SORTED_REFERENCE_XML_HH

#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include "oligo/Kmer.hh"

namespace isaac
{
namespace reference
{

class SortedReferenceXml : boost::property_tree::ptree
{
public:
    static const unsigned currentReferenceFormatVersion_ = 2;

    SortedReferenceXml()
    {
        const unsigned formatVersion = currentReferenceFormatVersion_;
        put("SortedReference.FormatVersion", formatVersion);
    }
    struct Contig
    {
        Contig() : index_(0), karyotypeIndex_(0), offset_(0), genomicPosition_(0), totalBases_(0), acgtBases_(0){}
        Contig(const unsigned int index, const unsigned karyotypeIndex, const std::string &name,
               const boost::filesystem::path &filePath, const unsigned long offset,
               const unsigned long genomicPosition, const unsigned long totalBases, const unsigned long acgtBases,
               const std::string &bamSqAs, const std::string &bamSqUr, const std::string &bamM5) :
                   index_(index), karyotypeIndex_(karyotypeIndex), name_(name), filePath_(filePath),
                   offset_(offset),
                   genomicPosition_(genomicPosition), totalBases_(totalBases), acgtBases_(acgtBases),
                   bamSqAs_(bamSqAs), bamSqUr_(bamSqUr), bamM5_(bamM5) {}
        unsigned int index_;
        unsigned karyotypeIndex_;
        std::string name_;
        boost::filesystem::path filePath_;
        unsigned long offset_;
        unsigned long genomicPosition_;
        unsigned long totalBases_;
        unsigned long acgtBases_;
        std::string bamSqAs_;
        std::string bamSqUr_;
        std::string bamM5_;
    };

    typedef std::vector<Contig> Contigs;

    void putContig(const unsigned long genomicOffset,
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
                   const std::string &bamM5);
    void addMaskFile(const std::string &permutationName, const unsigned int maskWidth,
                     const isaac::oligo::Kmer mask, const boost::filesystem::path &filePath,
                     const size_t kmers, unsigned maxPrefixRangeCount);

    unsigned int getDefaultMaskWidth() const;
    const boost::filesystem::path getMaskFile(const std::string &permutationName, const unsigned int maskWidth,
                                              const isaac::oligo::Kmer mask) const;
    /**
     ** \brief Return a list of contig where each contig is at the corresponding index.
     **
     ** Precondition: the contigs in the current instance are sequentially
     ** indexed from 0 and there are no duplicates.
     **/
    Contigs getContigs() const;

    size_t getContigsCount() const {return getContigs().size();}

    /**
     ** \brief Return a list of contig where each contig is at the corresponding karyotype index.
     **
     ** Precondition: the contigs in the current instance are sequentially
     ** indexed from 0 and there are no duplicates.
     **/
    Contigs getKaryotypeOrderedContigs() const;
    /**
     ** \return Total number of kmers in all mask files
     **/
    unsigned long getTotalKmers() const;

    /**
     ** \brief Return the max of all MaxPrefixRangeCount across all masks
     **/
    unsigned getMaxPrefixRangeCount() const;

    struct MaskFile
    {
        boost::filesystem::path path;
        unsigned maskWidth;
        unsigned mask;
        size_t kmers;
        unsigned maxPrefixRangeCount;
    };
    std::vector<MaskFile> getMaskFileList(const std::string &permutationName) const;

    friend std::ostream &operator << (std::ostream &os, const SortedReferenceXml &tree);
    friend std::istream &operator >> (std::istream &is, SortedReferenceXml &indexedTree);
};

typedef std::vector<SortedReferenceXml> SortedReferenceXmlList;

std::ostream &operator << (std::ostream &os, const SortedReferenceXml &tree);
std::istream &operator >> (std::istream &is, SortedReferenceXml &indexedTree);

reference::SortedReferenceXml loadSortedReferenceXml(
    const boost::filesystem::path &xmlPath);

inline size_t genomeLength(const SortedReferenceXml::Contigs &contigList)
{
    return std::accumulate(
        contigList.begin(), contigList.end(),
        size_t(0), boost::bind<size_t>(std::plus<size_t>(), _1, boost::bind(&SortedReferenceXml::Contig::totalBases_, _2)));
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_SORTED_REFERENCE_XML_HH
