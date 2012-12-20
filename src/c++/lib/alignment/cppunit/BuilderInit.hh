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
 **/

#include <vector>
#include <string>
#include <boost/foreach.hpp>

#include "flowcell/ReadMetadata.hh"
#include "alignment/SeedMetadata.hh"
#include "reference/Contig.hh"
#include "oligo/Nucleotides.hh"

inline isaac::flowcell::ReadMetadataList getReadMetadataList(const unsigned l0 = 100, const unsigned l1 = 100)
{
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, l0, 0, 0))
            (isaac::flowcell::ReadMetadata(l0 + 1, l0 + l1, 1, l0))
            ;
    return ret;
}

inline isaac::alignment::SeedMetadataList getSeedMetadataList()
{
    std::vector<isaac::alignment::SeedMetadata> ret =
        boost::assign::list_of
        (isaac::alignment::SeedMetadata( 0, 32, 0, 0))
        (isaac::alignment::SeedMetadata(32, 32, 0, 1))
        (isaac::alignment::SeedMetadata(64, 32, 0, 2))
        (isaac::alignment::SeedMetadata( 0, 32, 1, 3))
        (isaac::alignment::SeedMetadata(32, 32, 1, 4))
        (isaac::alignment::SeedMetadata(64, 32, 1, 5))
        ;
    return ret;
}

inline isaac::reference::Contig getContig(const std::string name, const unsigned length)
{
    char bases[] = {'A', 'C', 'G', 'T'};
    isaac::reference::Contig contig(0, name);
    contig.forward_.resize(length);
    BOOST_FOREACH(char &base, contig.forward_)
    {
        base = bases[rand() % 4];
    }
    return contig;
}

inline void show(const std::vector<char> &s)
{
    BOOST_FOREACH(char c, s) std::cerr << c;
}

inline std::vector<char> reverseComplement(const std::vector<char> &forward)
{
    std::vector<char> tr(256, 'N');
    tr['A'] = 'T';
    tr['C'] = 'G';
    tr['G'] = 'C';
    tr['T'] = 'A';
    std::vector<char> reverse;
    reverse.reserve(forward.size());
    for (std::vector<char>::const_reverse_iterator b = forward.rbegin(); forward.rend() != b; ++b)
    {
        reverse.push_back(tr[*b]);
    }
    return reverse;
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

inline std::string substr(const std::vector<char> &from, std::string::size_type __pos = 0,
                   std::string::size_type __n = std::string::npos)
{
    return std::string(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

inline std::vector<char> subv(const std::string &from, std::string::size_type __pos = 0,
                       std::string::size_type __n = std::string::npos)
{
    return std::vector<char>(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

inline std::vector<char> subv(const std::vector<char> &from, std::string::size_type __pos,
                       std::string::size_type __n)
{
    return std::vector<char>(from.begin() + __pos, from.begin() + __pos + __n);
}

inline std::vector<char> subv(const std::vector<char> &from, std::string::size_type __pos)
{
    return std::vector<char>(from.begin() + __pos, from.end());
}

inline std::vector<char> operator +(const std::vector<char> &right, const std::vector<char> &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}

inline std::vector<char> operator +(const std::vector<char> &right, const std::string &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}

inline std::vector<isaac::reference::Contig> getContigList(const unsigned l0 = 210, const unsigned l1 = 220, const unsigned l4 = 60)
{
    const isaac::reference::Contig c2 = getContig("c2", 230);
    isaac::reference::Contig c3(0, "c3");
    c3.forward_ = vectorFromString(std::string("AAAAA")) + c2.forward_;
    return boost::assign::list_of
        (getContig("c0", l0))
        (getContig("c1", l1))
        (c2)
        (c3)
        (getContig("c4", l4));
}

template <typename ContainerT>
std::vector<char> getBcl(const ContainerT &bases)
{
    std::vector<char> bcl;
    bcl.reserve(bases.size());
    BOOST_FOREACH(char b, bases)
    {
        using isaac::oligo::getValue;
        bcl.push_back((40 << 2) | getValue(b));
    }
    return bcl;
}

inline std::vector<char> getBcl(
    const std::vector<isaac::flowcell::ReadMetadata> &readMetadataList,
    const std::vector<isaac::reference::Contig> &contigList,
    const unsigned contigId, const int offset0, const int offset1,
    const bool reverse0 = false, const bool reverse1 = true)
{
    const isaac::reference::Contig &contig = contigList[contigId];
    const unsigned length0 = readMetadataList[0].getLength();
    const unsigned length1 = readMetadataList[1].getLength();

    const std::vector<char> &forward = contig.forward_;
    std::vector<char> reverse;
    reverse.reserve(forward.size());
    BOOST_REVERSE_FOREACH(const char base, forward)
    {
        using isaac::oligo::getReverseBase;
        using isaac::oligo::getValue;
        reverse.push_back(getReverseBase(getValue(base)));
    }

    const std::vector<char> &s0 = reverse0 ? reverse : forward;
    const std::vector<char> &s1 = reverse1 ? reverse : forward;
    std::string bases(s0.begin() + offset0, s0.begin() + offset0 + length0);
    bases += std::string(s1.begin() + offset1, s1.begin() + offset1 + length1);
    return getBcl(bases);
}
