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
 ** \file Contig.hh
 **
 ** \brief Definition of a contig
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_CONTIG_HH
#define iSAAC_REFERENCE_CONTIG_HH

#include <string>
#include <vector>

namespace isaac
{
namespace reference
{

struct Contig
{
    unsigned index_;
    std::string name_;
    std::vector<char> forward_;

    Contig(const unsigned index, const std::string &name) : index_(index), name_(name){;}
    size_t getLength() const {return forward_.size();}
};

typedef std::vector<reference::Contig> ContigList;

/// Total length of all the contigs of a genome
size_t genomeLength(const std::vector<Contig> &contigList);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIG_HH
