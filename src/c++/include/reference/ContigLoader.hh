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
 ** \file ContigsLoader.hh
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_CONTIGS_PRINTER_HH
#define iSAAC_REFERENCE_CONTIGS_PRINTER_HH

#include "common/Threads.hpp"
#include "reference/Contig.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

std::vector<std::vector<reference::Contig> > loadContigs(
    const reference::SortedReferenceXmlList &sortedReferenceXmlList,
    common::ThreadVector &loadThreads);

std::vector<reference::Contig> loadContigs(
    const reference::SortedReferenceXml &sortedReferenceXml,
    common::ThreadVector &loadThreads);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIGS_PRINTER_HH
