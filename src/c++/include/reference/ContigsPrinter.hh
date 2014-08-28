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
 ** \file ContigsPrinter.hh
 **
 ** Top level component that produces contigs metadata.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_CONTIGS_PRINTER_HH
#define iSAAC_REFERENCE_CONTIGS_PRINTER_HH

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

namespace isaac
{
namespace reference
{

class ContigsPrinter: boost::noncopyable
{
public:
    ContigsPrinter(
        const boost::filesystem::path &originalSortedReferenceXml,
        const boost::filesystem::path &genomeFile
    );
    void run();
private:
    const boost::filesystem::path originalSortedReferenceXml_;
    const boost::filesystem::path genomeFile_;
};

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIGS_PRINTER_HH
