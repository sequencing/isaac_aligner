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
