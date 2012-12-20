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
 ** \file DirectFragmentStorage.hh
 **
 ** \brief Immediate storing of aligned data in the output bin files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_STORAGE_HH

#include "alignment/BamTemplate.hh"
#include "alignment/BinMetadata.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

/**
 * \brief interface for various fragment storage implementations
 *        TODO: remove this when the implementation settles
 */
class FragmentStorage
{
protected:
    // prevent destruction of children in this base
    virtual ~FragmentStorage(){}
public:
    virtual const alignment::BinMetadataList &getBinPathList() const = 0;

    virtual void add(const BamTemplate &bamTemplate, const unsigned barcodeIdx) = 0;

    virtual void prepareFlush() = 0;
    virtual void flush() = 0;
    virtual void resize(const unsigned long clusters) = 0;
    virtual void unreserve() = 0;
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_STORAGE_HH
