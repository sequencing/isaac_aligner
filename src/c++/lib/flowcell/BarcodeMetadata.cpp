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
 ** \file BarcodeMetadata.cpp
 **
 ** Packaging of the metadata associated to a barcode.
 **
 ** \author Roman Petrovski
 **/

#include "flowcell/BarcodeMetadata.hh"

namespace isaac
{
namespace flowcell
{

const std::string BarcodeMetadata::NO_INDEX_BARCODE("none");
const std::string BarcodeMetadata::UNKNOWN_BARCODE("unknown");
const std::string BarcodeMetadata::UNKNOWN_SAMPLE("unknown");
const std::string BarcodeMetadata::DEFAULT_PROJECT("default");

} // namespace flowcell
} // namespace isaac
