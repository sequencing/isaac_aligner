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
