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
 ** \file SequencingAdapterMetadata.cpp
 **
 ** Packaging of the metadata associated to a barcode.
 **
 ** \author Roman Petrovski
 **/

#include <boost/assign.hpp>

#include "flowcell/SequencingAdapterMetadata.hh"

namespace isaac
{
namespace flowcell
{

const SequencingAdapterMetadataList STANDARD_ADAPTERS =
    boost::assign::list_of(SequencingAdapterMetadata("AGATCGGAAGAGC", false, 0))
                          (SequencingAdapterMetadata("GCTCTTCCGATCT",true, 0));

const SequencingAdapterMetadataList NEXTERA_STANDARD_ADAPTERS =
    boost::assign::list_of(SequencingAdapterMetadata("CTGTCTCTTATACACATCT", false, 0))
                          (SequencingAdapterMetadata("AGATGTGTATAAGAGACAG",true, 0));

const SequencingAdapterMetadataList NEXTERA_MATEPAIR_ADAPTERS =
    boost::assign::list_of(SequencingAdapterMetadata("CTGTCTCTTATACACATCT", false))
                          (SequencingAdapterMetadata("AGATGTGTATAAGAGACAG",false));

} // namespace flowcell
} // namespace isaac
