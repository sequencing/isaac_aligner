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

const SequencingAdapterMetadataList NEXTERA_STANDARD_ADAPTERS =
    boost::assign::list_of(SequencingAdapterMetadata("CTGTCTCTTATACACATCT", false, 0))
                          (SequencingAdapterMetadata("AGATGTGTATAAGAGACAG",true, 0));

const SequencingAdapterMetadataList NEXTERA_MATEPAIR_ADAPTERS =
    boost::assign::list_of(SequencingAdapterMetadata("CTGTCTCTTATACACATCT", false))
                          (SequencingAdapterMetadata("AGATGTGTATAAGAGACAG",false));

} // namespace flowcell
} // namespace isaac
