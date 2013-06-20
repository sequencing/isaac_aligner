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
 ** \file FileBufCache.hh
 **
 ** Vector of file buffers which are kept open to reduce the cost of closing/opening files.
 **
 ** \author Roman Petrovski
 **/

#include "io/FileBufCache.hh"

namespace isaac
{
namespace io
{

//const boost::filesystem::path FileBufWithReopenHolderBase::emptyPath_;

} // namespace io
} // namespace isaac
