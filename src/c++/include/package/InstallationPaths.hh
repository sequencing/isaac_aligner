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
 ** \file InstallationPaths.hh
 **
 ** \brief Path resolution for installed components
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_PACKAGE_INSTALLATION_PATHS_HH
#define iSAAC_PACKAGE_INSTALLATION_PATHS_HH

#include <boost/filesystem.hpp>

namespace isaac
{
namespace package
{

void initialize(const boost::filesystem::path &modulePath, const char *homeOffset);
boost::filesystem::path expandPath(const char *path);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_PACKAGE_INSTALLATION_PATHS_HH
