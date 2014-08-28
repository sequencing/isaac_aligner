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
 ** \file InstallationPaths.cpp
 **
 ** \brief Path resolution for installed components
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "package/InstallationPaths.hh"

namespace isaac
{
namespace package
{

static boost::filesystem::path installationRoot_;

void initialize(const boost::filesystem::path &modulePath, const char *homeOffset)
{
    ISAAC_ASSERT_MSG(installationRoot_.empty(), "Installation root is already set: " << installationRoot_);
    if (0 != *homeOffset)
    {
        installationRoot_ = modulePath.parent_path() / homeOffset;
    }
}

boost::filesystem::path expandPath(const char *path)
{
    if(installationRoot_.empty())
    {
        return path;
    }
    return installationRoot_ / path;
}

} // namespace package
} // namespace isaac
