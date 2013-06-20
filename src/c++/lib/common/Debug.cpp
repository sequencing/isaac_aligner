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
 ** \file Debug.cpp
 **
 ** \brief Various debugging-related helpers
 **
 ** \author Roman Petrovski
 **/

#include "common/Debug.hh"
#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace common
{
namespace detail
{

boost::recursive_mutex CerrLocker::cerrMutex_;


static bool malloc_warning_hook(size_t size, const void *caller)
{
    ISAAC_THREAD_CERR << "WARNING: blocked allocation of " << size << " bytes requested.\n";
    return true;
}

static bool malloc_strict_hook(size_t size, const void *caller)
{
    ISAAC_THREAD_CERR << "ERROR: blocked allocation of " << size << " bytes requested.\n";
    return false;
}

//static bool malloc_ignore_hook(size_t size, const void *caller)
//{
//    return true;
//}

} // namespace detail

ScoopedMallocBlock::ScoopedMallocBlock(const ScoopedMallocBlock::Mode mode) :
    mode_(mode)
{
    block();
}

void ScoopedMallocBlock::block()
{
    switch(mode_)
    {
    case Off:
//        hookMalloc(detail::malloc_ignore_hook);
        break;
    case Warning:
        hookMalloc(detail::malloc_warning_hook);
        break;
    case Strict:
        hookMalloc(detail::malloc_strict_hook);
        break;
    default:
        ISAAC_ASSERT_MSG(false, "invalid malloc block mode specified");
        break;
    }
}

ScoopedMallocBlock::~ScoopedMallocBlock()
{
    unblock();
}

void ScoopedMallocBlock::unblock()
{
    switch(mode_)
    {
    case Off:
//        unhookMalloc(detail::malloc_ignore_hook);
        break;
    case Warning:
        unhookMalloc(detail::malloc_warning_hook);
        break;
    case Strict:
        unhookMalloc(detail::malloc_strict_hook);
        break;
    default:
        ISAAC_ASSERT_MSG(false, "invalid malloc block mode specified");
        break;
    }
}

ScoopedMallocBlockUnblock::ScoopedMallocBlockUnblock(ScoopedMallocBlock &block) :
        block_(block)
{
    block_.unblock();
}

ScoopedMallocBlockUnblock::~ScoopedMallocBlockUnblock()
{
    block_.block();
}

} // namespace common
} // namespace isaac
