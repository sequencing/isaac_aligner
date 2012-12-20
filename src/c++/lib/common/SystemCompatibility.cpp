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
 ** \file SystemCompatibility.cpp
 **
 ** \brief See SystemCompatibility.hh
 ** 
 ** \author Come Raczy
 **/
#include <stdio.h>

#include <new>
#include <iostream>

#include <boost/thread.hpp>

#include "common/SystemCompatibility.hh"
#include "common/config.h"

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#else // #ifdef HAVE_SIGNAL_H
#error Only POSIX systems are supported. The header <signal.h> is required.
#endif // #ifdef HAVE_SIGNAL_H

/*
 * TODO: separate the implementation into different files, depending on the system
 */
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#else // #ifdef HAVE_MALLOC_H
#error Only POSIX systems are supported. The header <malloc.h> is required.
#endif // #ifdef HAVE_MALLOC_H

#ifdef HAVE_MCHECK_H
#include <mcheck.h>
#else // #ifdef HAVE_MCHECK_H
#error Only POSIX systems are supported. The header <mcheck.h> is required.
#endif // #ifdef HAVE_MCHECK_H

#include <sys/resource.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else // #ifdef HAVE_UNISTD_H
#error Only POSIX systems are supported. The header <unistd.h> is required.
#endif // #ifdef HAVE_UNISTD_H

#ifdef HAVE_TIME_H
#include <time.h>
#else // #ifdef HAVE_TIME_H
#error The header <time.h> is required.
#endif // #ifdef HAVE_TIME_H

#include <cassert>

namespace isaac
{
namespace common
{

unsigned int getMaxOpenFiles()
{
#ifdef HAVE_SYSCONF
    assert(0 < sysconf(_SC_OPEN_MAX));
    return sysconf(_SC_OPEN_MAX);
#else
#error 'sysconf' is required
#endif
}

long clock()
{
#ifdef HAVE_CLOCK
    return ::clock();
#else
#error 'clock' is required (from <time.h>)
#endif
}

bool isLittleEndian()
{
    const unsigned long v = 0x0706050403020100;
    const unsigned char * const p = reinterpret_cast<const unsigned char *>(&v);
    for (unsigned i = 0; i < sizeof(v); ++i)
    {
        if (p[i] != i)
        {
            return false;
        }
    }
    return true;
}

bool ulimitV(unsigned long availableMemory)
{
    const rlimit rl = {availableMemory, availableMemory};
    return !setrlimit(RLIMIT_AS, &rl);
}

bool ulimitV(unsigned long *pLimit)
{
    rlimit rl = {0, 0};
    if (-1 == getrlimit(RLIMIT_AS, &rl))
    {
        return false;
    }
    *pLimit = RLIM_INFINITY == rl.rlim_cur ? -1 : rl.rlim_cur;
    return true;
}

void mtrace()
{
    ::mtrace();
}

void muntrace()
{
    ::muntrace();
}

static boost::mutex block_malloc_hook_mutex_;

static void* (*old_malloc_hook_)(size_t, const void*) = 0;
static bool (*user_hook_)(size_t size, const void *caller) = 0;

unsigned mallocCount_(0);
static void * malloc_hook(size_t size, const void *caller)
{
    boost::unique_lock<boost::mutex> lock(block_malloc_hook_mutex_);
    ++mallocCount_;

    __malloc_hook = old_malloc_hook_;

    assert(user_hook_);
    if (!user_hook_(size, caller))
    {
        terminateWithCoreDump();
        // in case termination did not terminate...
        std::cerr << "ERROR: blocked allocation of " << size << " bytes. Returning 0 \n";
        __malloc_hook = malloc_hook;
        return 0;
    }
    else
    {
        void *result  = malloc(size);
        __malloc_hook = malloc_hook;
        return result;
    }
}

void hookMalloc(bool (*hook)(size_t size, const void *caller))
{
    boost::unique_lock<boost::mutex> lock(block_malloc_hook_mutex_);

    assert(__malloc_hook != malloc_hook);
    assert(!user_hook_);
    old_malloc_hook_ = __malloc_hook;
    user_hook_ = hook;
    mallocCount_ = 0;
    __malloc_hook = malloc_hook;

    //std::cerr << "malloc blocked\n";
}

unsigned unhookMalloc(bool (*hook)(size_t size, const void *caller))
{
    boost::unique_lock<boost::mutex> lock(block_malloc_hook_mutex_);

    assert(__malloc_hook == malloc_hook);
    assert(user_hook_ == hook);

    __malloc_hook = old_malloc_hook_;
    user_hook_ = 0;

//    std::cerr << "malloc unblocked\n";
    return mallocCount_;
}

void terminateWithCoreDump()
{
    raise(SIGSEGV);
}
} // namespace common
} // namespace isaac

