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
 ** \file SystemCompatibility.cpp
 **
 ** \brief See SystemCompatibility.hh
 ** 
 ** \author Come Raczy
 **/
#include <stdio.h>

#include <new>
#include <iostream>

#include <boost/format.hpp>
#include <boost/thread.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
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

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#else // #ifdef HAVE_SYS_STAT_H
#error The header <sys/stat.h> is required.
#endif // #ifdef HAVE_SYS_STAT_H


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


void terminateWithCoreDump()
{
    raise(SIGSEGV);
}

unsigned long getFileSize(const char *filePath)
{
#ifdef HAVE_STAT
    struct stat s;
    if (0 != stat(filePath, &s))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, std::string("Failed to stat file ") + filePath));
    }
    return s.st_size;
#else
#error 'stat' is required
#endif
}

boost::filesystem::path getModuleFileName()
{
    char szBuffer[10240];
    ISAAC_ASSERT_MSG(-1 != readlink("/proc/self/exe", szBuffer, sizeof(szBuffer)), "TODO: handle the readlink error: " << errno);
    return szBuffer;
}

} // namespace common
} // namespace isaac

#ifndef ISAAC_CYGWIN

namespace isaac
{
namespace common
{


void configureMemoryManagement(
    const bool disableMultipleArenas,
    const bool disableFastbins)
{
    if (disableMultipleArenas)
    {
        // By default linux creates  ((NUMBER_OF_CPU_CORES) * (sizeof(long) == 4 ? 2 : 8)) arenas to allow for
        // (I guess) faster concurrent memory management. As iSAAC pre-allocates all the memory and spends majority of
        // time doing something other than memory allocation/deallocation, having more than one arene is purely a waste
        // of virtual memory. Virtual memory is important for iSAAC as it normally tries to use all of the virtual memory
        // that is allowed to the process. So, this trick should save in the range of 2 gigabytes.
        // TODO: Interestingly enough it actually saves more and the value seems to depend on the number of threads.
        // Would be nice to figure out why...
        mallopt(M_ARENA_MAX, 1);
    }

    if (disableFastbins)
    {
        // in Linux, fastbins are enabled by default in order to improve performance of
        // applications that allocate and free small chunks of memory on multiple threads often.
        // This causes unnecessary memory overhead and fragmentation which easily amounts to
        // a loss of 10-20 gigabytes of RAM in 400-plex bam generation.
        // As isaac, where it matters, allocates memory up-front and controls concurrency,
        // fastbins are disabled.
        mallopt(M_MXFAST, 0);
    }
}

bool ulimitV(unsigned long availableMemory)
{
    const rlimit rl = {availableMemory, availableMemory};
    if(setrlimit(RLIMIT_AS, &rl))
    {
        BOOST_THROW_EXCEPTION(isaac::common::ResourceException(
            errno, (boost::format("Failed to set the memory consumption limit to: %d bytes") % availableMemory).str()));
    }
    return true;
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

    std::cerr << "malloc unblocked\n";
    return mallocCount_;
}

} // namespace common
} // namespace isaac

#else //ISAAC_CYGWIN

#include <w32api/windows.h>

namespace isaac
{
namespace common
{

void configureMemoryManagement(
    const bool disableMultipleArenas,
    const bool disableFastbins)
{
    if (disableMultipleArenas)
    {
        // By default linux creates  ((NUMBER_OF_CPU_CORES) * (sizeof(long) == 4 ? 2 : 8)) arenas to allow for
        // (I guess) faster concurrent memory management. As iSAAC pre-allocates all the memory and spends majority of
        // time doing something other than memory allocation/deallocation, having more than one arene is purely a waste
        // of virtual memory. Virtual memory is important for iSAAC as it normally tries to use all of the virtual memory
        // that is allowed to the process. So, this trick should save in the range of 2 gigabytes.
        // TODO: Interestingly enough it actually saves more and the value seems to depend on the number of threads.
        // Would be nice to figure out why...
        //mallopt(M_ARENA_MAX, 1);
        ISAAC_THREAD_CERR << "WARNING: M_ARENA_MAX is not available in CYGWIN." << std::endl;
    }

    if (disableFastbins)
    {
        // in Linux, fastbins are enabled by default in order to improve performance of
        // applications that allocate and free small chunks of memory on multiple threads often.
        // This causes unnecessary memory overhead and fragmentation which easily amounts to
        // a loss of 10-20 gigabytes of RAM in 400-plex bam generation.
        // As isaac, where it matters, allocates memory up-front and controls concurrency,
        // fastbins are disabled.
        mallopt(M_MXFAST, 0);
    }
}

/**$
 * \brief Exception thrown when there is insufficient resources to perform an operation. For example$
 *        if the adjusting the soft ulimit fails due to a set hard limit$
 */
class Win32Exception: public std::exception, public ExceptionData
{
public:
    Win32Exception(DWORD dwLastError, const std::string &message):
        ExceptionData(errno, message + ": " + (boost::format("0x%x") % dwLastError).str() + ": "+ getLastErrorText(dwLastError))
    {

    }
    
    static std::string getLastErrorText(DWORD nErrorCode)
    {
        char* msg;
        // Ask Windows to prepare a standard message for a GetLastError() code:
        if (!FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, 
                       NULL, nErrorCode, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&msg, 0, NULL))
        {
            return (boost::format("FormatMessage failed for error code: 0x%x") % nErrorCode).str();
        }
        else
        {
            return msg;
        }
    }
};

bool ulimitV(unsigned long availableMemory)
{
    const std::string eventName = "Local\\iSAACForkResumeEvent";
    HANDLE hWaitEvent = CreateEvent(NULL, FALSE, FALSE, eventName.c_str());

    pid_t childPid = fork();
    if (!childPid)
    {
        // Child needs to reopen event as otherwise the hWaitEvent is inaccessible
        hWaitEvent = OpenEvent(EVENT_ALL_ACCESS, true, eventName.c_str());
        if (INVALID_HANDLE_VALUE == hWaitEvent)
        {
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to open fork sycnrhonization event handle")).str()));
        }

        DWORD dwWaitResult = WaitForSingleObject(hWaitEvent, INFINITE);
        if (WAIT_OBJECT_0 != dwWaitResult)
        {
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to WaitForSingleObject on fork sycnrhonization event. Wait result: 0x%x") % dwWaitResult).str()));
        }
        CloseHandle(hWaitEvent);
        // Parent has finished attaching child to the job. Child is free to run...
        return true;
    }
    else
    {
        HANDLE hProcess = OpenProcess(PROCESS_ALL_ACCESS, TRUE, childPid);
        if (INVALID_HANDLE_VALUE == hProcess)
        {
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to open child process handle for process id %x") % childPid).str()));
        }

        HANDLE hJobObject = CreateJobObject(NULL, NULL);
        if (INVALID_HANDLE_VALUE == hJobObject)
        {
            CloseHandle(hProcess);
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to CreateJobObject")).str()));
        }

        JOBOBJECT_EXTENDED_LIMIT_INFORMATION extendedLimits = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        extendedLimits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_JOB_MEMORY;
        extendedLimits.JobMemoryLimit = availableMemory;
        BOOL res = SetInformationJobObject(hJobObject, JobObjectExtendedLimitInformation, &extendedLimits, sizeof(extendedLimits));

        if (!res)
        {
            CloseHandle(hJobObject);
            CloseHandle(hProcess);
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to SetInformationJobObject")).str()));
        }

        res = AssignProcessToJobObject(hJobObject, hProcess);
        if (!res)
        {
            CloseHandle(hJobObject);
            CloseHandle(hProcess);
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("AssignProcessToJobObject failed")).str()));
        }

        CloseHandle (hJobObject);


        if (!SetEvent(hWaitEvent))
        {
            CloseHandle(hProcess);
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("SetEvent on fork synchronizaqtion event failed")).str()));
        }

        if (WAIT_OBJECT_0 != WaitForSingleObject(hProcess, INFINITE))
        {
            CloseHandle(hProcess);
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed wile waiting for child process termination")).str()));
        }

        DWORD dwExitCode = 0;
        if (!GetExitCodeProcess(hProcess, &dwExitCode))
        {
            CloseHandle(hProcess);
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to get child process exit code")).str()));
        }

        CloseHandle(hProcess);
        if (dwExitCode)
        {
            CloseHandle(hWaitEvent);
            BOOST_THROW_EXCEPTION(Win32Exception(NOERROR, (boost::format("Child process exited with non-zero exit code: 0x%x") % dwExitCode).str()));
        }

        CloseHandle(hWaitEvent);


        // signal the caller that everything completed successfully and the parent process must exit now
        return false;
    }
}

bool ulimitV(unsigned long *pLimit)
{
    *pLimit = 0;
    return true;
}

void hookMalloc(bool (*hook)(size_t size, const void *caller))
{
    // memory control is not supported under cygwin
}

unsigned unhookMalloc(bool (*hook)(size_t size, const void *caller))
{
    return 0;
}

} // namespace common
} // namespace isaac

#endif //ISAAC_CYGWIN


