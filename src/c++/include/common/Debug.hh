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
 ** \file Debug.hh
 **
 ** \brief Various debugging-related helpers
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_LOG_THREAD_TIMESTAMP_HH
#define iSAAC_LOG_THREAD_TIMESTAMP_HH

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/thread.hpp>

#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace common
{

//TODO: check why simple CerrLocker(std::cerr) << ... << is not good enough
/**
 * \brief helper macro to simplify the thread-guarded logging. All elements on a single << line are serialized
 * under one CerrLocker
 */
#define ISAAC_THREAD_CERR if(::isaac::common::detail::CerrLocker isaac_cerr_lock = ::isaac::common::detail::CerrLocker()); else std::cerr << isaac::common::detail::ThreadTimestamp()

/**
 * \brief Evaluates expression always (even if NDEBUG is set and so on). Also uses ostream serialization which,
 *        unlike the standard assert, has shown not to allocate the dynamic memory at the time when you least
 *        expect this to happen.
 */
#define ISAAC_ASSERT_MSG(expr, msg) {if (expr) {} else \
{ ISAAC_THREAD_CERR << "ERROR: ***** Internal Program Error - assertion (" << #expr << ") failed in " \
    << (BOOST_CURRENT_FUNCTION) << ":" << __FILE__ << '(' << __LINE__ << "): " << msg << std::endl; \
    ::isaac::common::terminateWithCoreDump();}}

inline std::string parseStat(const std::string &stat)
{
    std::vector<std::string> statFields;
    boost::algorithm::split(statFields, stat,  boost::algorithm::is_any_of(" "));
    return std::string(statFields.at(22) + "vm " + statFields.at(23) + "res");
}

#define ISAAC_TRACE_STAT(prefix) {\
    std::string statm; std::ifstream ifs("/proc/self/stat"); \
    std::getline(ifs, statm); \
    ISAAC_THREAD_CERR << "STAT: " << prefix << isaac::common::parseStat(statm) << std::endl;\
    }

class ScoopedMallocBlock : boost::noncopyable
{
public:
    enum Mode
    {
        Invalid = 0,
        Off,
        Warning,
        Strict
    };

    ScoopedMallocBlock(const Mode mode);
    ~ScoopedMallocBlock();
private:
    const Mode mode_;

    friend class ScoopedMallocBlockUnblock;
    void block();
    void unblock();
};

class ScoopedMallocBlockUnblock : boost::noncopyable
{
    ScoopedMallocBlock &block_;
public:
    ScoopedMallocBlockUnblock(ScoopedMallocBlock & block);
    ~ScoopedMallocBlockUnblock();
};

namespace detail
{
class ThreadTimestamp
{
public:
};

/**
 * \brief formats time stamp and thread id to simplify threaded logging
 */
inline std::ostream & operator << (std::ostream &os, const ThreadTimestamp &) {

    const boost::posix_time::ptime currentTime = boost::posix_time::second_clock::local_time();

    // IMPORTANT: this is the way to serialize date without causing any dynamic memory operations to occur
    const std::tm t = boost::posix_time::to_tm(currentTime.time_of_day());
    const boost::posix_time::ptime::date_type d = currentTime.date();
    os << d.year() << '-' <<
        std::setfill('0') << std::setw(2) << d.month().as_number() << '-'  <<
        std::setfill('0') << std::setw(2) << d.day() << ' '  <<

        std::setfill('0') << std::setw(2) << t.tm_hour << ':' <<
        std::setfill('0') << std::setw(2) << t.tm_min << ':' <<
        std::setfill('0') << std::setw(2) << t.tm_sec << ' ' <<
        "\t[" << boost::this_thread::get_id() << "]\t";
    return os;
}

/**
 * \brief Guards std::cerr for the duration of CerrLocker existance
 *        Restores any changes made to ios::base
 */
class CerrLocker
{
	// some people allocate memory from under their trace code. For example by using boost::format.
    // if memory control is on, we don't want them to be dead-locked on their own thread cerrMutex_.
    static boost::recursive_mutex cerrMutex_;
    boost::lock_guard<boost::recursive_mutex> lock_;
    boost::io::ios_base_all_saver ias_;

public:

    CerrLocker(const CerrLocker &that) : lock_(cerrMutex_), ias_(std::cerr){
    }
    CerrLocker() : lock_(cerrMutex_), ias_(std::cerr) {
    }
    operator bool () const {
        return false;
    }
};

inline void assertion_failed_msg(char const * expr, char const * msg, char const * function,
                                 char const * file, long line)
{
    ISAAC_THREAD_CERR
    << "ERROR: ***** Internal Program Error - assertion (" << expr << ") failed in "
    << function << ":" << file << '(' << line << "): " << msg << std::endl;

    common::terminateWithCoreDump();
}

} // namespace detail

/**
 ** \brief Provide a mechanism for detailed level of debugging
 **/
#ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    #define ISAAC_THREAD_CERR_DEV_TRACE(trace) {ISAAC_THREAD_CERR << trace << std::endl;}
    #define ISAAC_DEV_TRACE_BLOCK(block) block
    #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, trace) ISAAC_THREAD_CERR_DEV_TRACE(trace)
#else
    #define ISAAC_THREAD_CERR_DEV_TRACE(blah)
    #ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED_CLUSTER_ID
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, trace) {if(ISAAC_THREAD_CERR_DEV_TRACE_ENABLED_CLUSTER_ID == (clusterId)) {ISAAC_THREAD_CERR << trace << std::endl;}}
        #define ISAAC_DEV_TRACE_BLOCK(block) block
    #else
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, blah)
        #define ISAAC_DEV_TRACE_BLOCK(block)
    #endif
#endif

struct TimeSpec : public timespec
{

};

inline TimeSpec tsdiff(TimeSpec start, TimeSpec end)
{
    TimeSpec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

static const long NS_IN_SEC = 1000000000;

inline TimeSpec tsadd(
    const struct TimeSpec & t1,
    const struct TimeSpec & t2)
{
    struct TimeSpec sum = t1;

    sum.tv_sec  += t2.tv_sec;
    sum.tv_nsec += t2.tv_nsec;

    if (sum.tv_nsec >= NS_IN_SEC)
    {
        sum.tv_sec++;
        sum.tv_nsec -= NS_IN_SEC;
    }

    return sum;
}

inline std::ostream &operator <<(std::ostream & os, const timespec &ts)
{
    return os << ts.tv_sec << "." << std::setw(3) << std::setfill('0') << ts.tv_nsec / 1000000;
}
} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_LOG_THREAD_TIMESTAMP_HH
