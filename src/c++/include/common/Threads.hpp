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
 ** \file Threads.hpp
 **
 ** Helpers for thread management.
 **
 ** \author Roman Petorvski
 **/

#ifndef iSAAC_COMMON_THREADS_HPP
#define iSAAC_COMMON_THREADS_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/thread.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac {
namespace common {

namespace detail
{
/**
 * \brief Base non-template class for ScopeEndCall
 */
class ScopeEndCallBase
{
public:
    virtual ~ScopeEndCallBase(){

    }

    operator bool () const
    {
        return false;
    }
};

/**
 * \brief Holds a functor which is called at the destruction time
 */
template <typename FuncT>
class ScopeEndCall : public ScopeEndCallBase
{
    FuncT f_;
    ScopeEndCall();
public:
    ScopeEndCall(FuncT f) : f_(f)
    {
    }

    ~ScopeEndCall()
    {
        f_();
    }
};


/**
 * \brief Helper to create ScopeEndCallHolder
 */
template <typename FuncT>
const detail::ScopeEndCall<FuncT> makeScopeEndCallHolder(FuncT f)
{
    return detail::ScopeEndCall<FuncT>(f);
}


}

/**
 * \brief ensures f is called during the stack unwind of the scope following the macro
 */
#define ISAAC_BLOCK_WITH_CLENAUP(f) \
    if(const common::detail::ScopeEndCallBase &b = common::detail::makeScopeEndCallHolder(f)) {(void)b;} else


/**
 * \brief Inversion of the boost::lock_guard
 */
template<typename Mutex>
class unlock_guard
{
private:
    Mutex& m;

    explicit unlock_guard(unlock_guard&);
    unlock_guard& operator=(unlock_guard&);
public:
    explicit unlock_guard(Mutex& m_):
        m(m_)
    {
        m.unlock();
    }

    ~unlock_guard()
    {
        m.lock();
    }
};


/**
 * \brief Use ThreadVector to execute a parallel operation using a vector of pre-allocated threads.
 */
template <bool crashOnExceptions>
class BasicThreadVector : boost::noncopyable, boost::ptr_vector<boost::thread>
{
    struct Executor
    {
        virtual void execute(size_type threadNum) = 0;
        virtual ~Executor() {}
    } *executor_;
    boost::mutex mutex_;
    boost::condition_variable stateChangedCondition_;

    // number of threads currently processing the request
    unsigned busyThreads_;

    // number of threads still required to process the request
    unsigned neededThreads_;

    // when executing with less threads than available, this prevents the higher number threads to carry out the request
    unsigned lowestBlockedThreadNumber_;

    // true when the whole thing goes down
    bool terminateRequested_;

    // constantly-incrementing number to make sure each thread processes one master call only once
    unsigned currentRequest_;

    boost::exception_ptr firstThreadException_;

    typedef boost::ptr_vector<boost::thread> base_type;
    typedef base_type::size_type size_type;
public:
    using base_type::size;
    /**
     * \brief Constructs a vector of size threads. All memory allocations that are required happen
     *        at this point.
     */
    BasicThreadVector(size_type size) : executor_(0), busyThreads_(size), neededThreads_(0),
        lowestBlockedThreadNumber_(0), terminateRequested_(false), currentRequest_(0)
    {
        boost::unique_lock<boost::mutex> lock(mutex_);
        ISAAC_ASSERT_MSG(0 < size, "Inadequate pool size");
        while (size--)
        {
            push_back(new boost::thread(boost::bind(&BasicThreadVector::threadFunc, this, size)));
        }
        // Initial wait for all threads to initialize their 'processedRequest'
        waitAll(lock);
    }

    /**
     * \brief Tells all threads to terminate and releases them.
     */
    ~BasicThreadVector()
    {
        {
            // TODO: locking is not really needed here, but it helps the assert to be a bit more definite
            boost::unique_lock<boost::mutex> lock(mutex_);
            ISAAC_ASSERT_MSG(!busyThreads_, "Workers must not be running at this point");

            terminateRequested_ = true;
            stateChangedCondition_.notify_all();
        }
        std::for_each(begin(), end(), boost::bind(&boost::thread::join, _1));
    }

    /**
     * \brief Executes func on requested number of threads. Only one execute call at a time is allowed. When
     *        execute returns, threads are guarranteed to have performed func once and only once,
     *        Multiple calls to execute are allowed.
     *
     * \threads number of threads to use. This must be less or equal than size()
     * \param func Unary function to execute. ThreadVector will supply a unique number in
     *             range [0, size()) to each invocation
     */

    template <typename F> void execute(F func, unsigned const threads)
    {
        ISAAC_ASSERT_MSG(threads <= size(), "Request must not exceed the amount of threads available");
        ISAAC_ASSERT_MSG(!executor_, "Queueing is not supported");
        struct FuncExecutor : public Executor
        {
            F &func_;
            FuncExecutor(F &func) : func_(func){}
            virtual void execute(size_type threadNum)
            {
                func_(threadNum);
            }
        }executor(func);

        executor_ = &executor;

        cycle(threads);

        executor_ = 0;
    }

    /**
     * \brief Executes func on requested number of size() threads.
     **/
    template <typename F> void execute(F func)
    {
        execute(func, size());
    }

private:
    void cycle(const unsigned threads)
    {
        boost::unique_lock<boost::mutex> lock(mutex_);
        ISAAC_ASSERT_MSG(!busyThreads_, "Only one at a time outstanding request is allowed");

        firstThreadException_ = boost::exception_ptr();
        if (1 == threads)
        {
            // Special case for one to simplify debugging. Just do it on the calling thread.
            executor_->execute(0);
        }
        else
        {
            lowestBlockedThreadNumber_ = threads;
            neededThreads_ = threads;
            ++currentRequest_;
            stateChangedCondition_.notify_all();
            waitAll(lock);
        }
        if (firstThreadException_)
        {
            ISAAC_THREAD_CERR << "WARNING: rethrowing a thread exception " << std::endl;
            boost::rethrow_exception(firstThreadException_);
        }
    }

    void waitAll(boost::unique_lock<boost::mutex> &lock)
    {
        while(busyThreads_ || neededThreads_){stateChangedCondition_.wait(lock);}
    }

    /**
     * \brief executes and allows the exception to escape
     */
    void unsafeExecute(size_type threadNum, boost::unique_lock<boost::mutex> &lock)
    {
        unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
        executor_->execute(threadNum);
    }

    /**
     * \brief executes and stores the exception information so that it can be rethrown on the main thread
     */
    void safeExecute(size_type threadNum, boost::unique_lock<boost::mutex> &lock)
    {
        try
        {
            unsafeExecute(threadNum, lock);
        }
        catch (...)
        {
            if (!firstThreadException_)
            {
                firstThreadException_ = boost::current_exception();
                ISAAC_THREAD_CERR << "ERROR: Thread: " << threadNum << " caught an exception first" << std::endl;
            }
            else
            {
                ISAAC_THREAD_CERR << "ERROR: Thread: " << threadNum << " also caught an exception" << std::endl;
            }
        }
    }

    void threadFunc(size_type threadNum)
    {
//        ISAAC_THREAD_CERR << "thread " << threadNum << " created\n";
        boost::unique_lock<boost::mutex> lock(mutex_);
        while(!terminateRequested_)
        {
            ISAAC_ASSERT_MSG(busyThreads_, "Thread is not accounted for!!!");
            --busyThreads_;
            const unsigned processedRequest = currentRequest_;
//            ISAAC_THREAD_CERR << "thread " << threadNum << " waiting for new request\n";
            stateChangedCondition_.notify_all();
            while(!terminateRequested_ && processedRequest == currentRequest_ ){stateChangedCondition_.wait(lock);}

            ++busyThreads_;

            if (!terminateRequested_)
            {
//                ISAAC_THREAD_CERR << "thread " << threadNum << " unblocked\n";
                if (lowestBlockedThreadNumber_ > threadNum)
                {
                    ISAAC_ASSERT_MSG(neededThreads_, "If thread is allowed to run, there must be a need for it!");
                    --neededThreads_;

                    if (crashOnExceptions)
                    {
                        unsafeExecute(threadNum, lock);
                    }
                    else
                    {
                        safeExecute(threadNum, lock);
                    }
                    // we're back under lock
                }
            }
        }

//        ISAAC_THREAD_CERR << "thread " << threadNum << " terminated\n";
    }

};

typedef BasicThreadVector<false> SafeThreadVector;
typedef BasicThreadVector<true> UnsafeThreadVector;
typedef UnsafeThreadVector ThreadVector;


} // namespacecommon
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_THREADS_HPP
