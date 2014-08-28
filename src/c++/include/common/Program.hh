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
 ** \file Program.hh
 **
 ** Declaration of the skeleton of all c++ programs.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_PROGRAM_HH
#define iSAAC_COMMON_PROGRAM_HH

#include <string>
#include <iostream>
#include <cstdlib>
#include <boost/program_options.hpp>
#include <boost/noncopyable.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "config.h"
#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace common
{

namespace bpo = boost::program_options;

/**
 ** Encapsulation of the pocessing of the command line options.
 **
 ** TODO: add config file and environment options
 **/
class Options: boost::noncopyable
{
public:
    enum Action
    {
        RUN, HELP, VERSION, ABORT
    };
    Options();
    virtual ~Options()
    {
    }
    Action parse(int argc, const char * const argv[]);
    std::string usage() const;
protected:
    bpo::options_description namedOptions_;
    bpo::options_description unnamedOptions_;
    bpo::positional_options_description positionalOptions_;

    typedef boost::shared_ptr<boost::program_options::option_description> OptionDescriptionPtr;
    typedef std::vector<OptionDescriptionPtr > OptionDescriptionPtrs;
    std::string helpDefaults(const OptionDescriptionPtrs &options) const;
    std::string help(const OptionDescriptionPtrs &options, const bool markdown) const;
private:
    virtual std::string usagePrefix() const = 0;
    virtual std::string usageSuffix() const
    {
        return "";
    }
    virtual void postProcess(bpo::variables_map &)
    {
    }

    static const unsigned MARKDOWN_LINE_LENGTH = 120;
    // "parse" will store the state in vm_ so that "usage" can access the details of parsed command line
    bpo::variables_map vm_;
};

/**
 ** Unified behavior of all programs.
 **/
template<class O>
void run(void(*callback)(const O &), int argc, char *argv[])
{
    // when statically-linked, the centos-originated binaries fail to properly discover
    // locale files on ubuntu. See if being conservative helps.
    setenv("LC_ALL", "C", 1);
    ISAAC_THREAD_CERR << "Forcing LC_ALL to C" << std::endl;
    try
    {
        O options;
        const typename O::Action action = options.parse(argc, argv);
        if (O::RUN == action)
        {
            callback(options);
        }
        else if (O::HELP == action)
        {
            std::cout << options.usage() << std::endl;
        }
        else if (O::VERSION == action)
        {
            std::cout << iSAAC_VERSION_FULL << std::endl;
        }
        else
        {
//            std::clog << options.usage() << std::endl;
            exit(1);
        }
    }
    catch (const isaac::common::ExceptionData &exception)
    {
        std::clog << "Error: " << exception.getContext() << ": " << exception.getMessage() << std::endl;
        exit (1);
    }
    catch (const boost::exception &e)
    {
        std::clog << "Error: boost::exception: " << boost::diagnostic_information(e) << std::endl;
        exit (1);
    }
    catch (const std::runtime_error &e)
    {
        std::clog << "runtime error: " << e.what() << std::endl;
        exit(2);
    }
    catch (const std::logic_error &e)
    {
        std::clog << "logic error: " << e.what() << std::endl;
        exit(3);
    }
}

} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_PROGRAM_HH
