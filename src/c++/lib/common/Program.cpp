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
 ** \file Program.cpp
 **
 ** Implementation of the skeleton of all c++ programs.
 **
 ** \author Come Raczy
 **/

#include <sys/ioctl.h>

#include "common/Program.hh"

namespace isaac
{
namespace common
{

winsize ioctlSTDERR_TIOCGWINSZ()
{
    winsize ret = {0,0,0,0};
    if (-1 == ioctl(STDERR_FILENO, TIOCGWINSZ, &ret))
    {
        ret.ws_col = bpo::options_description::m_default_line_length;
    }
    return ret;
}

Options::Options() :
    namedOptions_("Command line options", ioctlSTDERR_TIOCGWINSZ().ws_col, ioctlSTDERR_TIOCGWINSZ().ws_col - 50)
{
    namedOptions_.add_options()("help,h", "produce help message and exit");
    namedOptions_.add_options()("version,v", "print program version information");
}

Options::Action Options::parse(int argc, const char * const argv[])
{
    try
    {
        bpo::options_description allOptions("Allowed options");
        allOptions.add(namedOptions_).add(unnamedOptions_);
        bpo::variables_map vm;
        bpo::store(
                bpo::command_line_parser(argc, argv). options(allOptions).positional(
                        positionalOptions_).run(), vm);
        bpo::notify(vm);
        postProcess(vm);
        if (vm.count("help"))
        {
            return HELP;
        }
        else if (vm.count("version"))
        {
            return VERSION;
        }
    }
    catch (const boost::program_options::multiple_values &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
        return ABORT;
    }
    catch (const boost::program_options::multiple_occurrences &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
        return ABORT;
    }
    catch (const boost::program_options::required_option &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
        return ABORT;
    }
    catch (const std::exception &e)
    {
        std::clog << usage() << std::endl;
        std::clog << "Failed to parse the options: " << e.what() << std::endl;
        return ABORT;
    }
    return RUN;
}

std::string Options::usage() const
{
    std::ostringstream os;
    os << this->usagePrefix() << std::endl << std::endl;
    os << namedOptions_ << std::endl;
    return os.str();
}

} // namespace common
} // namespace isaac
