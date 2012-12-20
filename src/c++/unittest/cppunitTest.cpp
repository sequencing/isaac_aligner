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
 ** \file cppunitTest.cpp
 **
 ** Main program used for all the cppunit tests
 **
 ** \author Come Raczy
 **/

#include <sstream>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "RegistryName.hh"

int main()
{

  CppUnit::TextUi::TestRunner runner;
  // First add the tests from the named registries in the right order
  // To add/remove/modify a registry name, or to change its sequence
  // number, edit RegistryName.cpp
  for (std::vector<std::string>::const_iterator name = getRegistryNameList().begin();
       getRegistryNameList().end() != name; ++name)
  {
    CppUnit::Test *namedSuite = CppUnit::TestFactoryRegistry::getRegistry(*name).makeTest();
    runner.addTest(namedSuite);
  }

  // Add the top level (unnamed) suite from the list of tests to run
  CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
  runner.addTest( suite );

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(),
                                                       std::cerr ) );
  // Run the tests.
  bool wasSucessful = runner.run();

  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}
