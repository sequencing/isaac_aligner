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
 ** \file RegistryName.hh
 **
 ** Management of the registry names for the cppunit tests.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_UNIT_TEST_REGISTRY_NAME
#define iSAAC_UNIT_TEST_REGISTRY_NAME

#include <stdexcept>
#include <string>
#include <vector>

const std::vector<std::string> &getRegistryNameList();
std::string registryName(const std::string &name) throw (std::invalid_argument);

#endif // #ifndef iSAAC_UNIT_TEST_REGISTRY_NAME
