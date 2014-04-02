################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2012 Illumina, Inc.
##
## This software is provided under the terms and conditions of the
## Illumina Open Source Software License 1.
##
## You should have received a copy of the Illumina Open Source
## Software License 1 along with this program. If not, see
## <https://github.com/sequencing/licenses/>.
##
## The distribution includes the code libraries listed below in the
## 'redist' sub-directory. These are distributed according to the
## licensing terms governing each library.
##
################################################################################
##
## file cppunit.cmake
##
## Configuration file for the cppunit subfolders
##
## author Come Raczy
##
################################################################################

# unforce static linking. On many systems cppunitTest needs to be a dynamic executable
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")

##
## the location of the cppunit shared libraries will be needed for
## LD_LIBRARY_PATH
##

if (HAVE_CPPUNIT)
    get_filename_component(CPPUNIT_LOCATION ${HAVE_CPPUNIT} PATH)
    include_directories(${CPPUNIT_INCLUDE_DIR})
else (HAVE_CPPUNIT)
    message(FATAL_ERROR "cppunit not found")
endif (HAVE_CPPUNIT)

##
## find all the source files
##

file (GLOB iSAAC_TEST_SOURCE_LIST "*.cpp")

##
## create the targets to build the tests
##

set(iSAAC_CPPUNIT_TEST_NAME cppunitTest)
add_executable(${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME} ${iSAAC_TEST_SOURCE_LIST})
set_target_properties(${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME} PROPERTIES OUTPUT_NAME ${iSAAC_CPPUNIT_TEST_NAME})
add_test(cppunit_${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME} "${CMAKE_CURRENT_BINARY_DIR}/${iSAAC_CPPUNIT_TEST_NAME}")

include_directories   (${iSAAC_CXX_ALL_INCLUDES} ${iSAAC_OPT_INC} "${CMAKE_SOURCE_DIR}/c++/unittest")

set(iSAAC_LINK_LIBRARIES "-lpthread")
if    (HAVE_ZLIB)
    set(iSAAC_LINK_LIBRARIES "${iSAAC_LINK_LIBRARIES} -lz")
endif (HAVE_ZLIB)
if    (NOT iSAAC_FORCE_STATIC_LINK)
    set(iSAAC_LINK_LIBRARIES "${iSAAC_LINK_LIBRARIES} -ldl")
endif (NOT iSAAC_FORCE_STATIC_LINK)

target_link_libraries (${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME}
                       isaac_${iSAAC_LIB_DIR} ${iSAAC_AVAILABLE_LIBRARIES} 
                       isaac_cppunit ${Boost_LIBRARIES}
                       ${iSAAC_LINK_LIBRARIES} ${iSAAC_DEP_LIB} ${CPPUNIT_LIBRARY})

##
## Run some sanity check on the source file
##
foreach(iSAAC_CPPUNIT_SOURCE_FILE ${iSAAC_TEST_SOURCE_LIST})
    get_filename_component(FILE_NAME ${iSAAC_CPPUNIT_SOURCE_FILE} NAME)
    set(iSAAC_CPPUNIT_BINARY_FILE "${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}")
    add_custom_command(TARGET ${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME} 
                       PRE_BUILD
                       COMMAND ${CMAKE_SOURCE_DIR}/c++/unittest/check-source.sh ARGS ${iSAAC_CPPUNIT_BINARY_FILE}.checked ${iSAAC_CPPUNIT_SOURCE_FILE}
                       COMMENT "Sanity check on ${iSAAC_CPPUNIT_SOURCE_FILE}")
endforeach(iSAAC_CPPUNIT_SOURCE_FILE)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/RegistryNames.txt 
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/RegistryNames.txt ${CMAKE_CURRENT_BINARY_DIR}
		   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/RegistryNames.txt
                   COMMENT "Copying RegistryNames.txt for ${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME}")

##
## create the targets to run the tests
##
add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${iSAAC_CPPUNIT_TEST_NAME}.passed 
		   COMMAND export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:${CPPUNIT_LOCATION} && ./${iSAAC_CPPUNIT_TEST_NAME}
	           COMMAND touch ${iSAAC_CPPUNIT_TEST_NAME}.passed  
		   DEPENDS ${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME} ${CMAKE_CURRENT_BINARY_DIR}/RegistryNames.txt
		   COMMENT "Running unit tests ${iSAAC_UNIQUE_PREFIX}${iSAAC_CPPUNIT_TEST_NAME}")
add_custom_target(${iSAAC_UNIQUE_PREFIX}passed ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${iSAAC_CPPUNIT_TEST_NAME}.passed)

##
## Copy the data directory from the source tree if available
##

find_path(${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR data PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
if (${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR)
message (STATUS "Adding the data subdirectory for the cppunits under ${iSAAC_LIB_DIR}")
    add_subdirectory (data)
    add_dependencies(${iSAAC_UNIQUE_PREFIX}passed ${iSAAC_UNIQUE_PREFIX}data)
endif (${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR)
