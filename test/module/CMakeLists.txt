# SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Module support needs dynamic dependencies, which are only available with CMake >= 3.28 and Ninja/Visual Studio.
cmake_minimum_required (VERSION 3.28...3.31)
project (hibf_test_module CXX)

# These tests `import seqan.hibf;`, so the module has to be part of the library.
set (HIBF_MODULE
     ON
     CACHE BOOL "Additionally build the `seqan.hibf` C++20 module." FORCE)

include (../hibf-test.cmake)
include (GoogleTest OPTIONAL)

# With HIBF_POST_INSTALL_TEST, the setting above is ignored and `seqan::hibf` comes from `find_package`. Fail with a
# clear message instead of "no such module" if that installation was built without the module.
get_target_property (hibf_definitions seqan::hibf INTERFACE_COMPILE_DEFINITIONS)
if (NOT "HIBF_HAS_MODULE=1" IN_LIST hibf_definitions)
    message (FATAL_ERROR "  The hibf being tested was built without the module. Reinstall it with -DHIBF_MODULE=ON.")
endif ()

CPMGetPackage (googletest)

macro (hibf_module_test module_test_cpp)
    file (RELATIVE_PATH module_test "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${module_test_cpp}")
    hibf_test_component (target "${module_test}" TARGET_NAME)
    hibf_test_component (test_name "${module_test}" TEST_NAME)

    add_executable (${target} ${module_test_cpp})
    # This is `hibf::test::unit` spelled out. `seqan::hibf` has to be a *direct* link dependency: with
    # HIBF_POST_INSTALL_TEST it is an IMPORTED target, and CMake (4.0) does not propagate the module usage requirements
    # of an imported target through an intermediate INTERFACE library -- `import seqan.hibf;` then fails with "no such
    # module". Going through `hibf::test::unit` as well would put libhibf.a on the link line twice, which Clang's
    # ThinLTO rejects.
    target_link_libraries (${target} GTest::gtest_main seqan::hibf)
    target_include_directories (${target} PRIVATE "${HIBF_TEST_INCLUDE_DIR}")
    add_test (NAME "${test_name}" COMMAND ${target})

    unset (module_test)
    unset (target)
    unset (test_name)
endmacro ()

hibf_module_test (hibf_module_test.cpp)
hibf_module_test (mixed_mode_test.cpp)
