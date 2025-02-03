# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)

# cmake ../test/iwyu/ -DCMAKE_CXX_COMPILER="clang++-17" -DHIBF_IWYU_DIR=...
# make -j2 2> out.txt
# python3 fix_includes.py --nosafe_headers --comments --noreorder --ignore_re ".*/include/hibf/contrib/.*" < out.txt

project (hibf_iwyu CXX)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message (FATAL_ERROR "Please use clang for iwyu")
endif ()

find_program (IWYU "include-what-you-use" REQUIRED PATHS "${HIBF_IWYU_DIR}")

set (CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--max_line_length=120")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--comment_style=short")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--update_comments")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--no_fwd_decls")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--cxx17ns")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--transitive_includes_only")
# This is how iwyu glob works: https://godbolt.org/z/EPhhvjG7n
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--check_also=**/include/hibf/**.hpp")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--no_default_mappings")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--mapping_file=/home/infri/iwyu/mappings/foo.imp")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "--stdlib=libc++")

option (INSTALL_HIBF "" OFF)
add_subdirectory (../unit/ unit)
add_subdirectory (../performance/ performance)
add_subdirectory (../snippet/ snippet)
