# SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (hibf_test)

option (INSTALL_HIBF "" OFF)
option (HIBF_MODULE "Additionally build the `hibf` C++20 module." OFF)
enable_testing ()
add_subdirectory (../unit/ unit)
add_subdirectory (../performance/ performance)
add_subdirectory (../snippet/ snippet)

# Only meaningful with the module built, and needs CMake >= 3.28 plus Ninja/Visual Studio.
if (HIBF_MODULE)
    add_subdirectory (../module/ module)
endif ()
