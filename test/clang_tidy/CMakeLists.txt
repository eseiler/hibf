# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.10)

project (hibf_tidy NONE)

find_program (CLANG_TIDY "clang-tidy-18" REQUIRED)
set (CMAKE_CXX_CLANG_TIDY
     "${CLANG_TIDY};-checks=bugprone*,performance*,-bugprone-easily-swappable-parameters,-bugprone-exception-escape,-bugprone-narrowing-conversions"
)

add_subdirectory (../unit/ unit)
