# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library-template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

# Adds a macro benchmark target and a test which executes that macro benchmark
#
# Example:
#
# library_template_macro_benchmark(
#     edit_distance_macro_benchmark.cpp
#     REPETITIONS 3 # optional, default value is 3
#     META_DATA # optional
#         foo bar # user-specified data
#         library_template_version 3.0.1 # automatically set
#         library_template_git_version a1829272d7bafea48580af6eaac2bfbb678c3861 # automatically set
#         sdsl_version 3.0.0 # automatically set
#         sdsl_git_version 9a0d5676fd09fb8b52af214eca2d5809c9a32dbe # automatically set
# )
#
# This produces an executable target with the name `edit_distance_macro_benchmark`.
#
# When calling ctest, the executable `edit_distance_macro_benchmark` will be executed and a json file (i.e.
# `edit_distance_macro_benchmark.json`) with the benchmark results will be generated.
#
# Options:
#
# * REPETITIONS: [default value: 3]
#
#   Repeats the benchmark `REPETITIONS`-many times.
#
#   Under the hood executes the target with `--benchmark_repetitions = ${REPETITIONS}`.
#
# * META_DATA:
#
#   Adds the meta data to the generated json file.
function (library_template_macro_benchmark macro_benchmark_cpp)
    cmake_parse_arguments (
        "MACRO_BENCHMARK" #
        "" #
        "REPETITIONS" #
        "META_DATA" #
        ${ARGN})

    if (NOT MACRO_BENCHMARK_REPETITIONS)
        set (MACRO_BENCHMARK_REPETITIONS "3")
    endif ()

    file (RELATIVE_PATH macro_benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${macro_benchmark_cpp}")

    include (library_template_test_component)
    library_template_test_component (target "${macro_benchmark}" TARGET_NAME)
    library_template_test_component (test_name "${macro_benchmark}" TEST_NAME)

    add_executable (${target} ${macro_benchmark_cpp})
    target_link_libraries (${target} library_template::test::performance)
    add_test (NAME "${test_name}"
              COMMAND ${target} #
                      --benchmark_repetitions=${MACRO_BENCHMARK_REPETITIONS} #
                      --benchmark_min_time=0 #
                      --benchmark_enable_random_interleaving=true #
                      --benchmark_out=$<TARGET_FILE_NAME:${target}>.json #
                      --benchmark_out_format=json)
endfunction ()
