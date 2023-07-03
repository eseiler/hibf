cmake_minimum_required (VERSION 3.10)

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
# This is how iwyu glob works: https://godbolt.org/z/EPhhvjG7n
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--check_also=**/include/hibf/**.hpp")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--keep=**hibf/platform.hpp")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--keep=**hibf/version.hpp")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--no_default_mappings")
list (APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu;--mapping_file=${CMAKE_CURRENT_LIST_DIR}/mappings/iwyu.imp")

add_subdirectory (../unit/ unit)

# Make ./include/hibf and ./test/include SYSTEM header, such that they are included with -isystem
# IWYU will always use quoted includes ("a.hpp") for -I includes, and angular includes (<a.hpp>) for -isystem includes
get_target_property (HIBF_INCLUDES hibf INCLUDE_DIRECTORIES)
target_include_directories (hibf SYSTEM PUBLIC "${HIBF_INCLUDES}")
get_target_property (HIBF_TEST_INCLUDES hibf_test INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories (hibf_test SYSTEM INTERFACE "${HIBF_TEST_INCLUDES}")