// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, TestInfo, TestPartResult, TEST, CmpHelperEQ

#include <algorithm>  // for __generate, generate
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t, uint8_t
#include <filesystem> // for path
#include <fstream>    // for basic_ofstream, basic_ifstream, basic_ios, ofstream, ifs...
#include <random>     // for uniform_int_distribution, mt19937_64
#include <ranges>     // for __fn, iota, iota_view, views
#include <stdexcept>  // for runtime_error, invalid_argument
#include <string>     // for basic_string
#include <vector>     // for vector

#include <hibf/contrib/robin_hood.hpp>              // for unordered_flat_set
#include <hibf/contrib/std/chunk_view.hpp>          // for chunk, chunk_fn, chunk_view
#include <hibf/contrib/std/detail/adaptor_base.hpp> // for operator|
#include <hibf/misc/divide_and_ceil.hpp>            // for divide_and_ceil
#include <hibf/sketch/hyperloglog.hpp>              // for hyperloglog
#include <hibf/test/sandboxed_path.hpp>             // for operator/, sandboxed_path
#include <hibf/test/tmp_directory.hpp>              // for tmp_directory

TEST(hyperloglog, bit_widths)
{
    for (uint8_t i : std::views::iota(0u, 5u))
        EXPECT_THROW(seqan::hibf::sketch::hyperloglog{i}, std::invalid_argument);

    EXPECT_NO_THROW(seqan::hibf::sketch::hyperloglog{5u});
}

TEST(hyperloglog, initialization)
{
    size_t const bits = 6;
    size_t const size = 1 << bits;

    seqan::hibf::sketch::hyperloglog sketch(bits);

    EXPECT_EQ(sketch.data_size(), size);

    // No elements were inserted, so the small values correction should be used.
    // Since there are only zeros in the register, the correction formula should be:
    // m * log(m / #zeros) = m * log(m/m) = m * log(1) = 0
    EXPECT_EQ(sketch.estimate(), 0.0);
}

TEST(hyperloglog, add_and_estimate_small)
{
    seqan::hibf::sketch::hyperloglog sketch{5u};

    // first 4 bits of hash: 0000, rank: 3
    sketch.add(255881241332063154ULL);
    // first 4 bits of hash: 1011, rank: 2
    sketch.add(13132817195163223578ULL);
    // first 4 bits of hash: 0100, rank: 2
    sketch.add(5120631300412165844ULL);
    // first 4 bits of hash: 1110, rank: 1
    sketch.add(16862690419523416066ULL);
    // first 4 bits of hash: 0000, rank: 3
    sketch.add(148518882728022940ULL);
    // first 4 bits of hash: 1101, rank: 1
    sketch.add(15892358469365346306ULL);
    // first 4 bits of hash: 1001, rank: 2
    sketch.add(10885195586503739779ULL);
    // first 4 bits of hash: 1000, rank: 2
    sketch.add(9563173945158404745ULL);

    // estimate = alpha * m  * m  / sum(2^(-data[j]))
    //          = 0.697 * 32 * 32 / (89/8)
    //          = 0.697 * 32 * 32 / 11.125 = 64.155...

    // this still is in the range of small value corrections (< 2.5 * 32 = 80)
    // m * ln(m / #zeros) = 32 * ln(32/25) = 7.899522493... (with calculator)

    EXPECT_NEAR(sketch.estimate(), 7.899522493, 0.0000001);
}

TEST(hyperloglog, reset)
{
    // Same as add_and_estimate_small
    seqan::hibf::sketch::hyperloglog sketch{};

    sketch.add(255881241332063154ULL);
    sketch.add(13132817195163223578ULL);
    sketch.add(5120631300412165844ULL);
    sketch.add(16862690419523416066ULL);
    sketch.add(148518882728022940ULL);
    sketch.add(15892358469365346306ULL);
    sketch.add(10885195586503739779ULL);
    sketch.add(9563173945158404745ULL);

    EXPECT_NEAR(sketch.estimate(), 7.899522493, 0.0000001);

    // Actual reset test
    sketch.reset();
    EXPECT_EQ(sketch.estimate(), 0.0);
}

std::vector<uint64_t> const input_values = []()
{
    auto generator = []()
    {
        std::uniform_int_distribution<uint64_t> distribution{};
        std::mt19937_64 engine{0u};
        return distribution(engine);
    };

    std::vector<uint64_t> result(1500);
    std::ranges::generate(result, generator);
    return result;
}();

TEST(hyperloglog, add_and_estimate_large)
{
    seqan::hibf::sketch::hyperloglog sketch{};

    robin_hood::unordered_flat_set<uint64_t> control;

    for (uint64_t value : input_values)
    {
        control.insert(value);
        sketch.add(value);
    }

    // the estimate is greater than 2.5 * 16, therefore this is the raw estimate
    // 1.04 / sqrt(m) is the usual relative error bound proved in the paper
    EXPECT_NEAR(sketch.estimate(), control.size(), control.size() * 1.04 / 4);
}

TEST(hyperloglog, merge)
{
    size_t const chunks{10u};
    size_t const chunk_size = seqan::hibf::divide_and_ceil(input_values.size(), chunks);

    seqan::hibf::sketch::hyperloglog full_sketch{};
    seqan::hibf::sketch::hyperloglog merge_sketch{};
    seqan::hibf::sketch::hyperloglog merge_estimate_sketch{};

    std::vector<seqan::hibf::sketch::hyperloglog> partial_sketches;

    // put every chunk into the full full_sketch
    // and add a disjointed sketch for every chunk to partial_sketches
    for (auto && chunk : input_values | seqan::stl::views::chunk(chunk_size))
    {
        partial_sketches.emplace_back();

        for (uint64_t value : chunk)
        {
            partial_sketches.back().add(value);
            full_sketch.add(value);
        }
    }

    // merge all partial sketches into merge_sketch
    for (auto & partial_sketch : partial_sketches)
    {
        merge_sketch.merge(partial_sketch);
        merge_estimate_sketch.merge_and_estimate(partial_sketch);
    }

    // now full_sketch and merged_sketch should be equal
    EXPECT_EQ(full_sketch.estimate(), merge_sketch.estimate());
    EXPECT_EQ(full_sketch.estimate(), merge_estimate_sketch.estimate());
}

TEST(hyperloglog, failed_store)
{
    seqan::hibf::sketch::hyperloglog sketch{};
    std::ofstream ostrm{"hibf_non_existent_outputfile"};
    ostrm.close();

    try
    {
        sketch.store(ostrm);
        FAIL();
    }
    catch (std::runtime_error const & exception)
    {
        EXPECT_STREQ(exception.what(), "[HyperLogLog] Failed to store a HyperLogLog sketch to a file.");
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(hyperloglog, failed_load)
{
    seqan::hibf::test::tmp_directory tmp_dir{};
    std::filesystem::path file_name{tmp_dir.path() / "sketch.hll"};
    {
        uint8_t b{5u};
        std::ofstream ostrm{file_name};
        ostrm.write((char *)&b, sizeof(b));
    }
    seqan::hibf::sketch::hyperloglog sketch{};
    std::ifstream istrm{file_name};

    try
    {
        sketch.load(istrm);
        FAIL();
    }
    catch (std::runtime_error const & exception)
    {
        EXPECT_STREQ(exception.what(), "[HyperLogLog] Failed to load a HyperLogLog sketch from a file: I/O error.");
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(hyperloglog, failed_load_bit_width)
{
    seqan::hibf::test::tmp_directory tmp_dir{};
    std::filesystem::path file_name{tmp_dir.path() / "wrong.hll"};
    {
        uint8_t b{4u};
        std::ofstream ostrm{file_name};
        ostrm.write((char *)&b, sizeof(b));
    }
    seqan::hibf::sketch::hyperloglog sketch{};
    std::ifstream istrm{file_name};

    try
    {
        sketch.load(istrm);
        FAIL();
    }
    catch (std::runtime_error const & exception)
    {
        EXPECT_STREQ(exception.what(),
                     "[HyperLogLog] Failed to load a HyperLogLog sketch from a file: Invalid bit_width.");
    }
    catch (...)
    {
        FAIL();
    }
}

TEST(hyperloglog, store_and_load)
{
    seqan::hibf::sketch::hyperloglog store_sketch{};
    seqan::hibf::sketch::hyperloglog load_sketch{};

    for (uint64_t value : input_values)
        store_sketch.add(value);

    // create temp file
    seqan::hibf::test::tmp_directory tmp_dir{};
    std::filesystem::path store_filename{tmp_dir.path() / "dump.hll"};

    // store sketch
    std::ofstream ostrm(store_filename, std::ios::binary);
    store_sketch.store(ostrm);

    // load sketch
    std::ifstream istrm(store_filename, std::ios::binary);
    load_sketch.load(istrm);

    // now store_sketch and load_sketch should be equal
    EXPECT_EQ(store_sketch.estimate(), load_sketch.estimate());
}
