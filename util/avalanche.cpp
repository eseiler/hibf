// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

// Measures how well seqan::hibf::interleaved_bloom_filter::hash_and_fit satisfies the strict avalanche criterion
// (SAC): for every input bit i and output bit j, flipping bit i of the input should flip bit j of the output with
// probability ~0.5. `hash_and_fit` is tested exactly as implemented (bin_count == 1, so the interleaving is a no-op
// and this is equivalent to a plain Bloom filter).

#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include <sharg/parser.hpp>

#include <hibf/interleaved_bloom_filter.hpp>

#include "inspect/inspector.hpp"

inline constexpr size_t bits = std::numeric_limits<size_t>::digits; // 64

struct config
{
    size_t bin_size{8192};
    size_t hash_index{0};
    bool all_seeds{false};
    size_t samples{200'000};
    uint64_t seed{};
    std::string matrix_output{};
};

void init_parser(sharg::parser & parser, config & cfg)
{
    parser.add_option(
        cfg.bin_size,
        sharg::config{.short_id = '\0',
                      .long_id = "bin-size",
                      .description = "Size of a single bin in bits, i.e., the range that "
                                     "hash_and_fit maps into. This is the `bin_size` you would pass "
                                     "to construct a 1-bin (== plain Bloom filter) IBF.",
                      .validator = sharg::arithmetic_range_validator{size_t{1u}, //
                                                                     std::numeric_limits<size_t>::max()}});
    parser.add_option(cfg.hash_index,
                      sharg::config{.short_id = '\0',
                                    .long_id = "hash",
                                    .description = "Which of the 5 internal hash seeds to test (0-4). "
                                                   "Ignored if --all-seeds is given.",
                                    .validator = sharg::arithmetic_range_validator{size_t{0u}, size_t{4u}}});
    parser.add_flag(cfg.all_seeds,
                    sharg::config{.short_id = '\0',
                                  .long_id = "all-seeds",
                                  .description = "Test all 5 hash seeds, one after another."});
    parser.add_option(
        cfg.samples,
        sharg::config{.short_id = '\0',
                      .long_id = "samples",
                      .description = "Number of random inputs tested per input bit. "
                                     "The standard error of each reported probability is "
                                     "roughly 0.5 / sqrt(samples).",
                      .validator = sharg::arithmetic_range_validator{size_t{1u}, //
                                                                     std::numeric_limits<size_t>::max()}});
    parser.add_option(
        cfg.seed,
        sharg::config{.short_id = '\0', .long_id = "seed", .description = "Seed for the random number generator."});
    parser.add_option(cfg.matrix_output,
                      sharg::config{.short_id = '\0',
                                    .long_id = "matrix-output",
                                    .description = "Optional path prefix to additionally write the full "
                                                   "64x64 input-bit x output-bit flip-probability matrix as CSV. "
                                                   "One file `<prefix>.hash<seed_index>.csv` is written per tested "
                                                   "seed."});
}

// P(output bit j flips | input bit i flipped), for all 64x64 (i, j) pairs, estimated from `samples` random inputs.
std::vector<std::array<double, bits>> avalanche_matrix(seqan::hibf::interleaved_bloom_filter const & ibf,
                                                       size_t const hash_seed_value,
                                                       size_t const samples,
                                                       uint64_t const rng_seed)
{
    using inspector = seqan::hibf::inspector;

    std::vector<std::array<size_t, bits>> flips(bits); // flips[i][j] = number of samples where bit j flipped
    for (auto & row : flips)
        row.fill(0u);

    std::mt19937_64 gen{rng_seed};
    std::uniform_int_distribution<size_t> distrib{0u, std::numeric_limits<size_t>::max()};

    for (size_t sample = 0u; sample < samples; ++sample)
    {
        size_t const h = distrib(gen);
        size_t const baseline = inspector::hash_and_fit(ibf, h, hash_seed_value);

        for (size_t i = 0u; i < bits; ++i)
        {
            size_t const flipped = inspector::hash_and_fit(ibf, h ^ (size_t{1u} << i), hash_seed_value);
            size_t diff = baseline ^ flipped;

            while (diff != 0u)
            {
                size_t const j = std::countr_zero(diff);
                ++flips[i][j];
                diff &= diff - 1u; // clear lowest set bit
            }
        }
    }

    std::vector<std::array<double, bits>> probabilities(bits);
    for (size_t i = 0u; i < bits; ++i)
        for (size_t j = 0u; j < bits; ++j)
            probabilities[i][j] = static_cast<double>(flips[i][j]) / static_cast<double>(samples);

    return probabilities;
}

void write_matrix_csv(std::vector<std::array<double, bits>> const & matrix, std::string const & path)
{
    std::ofstream out{path};

    out << "input_bit\\output_bit";
    for (size_t j = 0u; j < bits; ++j)
        out << ',' << j;
    out << '\n';

    for (size_t i = 0u; i < bits; ++i)
    {
        out << i;
        for (size_t j = 0u; j < bits; ++j)
            out << ',' << matrix[i][j];
        out << '\n';
    }
}

void report(std::vector<std::array<double, bits>> const & matrix, size_t const bin_size, size_t const technical_bins)
{
    // hash_and_fit returns `position * technical_bins`, where `position` is in [0, bin_size) and technical_bins is a
    // power of two (64, for the 1-bin case this tool always uses). So the result is exactly `position` shifted left
    // by log2(technical_bins) bits: the low `shift` bits are structurally always 0, then `log2(bin_size)` bits carry
    // the actual entropy, then everything above that is 0 too (the output is bounded to a range smaller than 2^64).
    // Those always-0 bits trivially have a 0% flip probability; that's expected and not a weakness of the mixing
    // step, so they are excluded from the summary below.
    size_t const shift = static_cast<size_t>(std::countr_zero(technical_bins));
    size_t const live_output_bits = std::bit_width(bin_size == 0u ? 0u : bin_size - 1u);

    double sum_abs_deviation = 0.0;
    double max_abs_deviation = 0.0;
    size_t max_i = 0u;
    size_t max_j = 0u;
    size_t considered = 0u;

    for (size_t i = 0u; i < bits; ++i)
    {
        for (size_t j = shift; j < shift + live_output_bits; ++j)
        {
            double const deviation = std::abs(matrix[i][j] - 0.5);
            sum_abs_deviation += deviation;
            ++considered;
            if (deviation > max_abs_deviation)
            {
                max_abs_deviation = deviation;
                max_i = i;
                max_j = j;
            }
        }
    }

    double const mean_abs_deviation = considered == 0u ? 0.0 : sum_abs_deviation / static_cast<double>(considered);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Live output bits:                        [" << shift << ", " << (shift + live_output_bits - 1u)
              << "] (" << live_output_bits << " / " << bits << ")\n";
    std::cout << "Mean |P(flip) - 0.5| over live bits:     " << mean_abs_deviation
              << "  (0 = perfect avalanche, ideal is 0)\n";
    std::cout << "Max  |P(flip) - 0.5| over live bits:     " << max_abs_deviation << "  (input bit " << max_i
              << " -> output bit " << max_j << ", P(flip) = " << matrix[max_i][max_j] << ")\n";

    std::cout << "\nPer input-bit average P(flip) over live output bits:\n";
    std::cout << std::setprecision(4);
    for (size_t i = 0u; i < bits; ++i)
    {
        double sum = 0.0;
        for (size_t j = shift; j < shift + live_output_bits; ++j)
            sum += matrix[i][j];
        double const avg = live_output_bits == 0u ? 0.0 : sum / static_cast<double>(live_output_bits);

        std::cout << "  bit " << std::setw(2) << i << ": " << avg;
        if (i % 4u == 3u)
            std::cout << '\n';
        else
            std::cout << "    ";
    }
    if (bits % 4u != 0u)
        std::cout << '\n';
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"avalanche", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.short_description = "Tests the avalanche effect (strict avalanche criterion) of "
                                    "interleaved_bloom_filter::hash_and_fit.";
    config cfg{};
    init_parser(parser, cfg);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        return -1;
    }

    if (!parser.is_option_set("seed"))
    {
        std::random_device rd;
        cfg.seed = (static_cast<uint64_t>(rd()) << 32) | rd();
    }

    seqan::hibf::interleaved_bloom_filter const ibf{seqan::hibf::bin_count{1u},
                                                    seqan::hibf::bin_size{cfg.bin_size},
                                                    seqan::hibf::hash_function_count{1u}};

    std::vector<size_t> hash_indices =
        cfg.all_seeds ? std::vector<size_t>{0u, 1u, 2u, 3u, 4u} : std::vector<size_t>{cfg.hash_index};

    std::cout << "bin_size: " << cfg.bin_size << '\n';
    std::cout << "samples per input bit: " << cfg.samples << '\n';
    std::cout << "rng seed: " << cfg.seed << "\n\n";

    auto const & hash_seeds = seqan::hibf::inspector::hash_seeds(ibf);

    for (size_t const hash_index : hash_indices)
    {
        std::cout << "=== hash seed index " << hash_index << " (seed = " << hash_seeds[hash_index] << ") ===\n";
        auto const matrix = avalanche_matrix(ibf, hash_seeds[hash_index], cfg.samples, cfg.seed);
        report(matrix, cfg.bin_size, seqan::hibf::inspector::technical_bins(ibf));

        if (!cfg.matrix_output.empty())
        {
            std::string const path = cfg.matrix_output + ".hash" + std::to_string(hash_index) + ".csv";
            write_matrix_csv(matrix, path);
            std::cout << "\nFull matrix written to " << path << '\n';
        }
        std::cout << '\n';
    }
}
