// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>
#include <random>

#include <hibf/contrib/robin_hood.hpp>

#include <sharg/parser.hpp>

static constexpr size_t sequence_length{10'000'000ULL};

struct config
{
    // Set by user
    size_t kmer_size{12};

    // Internal
    size_t max_kmer_value{};

    void generate_histogram(robin_hood::unordered_map<size_t, size_t> const & counts) const
    {
        robin_hood::unordered_map<size_t, size_t> histogram{};
        for (auto const & [key, value] : counts)
            ++histogram[value];

        size_t total_kmers{};
        size_t const max_count = std::ranges::max_element(histogram)->first;
        for (size_t i = max_count; i > 0u; --i)
        {
            std::cout << i << '\t' << histogram[i] << '\n';
            total_kmers += histogram[i];
        }
        std::cout << "0\t" << (max_kmer_value - total_kmers + 1u) << '\n';
    }
};

void init_parser(sharg::parser & parser, config & cfg)
{
    parser.add_option(cfg.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size.",
                                    .validator = sharg::arithmetic_range_validator{1, 31}});
}

void random_kmers(config const & cfg)
{
    std::mt19937_64 gen(42u);
    std::uniform_int_distribution<uint64_t> distrib(0ULL, cfg.max_kmer_value);
    robin_hood::unordered_map<size_t, size_t> counts{};
    for (size_t i = 0; i < sequence_length; ++i)
        ++counts[distrib(gen)];

    cfg.generate_histogram(counts);
}

void generate_via_int(config const & cfg)
{
    std::mt19937_64 gen(42u);
    std::uniform_int_distribution<uint64_t> distrib(0ULL, 3ULL);

    size_t kmer_value{};
    size_t const kmer_mask{cfg.max_kmer_value};

    auto next = [&]()
    {
        kmer_value <<= 2;
        kmer_value |= distrib(gen);
        kmer_value &= kmer_mask;
    };

    for (size_t i = 0; i < cfg.kmer_size - 1u; ++i)
        next();

    robin_hood::unordered_map<size_t, size_t> counts{};
    for (size_t i = 0; i < sequence_length; ++i)
    {
        next();
        ++counts[kmer_value];
    }

    cfg.generate_histogram(counts);
}

void generate_via_double(config const & cfg)
{
    std::mt19937_64 gen(42u);
    std::uniform_real_distribution<double> distrib(0.0, 1.0);

    auto get_letter = [&]() -> size_t
    {
        double const value = distrib(gen);
        if (value < 0.25)
            return 0u;
        else if (value < 0.5)
            return 1u;
        else if (value < 0.75)
            return 2u;
        else
            return 3u;
    };

    size_t kmer_value{};
    size_t const kmer_mask{cfg.max_kmer_value};

    auto next = [&]()
    {
        kmer_value <<= 2;
        kmer_value |= get_letter();
        kmer_value &= kmer_mask;
    };

    for (size_t i = 0; i < cfg.kmer_size - 1u; ++i)
        next();

    robin_hood::unordered_map<size_t, size_t> counts{};
    for (size_t i = 0; i < sequence_length; ++i)
    {
        next();
        ++counts[kmer_value];
    }

    cfg.generate_histogram(counts);
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"fpr_correction_check", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.short_description = "Inserts a given amount of k-mers into an IBF and queries all possible k-mers. "
                                    "Reports the resulting FPR for both single and split bins.";
    config cfg{};
    init_parser(parser, cfg);
    parser.parse();

    cfg.max_kmer_value = (1ULL << (2 * cfg.kmer_size)) - 1u;
    random_kmers(cfg);
    std::cout << '\n';
    generate_via_double(cfg);
    std::cout << '\n';
    generate_via_int(cfg);
}
