// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Checks that headers and `import hibf;` can be mixed within one translation unit.
 * \details
 * Because src/hibf.cppm includes the headers in its global module fragment, an imported entity and an included entity
 * are the *same* entity. A project can therefore migrate to the module file by file instead of all at once.
 */

#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/misc/bit_vector.hpp>

import hibf;

TEST(mixed_mode_test, same_entity)
{
    // Named through the header, and through the import: one and the same type.
    static_assert(
        std::same_as<decltype(seqan::hibf::interleaved_bloom_filter{}), seqan::hibf::interleaved_bloom_filter>);

    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u}, seqan::hibf::bin_size{1024u}};
    ibf.emplace(126u, seqan::hibf::bin_index{0u});

    // `containment_agent()` is defined in the header; `bulk_contains` is compiled into libhibf.a. Both are reached
    // from a translation unit that also imported the module.
    auto agent = ibf.containment_agent();
    seqan::hibf::bit_vector const & result = agent.bulk_contains(126u);
    EXPECT_TRUE(result[0u]);
}

TEST(mixed_mode_test, entities_not_exported_are_still_usable_via_the_header)
{
    // seqan::stl is not part of the module's export list, but the header still provides it.
    std::vector<size_t> const values{4u, 5u, 6u};
    size_t sum{};
    for (auto const [index, value] : seqan::stl::views::enumerate(values))
        sum += index * value;

    EXPECT_EQ(sum, 0u * 4u + 1u * 5u + 2u * 6u);
}
