// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t

#include <hibf/build/build_data.hpp>         // for build_data
#include <hibf/contrib/robin_hood.hpp>       // for unordered_flat_set
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter
#include <hibf/layout/layout.hpp>            // for layout
#include <hibf/misc/timer.hpp>               // for concurrent_timer

namespace seqan::hibf::build
{

/*!\brief Inserts values into an IBF.
 * \ingroup hibf_build
 * \details
 * Automatically does naive splitting if number_of_bins > 1.
 */
void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan::hibf::interleaved_bloom_filter & ibf,
                     concurrent_timer & fill_ibf_timer);

//!\overload
void insert_into_ibf(build_data const & data,
                     layout::layout::user_bin const & record,
                     seqan::hibf::interleaved_bloom_filter & ibf);

} // namespace seqan::hibf::build
