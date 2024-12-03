// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Returns the number of technical bins available for use.
 * \param[in] tmax     The total number of bins.
 * \param[in] fraction The fraction of the total number of bins that should be empty.
 * \ingroup hibf
 * \sa https://godbolt.org/z/cMjbM39vj
 */
[[nodiscard]] constexpr size_t tmax_without_empty_bins(size_t const tmax, double const fraction) noexcept
{
    assert(fraction >= 0.0);
    assert(fraction < 1.0);
    size_t const number_of_empty_bins = std::clamp<size_t>(tmax * fraction, 1, tmax - 2) - (fraction == 0.0);
    assert(number_of_empty_bins == 0u || tmax > number_of_empty_bins);
    return tmax - number_of_empty_bins;
}

} // namespace seqan::hibf
