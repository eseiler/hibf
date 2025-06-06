// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <type_traits> // for is_base_of_v

#include <cereal/details/helpers.hpp> // for OutputArchiveBase, InputArchiveBase, BinaryInputArchive, BinaryOutputA...
#include <cereal/details/traits.hpp>  // for TextArchive

#include <hibf/platform.hpp>

namespace seqan::hibf
{

template <typename t>
concept cereal_output_archive = std::is_base_of_v<cereal::detail::OutputArchiveBase, t>;

template <typename t>
concept cereal_input_archive = std::is_base_of_v<cereal::detail::InputArchiveBase, t>;

template <typename t>
concept cereal_archive = cereal_output_archive<t> || cereal_input_archive<t>;

template <typename t>
concept cereal_text_archive = std::is_base_of_v<cereal::traits::TextArchive, t>;

} // namespace seqan::hibf
