// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <sharg/validators.hpp>

template <typename option_value_t>
    requires std::is_arithmetic_v<option_value_t>
struct positive_integer_validator
{
    using option_value_type = option_value_t;

    void operator()(option_value_type const & value) const
    {
        if (value < 1)
            throw sharg::validation_error{"The value must greater than 0."};
    }

    std::string get_help_page_message() const
    {
        return "Value must be greater than 0.";
    }
};
