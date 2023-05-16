// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <chopper/layout/layout.hpp>

namespace hibf
{

void update_user_bins(std::vector<int64_t> & filename_indices, chopper::layout::layout::user_bin const & record)
{
    std::fill_n(filename_indices.begin() + record.storage_TB_id, record.number_of_technical_bins, record.idx);
}

} // namespace hibf