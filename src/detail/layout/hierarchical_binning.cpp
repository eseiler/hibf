#include <hibf/detail/layout/hierarchical_binning.hpp>

namespace hibf::layout
{

size_t hierarchical_binning::execute()
{
    assert(data != nullptr);
    assert(data->positions.size() <= data->kmer_counts.size());

    static constexpr size_t max_size_t{std::numeric_limits<size_t>::max()};

    if (!data->user_bins_arranged)
    {
        sketch::toolbox sketch_toolbox{data->kmer_counts, data->sketches, data->positions};
        sketch_toolbox.sort_by_cardinalities();

        if (!config.disable_estimate_union && !config.disable_rearrangement)
            sketch_toolbox.rearrange_bins(config.max_rearrangement_ratio, config.threads);

        data->user_bins_arranged = true;
    }

    // technical bins (outer) = rows; user bins (inner) = columns
    std::vector<std::vector<size_t>> matrix(num_technical_bins, std::vector<size_t>(num_user_bins, max_size_t));

    // technical bins (outer) = rows; user bins (inner) = columns
    std::vector<std::vector<size_t>> ll_matrix(num_technical_bins, std::vector<size_t>(num_user_bins, 0u));

    // technical bins (outer) = rows; user bins (inner) = columns
    std::vector<std::vector<std::pair<size_t, size_t>>> trace(
        num_technical_bins,
        std::vector<std::pair<size_t, size_t>>(num_user_bins, {max_size_t, max_size_t}));

    initialization(matrix, ll_matrix, trace);

    recursion(matrix, ll_matrix, trace);

    // print_matrix(matrix, num_technical_bins, num_user_bins, max_size_t);
    // print_matrix(ll_matrix, num_technical_bins, num_user_bins, max_size_t);
    // print_matrix(trace, num_technical_bins, num_user_bins, std::make_pair(max_size_t, max_size_t));

    return backtracking(trace);
}

void hierarchical_binning::initialization(std::vector<std::vector<size_t>> & matrix,
                                          std::vector<std::vector<size_t>> & ll_matrix,
                                          std::vector<std::vector<std::pair<size_t, size_t>>> & trace)
{
    assert(data != nullptr);

    // initialize first column
    double const ub_cardinality = static_cast<double>(data->kmer_counts[data->positions[0]]);
    for (size_t i = 0; i < num_technical_bins; ++i)
    {
        size_t const corrected_ub_cardinality = static_cast<size_t>(ub_cardinality * data->fp_correction[i + 1]);
        matrix[i][0] = corrected_ub_cardinality / (i + 1);
        trace[i][0] = {0u, 0u}; // unnecessary?
    }

    // initialize first row
    size_t sum = data->kmer_counts[data->positions[0]];
    if (!config.disable_estimate_union)
    {
        sketch::toolbox::precompute_initial_union_estimates(data->union_estimates,
                                                            data->sketches,
                                                            data->kmer_counts,
                                                            data->positions);

        for (size_t j = 1; j < num_user_bins; ++j)
        {
            sum += data->kmer_counts[data->positions[j]];
            matrix[0][j] = data->union_estimates[j];
            ll_matrix[0][j] = max_merge_levels(j + 1) * sum;
            trace[0][j] = {0u, j - 1}; // unnecessary?
        }
    }
    else
    {
        for (size_t j = 1; j < num_user_bins; ++j)
        {
            assert(j < data->positions.size());
            assert(data->positions[j] < data->kmer_counts.size());
            sum += data->kmer_counts[data->positions[j]];
            matrix[0][j] = sum;
            ll_matrix[0][j] = max_merge_levels(j + 1) * sum;
            trace[0][j] = {0u, j - 1}; // unnecessary?
        }
    }
}

void hierarchical_binning::recursion(std::vector<std::vector<size_t>> & matrix,
                                     std::vector<std::vector<size_t>> & ll_matrix,
                                     std::vector<std::vector<std::pair<size_t, size_t>>> & trace)
{
    assert(data != nullptr);

    // we must iterate column wise
    for (size_t j = 1; j < num_user_bins; ++j)
    {
        size_t const current_weight = data->kmer_counts[data->positions[j]];
        double const ub_cardinality = static_cast<double>(current_weight);

        if (!config.disable_estimate_union)
            sketch::toolbox::precompute_union_estimates_for(data->union_estimates,
                                                            data->sketches,
                                                            data->kmer_counts,
                                                            data->positions,
                                                            j);

        for (size_t i = 1; i < num_technical_bins; ++i)
        {
            size_t minimum{std::numeric_limits<size_t>::max()};
            size_t full_minimum{std::numeric_limits<size_t>::max()};

            // check vertical cells
            for (size_t i_prime = 0; i_prime < i; ++i_prime)
            {
                // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                size_t const corrected_ub_cardinality =
                    static_cast<size_t>(ub_cardinality * data->fp_correction[(i - i_prime)]);
                size_t score = std::max<size_t>(corrected_ub_cardinality / (i - i_prime), matrix[i_prime][j - 1]);
                size_t full_score = score * (i + 1) /*#TBs*/ + config.alpha * ll_matrix[i_prime][j - 1];

                // std::cout << " ++ j:" << j << " i:" << i << " i':" << i_prime << " score:" << score << std::endl;

                if (full_score < full_minimum)
                {
                    minimum = score;
                    full_minimum = full_score;
                    trace[i][j] = {i_prime, j - 1};
                    ll_matrix[i][j] = ll_matrix[i_prime][j - 1];
                }
            }

            // seqan3::debug_stream << "current vertical minimum of " << "j:" << j << " i:" << i
            //                      << " -> score:" << full_minimum << " (M_ij=" << minimum << ")"
            //                      << " trace:" << trace[i][j]
            //                      << std::endl;

            // check horizontal cells
            size_t j_prime{j - 1};
            size_t weight{current_weight};

            auto get_weight = [&]()
            {
                // if we use the union estimate we plug in that value instead of the sum (weight)
                // union_estimates[j_prime] is the union of {j_prime, ..., j}
                // the + 1 is necessary because j_prime is decremented directly after weight is updated
                return config.disable_estimate_union ? weight : data->union_estimates[j_prime + 1];
            };

            // if the user bin j-1 was not split into multiple technical bins!
            // I may merge the current user bin j into the former
            while (j_prime != 0 && ((i - trace[i][j_prime].first) < 2) && get_weight() < minimum)
            {
                weight += data->kmer_counts[data->positions[j_prime]];
                --j_prime;

                // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                // ll_kmers: estimate for the number of k-mers that have to be resolved on lower levels
                // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                size_t const score = std::max<size_t>(matrix[i - 1][j_prime], get_weight());
                size_t const ll_kmers = ll_matrix[i - 1][j_prime] + max_merge_levels(j - j_prime) * weight;
                size_t const full_score = score * (i + 1) /*#TBs*/ + config.alpha * ll_kmers;

                // seqan3::debug_stream << " -- " << "j_prime:" << j_prime
                //                      << " -> full_score:" << full_score << " (M_{i-1,j'}=" << score << ")"
                //                      << std::endl;

                if (full_score < full_minimum)
                {
                    minimum = score;
                    full_minimum = full_score;
                    trace[i][j] = {i - 1, j_prime};
                    ll_matrix[i][j] = ll_kmers;
                }
            }

            matrix[i][j] = minimum;
        }
    }
}

size_t hierarchical_binning::backtracking(std::vector<std::vector<std::pair<size_t, size_t>>> const & trace)
{
    assert(data != nullptr);

    // backtracking starts at the bottom right corner:
    size_t trace_i = num_technical_bins - 1;
    size_t trace_j = num_user_bins - 1;

    // while backtracking, keep trach of the following variables
    size_t high_level_max_id{};   // the id of the technical bin with maximal size
    size_t high_level_max_size{}; // the maximum technical bin size seen so far
    size_t bin_id{};              // the current bin that is processed, we start naming the bins here!

    // process the trace starting at the bottom right call until you arrive at the first row or column
    while (trace_j > 0u && trace_i > 0u)
    {
        // std::cout << "\t I am now at " << trace_i << "," << trace_j << std::endl;
        size_t next_i = trace[trace_i][trace_j].first;
        size_t next_j = trace[trace_i][trace_j].second;

        size_t kmer_count = data->kmer_counts[data->positions[trace_j]];
        size_t number_of_bins = (trace_i - next_i);

        if (number_of_bins == 1 && next_j != trace_j - 1u) // merged bin
        {
            auto libf_data = initialise_libf_data(trace_j);

            // std::cout << "merged [" << trace_j;
            --trace_j;
            while (trace_j != next_j)
            {
                kmer_count += data->kmer_counts[data->positions[trace_j]];
                libf_data.positions.push_back(data->positions[trace_j]);
                // std::cout << "," << trace_j;
                --trace_j;
            }
            trace_i = next_i;
            trace_j = next_j; // unneccessary?

            process_merged_bin(libf_data, bin_id);

            update_max_id(high_level_max_id, high_level_max_size, bin_id, kmer_count);
            // std::cout << "]: " << kmer_count << std::endl;
        }
        else // split bin
        {
            size_t const kmer_count_per_bin = (kmer_count + number_of_bins - 1) / number_of_bins; // round up

            data->hibf_layout->user_bins.emplace_back(data->positions[trace_j],
                                                      data->previous.bin_indices,
                                                      number_of_bins,
                                                      bin_id);

            // std::cout << "split " << trace_j << " into " << number_of_bins << ": " << kmer_count_per_bin << std::endl;

            update_max_id(high_level_max_id, high_level_max_size, bin_id, kmer_count_per_bin);

            trace_i = trace[trace_i][trace_j].first;
            --trace_j;
        }

        bin_id += number_of_bins;
    }

    // process the first row or first column at last
    assert(trace_i == 0 || trace_j == 0);
    if (trace_i == 0u && trace_j > 0u) // the last UBs get merged into the remaining TB
    {
        size_t kmer_count = data->kmer_counts[data->positions[trace_j]];
        auto libf_data = initialise_libf_data(trace_j);

        // std::cout << "merged [" << trace_j;
        while (trace_j > 0)
        {
            --trace_j;
            kmer_count += data->kmer_counts[data->positions[trace_j]];
            libf_data.positions.push_back(data->positions[trace_j]);
            // std::cout << "," << trace_j;
        }
        assert(trace_j == 0);

        process_merged_bin(libf_data, bin_id);

        update_max_id(high_level_max_id, high_level_max_size, bin_id, kmer_count);

        // std::cout << "]: " << kmer_count << std::endl;
        // std::cout << "\t I am now at " << trace_i << "," << trace_j << std::endl;
    }
    else if (trace_j == 0u) // the last UB is split into the remaining TBs
    {
        // we only arrive here if the first user bin (UB-0) wasn't merged with some before so it is safe to assume
        // that the bin was split (even if only into 1 bin).
        size_t const kmer_count = data->kmer_counts[data->positions[0]];
        size_t const number_of_tbs = trace_i + 1;
        size_t const average_bin_size = (kmer_count + number_of_tbs - 1) / number_of_tbs; // round up

        data->hibf_layout->user_bins.emplace_back(data->positions[0],
                                                  data->previous.bin_indices,
                                                  number_of_tbs,
                                                  bin_id);

        update_max_id(high_level_max_id, high_level_max_size, bin_id, average_bin_size);
        // std::cout << "split " << trace_j << " into " << trace_i << ": " << kmer_count / number_of_tbs << std::endl;
    }

    return high_level_max_id;
}

} // namespace hibf::layout