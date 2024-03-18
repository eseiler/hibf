#include <hibf/interleaved_bloom_filter.hpp>

// Shared parameters
uint32_t const num_bins{32u};
uint8_t const num_threads{8u};
uint8_t const num_hashes{3u};

struct bin_data
{
    std::vector<uint64_t> data{};
};
struct bins
{
    std::vector<bin_data> bins_data{};
};
bins my_bins{};

void with_config()
{
    auto input_lambda = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        for (auto const hash : my_bins.bins_data[user_bin_id].data)
            it = hash;
    };

    seqan::hibf::config config{.input_fn = input_lambda,
                               .number_of_user_bins = num_bins,
                               .number_of_hash_functions = num_hashes,
                               .maximum_fpr = 0.05,
                               .threads = num_threads};

    // This will do I/O twice:
    // 1. Determine the maximum number of elements in any bin. Used to compute the optimal bin size for the given FPR.
    // 2. Insert all elements into the IBF.
    seqan::hibf::interleaved_bloom_filter cpu_ibf{config};
}

void without_config()
{
    // 8192u * 8192u might lead to a very high FPR, or might waste a lot of space.
    seqan::hibf::interleaved_bloom_filter cpu_ibf{seqan::hibf::bin_count{num_bins},
                                                  seqan::hibf::bin_size{8192u * 8192u},
                                                  seqan::hibf::hash_function_count{num_hashes}};

    // There will be concurrent write access for emplace.
    // Giving each thread a chunk of bins to insert into the IBF will reduce the likelyhood of concurrent writes to the
    // same memory. 64 bins per thread should be completely save, but 8 bins per thread allows for more parallelism.
    // I've never seen problems with using 8 bins per thread. I have seen it without any chunking :)
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(num_bins / num_threads), 8u, 64u);

    // Parallelise over the input bins.
    // For your setup, parallelising over the input data itself might be better.
    // You could test both versions. Not sure what the overhead of the parallel for is.
    // We usually do over the input bins because there might be some I/O - not the case for you.
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(num_threads)
    for (size_t user_bin_id = 0u; user_bin_id < num_bins; ++user_bin_id)
    {
        for (uint64_t const hash : my_bins.bins_data[user_bin_id].data)
            cpu_ibf.emplace(hash, seqan::hibf::bin_index{user_bin_id});
    }
}

void query(seqan::hibf::interleaved_bloom_filter const & cpu_ibf, std::vector<uint64_t> const & queries)
{
    auto agent = cpu_ibf.membership_agent();

// Each thread needs its own agent. The agent is not thread-safe.
// The agent holds a bit_vector to store the results.
// `bulk_contains` just returns a reference (!) to that vector.
#pragma omp parallel for schedule(dynamic) num_threads(num_threads) private(agent)
    for (size_t i = 0u; i < queries.size(); ++i)
    {
        auto & result = agent.bulk_contains(queries[i]);
    }

    // Above is the same as:
    // #pragma omp parallel num_threads(num_threads)
    //     {
    //         auto agent = cpu_ibf.membership_agent();
    // #pragma omp for schedule(dynamic)
    //         for (size_t i = 0u; i < queries.size(); ++i)
    //         {
    //             auto & result = agent.bulk_contains(queries[i]);
    //         }
    //     }
}
