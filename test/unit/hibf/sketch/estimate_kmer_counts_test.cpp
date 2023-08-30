#include <gtest/gtest.h> // for Test, Message, TestPartResult, EXPECT_EQ, TestInfo

#include <cinttypes>   // for uint8_t
#include <cstddef>     // for size_t
#include <string>      // for basic_string, string
#include <string_view> // for string_view
#include <vector>      // for allocator, vector

#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>          // for hyperloglog

TEST(estimate_kmer_counts_test, small_example)
{
    std::vector<std::string> const input_sequences{
        {"ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"
         "AGCTATTTCAGACCTACACTATCTAGCTTATTCACAAATATTATAACGGCATACGTCTAGTGCTCATCGTGATCTAGCGA"
         "GCTAGCGATCTGATTCACGAGCGTACGTGACGTACGTATCGTACTACGTATCGTACTACATGCATCGATCGACGTAGCTA"
         "TCAGCGTAGCGTACGAGTCAGCTGACTGACGTCGTAGCATCGTACGTAGCGTAGCGATCGAGTCACTTATCGTAGCTAGT"
         "CGACTAGCGTACGTAGTCAGCTATTATGACGAGGCGACTTAGCGACTACGAGCTAGCGAGGAGGCGAGGCGAGCGGACTG"},
        {"ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"
         "AGCTATTTCAGACCTACACTATCTAGCTTATTCACAAATATTATAACGGCATACGTCTAGTGCTCATCGTGATCTAGCGA"
         "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"
         "GCTAGCGATCTGATTCACGAGCGTACGTGACGTACGTATCGTACTACGTATCGTACTACATGCATCGATCGACGTAGCTA"
         "TCAGCGTAGCGTACGAGTCAGCTGACTGACGTCGTAGCATCGTACGTAGCGTAGCGATCGAGTCACTTATCGTAGCTAGT"
         "CGACTAGCGTACGTAGTCAGCTATTATGACGAGGCGACTTAGCGACTACGAGCTAGCGAGGAGGCGAGGCGAGCGGACTG"},
        {"ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"
         "AGCTATTTCAGACCTACACTATCTAGCTTATTCACAAATATTATAACGGCATACGTCTAGTGCTCATCGTGATCTAGCGA"
         "GCTAGCGATCTGATTCACGAGCGTACGTGACGTACGTATCGTACTACGTATCGTACTACATGCATCGATCGACGTAGCTA"
         "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"
         "TCAGCGTAGCGTACGAGTCAGCTGACTGACGTCGTAGCATCGTACGTAGCGTAGCGATCGAGTCACTTATCGTAGCTAGT"
         "CGACTAGCGTACGTAGTCAGCTATTATGACGAGGCGACTTAGCGACTACGAGCTAGCGAGGAGGCGAGGCGAGCGGACTG"
         "G"}};

    uint8_t const kmer_size{19};
    size_t const b = 12;
    seqan::hibf::sketch::hyperloglog sketch(b);

    for (std::string_view seq : input_sequences)
    {
        // we have to go C-style here for the HyperLogLog Interface
        char const * c_seq_it = seq.begin();
        char const * end = c_seq_it + seq.size();

        while (c_seq_it + kmer_size <= end)
        {
            sketch.add(c_seq_it, kmer_size);
            ++c_seq_it;
        }
    }

    std::vector<seqan::hibf::sketch::hyperloglog> sketches{sketch, sketch};
    std::vector<size_t> kmer_counts;

    seqan::hibf::sketch::estimate_kmer_counts(sketches, kmer_counts);

    ASSERT_EQ(kmer_counts.size(), 2);
    EXPECT_EQ(kmer_counts[0], 581);
    EXPECT_EQ(kmer_counts[1], 581);
}