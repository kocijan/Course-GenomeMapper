#include "crimson_minimizer_engine.hpp"
#include <algorithm>
#include <bitset>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

class MinimizeTest : public ::testing::Test {
protected:
  struct singleTestArgs {
    std::string sequence;
    unsigned kmer_len;
    unsigned window_len;
    unsigned expectedSize;
  };
  std::vector<singleTestArgs> singleTests;

  struct mapTestArgs {
    std::vector<const char *> reference_sequences;
    unsigned kmer_len;
    unsigned window_len;
    double frequency;
    std::string query_sequence;
    unsigned expectedSize;
  };
  std::vector<mapTestArgs> mapTests;

  void SetUp() override {
    singleTests.push_back({"GTCATGCACGTTCAC", 3, 4, 4});
    singleTests.push_back({"GTCATGCACGTTCAC", 3, 3, 4});
    singleTests.push_back({"ACCCAAC", 3, 4, 2});
    singleTests.push_back({"AAAAAAAAAAAAA", 3, 3, 3});
    mapTests.push_back({{"GTCATGCACGTTCAC"}, 3, 3, 0.0, "GTCATGCACGTTCAC", 4});
    mapTests.push_back({{"GTCATGCACGTTCAC"}, 3, 3, 0.5, "GTCATGCACGTTCAC", 2});
    mapTests.push_back(
        {{"AAAAATATACG", "GCATTGAC"}, 3, 3, 0.0, "AAATGCTATACGA", 3});
    mapTests.push_back(
        {{"AAAAAAACCCCCCCC", "CCCCAAAAAAAAAAA"}, 3, 3, 0.0, "AAAAAATAAAAA", 2});
    mapTests.push_back(
        {{"AAAAAAATTTTTTTT", "TTTTAAAAAAAAAAA"}, 3, 3, 0.0, "AAAAAAAAAAAA", 3});
    mapTests.push_back({{"AAAAAAATTTTTTTTTTTTTTTTTTT", "TTTTAAAAAAAAAAA"},
                        3,
                        3,
                        0.0,
                        "AAAAAAAAAAAAAAAAAAAAAA",
                        6});
  };

  void TearDown() override { crimson::ResetData(); }

  std::string KmerString(unsigned int kmer, unsigned int kmer_len) {
    std::string ret;
    for (unsigned i = 0; i < kmer_len; ++i) {
      ret += (kmer & 3) + '0';
      kmer >>= 2;
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
  }
};

TEST_F(MinimizeTest, Single) {
  for (singleTestArgs i : singleTests) {
    auto minimizers =
        crimson::Minimize(i.sequence.c_str(), (unsigned int)i.sequence.size(),
                          i.kmer_len, i.window_len);
    fprintf(stderr, "%zd\n", minimizers.size());
    for (auto j : minimizers) {
      std::cerr << KmerString(std::get<0>(j), i.kmer_len);
      fprintf(stderr, "\t%u\t%d\n", std::get<1>(j), std::get<2>(j));
    }
    fprintf(stderr, "\n");
    EXPECT_EQ(minimizers.size(), i.expectedSize);
  }
}

TEST_F(MinimizeTest, Map) {
  for (mapTestArgs i : mapTests) {
    std::vector<unsigned int> ref_seq_sizes;
    for (const char *j : i.reference_sequences)
      ref_seq_sizes.push_back((unsigned int)strlen(j));

    crimson::Minimize(i.reference_sequences, ref_seq_sizes, i.kmer_len,
                      i.window_len);

    crimson::Filter(i.frequency);

    auto overlaps = crimson::Map(i.query_sequence.c_str(),
                                 (unsigned)i.query_sequence.length());

    fprintf(stderr, "%zd\n", overlaps.size());

    for (crimson::Overlap j : overlaps) {
      std::cerr << KmerString(j.kmer, i.kmer_len);
      fprintf(stderr, "\t%u\t%u\t%u\t%d\t%d\n", j.reference_index, j.query_pos,
              j.reference_pos, j.is_original_query, j.is_original_reference);
    }

    fprintf(stderr, "\n");
    EXPECT_EQ(overlaps.size(), i.expectedSize);

    crimson::ResetData();
  }
}