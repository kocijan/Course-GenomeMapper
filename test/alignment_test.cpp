#include "crimson_alignment_engine.hpp"
#include <gtest/gtest.h>

class AlignTest : public ::testing::Test {
protected:
  unsigned globalTestCnt = 0;
  unsigned localTestCnt = 0;
  struct testArgs {
    std::string query;
    std::string target;
    int match;
    int mismatch;
    int gap;
    int gap_open;
    int gap_extend;
    int expected;
  };
  std::vector<testArgs> globalTests, localTests;

  void SetUp() override {
    globalTests.push_back({"GATTACA", "GCATGCU", 1, -1, -1, 0, 0, 0});

    globalTestCnt = (unsigned)globalTests.size();

    localTests.push_back({"GGTTGACTA", "TGTTACGG", 3, -3, -2, 0, 0, 13});
    localTests.push_back(
        {"TACGGGCCCGCTAC", "TAGCCCTATCGGTCA", 5, -4, -1, 0, 0, 39});
    localTests.push_back(
        {"TACGGGCCCGCTAC", "TAGCCCTATCGGTCA", 5, -4, -1, 0, -1, 39});
    localTests.push_back(
        {"TACGGGCCCGCTAC", "TAGCCCTATCGGTCA", 5, -4, -1, -5, -1, 25});

    localTestCnt = (unsigned)localTests.size();
  }

  // void TearDown() override {
  //   globalTestQuery.clear();
  //   globalTestTarget.clear();
  //   globalTestMatch.clear();
  //   globalTestMismatch.clear();
  //   globalTestGap.clear();
  // }
};

TEST_F(AlignTest, Global) {
  for (testArgs i : globalTests) {
    std::string cigar;
    unsigned int target_begin;
    int retAlign =
        crimson::Align(i.query.c_str(), (unsigned int)i.query.length(),
                       i.target.c_str(), (unsigned int)i.target.length(),
                       crimson::AlignmentType::global, i.match, i.mismatch,
                       i.gap, &cigar, &target_begin, i.gap_open, i.gap_extend);
    EXPECT_EQ(retAlign, i.expected);
    fprintf(stderr, "%d\n", retAlign);
    std::cerr << cigar << '\n';
    fprintf(stderr, "%d\n", target_begin);
  }
}

TEST_F(AlignTest, Local) {
  for (testArgs i : localTests) {
    std::string cigar;
    unsigned int target_begin;
    int retAlign = crimson::Align(
        i.query.c_str(), (unsigned int)i.query.length(), i.target.c_str(),
        (unsigned int)i.target.length(), crimson::AlignmentType::local, i.match,
        i.mismatch, i.gap, &cigar, &target_begin, i.gap_open, i.gap_extend);
    EXPECT_EQ(retAlign, i.expected);
    fprintf(stderr, "%d\n", retAlign);
    std::cerr << cigar << '\n';
    fprintf(stderr, "%d\n", target_begin);
  }
}