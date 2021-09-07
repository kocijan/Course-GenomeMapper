#include "crimson_alignment_engine.hpp"
#include <algorithm>
#include <string>
#include <vector>

namespace crimson {

int Align(const char *query, unsigned int query_len, const char *target,
          unsigned int target_len, AlignmentType type, int match, int mismatch,
          int gap, std::string *cigar, unsigned int *target_begin, int gap_open,
          int gap_extend) {

  using std::string;
  using std::vector;
  using uint = unsigned int;

  bool isAffineGap = gap_open != 0 && gap_extend != 0;
  bool needsCigar = cigar != nullptr && target_begin != nullptr;

  vector<vector<int>> dp(query_len + 1, vector<int>(target_len + 1));
  vector<vector<int>> dpD(query_len + 1, vector<int>(target_len + 1));
  vector<vector<int>> dpI(query_len + 1, vector<int>(target_len + 1));
  vector<vector<char>> prev(query_len + 1, vector<char>(target_len + 1));

  for (uint i = 1; i <= query_len; i++) {
    if (type == AlignmentType::local) {
      dp[i][0] = 0;
    } else {
      dp[i][0] = dp[i - 1][0] + gap;
    }
  }

  for (uint i = 1; i <= target_len; i++) {
    if (type == AlignmentType::local || type == AlignmentType::semiglobal) {
      dp[0][i] = 0;
    } else {
      dp[0][i] = dp[0][i - 1] + gap;
    }
  }

  int maxDpVal = 0;
  uint maxDpI = 0, maxDpJ = 0;

  for (uint i = 1; i <= query_len; i++) {
    for (uint j = 1; j <= target_len; j++) {
      int mscore =
          dp[i - 1][j - 1] + (query[i - 1] == target[j - 1] ? match : mismatch);
      dpI[i][j] =
          isAffineGap
              ? std::max(dp[i - 1][j] + gap_open, dpI[i - 1][j]) + gap_extend
              : dp[i - 1][j] + gap;
      dpD[i][j] =
          isAffineGap
              ? std::max(dp[i][j - 1] + gap_open, dpD[i][j - 1]) + gap_extend
              : dp[i][j - 1] + gap;

      dp[i][j] = std::max({mscore, dpI[i][j], dpD[i][j], 0});

      if (dp[i][j] == mscore) {
        prev[i][j] = 'M';
      } else if (dp[i][j] == dpI[i][j]) {
        prev[i][j] = 'I';
      } else if (dp[i][j] == dpD[i][j]) {
        prev[i][j] = 'D';
      }

      if (maxDpVal < dp[i][j]) {
        maxDpVal = dp[i][j];
        maxDpI = i;
        maxDpJ = j;
      }
    }
  }

  // for (uint i = 0; i <= query_len; ++i) {
  //   for (uint j = 0; j <= target_len; ++j) {
  //     printf("%d\t", dp[i][j]);
  //   }
  //   printf("\n");
  // }

  uint startI = query_len, startJ = target_len;

  int ret = dp[startI][startJ];

  if (type == AlignmentType::local) {
    startI = maxDpI;
    startJ = maxDpJ;
    ret = maxDpVal;
  }

  if (needsCigar) {
    uint i = startI;
    uint j = startJ;

    vector<char> longCigar;

    while (
        (type == AlignmentType::local && dp[i][j] > 0) ||
        ((type == AlignmentType::global || type == AlignmentType::semiglobal) &&
         i + j > 0)) {
      if (i == 0) {
        // D
        j--;
        longCigar.push_back('D');
        continue;
      }
      if (j == 0) {
        // I
        i--;
        longCigar.push_back('I');
        continue;
      }

      if (prev[i][j] == 'M') {
        // M
        i--;
        j--;
        longCigar.push_back('M');
        continue;
      } else if (prev[i][j] == 'I') {
        // I
        i--;
        longCigar.push_back('I');
        continue;
      } else if (prev[i][j] == 'D') {
        // D
        j--;
        longCigar.push_back('D');
        continue;
      }
    }

    std::reverse(longCigar.begin(), longCigar.end());

    int cigarNum = 0;
    for (uint i = 0; i < longCigar.size(); i++) {
      ++cigarNum;
      if (i + 1 == longCigar.size() || longCigar[i] != longCigar[i + 1]) {
        *cigar += std::to_string(cigarNum);
        *cigar += longCigar[i];
        cigarNum = 0;
      }
    }

    *target_begin = j;
  }

  return ret;
}

} // namespace crimson
