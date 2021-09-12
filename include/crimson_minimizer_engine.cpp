#include "crimson_minimizer_engine.hpp"
#include <algorithm>
#include <bitset>
#include <iostream>
#include <iterator>
#include <optional>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace crimson {

unsigned int CompressBase(char base, bool reverse = false) {
  if (reverse) {
    if (base == 'A')
      base = 'T';
    else if (base == 'C')
      base = 'G';
    else if (base == 'G')
      base = 'C';
    else if (base == 'T')
      base = 'A';
  }
  if (base == 'A')
    return 0;
  else if (base == 'C')
    return 1;
  else if (base == 'G')
    return 2;
  else // if (base == 'T')
    return 3;
}

bool MinKmerCmp(std::tuple<unsigned int, unsigned int, bool> a,
                std::tuple<unsigned int, unsigned int, bool> b) {
  using std::get;
  return std::tuple(get<0>(a), -(int)get<1>(a), get<2>(a)) <
         std::tuple(get<0>(b), -(int)get<1>(b), get<2>(b));
}

std::vector<std::tuple<unsigned int, unsigned int, bool>>
Minimize(const char *sequence, unsigned int sequence_len, unsigned int kmer_len,
         unsigned int window_len) {
  using std::get;
  using std::multiset;
  using std::queue;
  using std::tuple;
  using std::vector;
  typedef tuple<unsigned int, unsigned int, bool> uub;

  vector<uub> ret;

  multiset<uub, decltype(&MinKmerCmp)> windowKmers(&MinKmerCmp);
  queue<uub> kmerQueue;
  unsigned curKmer = 0, curKmerRev = 0;
  unsigned curMin, curMinI;
  bool curMinOrigin;

  std::cout << '\n';

  for (unsigned i = 0; i < sequence_len; ++i) {
    curKmer *= 4u;
    curKmer += CompressBase(sequence[i]);
    curKmer &= (1 << (kmer_len * 2)) - 1;
    curKmerRev *= 4u;
    curKmerRev += CompressBase(sequence[i], true);
    curKmerRev &= (1 << (kmer_len * 2)) - 1;

    if (i < kmer_len - 1)
      continue;

    unsigned minCurKmer;
    bool isOriginal;

    if (curKmer < curKmerRev) {
      minCurKmer = curKmer;
      isOriginal = true;
    } else {
      minCurKmer = curKmerRev;
      isOriginal = false;
    }

    uub curKmerData = {minCurKmer, i - kmer_len + 1, isOriginal};

    windowKmers.insert(curKmerData);
    kmerQueue.push(curKmerData);

    if (i < window_len + kmer_len - 2)
      continue;

    uub setMin = *windowKmers.begin();

    if (i == window_len + kmer_len - 2 ||
        i - kmer_len - curMinI + 1 >= window_len || get<0>(setMin) != curMin ||
        get<2>(setMin) != curMinOrigin) {
      curMin = get<0>(setMin);
      curMinI = get<1>(setMin);
      curMinOrigin = get<2>(setMin);
      ret.push_back(setMin);
    }

    windowKmers.erase(kmerQueue.front());
    kmerQueue.pop();
  }

  return ret;
}

unsigned int kmer_len_last, window_len_last;

size_t refsTotal;

std::unordered_map<unsigned int,
                   std::vector<std::tuple<unsigned int, unsigned int, bool>>>
    minLookup;

void ResetData() { minLookup.clear(); }

void Minimize(std::vector<const char *> sequence,
              std::vector<unsigned int> sequence_len, unsigned int kmer_len,
              unsigned int window_len) {
  using std::get;
  using std::multiset;
  using std::pair;
  using std::tuple;
  using std::unordered_map;
  using std::vector;
  typedef tuple<unsigned int, unsigned int, bool> uub;
  typedef vector<uub> vuub;

  kmer_len_last = kmer_len;
  window_len_last = window_len;
  refsTotal = sequence.size();

  for (unsigned i = 0; i < sequence.size(); ++i) {
    vuub mins = Minimize(sequence[i], sequence_len[i], kmer_len, window_len);

    std::reverse(mins.begin(), mins.end());

    for (uub j : mins) {
      minLookup[get<0>(j)].push_back({i, get<1>(j), get<2>(j)});
    }
  }
}

void Filter(double frequency) {
  std::vector<std::pair<unsigned, unsigned>> sizes;
  for (auto i : minLookup) {
    sizes.push_back({i.second.size(), i.first});
  }
  std::sort(sizes.begin(), sizes.end(),
            std::greater<std::pair<unsigned, unsigned>>());
  for (unsigned i = 0; i < frequency * double(minLookup.size()); ++i) {
    minLookup.erase(sizes[i].second);
  }
}

bool operator<(const Overlap &x, const Overlap &y) {
  if (x.reference_pos + kmer_len_last <= y.reference_pos &&
      x.query_pos + kmer_len_last <= y.query_pos)
    return true;

  if (x.reference_pos < y.reference_pos && x.query_pos < y.query_pos &&
      y.reference_pos - x.reference_pos == y.query_pos - x.query_pos)
    return true;

  return false;
}

std::vector<Overlap> Map(const char *sequence, unsigned int sequence_len) {
  using std::get;
  using std::multiset;
  using std::pair;
  using std::tuple;
  using std::unordered_map;
  using std::vector;
  typedef tuple<unsigned int, unsigned int, bool> uub;

  unsigned kmer_len = kmer_len_last;
  unsigned window_len = window_len_last;

  auto queryMins = Minimize(sequence, sequence_len, kmer_len, window_len);

  vector<vector<Overlap>> lis(refsTotal), overlaps(refsTotal);

  vector<vector<int>> ind(refsTotal), prev(refsTotal);

  vector<unsigned> jInd(refsTotal);

  for (uub i : queryMins) {
    for (uub j : minLookup[get<0>(i)]) {
      // if (get<2>(i) != get<2>(j))
      //   continue;
      Overlap curOverlap = {get<0>(i), get<0>(j), get<1>(i),
                            get<1>(j), get<2>(i), get<2>(j)};
      overlaps[get<0>(j)].push_back(curOverlap);
      auto lisPos = std::lower_bound(lis[get<0>(j)].begin(),
                                     lis[get<0>(j)].end(), curOverlap);
      unsigned lisPosInd =
          (unsigned)std::distance(lis[get<0>(j)].begin(), lisPos);
      if (lisPos != lis[get<0>(j)].end()) {
        *lisPos = curOverlap;
        ind[get<0>(j)][lisPosInd] = (int)jInd[get<0>(j)];
        prev[get<0>(j)].push_back(lisPosInd ? ind[get<0>(j)][lisPosInd - 1]
                                            : -1);
      } else {
        lis[get<0>(j)].push_back(curOverlap);
        ind[get<0>(j)].push_back((int)jInd[get<0>(j)]);
        prev[get<0>(j)].push_back(lisPosInd ? ind[get<0>(j)][lisPosInd - 1]
                                            : -1);
      }
      ++jInd[get<0>(j)];
    }
  }

  size_t maxLis = 0;
  unsigned maxLisI;

  for (unsigned i = 0; i < refsTotal; ++i) {
    if (maxLis < lis[i].size()) {
      maxLis = lis[i].size();
      maxLisI = i;
    }
  }

  vector<Overlap> ret;

  for (int i = (int)ind[maxLisI].back(); i != -1;
       i = prev[maxLisI][(unsigned)i]) {
    ret.push_back(overlaps[maxLisI][(unsigned)i]);
  }

  std::reverse(ret.begin(), ret.end());

  return ret;
}

} // namespace crimson