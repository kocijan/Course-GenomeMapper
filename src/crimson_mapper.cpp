#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "crimson_alignment_engine.hpp"
#include "crimson_minimizer_engine.hpp"
#include "include/crimson_mapperConfig.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

void help() {
  printf(R"((--version - show version
-h - show help
-c - calculate alignment (default: false)
-a <str> - alignment type (default: global)
-m <int> - match cost (default: 3)
-n <int> - mismatch cost (default: -5)
-g <int> - gap cost (default: -4)
-k <int> - k-mer size (default: 15)
-w <int> - window size (default: 10)
-f <int> - k-mer frequency threshold (default: 0.001)
))");
}

void version() {
  printf("v%d.%d.%d\n", crimson_mapper_VERSION_MAJOR,
         crimson_mapper_VERSION_MINOR, crimson_mapper_VERSION_PATCH);
}

int main(int argc, char **argv) {
  using std::cerr;
  using std::cout;
  using std::string;
  using std::unique_ptr;
  using std::vector;
  using seqsize_t = std::uint32_t;
  using namespace crimson;

  int opt;
  const struct option longOptions[] = {{"help", no_argument, 0, 0},
                                       {"version", no_argument, 0, 0}};
  int optionIndex;

  bool calcAlignment = false;
  AlignmentType alignType = AlignmentType::global;
  int matchCost = 3;
  int mismatchCost = -5;
  int gapCost = -4;
  unsigned int KmerSize = 15;
  unsigned int windowSize = 10;
  double freqThreshold = 0.001;

  while ((opt = getopt_long(argc, argv, "hca:m:n:g:k:w:f:", longOptions,
                            &optionIndex)) != -1) {
    if (opt == 0) {
      string curLongOpt = longOptions[optionIndex].name;
      if (curLongOpt == "help") {
        help();
        return 0;
      } else if (curLongOpt == "version") {
        version();
        return 0;
      }
    } else if (opt == 'h') {
      help();
      return 0;
    } else if (opt == 'c') {
      calcAlignment = true;
    } else if (opt == 'a') {
      cout << "[" << optarg << "]\n";
      if (std::strcmp(optarg, "global") == 0) {
        alignType = AlignmentType::global;
      } else if (std::strcmp(optarg, "local") == 0) {
        alignType = AlignmentType::local;
      } else if (std::strcmp(optarg, "semiglobal") == 0) {
        alignType = AlignmentType::semiglobal;
      }
    } else if (opt == 'm') {
      matchCost = std::stoi(optarg);
    } else if (opt == 'n') {
      mismatchCost = std::stoi(optarg);
    } else if (opt == 'g') {
      gapCost = std::stoi(optarg);
    } else if (opt == 'k') {
      KmerSize = (unsigned)std::stoi(optarg);
    } else if (opt == 'w') {
      windowSize = (unsigned)std::stoi(optarg);
    } else if (opt == 'f') {
      freqThreshold = std::stod(optarg);
    }
  }

  if (optind >= argc - 1) {
    fprintf(stderr,
            "[crimson_mapper] error: mandatory file arguments not provided");
    return 0;
  }

  int optindI = optind;
  string refFilename = argv[optindI];
  ++optindI;
  vector<string> fragFilenames(argv + optindI, argv + argc);

  struct Sequence {
    const char *name;
    seqsize_t nameLen;
    const char *data;
    seqsize_t dataLen;
    Sequence(const char *name, seqsize_t nameLen, const char *data,
             seqsize_t dataLen)
        : name(name), nameLen(nameLen), data(data), dataLen(dataLen) {}
  };

  auto refParser =
      bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(refFilename);

  const vector<unique_ptr<Sequence>> parsedRef =
      refParser->Parse(UINT64_MAX, false);

  // TODO support FASTQ

  vector<unique_ptr<Sequence>> parsedFrags;
  for (const string &fragFilename : fragFilenames) {
    unique_ptr<bioparser::Parser<Sequence>> fragParser =
        bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
            fragFilename);
    vector<unique_ptr<Sequence>> parsedFrag =
        fragParser->Parse(UINT64_MAX, false);
    parsedFrags.insert(parsedFrags.end(),
                       std::make_move_iterator(parsedFrag.begin()),
                       std::make_move_iterator(parsedFrag.end()));
  }

  const size_t refCnt = parsedRef.size();
  const size_t fragCnt = parsedFrags.size();

  vector<seqsize_t> fragLens(fragCnt);

  for (size_t i = 0; i < fragCnt; ++i) {
    fragLens[i] = parsedFrags[i]->dataLen;
  }

  const seqsize_t fragTotalLen =
      std::accumulate(fragLens.begin(), fragLens.end(), 0u);
  const seqsize_t fragMinLen =
      *std::min_element(fragLens.begin(), fragLens.end());
  const seqsize_t fragMaxLen =
      *std::max_element(fragLens.begin(), fragLens.end());
  double fragAvgLen = double(fragTotalLen) / double(fragCnt);
  seqsize_t fragN50 = 0;

  sort(fragLens.begin(), fragLens.end(), std::greater<seqsize_t>());
  seqsize_t fragCummLen = 0;
  for (size_t i = 0; i < fragLens.size(); ++i) {
    fragCummLen += fragLens[i];
    if (fragCummLen >= ceil(fragTotalLen / 2.0)) {
      fragN50 = fragLens[i];
      break;
    }
  }

  fprintf(stderr, "Reference genome statistics\n");
  fprintf(stderr, "Name: %.*s\n", parsedRef[0]->nameLen, parsedRef[0]->name);
  cerr << "Length: " << parsedRef[0]->dataLen << "\n\n";

  fprintf(stderr, "Fragment statistics\n");
  fprintf(stderr, "Number of fragments: %zd\n", fragCnt);
  cerr << "Total length: " << fragTotalLen << "\n";
  fprintf(stderr, "Average length: %f\n", fragAvgLen);
  cerr << "N50 length: " << fragN50 << "\n";
  cerr << "Minimal length: " << fragMinLen << "\n";
  cerr << "Maximal length: " << fragMaxLen << "\n\n";

  vector<const char *> refSequences;
  vector<unsigned int> refSeqLens;
  for (size_t i = 0; i < refCnt; ++i) {
    refSequences.push_back(parsedRef[i]->data);
    refSeqLens.push_back(parsedRef[i]->dataLen);
  }

  Minimize(refSequences, refSeqLens, KmerSize, windowSize);
  Filter(freqThreshold);

  // TODO parallelize

  for (size_t i = 0; i < fragCnt; i++) {
    vector<Overlap> overlaps =
        Map(parsedFrags[i]->data, parsedFrags[i]->dataLen);
    Overlap firstOverlap = overlaps[0];
    Overlap lastOverlap = overlaps.back();
    unsigned int j = firstOverlap.reference_index;

    unsigned int q_begin = firstOverlap.query_pos;
    unsigned int q_end = lastOverlap.query_pos + KmerSize;
    unsigned int t_begin = firstOverlap.reference_pos;
    unsigned int t_end = lastOverlap.reference_pos + KmerSize;

    printf("%.*s\t%d\t%d\t%d\t%c\t%.*s\t%d\t%d\t%d", parsedFrags[i]->nameLen,
           parsedFrags[i]->name, parsedFrags[i]->dataLen, q_begin, q_end, '+',
           parsedRef[j]->nameLen, parsedRef[j]->name, parsedRef[j]->dataLen,
           t_begin, t_end);

    if (calcAlignment) {
      string cigar;
      unsigned int target_begin;

      Align(parsedFrags[i]->data + q_begin, q_end - q_begin,
            parsedRef[j]->data + t_begin, t_end - t_begin, alignType, matchCost,
            mismatchCost, gapCost, &cigar, &target_begin);

      int curSum = 0, mSum = 0, totalSum = 0;
      for (unsigned i = 0; i < cigar.size(); ++i) {
        if (cigar[i] == 'M') {
          mSum += curSum;
          totalSum += curSum;
          curSum = 0;
        } else if (std::isdigit(cigar[i])) {
          curSum *= 10;
          curSum += cigar[i] - '0';
        } else {
          totalSum += curSum;
          curSum = 0;
        }
      }

      printf("\t%d\t%d\t%d\tcg:Z:", mSum, totalSum, 255);
      cout << cigar;
    } else {
      unsigned int lenQ = q_end - q_begin;
      unsigned int lenT = t_end - t_begin;
      unsigned int minLen = std::min(lenQ, lenT);
      printf("\t%d\t%d\t%d", minLen / 2, lenQ + lenT - minLen / 2, 255);
    }
    printf("\n");
  }

  return 0;
}
