#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "crimson_alignment_engine.hpp"
#include "include/crimson_mapperConfig.h"
#include <algorithm>
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

void help() { printf("Help\n"); } // TODO expand

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

  int opt;
  const struct option longOptions[] = {{"help", no_argument, 0, 0},
                                       {"version", no_argument, 0, 0}};
  int optionIndex;

  while ((opt = getopt_long(argc, argv, "h", longOptions, &optionIndex)) !=
         -1) {
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

  return 0;
}
