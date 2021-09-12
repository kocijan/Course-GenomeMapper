#ifndef CRIMSON_MINIMIZER_ENGINE_HPP_
#define CRIMSON_MINIMIZER_ENGINE_HPP_

#include <tuple>
#include <unordered_map>
#include <vector>
namespace crimson {

std::vector<std::tuple<unsigned int, unsigned int, bool>>
Minimize(const char *sequence, unsigned int sequence_len, unsigned int kmer_len,
         unsigned int window_len);

void Minimize(std::vector<const char *> sequence,
              std::vector<unsigned int> sequence_len, unsigned int kmer_len,
              unsigned int window_len);

void Filter(double frequency);

struct Overlap {
  unsigned kmer;
  unsigned reference_index;
  unsigned query_pos;
  unsigned reference_pos;
  bool is_original_query;
  bool is_original_reference;
};

void ResetData();

std::vector<Overlap> Map(const char *sequence, unsigned int sequence_len);

} // namespace crimson

#endif // CRIMSON_MINIMIZER_ENGINE_HPP_