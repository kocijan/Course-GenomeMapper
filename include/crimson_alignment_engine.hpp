#ifndef CRIMSON_ALIGNMENT_ENGINE_HPP_
#define CRIMSON_ALIGNMENT_ENGINE_HPP_

#include <string>

namespace crimson {

enum class AlignmentType { global, local, semiglobal };

int Align(const char *query, unsigned int query_len, const char *target,
          unsigned int target_len, AlignmentType type, int match, int mismatch,
          int gap, std::string *cigar = nullptr,
          unsigned int *target_begin = nullptr, int gap_open = 0,
          int gap_extend = 0);

} // namespace crimson

#endif // CRIMSON_ALIGNMENT_ENGINE_HPP_