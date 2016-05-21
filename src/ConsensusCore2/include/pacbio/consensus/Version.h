#include <string>

#define CC2_GIT_SHA1 "b156ec1*"
#define CC2_VERSION "0.11.0"

namespace PacBio {
namespace Consensus {

struct Version {
  static const size_t Major = 0;
  static const size_t Minor = 11;
  static const size_t Patch = 0;

  static std::string ToString() { return CC2_VERSION; }

  static std::string GitSha1() { return CC2_GIT_SHA1; }
};

} // namespace Consensus
} // namespace PacBio
