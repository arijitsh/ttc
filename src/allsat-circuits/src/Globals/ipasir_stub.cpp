#include "Globals/ipasir.h"

#include <stdexcept>

namespace
{
[[noreturn]] void ipasir_unavailable()
{
  throw std::runtime_error("IPASIR solver is not available in this build");
}
}  // namespace

extern "C" {

const char* ipasir_signature()
{
  return "ttc-ipasir-stub";
}

void* ipasir_init()
{
  ipasir_unavailable();
}

void ipasir_release(void*) {}

void ipasir_add(void*, int)
{
  ipasir_unavailable();
}

void ipasir_assume(void*, int)
{
  ipasir_unavailable();
}

int ipasir_solve(void*)
{
  ipasir_unavailable();
}

int ipasir_val(void*, int)
{
  ipasir_unavailable();
}

int ipasir_failed(void*, int)
{
  ipasir_unavailable();
}

void ipasir_set_terminate(void*, void*, int (*)(void*)) {}

void ipasir_set_learn(void*, void*, int, void (*)(void*, int*)) {}

}  // extern "C"

