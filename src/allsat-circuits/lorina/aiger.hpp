#pragma once

#include <cstdint>
#include <string>

namespace lorina
{

enum class return_code
{
  success,
  parse_error
};

class aiger_reader
{
 public:
  virtual ~aiger_reader() = default;

  virtual void on_header(uint64_t, uint64_t, uint64_t, uint64_t, uint64_t) const {}
  virtual void on_header(uint64_t, uint64_t, uint64_t, uint64_t, uint64_t,
                         uint64_t, uint64_t, uint64_t, uint64_t) const {}
  virtual void on_input(uint32_t, uint32_t) const {}
  virtual void on_output(uint32_t, uint32_t) const {}
  virtual void on_and(uint32_t, uint32_t, uint32_t) const {}
};

inline return_code read_aiger(const std::string&, const aiger_reader&)
{
  return return_code::parse_error;
}

inline return_code read_ascii_aiger(const std::string&, const aiger_reader&)
{
  return return_code::parse_error;
}

}  // namespace lorina

