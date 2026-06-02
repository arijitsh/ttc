#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <algorithm>
#include <chrono>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

class Logger {
  int verbosity = 0;
  std::chrono::steady_clock::time_point start_time =
      std::chrono::steady_clock::now();

public:
  class LogStream {
    int level;
    int verbosity;
    const Logger& logger;
    bool newLine = true;

  public:
    LogStream(int lvl, int v, const Logger& log)
        : level(lvl), verbosity(v), logger(log) {}

    template <typename T>
    LogStream& operator<<(const T& msg)
    {
      if (verbosity >= level)
      {
        if (newLine)
        {
          // logger.printPrefix();
          newLine = false;
        }
        std::cerr << msg;
      }
      return *this;
    }
    LogStream& operator<<(std::ostream& (*manip)(std::ostream&))
    {
      if (verbosity >= level)
      {
        if (manip == static_cast<std::ostream& (*)(std::ostream&)>(std::endl))
        {
          std::cerr << manip;
          newLine = true;
        }
        else
        {
          std::cerr << manip;
        }
      }
      return *this;
    }
  };

  void setVerbosity(int v) { verbosity = v; }
  int getVerbosity() const { return verbosity; }
  LogStream operator()(int level) const { return LogStream(level, verbosity, *this); }

  double elapsed() const
  {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - start_time).count();
  }

  void printPrefix() const
  {
    auto t = elapsed();
    auto f = std::cerr.flags();
    auto p = std::cerr.precision();
    std::cerr << "c [ttc] [" << std::fixed << std::setprecision(3) << t
              << "] ";
    std::cerr.flags(f);
    std::cerr.precision(p);
  }
};

inline Logger Log;

class Tracer {
  using TraceSet = std::unordered_set<std::string>;

  static std::string normalize(std::string_view name)
  {
    std::string lowered(name);
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return lowered;
  }

  TraceSet enabled;

public:
  class TraceStream {
    bool active;
    std::string component;
    const Logger* logger;
    bool newLine = true;

    void printPrefix()
    {
      std::cerr << "c [" << component << "] ";
    }

  public:
    TraceStream(bool enabled, std::string name, const Logger* log)
        : active(enabled), component(std::move(name)), logger(log)
    {
    }

    template <typename T>
    TraceStream& operator<<(const T& msg)
    {
      if (active)
      {
        if (newLine)
        {
          printPrefix();
          newLine = false;
        }
        std::cerr << msg;
      }
      return *this;
    }

    TraceStream& operator<<(std::ostream& (*manip)(std::ostream&))
    {
      if (active)
      {
        if (manip == static_cast<std::ostream& (*)(std::ostream&)>(std::endl))
        {
          std::cerr << manip;
          newLine = true;
        }
        else
        {
          std::cerr << manip;
        }
      }
      return *this;
    }
  };

  void clear() { enabled.clear(); }

  void enable(std::string_view name) { enabled.insert(normalize(name)); }

  void setEnabled(const std::vector<std::string>& names)
  {
    enabled.clear();
    for (const auto& name : names)
    {
      enable(name);
    }
  }

  bool isEnabled(std::string_view name) const
  {
    return enabled.find(normalize(name)) != enabled.end();
  }

  TraceStream operator()(std::string_view name) const
  {
    return TraceStream(isEnabled(name), std::string(name), &Log);
  }
};

inline Tracer Trace;

#endif // LOGGER_HPP
