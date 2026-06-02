#ifndef PROFILER_HPP
#define PROFILER_HPP

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdint>
#include <sys/resource.h>
#include <sys/utsname.h>

#include "logger.hpp"

#define TTC_SECTION_WIDTH 70

inline std::string& banner_solver_info()
{
  static std::string info;
  return info;
}

inline void set_banner_solver_info(std::string info)
{
  banner_solver_info() = std::move(info);
}

inline void print_section(const std::string& name, bool pad = true)
{
  std::string line = "c --- [ " + name + " ] ";
  if (line.size() < TTC_SECTION_WIDTH)
    line += std::string(TTC_SECTION_WIDTH - line.size(), '-');
  std::cout << line << std::endl;
  if (pad) std::cout << "c" << std::endl;
}

inline void print_banner()
{
  print_section("banner");
  std::cout << "c ttc: Toolbox for Theory Counting" << std::endl;
  std::cout << "c" << std::endl;
  std::cout << "c Version 0.1.0" << std::endl;
  std::cout << "c g++ " << __VERSION__ << std::endl;
  std::time_t now = std::time(nullptr);
  char timebuf[64];
  std::strftime(timebuf, sizeof(timebuf), "%c", std::localtime(&now));
  struct utsname uts;
  uname(&uts);
  std::cout << "c " << timebuf << " " << uts.sysname << " " << uts.release
            << " " << uts.machine << std::endl;
  const std::string& solverInfo = banner_solver_info();
  if (!solverInfo.empty())
  {
    std::cout << "c " << solverInfo << std::endl;
  }
  std::cout << "c" << std::endl;
}
class Profiler {
  double parse_time = 0.0;
  double search_time = 0.0;
  double rewrite_time = 0.0;
  double preprocess_time = 0.0;
  double tree_decomp_time = 0.0;
  double component_time = 0.0;
  double decompose_time = 0.0;
  double residual_substitute_time = 0.0;
  double residual_simplify_time = 0.0;
  double dnfization_time = 0.0;
  double polytope_volume_time = 0.0;
  double sampling_time = 0.0;
  std::uint64_t cache_hits = 0;
  std::uint64_t cache_lookups = 0;
  double cache_time = 0.0;
  std::size_t cache_memory = 0;
  std::uint64_t decompositions = 0;
  std::uint64_t stopped_decompositions = 0;
  std::uint64_t monotone_clauses = 0;
  std::uint64_t checksat_calls = 0;
  std::uint64_t checksat_unsat = 0;
  std::uint64_t sat_calls = 0;
  std::uint64_t sat_unsat = 0;
  std::uint64_t decisions = 0;
  std::size_t samples_generated = 0;
  std::size_t samples_deleted = 0;

public:
  void addParse(double t) { parse_time += t; }
  void addSearch(double t) { search_time += t; }
  void addRewrite(double t) { rewrite_time += t; }
  void addPreprocess(double t) { preprocess_time += t; }
  void addTreeDecomp(double t) { tree_decomp_time += t; }
  void addComponent(double t) { component_time += t; }
  void addDecompose(double t) { decompose_time += t; }
  void addDnfization(double t) { dnfization_time += t; }
  void addPolytopeVolume(double t) { polytope_volume_time += t; }
  void addSampling(double t) { sampling_time += t; }
  void addCacheHit() { ++cache_hits; }
  void addCacheLookup() { ++cache_lookups; }
  void addCacheTime(double t) { cache_time += t; }
  void setCacheMemory(std::size_t m) { cache_memory = m; }
  void setCacheStats(std::uint64_t lookups, std::uint64_t hits)
  {
    cache_lookups = lookups;
    cache_hits = hits;
  }
  void addDecomposition() { ++decompositions; }
  void addStoppedDecomposition() { ++stopped_decompositions; }
  void addMonotoneClause() { ++monotone_clauses; }
  void addCheckSatCall() { ++checksat_calls; }
  void addUnsatCheckSatCall() { ++checksat_unsat; }
  void addSatCall() { ++sat_calls; }
  void addUnsatSatCall() { ++sat_unsat; }
  void addDecision() { ++decisions; }
  void addResidualSubstitute(double t) { residual_substitute_time += t; }
  void addResidualSimplify(double t) { residual_simplify_time += t; }
  void setSampleStats(std::size_t generated, std::size_t deleted)
  {
    samples_generated = generated;
    samples_deleted = deleted;
  }
  std::uint64_t getDecompositions() const { return decompositions; }
  std::uint64_t getStoppedDecompositions() const {
    return stopped_decompositions;
  }
  std::uint64_t getMonotoneClauses() const { return monotone_clauses; }
  std::uint64_t getCheckSatCalls() const { return checksat_calls; }
  std::uint64_t getUnsatCheckSatCalls() const { return checksat_unsat; }
  std::uint64_t getSatCalls() const { return sat_calls; }
  std::uint64_t getUnsatSatCalls() const { return sat_unsat; }
  std::uint64_t getDecisions() const { return decisions; }
  void resetDecompositions() { decompositions = 0; }
  void resetSearchStats() {
    decompositions = 0;
    stopped_decompositions = 0;
    monotone_clauses = 0;
    checksat_calls = 0;
    checksat_unsat = 0;
    sat_calls = 0;
    sat_unsat = 0;
    decisions = 0;
    component_time = 0.0;
    decompose_time = 0.0;
    residual_substitute_time = 0.0;
    residual_simplify_time = 0.0;
    dnfization_time = 0.0;
    polytope_volume_time = 0.0;
    sampling_time = 0.0;
    samples_generated = 0;
    samples_deleted = 0;
  }

  double getParse() const { return parse_time; }
  double getSearch() const { return search_time; }
  double getRewrite() const { return rewrite_time; }
  double getPreprocess() const { return preprocess_time; }
  double getTreeDecomp() const { return tree_decomp_time; }
  double getComponent() const { return component_time; }
  double getDecompose() const { return decompose_time; }
  double getResidualSubstitute() const { return residual_substitute_time; }
  double getResidualSimplify() const { return residual_simplify_time; }
  double getDnfization() const { return dnfization_time; }
  double getPolytopeVolume() const { return polytope_volume_time; }
  double getSampling() const { return sampling_time; }
  std::uint64_t getCacheHits() const { return cache_hits; }
  std::uint64_t getCacheLookups() const { return cache_lookups; }
  double getCacheTime() const { return cache_time; }
  std::size_t getCacheMemory() const { return cache_memory; }
  std::size_t getSamplesGenerated() const { return samples_generated; }
  std::size_t getSamplesDeleted() const { return samples_deleted; }

  void print(double total_process) const
  {
    (void)total_process;
    double solve_time =
        parse_time + dnfization_time + polytope_volume_time + sampling_time;
    if (solve_time <= 0) solve_time = 1e-9; // avoid division by zero

    print_section("run-time profiling");
    std::cout << "c process time taken by individual solving procedures" << std::endl;
    std::cout << "c (percentage relative to process time for solving)" << std::endl;
    std::cout << "c" << std::endl;

    auto printLine = [&](double t, double pct, const char* name) {
      std::cout << "c " << std::setw(10) << std::fixed << std::setprecision(2) << t
                << "   " << std::setw(6) << std::fixed << std::setprecision(2) << pct
                << "% " << name << std::endl;
    };

    printLine(parse_time, (parse_time / solve_time) * 100.0, "parse");
    printLine(dnfization_time, (dnfization_time / solve_time) * 100.0, "DNFization");
    printLine(polytope_volume_time,
              (polytope_volume_time / solve_time) * 100.0,
              "polytope volume");
    printLine(sampling_time, (sampling_time / solve_time) * 100.0, "sampling");
    std::cout << "c =================================" << std::endl;
    printLine(solve_time, (solve_time / total_process) * 100.0, "solve");
    std::cout << "c" << std::endl;
    print_section("statistics");
    std::cout << "c total samples generated: " << samples_generated << std::endl;
    std::cout << "c total samples deleted: " << samples_deleted << std::endl;
    std::cout << "c" << std::endl;
  }
};

inline Profiler Profile;

inline double progress_start_time = 0.0;
inline double progress_next_percent = 0.01;
inline double progress_last_time = 0.0;
inline double progress_interval = 5.0;

inline void start_progress(double start)
{
  progress_start_time = start;
  progress_next_percent = 0.01;
  progress_last_time = 0.0;
}

inline void print_progress(double percent, std::size_t cache_size = 0)
{
  if (Log.getVerbosity() == 0) return;
  double elapsed = Log.elapsed() - progress_start_time;
  bool doPrint = false;
  if (elapsed - progress_last_time >= progress_interval)
  {
    progress_last_time = elapsed;
    doPrint = true;
  }
  if (!doPrint) return;
  std::uint64_t lookups = Profile.getCacheLookups();
  std::uint64_t hits = Profile.getCacheHits();
  std::uint64_t decomps = Profile.getDecompositions();
  std::uint64_t mono = Profile.getMonotoneClauses();
  std::uint64_t checks = Profile.getCheckSatCalls();
  std::uint64_t unsat_checks = Profile.getUnsatCheckSatCalls();
  std::uint64_t decs = Profile.getDecisions();
  double hitPct = lookups ? (static_cast<double>(hits) / lookups) * 100.0 : 0.0;
  std::cout << "c " << std::fixed << std::setprecision(3) << elapsed
            << "  " << std::fixed << std::setprecision(5) << percent
            << "  " << cache_size
            << "  " << lookups
            << "  " << std::fixed << std::setprecision(2) << hitPct
            << "  " << decomps
            << "  " << mono
            << "  " << checks << "(" << unsat_checks << ")"
            << "  " << decs
            << std::endl;
}

inline void print_progress_d4(double percent)
{
  if (Log.getVerbosity() == 0) return;
  double elapsed = Log.elapsed() - progress_start_time;
  bool doPrint = false;
  if (elapsed - progress_last_time >= progress_interval)
  {
    progress_last_time = elapsed;
    doPrint = true;
  }
  if (!doPrint) return;
  (void)percent;
  std::uint64_t hits = Profile.getCacheHits();
  std::uint64_t lookups = Profile.getCacheLookups();
  std::uint64_t miss = lookups > hits ? lookups - hits : 0;
  std::uint64_t smt = Profile.getCheckSatCalls();
  std::uint64_t smt_unsat = Profile.getUnsatCheckSatCalls();
  std::uint64_t decs = Profile.getDecisions();
  std::uint64_t decomps = Profile.getDecompositions();
  std::uint64_t stopped = Profile.getStoppedDecompositions();
  std::cout << "c " << std::fixed << std::setprecision(3) << elapsed
            << "  " << hits << "/" << miss
            << "  " << smt << "(" << smt_unsat << ")"
            << "  " << decs
            << "  " << decomps << "(" << stopped << ")"
            << std::endl;
}

inline void print_resources(double total_process)
{
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  double maxrss = usage.ru_maxrss / 1024.0;
  print_section("resources");
  std::cout << "c total process time since initialization: "
            << std::fixed << std::setprecision(2) << total_process
            << "    seconds" << std::endl;
  std::cout << "c total real time since initialization: "
            << std::fixed << std::setprecision(2) << total_process
            << "    seconds" << std::endl;
  std::cout << "c maximum resident set size of process: "
            << std::fixed << std::setprecision(2) << maxrss
            << "    MB" << std::endl;
  std::cout << "c" << std::endl;
}

#endif // PROFILER_HPP
