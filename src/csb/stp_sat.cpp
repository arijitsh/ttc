#include "stp_sat.hpp"

#include <optional>
#include <string>

#ifdef TTC_HAS_STP
#include <stp/cpp_interface.h>
#include <stp/Parser/parser.h>
#include <stp/Globals/Globals.h>

namespace csb {
std::optional<bool> stpCheckSat(const std::string& smtFormula) {
    using namespace stp;
    STPMgr mgr;
    GlobalParserBM = &mgr;
    Cpp_interface iface(mgr);
    iface.startup();
    GlobalParserInterface = &iface;
    STP stp(&mgr);
    GlobalSTP = &stp;
    mgr.UserFlags.smtlib2_parser_flag = true;
    mgr.UserFlags.print_output_flag = false;
    SMT2ScanString(smtFormula.c_str());
    if (SMT2Parse() != 0) {
        GlobalSTP = nullptr;
        GlobalParserInterface = nullptr;
        GlobalParserBM = nullptr;
        return std::nullopt;
    }
    bool isSat = !mgr.ValidFlag;
    GlobalSTP = nullptr;
    GlobalParserInterface = nullptr;
    GlobalParserBM = nullptr;
    return isSat;
}
} // namespace csb

#else
namespace csb {
std::optional<bool> stpCheckSat(const std::string&) {
    return std::nullopt;
}
} // namespace csb
#endif
