#pragma once

#include <vector>
#include "array_id_func.hpp"

namespace TWD {
class IFlowCutter {
public:
    IFlowCutter();
    void importGraph(int nodeCount, const std::vector<int>& tails, const std::vector<int>& heads);
    ArrayIDIDFunc constructTD();
private:
    int m_nodeCount;
    ArrayIDIDFunc m_tail;
    ArrayIDIDFunc m_head;
};
}

