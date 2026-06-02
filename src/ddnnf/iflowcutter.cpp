#include "iflowcutter.hpp"
#include "greedy_order.hpp"

namespace TWD {

IFlowCutter::IFlowCutter() : m_nodeCount(0) {}

void IFlowCutter::importGraph(int nodeCount, const std::vector<int>& tails, const std::vector<int>& heads)
{
    m_nodeCount = nodeCount;
    size_t arcCount = tails.size();
    m_tail = ArrayIDIDFunc(static_cast<int>(arcCount), m_nodeCount);
    m_head = ArrayIDIDFunc(static_cast<int>(arcCount), m_nodeCount);
    for (size_t i = 0; i < arcCount; ++i)
    {
        m_tail.set(static_cast<int>(i), tails[i]);
        m_head.set(static_cast<int>(i), heads[i]);
    }
}

ArrayIDIDFunc IFlowCutter::constructTD()
{
    return compute_greedy_min_degree_order(m_tail, m_head);
}

} // namespace TWD

