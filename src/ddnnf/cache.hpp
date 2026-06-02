#ifndef CACHE_HPP
#define CACHE_HPP

#include <cvc5/cvc5.h>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>


struct DDNNFCacheKey
{
    cvc5::Term formula;
    std::vector<cvc5::Term> vars;
    std::vector<std::pair<std::size_t, bool>> assignments;
    std::vector<std::size_t> constraints;
    bool operator==(const DDNNFCacheKey& other) const
    {
        return formula == other.formula && vars == other.vars &&
               assignments == other.assignments && constraints == other.constraints;
    }
};

struct DDNNFCacheKeyHash
{
    std::size_t operator()(const DDNNFCacheKey& k) const
    {
        std::size_t h = std::hash<cvc5::Term>{}(k.formula);
        for (const auto& v : k.vars)
        {
            h ^= std::hash<cvc5::Term>{}(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        for (const auto& a : k.assignments)
        {
            std::size_t ah = std::hash<std::size_t>{}(a.first);
            ah ^= std::hash<bool>{}(a.second) + 0x9e3779b9 + (ah << 6) + (ah >> 2);
            h ^= ah + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        for (const auto& c : k.constraints)
        {
            h ^= std::hash<std::size_t>{}(c) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

bool termsEquivalent(cvc5::Solver& solver, const cvc5::Term& a, const cvc5::Term& b);

template <typename Key, typename Hash>
class ComponentCache
{
public:
    using Equiv = std::function<bool(const Key&, const Key&)>;

    ComponentCache(bool useHash, Equiv equiv)
        : d_useHash(useHash), d_equiv(equiv)
    {}

    void clear()
    {
        if (d_useHash)
        {
            d_hashCache.clear();
        }
        else
        {
            d_structCache.clear();
        }
    }

    bool lookup(const Key& key, std::uint64_t& value)
    {
        if (d_useHash)
        {
            auto it = d_hashCache.find(key);
            if (it != d_hashCache.end())
            {
                value = it->second;
                return true;
            }
            return false;
        }
        for (const auto& entry : d_structCache)
        {
            if (d_equiv(entry.key, key))
            {
                value = entry.value;
                return true;
            }
        }
        return false;
    }

    void insert(const Key& key, std::uint64_t value)
    {
        if (d_useHash)
        {
            d_hashCache[key] = value;
        }
        else
        {
            d_structCache.push_back({key, value});
        }
    }

    std::size_t memoryUsage() const
    {
        if (d_useHash)
        {
            return d_hashCache.size() * (sizeof(Key) + sizeof(std::uint64_t));
        }
        return d_structCache.size() * sizeof(Entry);
    }

    std::size_t size() const
    {
        if (d_useHash)
        {
            return d_hashCache.size();
        }
        return d_structCache.size();
    }

private:
    struct Entry
    {
        Key key;
        std::uint64_t value;
    };

    bool d_useHash;
    Equiv d_equiv;
    std::vector<Entry> d_structCache;
    std::unordered_map<Key, std::uint64_t, Hash> d_hashCache;
};

#endif // CACHE_HPP
