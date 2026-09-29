#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <variant>
#include <vector>

#include "graph/graph.h"
#include "utils/types.h"

class ICTCache
{
public:
    class RawRange
    {
    public:
        using const_iterator = std::vector<VertexID>::const_iterator;

        RawRange(const_iterator begin, const_iterator end)
            : begin_(begin), end_(end)
        {}

        explicit RawRange(const std::vector<VertexID>& values)
            : begin_(values.begin()), end_(values.end())
        {}

        const_iterator begin() const { return begin_; }
        const_iterator end() const { return end_; }
        size_t size() const { return static_cast<size_t>(end_ - begin_); }
        bool empty() const { return begin_ == end_; }

    private:
        const_iterator begin_;
        const_iterator end_;
    };

    struct MemoryStats
    {
        size_t active_nodes = 0;
        size_t label_roots = 0;
        size_t pool_capacity_nodes = 0;
        size_t child_edges = 0;
        size_t linear_child_capacity = 0;
        size_t hashed_entries = 0;
        size_t hash_buckets = 0;
        size_t intersection_elements = 0;
        size_t intersection_capacity_elements = 0;

        size_t IntersectionBytes() const { return intersection_elements * sizeof(VertexID); }
        size_t IntersectionCapacityBytes() const { return intersection_capacity_elements * sizeof(VertexID); }

    private:
        friend class ICTCache;
        size_t node_size_bytes = 0;

    public:
        size_t ActiveNodeStorageBytes() const { return active_nodes * node_size_bytes; }
        size_t PoolReservedBytes() const { return pool_capacity_nodes * node_size_bytes; }
        size_t EstimatedHashBytes() const
        {
            return hash_buckets * sizeof(void*) +
                hashed_entries * (sizeof(VertexID) + sizeof(void*) + 2 * sizeof(void*));
        }
        size_t LinearChildCapacityBytes() const
        {
            return linear_child_capacity * (sizeof(VertexID) + sizeof(void*));
        }
        size_t EstimatedTotalReservedBytes() const
        {
            return PoolReservedBytes() + IntersectionCapacityBytes() +
                EstimatedHashBytes() + LinearChildCapacityBytes();
        }
    };

private:
    static constexpr VertexID INVALID_VERTEX = static_cast<VertexID>(-1);

    struct TreeNode;

    struct ChildStore
    {
        // Favor lookup speed: only very small fan-outs stay linear.
        static constexpr size_t HASH_THRESHOLD = 4;
        using Entry = std::pair<VertexID, TreeNode*>;

        std::vector<Entry> linear;
        std::unique_ptr<std::unordered_map<VertexID, TreeNode*>> hashed;

        TreeNode* Find(VertexID vid) const
        {
            if (hashed)
            {
                auto it = hashed->find(vid);
                return it == hashed->end() ? nullptr : it->second;
            }
            for (const Entry& entry : linear)
            {
                if (entry.first == vid) return entry.second;
            }
            return nullptr;
        }

        void Insert(VertexID vid, TreeNode* node)
        {
            if (hashed)
            {
                hashed->emplace(vid, node);
                return;
            }

            linear.emplace_back(vid, node);
            if (linear.size() <= HASH_THRESHOLD) return;

            auto table = std::make_unique<std::unordered_map<VertexID, TreeNode*>>();
            table->reserve(linear.size() * 2);
            for (const Entry& entry : linear)
            {
                table->emplace(entry.first, entry.second);
            }
            std::vector<Entry>().swap(linear);
            hashed = std::move(table);
        }

        void Clear()
        {
            linear.clear();
            hashed.reset();
        }

        bool Empty() const { return hashed ? hashed->empty() : linear.empty(); }
        size_t Size() const { return hashed ? hashed->size() : linear.size(); }
        size_t LinearCapacity() const { return hashed ? 0 : linear.capacity(); }
        size_t HashedEntries() const { return hashed ? hashed->size() : 0; }
        size_t HashBuckets() const { return hashed ? hashed->bucket_count() : 0; }
    };

    struct TreeNode
    {
        struct BorrowedRaw
        {
            RawRange::const_iterator begin;
            RawRange::const_iterator end;
        };

        VertexID vid = INVALID_VERTEX;
        bool built = false;

        // First-level nodes borrow the immutable label-neighbor range from the
        // data graph. Deeper nodes own their materialized intersection.
        std::variant<std::vector<VertexID>, BorrowedRaw> raw;

        // child key = next dependency data vertex.
        ChildStore children;

        void Reset(VertexID v)
        {
            vid = v;
            built = false;
            if (auto* owned = std::get_if<std::vector<VertexID>>(&raw))
                owned->clear();
            else
                raw = std::vector<VertexID>();
            children.Clear();
        }

        RawRange Raw() const
        {
            if (const auto* borrowed = std::get_if<BorrowedRaw>(&raw))
                return RawRange(borrowed->begin, borrowed->end);
            return RawRange(std::get<std::vector<VertexID>>(raw));
        }

        void SetBorrowedRaw(RawRange::const_iterator begin, RawRange::const_iterator end)
        {
            raw = BorrowedRaw{begin, end};
        }

        std::vector<VertexID>& OwnedRaw()
        {
            if (!std::holds_alternative<std::vector<VertexID>>(raw))
                raw = std::vector<VertexID>();
            return std::get<std::vector<VertexID>>(raw);
        }

        const std::vector<VertexID>* OwnedRawIfPresent() const
        {
            return std::get_if<std::vector<VertexID>>(&raw);
        }
    };

    struct NodePool
    {
        static constexpr size_t CHUNK_SIZE = 4096;

        std::vector<std::unique_ptr<TreeNode[]>> chunks;
        size_t chunk_id = 0;
        size_t offset = 0;

        TreeNode* Alloc(VertexID vid)
        {
            if (chunks.empty())
            {
                chunks.emplace_back(std::make_unique<TreeNode[]>(CHUNK_SIZE));
            }

            if (offset == CHUNK_SIZE)
            {
                ++chunk_id;
                offset = 0;

                if (chunk_id == chunks.size())
                {
                    chunks.emplace_back(std::make_unique<TreeNode[]>(CHUNK_SIZE));
                }
            }

            TreeNode* node = &chunks[chunk_id][offset++];
            node->Reset(vid);
            return node;
        }

        void Clear()
        {
            chunks.clear();
            chunk_id = 0;
            offset = 0;
        }
    };

private:
    // roots_[label] is the root of the ICT for this label.
    std::unordered_map<LabelID, TreeNode*> roots_;

    NodePool pool_;

private:
    template <typename LeftRange, typename RightRange>
    static void IntersectSorted(
        const LeftRange& left,
        const RightRange& right,
        std::vector<VertexID>& output)
    {
        if (left.empty() || right.empty()) return;

        const size_t left_size = left.size();
        const size_t right_size = right.size();
        const size_t max_output_size = std::min(left_size, right_size);
        auto append = [&](VertexID value)
        {
            if (output.empty())
            {
                // Do not allocate for disjoint inputs.  A modest first reserve
                // avoids both repeated growth and large unused capacities.
                output.reserve(std::min<size_t>(max_output_size, 1024));
            }
            output.push_back(value);
        };

        // Binary probing wins for highly skewed inputs and preserves the
        // ascending output order.  Balanced inputs use the linear merge.
        if (left_size < right_size && right_size / left_size >= 32)
        {
            for (VertexID value : left)
            {
                if (std::binary_search(right.begin(), right.end(), value))
                    append(value);
            }
        }
        else if (right_size < left_size && left_size / right_size >= 32)
        {
            for (VertexID value : right)
            {
                if (std::binary_search(left.begin(), left.end(), value))
                    append(value);
            }
        }
        else
        {
            auto left_it = left.begin();
            auto right_it = right.begin();
            while (left_it != left.end() && right_it != right.end())
            {
                if (*left_it < *right_it)
                {
                    ++left_it;
                }
                else if (*right_it < *left_it)
                {
                    ++right_it;
                }
                else
                {
                    append(*left_it);
                    ++left_it;
                    ++right_it;
                }
            }
        }
    }

    TreeNode* NewNode(VertexID vid)
    {
        ++m_newNodeCount;
        return pool_.Alloc(vid);
    }

    TreeNode* GetRoot(LabelID label)
    {
        auto it = roots_.find(label);
        if (it != roots_.end())
        {
            return it->second;
        }

        TreeNode* root = NewNode(INVALID_VERTEX);
        root->built = true;
        roots_.emplace(label, root);
        return root;
    }

public:
    size_t m_getCount = 0;
    size_t m_hitCount = 0;
    size_t m_reuseCount = 0;
    size_t m_reusedPrefixLength = 0;
    size_t m_newNodeCount = 0;
    size_t m_buildCount = 0;
    size_t m_intersectionElementCount = 0;
    size_t m_intersectionCapacityCount = 0;
    size_t m_childEdgeCount = 0;

public:
    ICTCache() = default;

    void Clear()
    {
        roots_.clear();
        pool_.Clear();

        m_getCount = 0;
        m_hitCount = 0;
        m_reuseCount = 0;
        m_reusedPrefixLength = 0;
        m_newNodeCount = 0;
        m_buildCount = 0;
        m_intersectionElementCount = 0;
        m_intersectionCapacityCount = 0;
        m_childEdgeCount = 0;
    }

    RawRange GetRaw(
        LabelID label,
        std::vector<VertexID> depVertices,
        const Graph& dataGraph)
    {
        std::sort(depVertices.begin(), depVertices.end());
        depVertices.erase(std::unique(depVertices.begin(), depVertices.end()), depVertices.end());
        return GetRawCanonical(label, depVertices.data(), depVertices.size(), dataGraph);
    }

    // Fast path for the matcher: keys must be strictly increasing.
    RawRange GetRawCanonical(
        LabelID label,
        const VertexID* keys,
        size_t key_count,
        const Graph& dataGraph)
    {
        ++m_getCount;

        if (key_count == 0)
        {
            return RawRange(dataGraph.GetVerticesByLabel(label));
        }

#ifndef NDEBUG
        for (size_t i = 1; i < key_count; ++i)
        {
            assert(keys[i - 1] < keys[i]);
        }
#endif

        bool full_hit = true;
        size_t reused_prefix_length = 0;

        TreeNode* cur = GetRoot(label);

        for (size_t key_index = 0; key_index < key_count; ++key_index)
        {
            const VertexID vid = keys[key_index];
            TreeNode* nxt = nullptr;

            nxt = cur->children.Find(vid);
            if (nxt != nullptr)
            {
                ++reused_prefix_length;
            }
            else
            {
                full_hit = false;
                nxt = NewNode(vid);
                cur->children.Insert(vid, nxt);
                ++m_childEdgeCount;
            }

            if (!nxt->built)
            {
                full_hit = false;
                ++m_buildCount;

                if (cur->vid == INVALID_VERTEX)
                {
                    const Graph::NeighborRange labelNbrs =
                        dataGraph.GetNeighborsByLabel(vid, label);
                    nxt->SetBorrowedRaw(labelNbrs.begin(), labelNbrs.end());
                }
                else
                {
                    const RawRange parent_raw = cur->Raw();
                    if (!parent_raw.empty())
                    {
                        const Graph::NeighborRange labelNbrs =
                            dataGraph.GetNeighborsByLabel(vid, label);
                        IntersectSorted(
                            parent_raw,
                            labelNbrs,
                            nxt->OwnedRaw()
                        );
                    }
                }

                nxt->built = true;
                if (const auto* owned = nxt->OwnedRawIfPresent())
                {
                    m_intersectionElementCount += owned->size();
                    m_intersectionCapacityCount += owned->capacity();
                }
            }

            cur = nxt;
        }

        if (full_hit)
        {
            ++m_hitCount;
        }
        if (reused_prefix_length != 0)
        {
            ++m_reuseCount;
            m_reusedPrefixLength += reused_prefix_length;
        }

        return cur->Raw();
    }

    MemoryStats GetMemoryStats() const
    {
        MemoryStats stats;
        stats.node_size_bytes = sizeof(TreeNode);
        stats.active_nodes = m_newNodeCount;
        stats.label_roots = roots_.size();
        stats.pool_capacity_nodes = pool_.chunks.size() * NodePool::CHUNK_SIZE;
        stats.child_edges = m_childEdgeCount;
        stats.intersection_elements = m_intersectionElementCount;
        stats.intersection_capacity_elements = m_intersectionCapacityCount;

        stats.hash_buckets = roots_.empty() ? 0 : roots_.bucket_count();
        stats.hashed_entries = roots_.size();
        size_t remaining = stats.active_nodes;
        for (const auto& chunk : pool_.chunks)
        {
            const size_t active_in_chunk = std::min(remaining, NodePool::CHUNK_SIZE);
            for (size_t i = 0; i < active_in_chunk; ++i)
            {
                const TreeNode* node = &chunk[i];
                stats.linear_child_capacity += node->children.LinearCapacity();
                stats.hashed_entries += node->children.HashedEntries();
                stats.hash_buckets += node->children.HashBuckets();
            }
            remaining -= active_in_chunk;
            if (remaining == 0) break;
        }
        return stats;
    }
};
