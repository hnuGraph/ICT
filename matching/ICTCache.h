#pragma once

#include <algorithm>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "utils/types.h"

class ICTCache
{
private:
    static constexpr VertexID INVALID_VERTEX = static_cast<VertexID>(-1);

    struct TreeNode
    {
        VertexID vid = INVALID_VERTEX;
        bool built = false;

        // The intersection result represented by the path from the label root to this node.
        std::vector<VertexID> raw_intersection;

        // child key = next dependency data vertex.
        std::unordered_map<VertexID, TreeNode*> children;

        void Reset(VertexID v)
        {
            vid = v;
            built = false;
            raw_intersection.clear();
            children.clear();
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

        void ResetAllocation()
        {
            chunk_id = 0;
            offset = 0;
        }
    };

private:
    // roots_[label] is the root of the ICT for this label.
    std::unordered_map<LabelID, TreeNode*> roots_;

    NodePool pool_;

private:
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
    size_t m_newNodeCount = 0;
    size_t m_buildCount = 0;

public:
    ICTCache() = default;

    void Clear()
    {
        roots_.clear();
        pool_.ResetAllocation();

        m_getCount = 0;
        m_hitCount = 0;
        m_newNodeCount = 0;
        m_buildCount = 0;
    }

    const std::vector<VertexID>& GetRaw(
        LabelID label,
        const std::vector<VertexID>& depVertices,
        const Graph& dataGraph)
    {
        ++m_getCount;

        // If there is no dependency vertex, the result is all vertices with this label.
        // DsqlMeta usually handles this case before calling ICT.
        if (depVertices.empty())
        {
            return dataGraph.GetVerticesByLabel(label);
        }

        std::vector<VertexID> keys = depVertices;
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

        bool full_hit = true;

        TreeNode* cur = GetRoot(label);

        for (VertexID vid : keys)
        {
            TreeNode* nxt = nullptr;

            auto it = cur->children.find(vid);
            if (it != cur->children.end())
            {
                nxt = it->second;
            }
            else
            {
                full_hit = false;
                nxt = NewNode(vid);
                cur->children.emplace(vid, nxt);
            }

            if (!nxt->built)
            {
                full_hit = false;
                ++m_buildCount;

                const std::vector<VertexID>& labelNbrs =
                    dataGraph.GetNeighborsByLabel(vid, label);

                if (cur->vid == INVALID_VERTEX)
                {
                    nxt->raw_intersection.assign(labelNbrs.begin(), labelNbrs.end());
                }
                else
                {
                    std::set_intersection(
                        cur->raw_intersection.begin(),
                        cur->raw_intersection.end(),
                        labelNbrs.begin(),
                        labelNbrs.end(),
                        std::back_inserter(nxt->raw_intersection)
                    );
                }

                nxt->built = true;
            }

            cur = nxt;
        }

        if (full_hit)
        {
            ++m_hitCount;
        }

        return cur->raw_intersection;
    }
};