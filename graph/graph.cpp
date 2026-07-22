#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>

#include "utils/types.h"
#include "utils/utils.h"
#include "graph/graph.h"

Graph::Graph()
    : edge_count_(0)
    , vlabel_count_(0)
    , elabel_count_(0)
    , max_degree_(0)
    , max_label_frequency_(0)
    , neighbors_{}
    , elabels_{}
    , verterbylabel_{}
    , neighbors_by_label_{}
    , vlabels_{}
{
    verterbylabel_.clear();
    neighbors_by_label_.clear();
}

void Graph::AddVertex(uint id, uint label)
{
    if (id >= vlabels_.size())
    {
        vlabels_.resize(id + 1, NOT_EXIST);
        vlabels_[id] = label;

        neighbors_.resize(id + 1);
        elabels_.resize(id + 1);
    }
    else if (vlabels_[id] == NOT_EXIST)
    {
        vlabels_[id] = label;
    }

    vlabel_count_ = std::max(vlabel_count_, label + 1);

    verterbylabel_[label].push_back(id);
    max_label_frequency_ = std::max(
        max_label_frequency_,
        static_cast<uint>(verterbylabel_[label].size())
    );
}

void Graph::RemoveVertex(uint id)
{
    vlabels_[id] = NOT_EXIST;
    neighbors_[id].clear();
    elabels_[id].clear();

    // Label-aware index becomes stale after graph updates.
    // Rebuild it manually if dynamic updates are used.
    neighbors_by_label_.clear();
}

void Graph::AddEdge(uint v1, uint v2, uint label)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return;

    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].insert(lower, v2);
    elabels_[v1].insert(elabels_[v1].begin() + dis, label);

    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].insert(lower, v1);
    elabels_[v2].insert(elabels_[v2].begin() + dis, label);

    edge_count_++;
    elabel_count_ = std::max(elabel_count_, label + 1);
    max_degree_ = std::max(max_degree_, std::max(neighbors_[v1].size(), neighbors_[v2].size()));

    // Label-aware index becomes stale after graph updates.
    // LoadFromFile() will call BuildLabelNeighborIndex() after all edges are inserted.
    neighbors_by_label_.clear();
}

void Graph::RemoveEdge(uint v1, uint v2)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower == neighbors_[v1].end() || *lower != v2)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }

    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].erase(lower);
    elabels_[v1].erase(elabels_[v1].begin() + dis);

    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    if (lower == neighbors_[v2].end() || *lower != v1)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }

    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].erase(lower);
    elabels_[v2].erase(elabels_[v2].begin() + dis);

    edge_count_--;

    // Label-aware index becomes stale after graph updates.
    // Rebuild it manually if dynamic updates are used.
    neighbors_by_label_.clear();
}

uint Graph::GetVertexLabel(uint u) const
{
    return vlabels_[u];
}

const std::vector<uint>& Graph::GetNeighbors(uint v) const
{
    return neighbors_[v];
}

const std::vector<uint>& Graph::GetNeighborLabels(uint v) const
{
    return elabels_[v];
}

const std::vector<uint>& Graph::GetNeighborsByLabel(uint v, uint label) const
{
    static const std::vector<uint> empty;

    if (v >= neighbors_by_label_.size()) return empty;
    if (label >= neighbors_by_label_[v].size()) return empty;

    return neighbors_by_label_[v][label];
}

const std::vector<uint>& Graph::GetVerticesByLabel(uint label) const
{
    static const std::vector<uint> empty;

    auto it = verterbylabel_.find(label);
    if (it == verterbylabel_.end()) return empty;

    return it->second;
}

uint Graph::getGraphMaxLabelFrequency() const
{
    return max_label_frequency_;
}

bool Graph::checkEdgeExistence(uint v1, uint v2) const
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return true;
    return false;
}

std::tuple<uint, uint, uint> Graph::GetEdgeLabel(uint v1, uint v2) const
{
    uint v1_label, v2_label, e_label;
    v1_label = GetVertexLabel(v1);
    v2_label = GetVertexLabel(v2);

    const std::vector<uint>* nbrs;
    const std::vector<uint>* elabel;
    uint other;

    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        elabel = &elabels_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        elabel = &elabels_[v2];
        other = v1;
    }

    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;

        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_label = elabel->at(mid);
            return {v1_label, v2_label, e_label};
        }
    }

    return {v1_label, v2_label, static_cast<uint>(-1)};
}

uint Graph::GetDegree(uint v) const
{
    return neighbors_[v].size();
}

uint Graph::GetDiameter() const
{
    uint diameter = 0;

    for (uint i = 0u; i < NumVertices(); i++)
    {
        if (GetVertexLabel(i) != NOT_EXIST)
        {
            std::queue<uint> bfs_queue;
            std::vector<bool> visited(NumVertices(), false);
            uint level = UINT_MAX;

            bfs_queue.push(i);
            visited[i] = true;

            while (!bfs_queue.empty())
            {
                level++;

                uint size = bfs_queue.size();
                for (uint j = 0u; j < size; j++)
                {
                    uint front = bfs_queue.front();
                    bfs_queue.pop();

                    const auto& nbrs = GetNeighbors(front);
                    for (const uint nbr : nbrs)
                    {
                        if (!visited[nbr])
                        {
                            bfs_queue.push(nbr);
                            visited[nbr] = true;
                        }
                    }
                }
            }

            if (level > diameter) diameter = level;
        }
    }

    return diameter;
}

void Graph::BuildLabelNeighborIndex()
{
    // Sort and deduplicate vertices grouped by label.
    for (auto& kv : verterbylabel_)
    {
        std::vector<uint>& vertices = kv.second;
        std::sort(vertices.begin(), vertices.end());
        vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    }

    max_label_frequency_ = 0;
    for (const auto& kv : verterbylabel_)
    {
        max_label_frequency_ = std::max(
            max_label_frequency_,
            static_cast<uint>(kv.second.size())
        );
    }

    neighbors_by_label_.clear();
    neighbors_by_label_.resize(NumVertices());

    for (uint v = 0; v < NumVertices(); ++v)
    {
        neighbors_by_label_[v].resize(vlabel_count_);

        const std::vector<uint>& nbrs = neighbors_[v];

        for (uint nbr : nbrs)
        {
            uint nbr_label = GetVertexLabel(nbr);

            if (nbr_label == NOT_EXIST) continue;
            if (nbr_label >= vlabel_count_) continue;

            neighbors_by_label_[v][nbr_label].push_back(nbr);
        }
    }

    // Because neighbors_[v] is sorted by vertex id, every
    // neighbors_by_label_[v][label] is also sorted by vertex id.
}

void Graph::LoadFromFile(const std::string& path)
{
    size_t last_slash_pos = path.find_last_of('/');
    size_t last_point_pos = path.find_last_of('.');

    std::string extracted_str = path.substr(
        last_slash_pos + 1,
        last_point_pos - last_slash_pos - 1
    );

    SetName(extracted_str);

    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }

    std::ifstream ifs(path);

    char type;
    while (ifs >> type)
    {
        if (type == 't')
        {
            char temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label, degree;
            ifs >> vertex_id >> label >> degree;
            AddVertex(vertex_id, label);
        }
        else if (type == 'e')
        {
            uint from_id, to_id;
            ifs >> from_id >> to_id;
            AddEdge(from_id, to_id, 0);
        }
        else
        {
            // Ignore unknown line type.
        }
    }

    ifs.close();

    BuildLabelNeighborIndex();
}

void Graph::PrintMetaData() const
{
    std::cout << "# vertices = " << NumVertices() << std::endl;
    std::cout << "# edges = " << NumEdges() << std::endl;
    std::cout << "# max_degree =  " << max_degree_ << std::endl;
}