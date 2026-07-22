#ifndef GRAPH_GRAPH
#define GRAPH_GRAPH

#include <queue>
#include <tuple>
#include <vector>
#include <string>
#include <unordered_map>
#include <bitset>

#include "utils/types.h"
#include "utils/utils.h"

class Graph
{
protected:
    std::string name_;
    uint edge_count_;
    uint vlabel_count_;
    uint elabel_count_;
    size_t max_degree_;
    uint max_label_frequency_;

    // neighbors_[v] = all neighbors of vertex v, sorted by vertex id.
    std::vector<std::vector<uint>> neighbors_;

    // elabels_[v][i] = edge label of edge (v, neighbors_[v][i]).
    std::vector<std::vector<uint>> elabels_;

    // verterbylabel_[label] = all vertices with this label.
    // Keep the original variable name to avoid changing too many places.
    std::unordered_map<uint, std::vector<uint>> verterbylabel_;

    // neighbors_by_label_[v][label] = neighbors of v whose vertex label is label.
    std::vector<std::vector<std::vector<uint>>> neighbors_by_label_;

public:
    std::vector<uint> vlabels_;

public:
    Graph();

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }

    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;
    uint GetMaxDegree() const { return max_degree_; }

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);

    void AddEdge(uint v1, uint v2, uint label);
    void RemoveEdge(uint v1, uint v2);

    void SetName(const std::string& name) { name_ = name; }
    std::string GetName() const { return name_; }

    uint GetVertexLabel(uint u) const;

    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;

    // New interface:
    // Return neighbors of v whose vertex label is label.
    const std::vector<uint>& GetNeighborsByLabel(uint v, uint label) const;

    uint GetDegree(uint v) const;

    // Return all vertices with the given label.
    const std::vector<uint>& GetVerticesByLabel(uint label) const;

    uint getGraphMaxLabelFrequency() const;

    bool checkEdgeExistence(uint v1, uint v2) const;

    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    // Build label-aware indexes after the graph is loaded.
    void BuildLabelNeighborIndex();

    void LoadFromFile(const std::string& path);
    void PrintMetaData() const;
};

#endif // GRAPH_GRAPH