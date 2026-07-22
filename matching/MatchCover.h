#ifndef MATCHING_MatCo
#define MATCHING_MatCo

#include <vector>

#include "utils/types.h"
#include "graph/graph.h"
#include "matching/matching.h"
#include <unordered_map>
#include "utils/timer.hpp"
#include <unordered_set>
#include <set>
#include <fstream>
#include <sstream>

class MatCo : public matching
{
private:

    /// <summary>
    /// 第一维下标是深度，代表深度为index的部分匹配gi
    /// 第二维下标也是深度，
    /// CandidateSets[i][j],表示部分匹配gi下，深度为j的顶点的候选集,相当于论文中的Cgi(j);
    /// </summary>
    std::vector<std::vector<std::vector<uint>>> CandidateSets; 
    std::vector<std::vector<bool>> CandidateSetFlag; 
    std::vector<std::vector<bool>> adjacency_matrix ; 
    std::vector<bool> KeyVertexSet;
    std::vector<uint> match_order;
    uint mutiexp_depth_;
    uint prune_depth_;

public:
    MatCo(Graph& query_graph, Graph& data_graph, size_t max_num_results,
            bool print_prep, bool print_enum, bool homo);
    ~MatCo() override {};

    void Preprocessing() override;
    void InitialMatching() override;
    void PrintKeyVertexSet();
    bool VerifyCorrectness(const std::string &kv_path);
    


private:
    void GenerateMatchingOrder();
    void BuildCP2LEOrder();
    void BuildAdjMatrix();

    void BuildCover();
    
    bool ComputeCand(uint depth,std::vector<uint> m);
    bool FullCoveragePrune(const uint depth,const std::vector<uint> &m);

    void MutiExpansion(std::vector<uint> m);
    bool MutiExpTest(int depth,const std::vector<uint> &label_same_index, std::vector<uint> &m);

    void CountRes(std::vector<uint> m);
    void FlushFlag(uint flush_depth);


    void PrintMatch(const std::vector<uint> &m);

    
    
    void FindMatCo(uint depth, std::vector<uint> m);
    
    
    
};

#endif 
