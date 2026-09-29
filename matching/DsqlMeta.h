#pragma once
#include "matching.h"
#include "unordered_set"

class ICTCache;

class DsqlMeta : public matching
{

public:
    // �Ż�
    struct OPT
    {
    public:
        OPT() = default;
        /* �Ż��� */
        std::vector<uint> labelRM;
        std::vector<uint> neighborRM;


        /* �Ż��� */
        // �±���depth��Ԫ�ر�����������Ƿ�ƥ��ʧ�ܣ����Ҳ���һ��ƥ�䶥��
        std::vector<bool> matchSuc; 
        // ��ά�±궼�ǲ�ѯ����ID��Ԫ�ر����ڶ�άID�������һάID�����Ƿ��Ƕ�̬��ͻ
        std::vector<std::vector<bool>> dynamicConflictTable;
        int failDepth = NOFAILNODE;
        size_t skipQueryCount = 0;

        
        /* �Ż��� */
        // �±���depth��badVertices[depth] �������depth+1�� "bad" vertices
        std::vector<std::unordered_set<uint>> badVertices;
        size_t skipDataVertexCount = 0;


    };
    OPT opt;

    void InitOptSESM();
    void SetOptSESM(uint depth,std::vector<uint>& cand);
    void InitOptSKNB(uint depth);
    void SetOptDynamicCT(uint depth, const std::vector<uint>& m, uint v);
    bool CheckOptSkipNode(uint curDepth);
    void ClearOptBadVertices(int depth);
    void ClearOptAllBadVertices();
    void SetOptBadVertices(uint depth, uint v);
    bool CheckOptSkipDataVertice(uint depth, uint v);


private:

    /// <summary>
    /// ��һά�±�����ȣ��������Ϊindex�Ĳ���ƥ��gi
    /// �ڶ�ά�±�Ҳ����ȣ�
    /// CandidateSets[i][j],��ʾ����ƥ��gi�£����Ϊj�Ķ���ĺ�ѡ��,�൱�������е�Cgi(j);
    /// </summary>
    std::vector<std::vector<std::vector<uint>>> CandidateSets;
    std::vector<std::vector<bool>> CandidateSetFlag;
    std::vector<std::vector<bool>> adjacency_matrix;
std::vector<bool> KeyVertexSet;
std::vector<uint> match_order;
uint mutiexp_depth_ = 0;
uint prune_depth_ = 0;

// Shared intersection cache tree.
// It is owned by the caller, not by DsqlMeta.
ICTCache* ict_cache_ = nullptr;
bool use_ict_cache_ = false;

// DFS only has one active call at each depth. Reusing one candidate buffer per
// depth avoids an allocation at every search node without changing traversal
// order or candidate contents.
std::vector<std::vector<VertexID>> candidate_buffers_;

// During DSQL-P2 this is the number of vertices covered by the P2 snapshot
// that are not covered by the current solution. It replaces a full data-graph
// scan in every early-termination check.
size_t p2_uncovered_snapshot_vertices_ = 0;
bool p2_coverage_tracking_active_ = false;

public:
    /*   ����DSQL�㷨  */
    /// �±������ݶ���ID��ֵΪ�����ݶ�����ֵĴ�����m_VT[i]>=1����ʾ���ݶ���i�Ѿ������ǣ��Ҹ��ǵĴ���Ϊm_VT[i]
    std::vector<uint> m_VT;
    // �洢������ͼ,std::vector<uint>����match
    std::vector<std::vector<uint>> m_T;
    // m_VM[i]��ʾ��ѯ����i�Ѿ�ƥ������ݶ��㡣
    std::vector<std::unordered_set<uint>> m_VM;


    std::vector <std::vector<bool>> qOverlapList;
    // �ɵ�match_order,����ReSort()
    std::vector<uint> m_oldSort;

    // �ɵ�VM��VT������DSQL-P2
    std::vector<uint> m_p2VT;
    std::vector<std::unordered_set<uint>> m_p2VM;
    double m_gama = 0.0; // ���Ʊ�
    double m_alpha = 0.0; // ÿ�ε�����alpha


    // ������Ϣ
    size_t m_stopLevel = 0;
    double m_stopRatio = 0.0;
    bool m_isEnterP2 = false;
    bool m_isEarlyTer = false;
    size_t m_swapCount = 0;
    size_t m_QSDCallCount = 0;

public:
    DsqlMeta(Graph& query_graph, Graph& data_graph, size_t max_num_results,
    bool print_prep, bool print_enum, bool homo);

void SetICTCache(ICTCache* cache)
{
    ict_cache_ = cache;
    use_ict_cache_ = (cache != nullptr);
}

~DsqlMeta() override {};

    void Preprocessing() override;
    void InitialMatching() override;
    void PrintKeyVertexSet();
    bool VerifyCorrectness(const std::string& kv_path);



private:
    void GenerateMatchingOrder();
    void BuildCP2LEOrder();
    void BuildAdjMatrix();

    void BuildCover();

    bool ComputeCand(uint depth, std::vector<uint> m);
    bool FullCoveragePrune(const uint depth, const std::vector<uint>& m);

    void MutiExpansion(std::vector<uint> m);
    bool MutiExpTest(int depth, const std::vector<uint>& label_same_index, std::vector<uint>& m);

    void CountRes(std::vector<uint> m);
    void FlushFlag(uint flush_depth);


    void PrintMatch(const std::vector<uint>& m);



    void FindMatCo(uint depth, std::vector<uint> m);



public:
    void StartDsql();
    size_t DSQLP1();
    size_t DSQLP2(uint leveli);

    void ReSort(const std::vector<bool>& qovp);
    void Q1iSearch(uint depth, std::vector<uint>& m, const std::vector<bool>& qovp);
    bool Q2Search(uint depth, std::vector<uint>& m, const std::vector<bool>& qOvp, uint levelj);
    bool QSearchD(uint depth, std::vector<uint>& inM, const std::vector<bool>& qOvp);
    bool QSearchDP2(uint depth, std::vector<uint>& inM, const std::vector<bool>& qOvp);
    void ComputeCand(uint depth, const std::vector<uint>& m, const std::vector<bool>& qOvp,
        bool phase2, std::vector<uint>& cand);
    void SetCandidates(uint depth, const std::vector<bool>& qOvp, std::vector<uint>& cand);
    void SetCandidatesP2(uint depth, const std::vector<bool>& qOvp,std::vector<uint>& cand);
    void AddMatchedSubgraph(const std::vector<uint>& m);
    void GetOverlapList(uint leveli);
    void SetGamaAndAlpha();
    void SwapSubgraph(const std::vector<uint>& h);
    bool CheckEarlyTerminated(uint level);
    void RemoveSubgraph(const std::vector<uint>& f);
    uint ComputeBenefit(const std::vector<uint>& h);
    uint ComputeLoss(const std::vector<uint>& f);





    size_t GetNumKeyVertices() const { return this->num_keyvertex_; }
};
