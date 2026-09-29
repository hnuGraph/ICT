#include "DsqlMeta.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph/graph.h"
#include "matching/MatchCover.h"
#include "matching/ICTCache.h"

#include <unordered_map>

void DsqlMeta::InitOptSESM()
{
	size_t verticesCount = query_.NumVertices();
	std::vector<bool> visited(verticesCount, false);
	std::unordered_map<LabelID, int> countMap;

	opt.labelRM.assign(verticesCount, 0);
	opt.neighborRM.assign(verticesCount, 0);

	// neighbor
	for (size_t k = 0; k < verticesCount; ++k)
	{
		VertexID curU = match_order[k];
		visited[curU] = true;

		LabelID label = query_.GetVertexLabel(curU);
		countMap[label]++;

		const auto& nbrs = query_.GetNeighbors(curU);
		const size_t nextNbrCount = nbrs.size();
		for (size_t i = 0; i < nextNbrCount; ++i)
		{
			VertexID nbrID = nbrs[i];
			if (!visited[nbrID])
			{
				++opt.neighborRM[k];
			}
		}
	}

	// label
	for (size_t k = 0; k < verticesCount; ++k)
	{
		VertexID curU = match_order[k];
		LabelID label = query_.GetVertexLabel(curU);
		opt.labelRM[k] = countMap[label];
	}
}

void DsqlMeta::SetOptSESM(uint depth, std::vector<uint>& cand)
{
	if (opt.neighborRM[depth] != 0) return;
	if (cand.size() > opt.labelRM[depth])
		cand.resize(opt.labelRM[depth]);

}

void DsqlMeta::InitOptSKNB(uint depth)
{
	opt.failDepth = NOFAILNODE;
	opt.matchSuc[depth] = false;
	opt.dynamicConflictTable[match_order[depth]].assign(query_.NumVertices(), false);
}

void DsqlMeta::SetOptDynamicCT(uint depth, const std::vector<uint>& inM,uint v)
{
	VertexID u = match_order[depth];
	for (uint p = 0; p < depth; ++p)
	{
		VertexID uid = match_order[p];
		if (inM[uid] == v) // 动态冲突
		{
			opt.dynamicConflictTable[u][uid] = true;
			break;
		}
	}
}

bool DsqlMeta::CheckOptSkipNode(uint curDepth)
{
	if (opt.failDepth == NOFAILNODE) return false;
	VertexID failVertice = match_order[opt.failDepth];
	VertexID curVertice = match_order[curDepth];
	if (adjacency_matrix[failVertice][curVertice]) return false;  // 静态冲突表
	if (opt.dynamicConflictTable[failVertice][curVertice]) return false; // 动态冲突表
	++opt.skipQueryCount;
	return true;
}

void DsqlMeta::ClearOptBadVertices(int depth)
{
	uint maxDepth = query_.NumVertices();
	if (depth == static_cast<int>(maxDepth) - 1) return;
	for(int i = depth + 1; i < static_cast<int>(maxDepth); ++i)
		opt.badVertices[i].clear();
}

void DsqlMeta::ClearOptAllBadVertices()
{
	for (uint i = 0; i < query_.NumVertices(); ++i)
		opt.badVertices[i].clear();
}

void DsqlMeta::SetOptBadVertices(uint depth,uint v)
{
	if (opt.failDepth == NOFAILNODE) return;
	if (depth == 0) return;
	VertexID failVertex = match_order[opt.failDepth];
	VertexID curVertex = match_order[depth];
	if (!adjacency_matrix[failVertex][curVertex]) return; // 冲突是因为动态冲突，而不是静态冲突。不可标记

	VertexID preVertex = match_order[depth-1];
	if (adjacency_matrix[failVertex][preVertex]) return;  // 静态冲突表
	if (opt.dynamicConflictTable[failVertex][preVertex]) return; // 动态冲突表
	opt.badVertices[depth - 1].insert(v);
}

bool DsqlMeta::CheckOptSkipDataVertice(uint depth,uint v)
{
	if (depth == 0) return false;
	if (opt.badVertices[depth - 1].find(v) == opt.badVertices[depth - 1].end()) return false;
	++opt.skipDataVertexCount;
	return true;
}

DsqlMeta::DsqlMeta(Graph& query_graph, Graph& data_graph,
	size_t max_num_results,
	bool print_prep,
	bool print_enum,
	bool homo)
	: matching(
		query_graph, data_graph, max_num_results,
		print_prep, print_enum, homo)
{
	opt.labelRM.resize(query_.NumVertices(), 0);
	opt.neighborRM.resize(query_.NumVertices(), 0);

	opt.matchSuc.resize(query_.NumVertices(), 0);
	opt.dynamicConflictTable.resize(query_.NumVertices(), std::vector<bool>(query_.NumVertices(), false));
	opt.badVertices.resize(query_.NumVertices(), std::unordered_set<uint>());
	candidate_buffers_.resize(query_.NumVertices());
}


void DsqlMeta::Preprocessing()
{

	GenerateMatchingOrder();
#ifdef CP2LE
	BuildCP2LEOrder();
#endif
	BuildAdjMatrix();
}


void DsqlMeta::BuildAdjMatrix() {

	//build adjacency matrix of query graph
	adjacency_matrix.reserve(query_.NumVertices());
	for (uint i = 0; i < query_.NumVertices(); i++) {
		adjacency_matrix.emplace_back();
		adjacency_matrix[i].reserve(query_.NumVertices());
		adjacency_matrix[i].resize(query_.NumVertices(), false);
	}

	for (uint i = 0; i < query_.NumVertices(); i++) {
		const auto& q_nbrs = query_.GetNeighbors(i);
		for (auto j : q_nbrs) {
			adjacency_matrix[i][j] = true;
			adjacency_matrix[j][i] = true;
		}
	}

	//compute prune_depth_ , 从prune_depth_开始，C(j)所有元素都填充了元素
	std::vector<bool> AllRN(query_.NumVertices(), false);
	for (uint i = 0; i < query_.NumVertices() - 1; i++)
	{
		std::fill(AllRN.begin(), AllRN.end(), false);
		for (uint s = 0; s <= i; s++)
		{
			uint uf = match_order[s];
			const auto& nbrs = query_.GetNeighbors(uf);
			for (auto j : nbrs)
			{
				AllRN[j] = true;
			}
		}
		bool flag = true;
		for (uint j = i + 1; j < query_.NumVertices(); j++)
		{
			uint ub = match_order[j];
			if (!AllRN[ub]) {
				flag = false;
				break;
			}
		}
		if (flag) {
			prune_depth_ = i;
			break;
		}
	}
	if (print_preprocessing_results_) std::cout << "prune_depth_: " << prune_depth_ << std::endl;
}
void DsqlMeta::BuildCP2LEOrder() {
	if (print_preprocessing_results_) std::cout << "BuildCP2LEOrder...." << std::endl;
	mutiexp_depth_ = UINT32_MAX;
	int sum = 0;
	const int query_size = static_cast<int>(query_.NumVertices());
	for (int i = query_size - 1; i >= 0; i--) {
		uint u = match_order[i];
		bool flag = true;
		for (int j = i + 1; j < query_size; j++) {
			uint v = match_order[j];
			if (query_.checkEdgeExistence(u, v)) {
				flag = false;
				break;
			}
		}
		if (flag) {
			sum++;
			for (int j = i; j < query_size - 1; j++) {
				match_order[j] = match_order[j + 1];
			}
			match_order[query_.NumVertices() - 1] = u;
		}
	}
	mutiexp_depth_ = query_.NumVertices() - sum;
	std::reverse(match_order.begin() + mutiexp_depth_, match_order.end());
	if (print_preprocessing_results_) {
		std::cout << "CP2LEOrder: " << std::endl;
		for (int i = 0; i < query_size; i++) {
			std::cout << match_order[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "CP2LE_depth_:" << mutiexp_depth_;
	}
	std::cout << std::endl;

}
void DsqlMeta::GenerateMatchingOrder()
{

	uint max_degree = 0u;
	uint max_degree_id = 0u;
	for (uint i = 0; i < query_.NumVertices(); i++)
	{
		if (query_.GetDegree(i) > max_degree)
		{
			max_degree = query_.GetDegree(i);
			max_degree_id = i;
		}
	}
	match_order.resize(query_.NumVertices());
	std::vector<bool> visited(query_.NumVertices(), false);
	match_order[0] = max_degree_id;
	visited[max_degree_id] = true;

	for (uint i = 1; i < query_.NumVertices(); ++i)
	{
		uint max_adjacent = 0;
		uint max_adjacent_u = NOT_EXIST;

		for (size_t j = 0; j < query_.NumVertices(); j++)
		{
			uint cur_adjacent = 0u;
			if (visited[j]) continue;

			auto& q_nbrs = query_.GetNeighbors(j);
			for (auto& other : q_nbrs)
				if (visited[other])
					cur_adjacent++;
			if (!cur_adjacent) continue;

			if (
				max_adjacent_u == NOT_EXIST ||
				(cur_adjacent == max_adjacent &&
					query_.GetDegree(j) > query_.GetDegree(max_adjacent_u)) ||
				cur_adjacent > max_adjacent
				) {
				max_adjacent = cur_adjacent;
				max_adjacent_u = j;
			}
		}

		match_order[i] = max_adjacent_u;
		visited[max_adjacent_u] = true;
	}

	//print order
	if (print_preprocessing_results_) {
		std::cout << "origin match order: " << std::endl;
		for (uint i = 0; i < query_.NumVertices(); i++)
		{
			std::cout << match_order[i] << " ";
		}
		std::cout << std::endl;
	}

}


void DsqlMeta::InitialMatching()
{
	StartDsql();
}

void DsqlMeta::BuildCover()
{
	uint max_degree_ = data_.GetMaxDegree();
	uint qsize_ = query_.NumVertices();
	CandidateSets.reserve(qsize_);
	for (uint i = 0; i < qsize_; ++i) {
		CandidateSets.emplace_back();
		CandidateSets[i].reserve(qsize_);
		for (uint j = 0; j < qsize_; ++j) {
			CandidateSets[i].emplace_back();
			CandidateSets[i][j].reserve(max_degree_);
		}
	}
	CandidateSetFlag.reserve(qsize_);
	for (uint i = 0; i < qsize_; i++) {
		CandidateSetFlag.emplace_back();
		CandidateSetFlag[i].resize(qsize_, false);
	}

	// 下标就是顶点ID
	std::vector<uint> m(query_.NumVertices(), UNMATCHED);

	KeyVertexSet = std::vector<bool>(data_.NumVertices(), false);
	LabelID firstLabel = query_.GetVertexLabel(match_order[0]);
const std::vector<uint>& firstCandidates = data_.GetVerticesByLabel(firstLabel);

for (uint i : firstCandidates)
{
    m[match_order[0]] = i;
    visited_[i] = true;

    FindMatCo(1, m);

    visited_[i] = false;
    m[match_order[0]] = UNMATCHED;
}
}

void DsqlMeta::FindMatCo(uint depth, std::vector<uint> m)
{
	if (reach_time_limit) return;
	if (num_initial_results_ >= max_num_results_) return;
	batches_num++;
	if (!ComputeCand(depth, m)) {
		empty_set_num++;
		return;
	}

#ifdef FullCoverage
	if (depth >= prune_depth_ && FullCoveragePrune(depth, m)) return;
#endif

#ifdef CP2LE
	if (depth == mutiexp_depth_) {
		muti_exp_num++;
		MutiExpansion(m);
		return;
	}
#endif
	uint u = match_order[depth];
	if (CandidateSets[depth][depth].size() == 0) {
    LabelID uLabel = query_.GetVertexLabel(u);
    const std::vector<uint>& vertices = data_.GetVerticesByLabel(uLabel);
    CandidateSets[depth][depth].assign(vertices.begin(), vertices.end());
}
	if (CandidateSets[depth][depth].size() == 0) return;

	for (uint i = 0; i < CandidateSets[depth][depth].size(); i++) {
		uint v = CandidateSets[depth][depth][i];
		if (!homomorphism_ && visited_[v]) continue;
		m[u] = v;
		visited_[v] = true;

		if (depth == query_.NumVertices() - 1)
		{
			num_initial_results_++;
			for (auto j : m) {
				if (KeyVertexSet[j]) continue;
				else {
					num_keyvertex_++;
					KeyVertexSet[j] = true;
				}
			}
			if (print_enumeration_results_) PrintMatch(m);
		}
		else
		{
			FindMatCo(depth + 1, m);
		}
		visited_[v] = false;
		m[u] = UNMATCHED;
		if (num_initial_results_ >= max_num_results_) return;
		if (reach_time_limit) return;
	}
}

void DsqlMeta::StartDsql()
{

	m_oldSort = match_order;
	m_VT = std::vector<uint>(data_.NumVertices(), 0);
	m_VM.resize(query_.NumVertices());
	p2_uncovered_snapshot_vertices_ = 0;
	p2_coverage_tracking_active_ = false;

	size_t level_i = DSQLP1();

	this->SetGamaAndAlpha();
	size_t matchCount = m_T.size();
	this->m_stopLevel = level_i;
	
	if (!reach_time_limit && level_i > 0 && matchCount == this->max_num_results_ && m_gama < 0.5)
	{
		// 第二阶段
		this->m_isEnterP2 = true;
		this->m_stopLevel = DSQLP2(level_i);
	}
	// 更新调试信息
	if ((level_i == 0) || (level_i == query_.NumVertices() - 1 && matchCount < max_num_results_))
	{
		this->m_stopRatio = 1.0;
	}
	else
	{
		this->SetGamaAndAlpha();
		this->m_stopRatio = this->m_gama;
	}

}

size_t DsqlMeta::DSQLP1()
{
	// 下标就是顶点ID
	uint verticesCount = query_.NumVertices();
	std::vector<uint> m(verticesCount, UNMATCHED);
	ReSort(std::vector<bool>());

#ifdef DSQL1
	InitOptSESM();
#endif // DSQL1
#ifdef DSQL3
	ClearOptAllBadVertices();
#endif // DSQL3
	Q1iSearch(0, m, std::vector<bool>());
	uint matchCount = m_T.size();
	if (matchCount >= max_num_results_)
	{
		return 0;
	}
	if (reach_time_limit) return 0;
	uint level_i = 1;
	for (; level_i <= verticesCount - 1; ++level_i)
	{
		GetOverlapList(level_i);
		for (uint over = 0; over < qOverlapList.size(); ++over)
		{
			const auto& Qovp = qOverlapList[over];
			ReSort(Qovp);
#ifdef DSQL1
			InitOptSESM();
#endif // DSQL1
#ifdef DSQL3
			ClearOptAllBadVertices();
#endif // DSQL3
			Q1iSearch(0, m, Qovp);
			matchCount = m_T.size();
			if (matchCount >= max_num_results_)
			{
				return level_i;
			}
			if (reach_time_limit) return level_i;
		}
	}
	return level_i - 1;
}

size_t DsqlMeta::DSQLP2(uint leveli)
{
	size_t verticesCount = query_.NumVertices();
	std::vector<uint> m(verticesCount, UNMATCHED);
	this->m_p2VT = this->m_VT;
	p2_uncovered_snapshot_vertices_ = 0;
	p2_coverage_tracking_active_ = true;

	// leveli > 0 才能进这个函数
	uint levelj = leveli;
	this->m_isEarlyTer = false;
	for (; levelj <= verticesCount - 1; ++levelj)
	{
		SetGamaAndAlpha();
		GetOverlapList(levelj);
		for (uint over = 0; over < qOverlapList.size(); ++over)
		{
			const auto& Qovp = qOverlapList[over];
			ReSort(Qovp);
#ifdef DSQL1
			InitOptSESM();
#endif // DSQL1
#ifdef DSQL3
			ClearOptBadVertices(-1);
#endif // DSQL3
			if (Q2Search(0, m, Qovp, levelj))
			{
				break; //提前终止
			}
		}
	}
	return levelj - 1;
}

void DsqlMeta::ReSort(const std::vector<bool>& qovp)
{

	uint verticesCount = query_.NumVertices();
	std::vector<bool> visited(verticesCount, false);
	std::vector<uint> newOrder;

	// 找到第一个在qovp的查询顶点；
	VertexID firstU = -1;
	if (qovp.empty())
	{
		return;
	}
	else
	{
		for (uint depth = 0; depth < verticesCount; ++depth)
		{
			if (qovp[this->m_oldSort[depth]])
			{
				firstU = this->m_oldSort[depth];
				break;
			}
		}
	}
	const auto& nbrs = query_.GetNeighbors(firstU);
	uint nbrsCount = nbrs.size();
	newOrder.push_back(firstU);
	visited[firstU] = true;
	for (uint i = 0; i < nbrsCount; ++i)
	{
		newOrder.push_back(nbrs[i]);
		visited[nbrs[i]] = true;
	}
	size_t k = 1;
	while (newOrder.size() < verticesCount)
	{
		const auto& nbrs = query_.GetNeighbors(newOrder[k]);
		nbrsCount = nbrs.size();
		for (uint i = 0; i < nbrsCount; ++i)
		{
			if (visited[nbrs[i]]) continue;
			newOrder.push_back(nbrs[i]);
			visited[nbrs[i]] = true;
		}
		++k;
	}
	match_order = std::move(newOrder);
}

void DsqlMeta::Q1iSearch(uint depth, std::vector<uint>& m, const std::vector<bool>& qOvp)
{
	if (reach_time_limit) return;
	std::vector<uint>& cand = candidate_buffers_[depth];
	ComputeCand(depth, m, qOvp, false, cand);

#ifdef DSQL2
	InitOptSKNB(depth);
#endif // DSQL2
#ifdef DSQL3
	InitOptSKNB(depth);
	ClearOptBadVertices(depth);
#endif // DSQL3


	VertexID u = match_order[depth];
	for (uint i = 0; i < cand.size(); i++) {
		uint v = cand[i];
		if (!homomorphism_ && visited_[v])
		{
#ifdef DSQL2
			SetOptDynamicCT(depth, m, v);
#endif // DSQL2
#ifdef DSQL3
			SetOptDynamicCT(depth, m, v);
#endif // DSQL3
			continue;
		}
		opt.matchSuc[depth] = true;

		// Recursive matches can change coverage after cand was constructed.
		if ((qOvp.empty() || !qOvp[u]) && (m_VT[v] != 0)) continue;

#ifdef DSQL3
		if (CheckOptSkipDataVertice(depth, v)) continue;
#endif // DSQL3

		m[u] = v;
		visited_[v] = true;
		if (depth == query_.NumVertices() - 1)
		{
			AddMatchedSubgraph(m);
			if (print_enumeration_results_) PrintMatch(m);
		}
		else
		{
			if (qOvp.empty() || !qOvp[u])
			{
				std::vector<uint> inM = m;
				if (QSearchD(depth + 1, inM, qOvp))
				{
					AddMatchedSubgraph(inM);
					if (print_enumeration_results_) PrintMatch(inM);
				}
			}
			else
			{
				Q1iSearch(depth + 1, m, qOvp);
			}
		}
		visited_[v] = false;
		m[u] = UNMATCHED;
		if (m_T.size() >= max_num_results_) return;
		if (reach_time_limit) return;
#ifdef DSQL2
		if (CheckOptSkipNode(depth)) break;
		opt.failDepth = NOFAILNODE;
#endif // DSQL2
#ifdef DSQL3
		ClearOptBadVertices(depth);
		if (CheckOptSkipNode(depth)) break;
		SetOptBadVertices(depth, v);
		opt.failDepth = NOFAILNODE;
#endif // DSQL3

	}
#ifdef DSQL2
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL2
#ifdef DSQL3
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL3
}

bool DsqlMeta::Q2Search(uint depth, std::vector<uint>& m, const std::vector<bool>& qOvp, uint levelj)
{
	if (reach_time_limit) return true;
	std::vector<uint>& cand = candidate_buffers_[depth];
	ComputeCand(depth, m, qOvp, true, cand);

#ifdef DSQL2
	InitOptSKNB(depth);
#endif // DSQL2
#ifdef DSQL3
	InitOptSKNB(depth);
	ClearOptBadVertices(depth);
#endif // DSQL3

	VertexID u = match_order[depth];
	for (uint i = 0; i < cand.size(); i++) {
		uint v = cand[i];
		if (!homomorphism_ && visited_[v])
		{
#ifdef DSQL2
			SetOptDynamicCT(depth, m, v);
#endif // DSQL2
#ifdef DSQL3
			SetOptDynamicCT(depth, m, v);
#endif // DSQL3
			continue;
		}
		opt.matchSuc[depth] = true;

		// Recursive swaps can change coverage after cand was constructed.
		if ((qOvp.empty() || !qOvp[u]) && (m_VT[v] != 0)) continue;

#ifdef DSQL3
		if (CheckOptSkipDataVertice(depth, v)) continue;
#endif // DSQL3


		m[u] = v;
		visited_[v] = true;
		if (depth == query_.NumVertices() - 1)
		{
			SwapSubgraph(m);
			if (print_enumeration_results_) PrintMatch(m);
			this->m_isEarlyTer = CheckEarlyTerminated(levelj);
		}
		else
		{
			if (qOvp.empty() || !qOvp[u])
			{
				std::vector<uint> inM = m;
				if (QSearchDP2(depth + 1, inM, qOvp))
				{
					SwapSubgraph(inM);
					if (print_enumeration_results_) PrintMatch(inM);
					this->m_isEarlyTer = CheckEarlyTerminated(levelj);
				}
			}
			else
			{
				Q2Search(depth + 1, m, qOvp, levelj);
			}
		}
		visited_[v] = false;
		m[u] = UNMATCHED;
		if (this->m_isEarlyTer) return true;
		if (reach_time_limit) return true;
#ifdef DSQL2
		if (CheckOptSkipNode(depth)) break;
		opt.failDepth = NOFAILNODE;
#endif // DSQL2
#ifdef DSQL3
		if (CheckOptSkipNode(depth)) break;
		SetOptBadVertices(depth, v);
		ClearOptBadVertices(depth);
		opt.failDepth = NOFAILNODE;
#endif // DSQL3
	}
#ifdef DSQL2
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL2
#ifdef DSQL3
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL2
	return false;
}

bool DsqlMeta::QSearchD(uint depth, std::vector<uint>& inM, const std::vector<bool>& qOvp)
{
	if (reach_time_limit) return false;
	std::vector<uint>& cand = candidate_buffers_[depth];
	ComputeCand(depth, inM, qOvp, false, cand);

#ifdef DSQL1
	SetOptSESM(depth, cand);
#endif // DSQL1
#ifdef DSQL2
	InitOptSKNB(depth);
#endif // DSQL2
#ifdef DSQL3
	InitOptSKNB(depth);
	ClearOptBadVertices(depth);
#endif // DSQL3


	VertexID u = match_order[depth];

	for (uint i = 0; i < cand.size(); i++) {
		uint v = cand[i];
		if (!homomorphism_ && visited_[v])
		{
#ifdef DSQL2
			SetOptDynamicCT(depth, inM, v);
#endif // DSQL2
#ifdef DSQL3
			SetOptDynamicCT(depth, inM, v);
#endif // DSQL3
			continue;
		}
		opt.matchSuc[depth] = true;

#ifdef DSQL3
		if (CheckOptSkipDataVertice(depth, v)) continue;
#endif // DSQL3
			
		inM[u] = v;
		visited_[v] = true;
		if (depth == query_.NumVertices() - 1)
		{
			visited_[v] = false;
			return true;
		}
		if (QSearchD(depth + 1, inM, qOvp)) // 找到了一个匹配
		{ 
			visited_[v] = false; // 当前depth的查询顶点匹配的数据顶点标记为未访问
			return true;
		}
		visited_[v] = false;
		inM[u] = UNMATCHED;
		if (reach_time_limit) return false;
#ifdef DSQL2
		if (CheckOptSkipNode(depth)) break;
		opt.failDepth = NOFAILNODE;
#endif // DSQL2
#ifdef DSQL3
		ClearOptBadVertices(depth);
		if (CheckOptSkipNode(depth)) break;
		SetOptBadVertices(depth,v);
		opt.failDepth = NOFAILNODE;
#endif // DSQL3
	}

#ifdef DSQL2
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL2
#ifdef DSQL3
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL3
	return false;
}

bool DsqlMeta::QSearchDP2(uint depth, std::vector<uint>& inM, const std::vector<bool>& qOvp)
{
	if (reach_time_limit) return false;
	std::vector<uint>& cand = candidate_buffers_[depth];
	ComputeCand(depth, inM, qOvp, true, cand);

#ifdef DSQL1
	SetOptSESM(depth, cand);
#endif // DSQL1
#ifdef DSQL2
	InitOptSKNB(depth);
#endif // DSQL2
#ifdef DSQL3
	InitOptSKNB(depth);
	ClearOptBadVertices(depth);
#endif // DSQL3

	VertexID u = match_order[depth];
	for (uint i = 0; i < cand.size(); i++) {
		uint v = cand[i];
		if (!homomorphism_ && visited_[v])
		{
#ifdef DSQL2
			SetOptDynamicCT(depth, inM, v);
#endif // DSQL2
#ifdef DSQL3
			SetOptDynamicCT(depth, inM, v);
#endif // DSQL3
			continue;
		}
		opt.matchSuc[depth] = true;


#ifdef DSQL3
		if (CheckOptSkipDataVertice(depth, v)) continue;
#endif // DSQL3
		
		inM[u] = v;
		visited_[v] = true;
		if (depth == query_.NumVertices() - 1)
		{
			visited_[v] = false;
			return true;
		}
		if (QSearchDP2(depth + 1, inM, qOvp)) // 找到了一个匹配
		{
			visited_[v] = false; // 当前depth的查询顶点匹配的数据顶点标记为未访问
			return true;
		}
		visited_[v] = false;
		inM[u] = UNMATCHED;
		if (reach_time_limit) return false;
#ifdef DSQL2
		if (CheckOptSkipNode(depth)) break;
		opt.failDepth = NOFAILNODE;
#endif // DSQL2
#ifdef DSQL3
		if (CheckOptSkipNode(depth)) break;
		SetOptBadVertices(depth, v);
		ClearOptBadVertices(depth);
		opt.failDepth = NOFAILNODE;
#endif // DSQL3
	}
#ifdef DSQL2
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL2
#ifdef DSQL3
	if (!opt.matchSuc[depth])  // 匹配失败。
	{
		opt.failDepth = depth;
	}
#endif // DSQL3
	return false;
}

void DsqlMeta::ComputeCand(
    uint depth,
    const std::vector<uint>& m,
    const std::vector<bool>& qOvp,
    bool phase2,
    std::vector<uint>& cand)
{
    cand.clear();

    VertexID curU = match_order[depth];
    LabelID curLabel = query_.GetVertexLabel(curU);
    const bool overlap_vertex = !qOvp.empty() && qOvp[curU];

    auto append_filtered = [&](const auto& source)
    {
        cand.reserve(source.size());
        for (VertexID dataId : source)
        {
            bool keep;
            if (overlap_vertex)
            {
                keep = phase2 ? (m_p2VT[dataId] != 0) : (m_VT[dataId] != 0);
            }
            else
            {
                keep = phase2
                    ? (m_VT[dataId] == 0 && m_p2VT[dataId] == 0)
                    : (m_VT[dataId] == 0);
            }
            if (keep) cand.push_back(dataId);
        }
    };

    if (depth == 0)
    {
        const std::vector<uint>& vertices = data_.GetVerticesByLabel(curLabel);
        append_filtered(vertices);
        return;
    }

    static constexpr size_t INLINE_DEPENDENCIES = 16;
    std::array<VertexID, INLINE_DEPENDENCIES> inline_dependencies;
    std::vector<VertexID> overflow_dependencies;
    const bool use_inline_dependencies = depth <= INLINE_DEPENDENCIES;
    if (!use_inline_dependencies) overflow_dependencies.reserve(depth);
    size_t dependency_count = 0;

    for (uint i = 0; i < depth; ++i)
    {
        VertexID preU = match_order[i];

        if (adjacency_matrix[curU][preU])
        {
            if (use_inline_dependencies)
                inline_dependencies[dependency_count] = m[preU];
            else
                overflow_dependencies.push_back(m[preU]);
            ++dependency_count;
        }
    }

    if (dependency_count == 0)
    {
        const std::vector<uint>& vertices = data_.GetVerticesByLabel(curLabel);
        append_filtered(vertices);
        return;
    }

    if (use_ict_cache_ && ict_cache_ != nullptr)
    {
        const VertexID* canonical_dependencies;
        if (use_inline_dependencies)
        {
            auto first = inline_dependencies.begin();
            auto last = first + dependency_count;
            std::sort(first, last);
            dependency_count = static_cast<size_t>(std::unique(first, last) - first);
            canonical_dependencies = inline_dependencies.data();
        }
        else
        {
            std::sort(overflow_dependencies.begin(), overflow_dependencies.end());
            overflow_dependencies.erase(
                std::unique(overflow_dependencies.begin(), overflow_dependencies.end()),
                overflow_dependencies.end()
            );
            dependency_count = overflow_dependencies.size();
            canonical_dependencies = overflow_dependencies.data();
        }

        const ICTCache::RawRange cached = ict_cache_->GetRawCanonical(
            curLabel, canonical_dependencies, dependency_count, data_);
        append_filtered(cached);
        return;
    }

    // Fallback: stage-2 label-aware intersection logic.
    std::vector<uint> temp;
    bool initialized = false;

    for (uint i = 0; i < depth; ++i)
    {
        VertexID preU = match_order[i];

        if (!adjacency_matrix[curU][preU])
        {
            continue;
        }

        VertexID preV = m[preU];
        const Graph::NeighborRange labelNbrs = data_.GetNeighborsByLabel(preV, curLabel);

        if (!initialized)
        {
            cand.assign(labelNbrs.begin(), labelNbrs.end());
            initialized = true;
        }
        else
        {
            temp.clear();

            std::set_intersection(
                cand.begin(), cand.end(),
                labelNbrs.begin(), labelNbrs.end(),
                std::back_inserter(temp)
            );

            cand.swap(temp);

            if (cand.empty())
            {
                break;
            }
        }
    }

    if (!initialized)
    {
        const std::vector<uint>& vertices = data_.GetVerticesByLabel(curLabel);
        cand.assign(vertices.begin(), vertices.end());
    }

    if (phase2) SetCandidatesP2(depth, qOvp, cand);
    else SetCandidates(depth, qOvp, cand);
}

void DsqlMeta::SetCandidates(uint depth, const std::vector<bool>& qOvp, std::vector<uint>& cand)
{
	VertexID queryId = match_order[depth];
	size_t write = 0;
	for (const uint dataId : cand)
	{
		bool keep = false;
		if (!qOvp.empty() && qOvp[queryId])
		{
			keep = (m_VT[dataId] != 0);
		}
		else
		{
			keep = (m_VT[dataId] == 0);
		}
		if (keep)
		{
			cand[write++] = dataId;
		}
	}
	cand.resize(write);
}

void DsqlMeta::SetCandidatesP2(uint depth, const std::vector<bool>& qOvp, std::vector<uint>& cand)
{
	VertexID queryId = match_order[depth];
	size_t write = 0;
	for (const uint dataId : cand)
	{
		bool keep = false;
		if (!qOvp.empty() && qOvp[queryId])
		{
			keep = (m_p2VT[dataId] != 0);
		}
		else
		{
			// 第二阶段的非重叠顶点也不能出现在m_p2VT中
			keep = (m_VT[dataId] == 0 && m_p2VT[dataId] == 0);
		}
		if (keep)
		{
			cand[write++] = dataId;
		}
	}
	cand.resize(write);
}

void DsqlMeta::AddMatchedSubgraph(const std::vector<uint>& m)
{
	for (uint i = 0; i < m.size(); ++i)
	{
		uint dataId = m[i];
		if (m_VT[dataId] == 0)
		{
			num_keyvertex_++;
			if (p2_coverage_tracking_active_ && m_p2VT[dataId] > 0)
			{
				assert(p2_uncovered_snapshot_vertices_ > 0);
				--p2_uncovered_snapshot_vertices_;
			}
		}
		++m_VT[dataId];
	}
	m_T.push_back(m);
}

void DsqlMeta::GetOverlapList(uint leveli)
{
	qOverlapList.clear();
	int n = query_.NumVertices();

	// 生成掩码：i 个 1 + (n - i) 个 0
	std::vector<int> bitmask(leveli, 1);
	bitmask.resize(n, 0); // 扩展成长度 n

	do {
		std::vector<bool> subset(n, false);
		for (int j = 0; j < n; ++j) {
			if (bitmask[j]) subset[this->m_oldSort[j]] = true;
		}
		qOverlapList.push_back(subset);
	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

void DsqlMeta::SetGamaAndAlpha()
{
	m_gama = num_keyvertex_ * 1.0 / (max_num_results_ * query_.NumVertices());
	m_alpha = 1.0 - 2.0 * m_gama;
}

void DsqlMeta::SwapSubgraph(const std::vector<uint>& h)
{
	uint B = ComputeBenefit(h);
	for (auto it = m_T.begin(); it != m_T.end(); ++it)
	{
		uint L = ComputeLoss(*it);
		// 满足交换条件
		if (B >= (1 + m_alpha) * L)
		{
			this->m_swapCount++;
			RemoveSubgraph(*it);
			m_T.erase(it);
			AddMatchedSubgraph(h);
			break;
		}
	}
}

bool DsqlMeta::CheckEarlyTerminated(uint level)
{
	// 条件1：C(T1) ∈ C(T);
	if (p2_uncovered_snapshot_vertices_ != 0) return false;

	// 条件2：L(f, T) ≥ (q − i)/(1 + α).
	double ra = (query_.NumVertices() - level) * 1.0 / (1 + m_alpha);
	for (const auto& f : m_T)
	{
		auto L = ComputeLoss(f);
		if (L < ra)
			return false;
	}

	return true;
}

void DsqlMeta::RemoveSubgraph(const std::vector<uint>& f)
{
	for (VertexID queryId = 0; queryId < f.size(); ++queryId)
	{
		VertexID dataId = f[queryId];
		if (--m_VT[dataId] == 0) {
			--num_keyvertex_;
			if (p2_coverage_tracking_active_ && m_p2VT[dataId] > 0)
			{
				++p2_uncovered_snapshot_vertices_;
			}
		}
	}
}

uint DsqlMeta::ComputeBenefit(const std::vector<uint>& h)
{
	uint benefit = 0;
	for (VertexID v : h) {
		if (m_VT[v] == 0) {  // 该顶点未被覆盖
			benefit++;
		}
	}
	return benefit;
}

uint DsqlMeta::ComputeLoss(const std::vector<uint>& f)
{
	uint loss = 0;
	for (VertexID v : f) {
		if (m_VT[v] == 1) {  // 仅属于 f 的顶点
			loss++;
		}
	}
	return loss;
}

void DsqlMeta::MutiExpansion(std::vector<uint> m) {

	std::vector<uint> label_same_index;
	std::vector<bool> label_check(query_.NumVertices(), false);
	for (uint i = mutiexp_depth_; i < query_.NumVertices(); i++) {
		if (label_check[i]) continue;
		uint uf = match_order[i];
		label_check[i] = true;
		label_same_index.clear();
		label_same_index.push_back(i);
		for (uint j = mutiexp_depth_; j < query_.NumVertices(); j++) {
			if (label_check[j]) continue;
			uint ub = match_order[j];
			if (query_.GetVertexLabel(uf) == query_.GetVertexLabel(ub)) {
				label_same_index.push_back(j);
				label_check[j] = true;
			}
		}
		bool flag = false;
		uint u_index = label_same_index[0];
		uint u = match_order[u_index];
		for (size_t s = 0; s < CandidateSets[mutiexp_depth_][u_index].size(); s++) {
			uint v = CandidateSets[mutiexp_depth_][u_index][s];
			if (visited_[v]) continue;
			visited_[v] = true;
			m[u] = v;
			if (MutiExpTest(1, label_same_index, m)) {
				visited_[v] = false;
				flag = true;
				break;
			}
			m[u] = UNMATCHED;
			visited_[v] = false;
		}
		if (!flag) return;
	}
	CountRes(m);
	FlushFlag(mutiexp_depth_ - 1);
}
void DsqlMeta::CountRes(std::vector<uint> m) {
	bool flag_all_cv = true;
	for (uint i = 0; i < query_.NumVertices(); i++) {
		uint u = match_order[i];
		if (KeyVertexSet[m[u]] == false) {
			KeyVertexSet[m[u]] = true;
			num_keyvertex_++;
			flag_all_cv = false;
		}
	}
	if (!flag_all_cv) {
		num_initial_results_++;
		if (print_enumeration_results_) PrintMatch(m);
	}
	for (uint i = mutiexp_depth_; i < query_.NumVertices(); i++) {
		for (auto cand : CandidateSets[mutiexp_depth_][i]) {
			if (KeyVertexSet[cand] == false) {
				KeyVertexSet[cand] = true;
				num_keyvertex_++;
				num_initial_results_++;
				if (print_enumeration_results_) {
					uint u = match_order[i];
					m[u] = cand;
					PrintMatch(m);
				}
			}
		}
	}
	for (uint i = mutiexp_depth_; i < query_.NumVertices(); i++) {
		uint u = match_order[i];
		m[u] = UNMATCHED;
	}
}

void DsqlMeta::FlushFlag(uint flush_depth) {
	uint u = match_order[flush_depth];
	for (uint i = flush_depth + 1; i < query_.NumVertices(); i++) {
		uint ub = match_order[i];
		if (adjacency_matrix[u][ub]) continue;
		else {
			CandidateSetFlag[flush_depth][i] = true;
		}
	}
}
void DsqlMeta::PrintMatch(const std::vector<uint>& m) {
	for (auto j : m)
	{
		std::cout << j << " ";
	}
	std::cout << std::endl;
}

void DsqlMeta::PrintKeyVertexSet() {
	std::cout << "Key vertices: ";
	for (uint i = 0; i < m_VT.size(); i++) {
		if (m_VT[i]) std::cout << i << " ";
	}
	std::cout << std::endl;
}

bool DsqlMeta::VerifyCorrectness(const std::string& kvPath) {
	std::ifstream inputFile(kvPath);
	if (!inputFile.is_open()) {
		std::cerr << "file path error: " << kvPath << std::endl;
		return false;
	}

	std::set<int> kvSet;
	for (size_t i = 0; i < m_VT.size(); i++) {
		if (m_VT[i]) kvSet.insert(i);
	}

	std::set<int> verifySet;
	std::string line;
	while (std::getline(inputFile, line)) {
		std::stringstream ss(line);
		std::string number;

		while (std::getline(ss, number, ',')) {
			try {
				verifySet.insert(std::stoi(number));
			}
			catch (const std::invalid_argument&) {
				std::cerr << "Invalid: " << number << std::endl;
				return false;
			}
		}
	}

	//VerifyCorrectness
	return kvSet == verifySet;
}


bool DsqlMeta::MutiExpTest(int depth, const std::vector<uint>& label_same_index, std::vector<uint>& m) {
	if (reach_time_limit) return false;
	if (depth == static_cast<int>(label_same_index.size())) return true;
	uint u_index = label_same_index[depth];
	uint u = match_order[u_index];
	for (size_t i = 0; i < CandidateSets[mutiexp_depth_][u_index].size(); i++) {
		uint v = CandidateSets[mutiexp_depth_][u_index][i];
		if (visited_[v]) continue;
		visited_[v] = true;
		m[u] = v;
		if (MutiExpTest(depth + 1, label_same_index, m)) {
			visited_[v] = false;
			return true;
		}
		m[u] = UNMATCHED;
		visited_[v] = false;
	}
	return false;

}

/// <summary>
/// 填充Cgi(j)数组，其中i=depth,j>=depth
/// </summary>
/// <param name="depth">当前查询深度</param>
/// <param name="m">当前匹配集合</param>
/// <returns>false表示Cgi(j)有空数组，否则没有</returns>
bool DsqlMeta::ComputeCand(uint depth, std::vector<uint> m) {
    uint uf = match_order[depth - 1];
    uint vf = m[uf];

    for (uint i = depth; i < query_.NumVertices(); i++) {
        uint ub = match_order[i];

        if (adjacency_matrix[uf][ub]) {
            CandidateSets[depth][i].clear();

            LabelID ubLabel = query_.GetVertexLabel(ub);
            const Graph::NeighborRange vf_label_nbrs = data_.GetNeighborsByLabel(vf, ubLabel);

            if (CandidateSets[depth - 1][i].size() == 0) {
                CandidateSets[depth][i].assign(vf_label_nbrs.begin(), vf_label_nbrs.end());
            }
            else {
                std::set_intersection(
                    CandidateSets[depth - 1][i].begin(),
                    CandidateSets[depth - 1][i].end(),
                    vf_label_nbrs.begin(),
                    vf_label_nbrs.end(),
                    std::back_inserter(CandidateSets[depth][i])
                );
            }

            if (CandidateSets[depth][i].size() == 0) {
                return false;
            }
        }
        else {
            // ub is not a right-neighbor of uf, so its candidate set is unchanged.
            CandidateSets[depth][i] = CandidateSets[depth - 1][i];
        }
    }

    return true;
}

bool DsqlMeta::FullCoveragePrune(uint depth, const std::vector<uint>& m) {

	full_coverage_test_num++;
	for (uint i = 0; i < depth; i++) {
		uint u = match_order[i];
		if (!KeyVertexSet[m[u]]) {
			for (uint i = depth; i < query_.NumVertices(); i++) {
				CandidateSetFlag[depth][i] = CandidateSetFlag[depth - 1][i];
			}
			return false;
		}
	}
	partial_match_is_all_kv_num++;
	for (uint i = depth; i < query_.NumVertices(); i++) {
		if (CandidateSets[depth][i].size() == 0) {
			for (uint j = i; j < query_.NumVertices(); j++) {
				CandidateSetFlag[depth][j] = CandidateSetFlag[depth - 1][j];
			}
			return false;
		}
		candidate_set_test_num++;
		if (CandidateSetFlag[depth - 1][i]) {
			candidate_set_flag_reuse_num++;
			CandidateSetFlag[depth][i] = true;
		}
		else {
			CandidateSetFlag[depth][i] = true;
			for (auto v : CandidateSets[depth][i]) {
				if (!KeyVertexSet[v]) {
					for (uint j = i; j < query_.NumVertices(); j++) {
						CandidateSetFlag[depth][j] = CandidateSetFlag[depth - 1][j];
					}
					return false;
				}
			}
		}
	}
	full_coverage_suc_num++;
	FlushFlag(depth - 1);
	return true;
}
