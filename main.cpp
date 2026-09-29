#include <iostream>
#include <filesystem>
#include "utils/timer.hpp"
#include "matching/MatchCover.h"
#include "matching/DsqlMeta.h"
#include "utils/globals.h"
#include "matching/ICTCache.h"
#include "utils/CLI11.hpp"
#include "utils/types.h"
#include <time.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <climits>
#include <memory>
namespace fs = std::filesystem;

static std::string s_thisProgramName = "";

std::vector<fs::path> getAllFilesInDirectory(const fs::path& directory) {
	std::vector<fs::path> filePaths;

	// �����ļ��в����ļ�·���洢�� vector ��
	for (const auto& entry : fs::directory_iterator(directory)) {
		if (fs::is_regular_file(entry)) {
			filePaths.push_back(entry.path());
		}
	}

	return filePaths;
}

void SaveRunningInfo(const std::string& dataGraphPath, const fs::path& queryGraphPath, double totalTime, const DsqlMeta* dsqlMeta, int numOfk,size_t memUsage)
{
	std::ofstream ofs;

	fs::path dataPath(dataGraphPath);
	std::string dataGraphName = dataPath.stem().string();
	std::string base = queryGraphPath.parent_path().string();
	std::string size = queryGraphPath.parent_path().filename().string();
	std::string path = base + "/" + s_thisProgramName + "_" + dataGraphName + "_" + size + "_" + std::to_string(numOfk);
	std::string subfix = "_DSQL0.csv";
#ifdef DSQL1
	subfix = "_DSQL1.csv";
#endif // DSQL1
#ifdef DSQL2
	subfix = "_DSQL2.csv";
#endif // DSQL2
#ifdef DSQL3
	subfix = "_DSQL3.csv";
#endif // DSQL3
#ifdef ALL
	subfix = "_ALL.csv";
#endif // ALL
	path += subfix;


	ofs.open(path, std::ios::app);
	ofs << std::fixed << std::setprecision(6) << totalTime << "," << dsqlMeta->m_T.size() << "," << 0 << "," << dsqlMeta->GetNumKeyVertices() << "," << dsqlMeta->m_stopLevel << ","
		<< dsqlMeta->m_stopRatio << "," << dsqlMeta->m_isEnterP2 << "," << dsqlMeta->m_swapCount << "," << dsqlMeta->m_isEarlyTer << "," << memUsage << "," << queryGraphPath.filename() << std::endl;
	ofs.close();

}

void SaveMemUsage(const std::string& dataGraphPath, const fs::path& queryGraphPath, double totalTime, const DsqlMeta* dsqlMeta, int numOfk,size_t memUsage)
{
	std::ofstream ofs;

	fs::path dataPath(dataGraphPath);
	std::string dataGraphName = dataPath.stem().string();
	std::string base = queryGraphPath.parent_path().string();
	std::string size = queryGraphPath.parent_path().filename().string();
	std::string path = base + "/" + s_thisProgramName + "_" + dataGraphName + "_" + size + "_" + std::to_string(numOfk);
	std::string subfix = "_Mem.csv";
	path += subfix;


	ofs.open(path, std::ios::app);
	ofs << std::fixed << std::setprecision(6) << totalTime << "," << dsqlMeta->m_T.size() << "," << 0 << "," << dsqlMeta->GetNumKeyVertices() << "," << dsqlMeta->m_stopLevel << ","
		<< dsqlMeta->m_stopRatio << "," << dsqlMeta->m_isEnterP2 << "," << dsqlMeta->m_swapCount << "," << dsqlMeta->m_isEarlyTer << "," << memUsage << "," << queryGraphPath.filename() << std::endl;
	ofs.close();
}


int main(int argc, char* argv[]) {


	std::string exec_name = argv[0];  // 获取 argv[0]
    
    // 找到最后一个 '/'，获取文件名部分
    size_t pos = exec_name.find_last_of("/\\");
    if (pos != std::string::npos) {
        exec_name = exec_name.substr(pos + 1);
    }
	s_thisProgramName = exec_name;

	CLI::App app{ "App description" };
	bool print_enum = false, print_prep = false, print_key_vertice = false;
	std::string data_graph_path;
	std::string query_graph_path;
	bool verify_correctness;

	size_t num_of_k = ULONG_MAX;
	bool save_run_info = false;

	app.add_option("-d,--data_graph_path", data_graph_path, "data_graph_path")->required();
	app.add_option("-q,--query_graph_path", query_graph_path, "query_graph_path")->required();
	app.add_option("-n,--num_of_k", num_of_k, "num_of_k")->default_val(ULONG_MAX);
	app.add_option("-e,--print_enum", print_enum, "print_enum")->default_val(false);
	app.add_option("-s,--save_run_info", save_run_info, "save_run_info")->default_val(false);
	app.add_option("-p,--print_prep", print_prep, "print_prep")->default_val(false);
	app.add_option("-k,--print_key_vertice", print_key_vertice, "print_key_vertice")->default_val(false);
	app.add_option("-v,--verify_correctness", verify_correctness, "verify_correctness")->default_val(false);

	CLI11_PARSE(app, argc, argv);

	std::cout << " n ===== " << num_of_k << std::endl;

	Graph data_graph;
	std::cout << "data graph loading..." << std::endl;
	data_graph.LoadFromFile(data_graph_path);
	data_graph.PrintMetaData();
	ICTCache sharedICT;


	fs::path directoryPath = query_graph_path;  // �뽫�˴��滻Ϊ����ļ���·��
	auto files = getAllFilesInDirectory(directoryPath);


	size_t totalKeyVertices = 0;
	size_t totalQueryGraph = 0;
	size_t totalSkipCount = 0;
	size_t totalSkipDataCount = 0;
	double totalTimeCost = 0.0;
	size_t totalMatching = 0;
	size_t totalMemUsage = 0;

	for (const auto& file : files)
	{
		if (file.extension() != ".graph")
			continue;
		++totalQueryGraph;

		Graph query_graph;
		std::cout << "the" << totalQueryGraph << " query graph" << std::endl;
		std::cout << "query graph " << file.string() << " loading..." << std::endl;

		query_graph.LoadFromFile(file.string());
		query_graph.PrintMetaData();

		int mem_start = mem::getMemUsage();

size_t ict_get_before = sharedICT.m_getCount;
size_t ict_hit_before = sharedICT.m_hitCount;
size_t ict_reuse_before = sharedICT.m_reuseCount;
size_t ict_reused_prefix_before = sharedICT.m_reusedPrefixLength;
size_t ict_new_node_before = sharedICT.m_newNodeCount;
size_t ict_build_before = sharedICT.m_buildCount;
size_t ict_intersection_elements_before = sharedICT.m_intersectionElementCount;
size_t ict_intersection_capacity_before = sharedICT.m_intersectionCapacityCount;

Timer t;
t.StartTimer();

auto dsqlMeta = std::make_unique<DsqlMeta>(
    query_graph, data_graph, num_of_k, print_prep, print_enum, false);
dsqlMeta->SetICTCache(&sharedICT);

std::cout << "DsqlMeta preprocessing..." << std::endl;

		dsqlMeta->Preprocessing();

		std::cout << "DsqlMeta initial matching..." << std::endl;

		auto InitialFun = [&dsqlMeta]()
			{
				dsqlMeta->InitialMatching();
			};
		execute_with_time_limit(InitialFun, 300, reach_time_limit);

		int mem_end = mem::getMemUsage();
		double time_cost = t.StopTimer_ms();

		std::cout << "DsqlMeta finished..." << std::endl;
		size_t DsqlMeata_num_results = 0ul;
		size_t DsqlMeata_num_kv = 0ul;
		dsqlMeta->GetNumKeyVertex(DsqlMeata_num_kv);
		dsqlMeta->GetNumInitialResults(DsqlMeata_num_results);

		const int mem_delta = (mem_start >= 0 && mem_end >= mem_start)
			? (mem_end - mem_start) : 0;
		size_t memUsage = static_cast<size_t>(mem_delta);
		std::cout << "Memory usage: " << mem_delta << " KB" << std::endl;
		std::cout << "DsqlMeta time: " << time_cost << " ms" << std::endl;
		std::cout << "DsqlMeta match size: " << dsqlMeta->m_T.size() << std::endl;
		std::cout << "DsqlMeta the number of key vertices: " << DsqlMeata_num_kv << std::endl;
		std::cout << "DsqlMeta stop level: " << dsqlMeta->m_stopLevel << std::endl;
		std::cout << "DsqlMeta the number of skipping query vertices: " << dsqlMeta->opt.skipQueryCount << std::endl;
		std::cout << "DsqlMeta the number of skipping data vertices: " << dsqlMeta->opt.skipDataVertexCount << std::endl;
		size_t ict_get_this_query = sharedICT.m_getCount - ict_get_before;
size_t ict_hit_this_query = sharedICT.m_hitCount - ict_hit_before;
size_t ict_reuse_this_query = sharedICT.m_reuseCount - ict_reuse_before;
size_t ict_reused_prefix_this_query = sharedICT.m_reusedPrefixLength - ict_reused_prefix_before;
size_t ict_new_node_this_query = sharedICT.m_newNodeCount - ict_new_node_before;
size_t ict_build_this_query = sharedICT.m_buildCount - ict_build_before;
size_t ict_intersection_elements_this_query =
    sharedICT.m_intersectionElementCount - ict_intersection_elements_before;
size_t ict_intersection_capacity_this_query =
    sharedICT.m_intersectionCapacityCount - ict_intersection_capacity_before;

std::cout << "ICT get count of this query: " << ict_get_this_query << std::endl;
std::cout << "ICT exact-hit count of this query: " << ict_hit_this_query << std::endl;
std::cout << "ICT reuse count of this query: " << ict_reuse_this_query << std::endl;
std::cout << "ICT new node count of this query: " << ict_new_node_this_query << std::endl;
std::cout << "ICT build count of this query: " << ict_build_this_query << std::endl;
std::cout << "ICT intersection elements added by this query: "
          << ict_intersection_elements_this_query << std::endl;
std::cout << "ICT intersection capacity bytes added by this query: "
          << ict_intersection_capacity_this_query * sizeof(VertexID) << std::endl;

if (ict_get_this_query > 0)
{
    std::cout << "ICT exact-hit ratio of this query: "
              << (ict_hit_this_query * 1.0 / ict_get_this_query)
              << std::endl;
    std::cout << "ICT reuse ratio of this query: "
              << (ict_reuse_this_query * 1.0 / ict_get_this_query)
              << std::endl;
    std::cout << "ICT average reused prefix length of this query: "
              << (ict_reuse_this_query == 0 ? 0.0 :
                  ict_reused_prefix_this_query * 1.0 / ict_reuse_this_query)
              << std::endl;
}
else
{
    std::cout << "ICT exact-hit ratio of this query: 0" << std::endl;
    std::cout << "ICT reuse ratio of this query: 0" << std::endl;
    std::cout << "ICT average reused prefix length of this query: 0" << std::endl;
}

		if (print_key_vertice) dsqlMeta->PrintKeyVertexSet();
		if (verify_correctness) {
			std::string base = file.parent_path().string();
			std::string stem = file.stem().string();
			std::string path = base + "/" + stem + "_kv.txt";
			if (dsqlMeta->VerifyCorrectness(path)) std::cout << "Correctness verified!" << std::endl;
			else
			{
				std::cout << "Correctness not verified!" << std::endl;
				exit(-1);
			}
		}
		totalKeyVertices += DsqlMeata_num_kv;
		totalSkipCount += dsqlMeta->opt.skipQueryCount;
		totalSkipDataCount += dsqlMeta->opt.skipDataVertexCount;
		totalTimeCost += time_cost;
		totalMatching += dsqlMeta->m_T.size();
		totalMemUsage += memUsage;

		// if (save_run_info) SaveRunningInfo(data_graph_path, file, time_cost, dsqlMeta, num_of_k,memUsage);
		if (save_run_info) SaveMemUsage(
			data_graph_path, file, time_cost, dsqlMeta.get(), num_of_k, memUsage);
	}
	std::cout << "##############################################" << std::endl;
	std::cout << "# total query graph : " << totalQueryGraph << std::endl;
	std::cout << "# total time cost : " << totalTimeCost << std::endl;
	std::cout << "# total matching : " << totalMatching << std::endl;
	std::cout << "# total key vertices : " << totalKeyVertices << std::endl;
	std::cout << "# total skip query node count : " << totalSkipCount << std::endl;
	std::cout << "# total skip data vertex count : " << totalSkipDataCount << std::endl;
	std::cout << "# total mem usage : " << totalMemUsage << std::endl;
	std::cout << "# ICT get count : " << sharedICT.m_getCount << std::endl;
	std::cout << "# ICT exact-hit count : " << sharedICT.m_hitCount << std::endl;
	std::cout << "# ICT reuse count : " << sharedICT.m_reuseCount << std::endl;
	std::cout << "# ICT new node count : " << sharedICT.m_newNodeCount << std::endl;
	std::cout << "# ICT build count : " << sharedICT.m_buildCount << std::endl;

	const ICTCache::MemoryStats ict_memory = sharedICT.GetMemoryStats();
	std::cout << "# ICT label root count : " << ict_memory.label_roots << std::endl;
	std::cout << "# ICT cache path node count : " << ict_memory.child_edges << std::endl;
	std::cout << "# ICT active node storage bytes : "
	          << ict_memory.ActiveNodeStorageBytes() << std::endl;
	std::cout << "# ICT node-pool reserved bytes : "
	          << ict_memory.PoolReservedBytes() << std::endl;
	std::cout << "# ICT intersection logical bytes : "
	          << ict_memory.IntersectionBytes() << std::endl;
	std::cout << "# ICT intersection capacity bytes : "
	          << ict_memory.IntersectionCapacityBytes() << std::endl;
	std::cout << "# ICT estimated hash bytes : "
	          << ict_memory.EstimatedHashBytes() << std::endl;
	std::cout << "# ICT estimated total reserved bytes : "
	          << ict_memory.EstimatedTotalReservedBytes() << std::endl;

if (sharedICT.m_getCount > 0)
{
    std::cout << "# ICT exact-hit ratio : "
              << (sharedICT.m_hitCount * 1.0 / sharedICT.m_getCount)
              << std::endl;
    std::cout << "# ICT reuse ratio : "
              << (sharedICT.m_reuseCount * 1.0 / sharedICT.m_getCount)
              << std::endl;
    std::cout << "# ICT average reused prefix length : "
              << (sharedICT.m_reuseCount == 0 ? 0.0 :
                  sharedICT.m_reusedPrefixLength * 1.0 / sharedICT.m_reuseCount)
              << std::endl;
}

	return 0;
}
