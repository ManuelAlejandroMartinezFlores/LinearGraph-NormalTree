#include <cwchar>
#include "build_tests.cpp"
#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <omp.h>
#include <eigen3/Eigen/Sparse>

using namespace std;
using json = nlohmann::json;


int main(int argc, char** argv) {
    // vector<int> sources = {0, 1, 1};
    // vector<int> targets = {1, 0, 0};
    // save_graph(sources, targets, "test10");
    test(argv[1]); 
    return 0;
}

// static void BM_test17(benchmark::State& state) {
//     omp_set_num_threads(8);
//     LinearGraph g = graph_from_json("../assets/test17.json");
//     g.build_normal();
//     for (auto _ : state) {
//         // Code to benchmark
//         Eigen::SparseMatrix<double> result = g.generate_state_eq("");
//         benchmark::DoNotOptimize(result);
//         benchmark::DoNotOptimize(g);
//         benchmark::ClobberMemory();
//     }
// }


// BENCHMARK(BM_test17)->MinTime(2.0)->Repetitions(20)->UseRealTime();
// BENCHMARK_MAIN();

