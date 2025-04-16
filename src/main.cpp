#include <cwchar>
#include "build_tests.cpp"
#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <omp.h>

using namespace std;
using json = nlohmann::json;


int main(int argc, char** argv) {

    // auto a = make_shared<Symbol>("A");
    // auto b = make_shared<Symbol>("B");
    // auto c = make_shared<Symbol>("C");

    // cout << to_string((-a) / (-b)) << endl;


    // // auto sum = make_shared<Sum>();
    // // auto binv = make_shared<Reciprocal>(b);
    // // sum->add_sumand(a);
    // // sum->add_sumand(b);
    // cout << to_string(2 * (a + b) / (a + b)) << endl;
    // cout << to_string(sum) << " * " << to_string(sum2) << endl;
    // cout << "=" << endl;
    // auto m = sum * sum2;
    // cout << to_string(m) << endl;
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

