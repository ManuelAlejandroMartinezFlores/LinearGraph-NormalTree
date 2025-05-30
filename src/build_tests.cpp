#include <iostream>
#include "normal_tree.cpp"
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <omp.h>

using namespace std;
using json = nlohmann::json;


void save_graph(vector<int> sources, vector<int> targets, string name) {
    /*
    Saves a graph with given sources and targets
    */
    vector<Edge> edges;
    if (sources.size() != targets.size()) {
        throw runtime_error("Different sizes");
    }
    for (int i=0; i<sources.size(); ++i) {
        Edge e;
        e.source_node = sources[i];
        e.target_node = targets[i];
        edges.push_back(e);
    }
    json j = edges;
    ofstream file("../assets/" + name + ".json");
    file << j.dump(4);
    file.close();
}


void test(string name){
    /*
    Builds the normal tree and generates state equations for a graph
    */
    LinearGraph g = graph_from_json(join_paths({"..", "assets", "tests", name, "edges.json"}));
    g.build_normal();
    g.generate_state_eq(join_paths({"..", "assets", "tests", name}), true);
    g.generate_state_eq_symbolic(join_paths({"..", "assets", "tests", name, "symb"}), true);
    g.save_to_json(join_paths({"..", "assets", "tests", name, "normal_tree.json"}));

}

void benchmark(string test_name) {
    /*
    Benchmarks the time performance of operations in different number of threads
    */
    json j;
    ifstream file("../assets/" + test_name + ".json");
    file >> j;
    file.close();

    for (int threads = 1; threads<9; threads++){
        omp_set_num_threads(threads);
        cout << "Threads: " << threads << endl;
        std::vector<double> timings;
        for (int run = 0; run < 100; run++) {
            auto start = std::chrono::high_resolution_clock::now();
            LinearGraph g(j.get<vector<Edge>>());
            g.build_normal();
            g.generate_state_eq("");
            auto end = std::chrono::high_resolution_clock::now();
            
            double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
            timings.push_back(elapsed);
        }

        double sum = std::accumulate(timings.begin(), timings.end(), 0.0);
        double mean = sum / timings.size();
        
        double sq_sum = std::inner_product(timings.begin(), timings.end(), 
                                            timings.begin(), 0.0,
                                            std::plus<>(), 
                                            [mean](double x, double y) {
                                                return (x - mean) * (y - mean);
                                            });
        double stdev = std::sqrt(sq_sum / timings.size());
        
        std::cout << "Average: " << mean << " ms, Stdev: " << stdev << " ms" << std::endl;
    }
}