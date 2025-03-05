#include <iostream>
#include "normal_tree.cpp"
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>

using namespace std;
using json = nlohmann::json;


void save_graph(vector<int> sources, vector<int> targets, string name) {
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
    LinearGraph g = graph_from_json("../assets/" + name + ".json");
    g.build_normal();
    g.generate_state_eq(name);
    // g.get_path(0, 2);
    g.save_to_json("../assets/" + name + "_normal_tree.json");

}