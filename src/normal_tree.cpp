#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/SparseLU>
#include <unordered_set>
#include <stack>
#include "state_eq.cpp"

using namespace std;
using json = nlohmann::json;
using namespace Eigen;


class Edge {
    public:
    int id;
    int source_node;
    int target_node;
    string type;
    string constant_type;
    double value;
    bool across_source;
    bool through_source;
    string across;
    string through;
    bool in_normal;

    Edge() {
        id = 0;
        source_node = 0; 
        target_node = 0;
        type = "A";
        constant_type = "";
        value = 0;
        across_source = false;
        through_source = false;
        across = "";
        through = "";
        in_normal = false;
    }

    friend std::ostream& operator<<(std::ostream& os, const Edge& e) {
        return os << e.source_node << " - " << e.target_node << " - " << e.type << " - normal: " << e.in_normal;
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Edge, source_node, target_node, type, constant_type, value, across_source, through_source, across, through, in_normal);

};

class LinearGraph
{
public:
    int num_nodes;
    int num_edges;
    int num_sources_across;
    int num_sources_through;
    unordered_map<int, int> rank; 
    unordered_map<int, int> parent;
    unordered_set<int> independent;
    unordered_set<int> sources;
    bool normal_tree_built;
    unordered_map<int, unordered_map<int, vector<int>>> edgelist;
    vector<Edge> edges;
    unordered_map<int, unordered_map<int,int>> normal_neighbors;

    LinearGraph(vector<Edge> edges){
        num_nodes = 0;
        num_edges = 0;
        num_sources_across = 0;
        num_sources_through = 0;
        normal_tree_built = false;
        for (Edge& edge : edges) {
            add_edge(edge);
        }
    }
    LinearGraph(){
        num_nodes = 0;
        num_edges = 0;
        normal_tree_built = false;
    }

    void add_edge(Edge& edge) {
        edge.id = num_edges;
        rank[edge.source_node] = 1;
        parent[edge.source_node] = edge.source_node;
        rank[edge.target_node] = 1;
        parent[edge.target_node] = edge.target_node;
        if (edgelist.find(edge.source_node) != edgelist.end()) {
            if (edgelist[edge.source_node].find(edge.target_node) == edgelist[edge.source_node].end()) {
                edgelist[edge.source_node].emplace(edge.target_node, vector<int>());
            }
        } else {
            edgelist.emplace(edge.source_node, unordered_map<int, vector<int>>()); 
            edgelist[edge.source_node].emplace(edge.target_node, vector<int>());
        }
        edgelist[edge.source_node][edge.target_node].push_back(edge.id);
        if (edge.across_source) {
            ++num_sources_across;
        }
        if (edge.through_source) {
            ++num_sources_through;
        }
        edges.push_back(edge);
        ++num_edges;
        normal_tree_built = false;
        num_nodes = rank.size();
        if (edge.across_source or edge.through_source) {
            sources.insert(edge.id);
        }
    }


    void save_to_json(string path) {
        json j = edges;
        ofstream file(path);
        file << j.dump(4);
        file.close();
    }

    friend std::ostream& operator<<(std::ostream& os, const LinearGraph& g) {
        return os << "nodes: " << g.num_nodes << " - branches: " << g.num_edges;
    }
    
    void add_neighbors(int n, int m, int id) {
        if (normal_neighbors.end() == normal_neighbors.find(n)) {
            normal_neighbors[n] = unordered_map<int, int>();
        }
        normal_neighbors[n][m] = id;
    }

    void build_normal() {
        int tree_branches = 0;
        // All source across
        for (Edge& edge : edges){
            if (edge.across_source) {
                if (!add_branch_to_normal(edge.source_node, edge.target_node)) {
                    throw runtime_error("Sources form a loop");
                };
                edge.in_normal = true;
                ++tree_branches;
                add_neighbors(edge.source_node, edge.target_node, edge.id);
                add_neighbors(edge.target_node, edge.source_node, edge.id);
            } 
        }
        for (Edge& edge : edges) {
            if (edge.type == "A") {
                if (add_branch_to_normal(edge.source_node, edge.target_node)) {
                    edge.in_normal = true;
                    ++tree_branches;
                    add_neighbors(edge.source_node, edge.target_node, edge.id);
                    add_neighbors(edge.target_node, edge.source_node, edge.id);
                    independent.insert(edge.id);
                }
            }
        }
        for (Edge& edge : edges) {
            if (tree_branches == num_nodes - 1) {
                break;
            }
            if (edge.type == "D") {
                if (add_branch_to_normal(edge.source_node, edge.target_node)) {
                    edge.in_normal = true;
                    ++tree_branches;
                    add_neighbors(edge.source_node, edge.target_node, edge.id);
                    add_neighbors(edge.target_node, edge.source_node, edge.id);
                }
            }
        }
        for (Edge& edge : edges) {
            if (edge.type == "T") {
                // cout << edge.id << "\n";
                if (add_branch_to_normal(edge.source_node, edge.target_node)) {
                    edge.in_normal = true;
                    ++tree_branches;
                    add_neighbors(edge.source_node, edge.target_node, edge.id);
                    add_neighbors(edge.target_node, edge.source_node, edge.id);
                } else {
                    independent.insert(edge.id);
                }
            }
        }
        for (Edge& edge : edges) {
            if (edge.through_source) {
                if (add_branch_to_normal(edge.source_node, edge.target_node)) {
                    throw runtime_error("Continuity violated: through sources can be added");
                }
            }
        }

        if (tree_branches != num_nodes - 1) {
            throw runtime_error("Could not build normal tree");
        }
        normal_tree_built = true;
    }

    void generate_state_eq(string path) {
        if (!normal_tree_built) {
            build_normal();
        }
        Eigen::SparseMatrix<double> matrix(3*num_edges-3*(num_sources_across+num_sources_through), 4*num_edges);
        vector<string> primary;
        vector<string> secondary;
        int row = 0;

        // Elemental equations
        ElementalEq eqs; 
        for (Edge& edge : edges) {
            if (!edge.across_source and !edge.through_source) {
                array<double,4> vals = eqs.equations[edge.constant_type] (edge.value);
                for (int i=0; i<4; ++i) {
                    matrix.insert(row, edge.id + i*num_edges) = vals[i];
                }
                ++row;
            }
        }
        cout << "Elemental: " << row << endl;
        // Continuity equations
        for (auto& [k, map] : normal_neighbors) {
            for (auto& [m, edgeid] : map) {
                if (m <= k) continue;
                // if (edges[edgeid].across_source) continue;
                int n = edges[edgeid].source_node;
                vector<int> vedge = edges_from_node(n);
                if (vedge.size() == 0) continue;   
                for (int edgeid : vedge) {
                    Edge edge = edges[edgeid];
                    if (edge.source_node == n){
                        matrix.insert(row, edge.id + 3*num_edges) = 1;
                        matrix.insert(row+1, edge.id + 2*num_edges) = 1;
                    } else {
                        matrix.insert(row, edge.id + 3*num_edges) = -1;
                        matrix.insert(row+1, edge.id + 2*num_edges) = -1;
                    }
                }
                ++row;
                ++row;
            }
        }
        cout << "Continuity: " << row << endl;
        // for (int n=0; n<num_nodes; ++n) {
        //     vector<Edge> edges = edges_from_node(n);
        //     if (edges.size() == 0) {
        //         continue;
        //     }
        //     for (Edge& edge: edges_from_node(n)) {
        //         if (edge.source_node == n){
        //             matrix.insert(row, edge.id + 3*num_edges) = 1;
        //         } else {
        //             matrix.insert(row, edge.id + 3*num_edges) = -1;
        //         }
        //     }
        //     ++row;
        // }

        // Compatibility
        for (Edge& edge : edges) {
            if (edge.in_normal or edge.through_source) continue;
            matrix.insert(row, edge.id + num_edges) = -1;
            matrix.insert(row+1, edge.id) = -1;
            for (auto& [edgeid, value] : get_path(edge.source_node, edge.target_node)) {
                matrix.insert(row, edgeid + num_edges) = value;
                matrix.insert(row+1, edgeid) = value;
            }
            ++row;
            ++row;
        }
        cout << "Compatibility: " << row << endl;


        saveSparseMatrixToCSV(matrix, "../assets/" + path + "_eqs.csv");
        matrix.makeCompressed();
        vector<int> priorityVars;
        for (Edge& edge : edges) {
            if ((independent.find(edge.id) == independent.end()) and  (sources.find(edge.id) == sources.end())) {
                int var = edge.id + 2 * edge.in_normal * num_edges;
                priorityVars.push_back(var);
                priorityVars.push_back(var + num_edges);
            }  
        }
        
        for (Edge& edge : edges) {
            if ((independent.find(edge.id) == independent.end()) and (sources.find(edge.id) == sources.end())) {
                int var = edge.id + 2 * !edge.in_normal * num_edges;
                priorityVars.push_back(var);
                priorityVars.push_back(var + num_edges);
            }  
        }
        
        for (int edgeid : independent) {
            int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            priorityVars.push_back(var + num_edges);
        }
        matrix = gaussianElimination(matrix, priorityVars);
        priorityVars = {};
        // for (int edgeid : independent) {
        //     int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
        //     priorityVars.push_back(var);
        //     priorityVars.push_back(var + num_edges);
        // }
        // for (int edgeid : sources) {
        //     int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
        //     priorityVars.push_back(var);
        //     priorityVars.push_back(var + num_edges);
        // }
        // for (int id : independent) {
        //     cout << id << "\n";
        // }
        // for (int id : priorityVars) {
        //     cout << id << "\n";
        // }
        // for (int edgeid : independent) {
        //     int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
        //     priorityVars.push_back(var);
        //     priorityVars.push_back(var + num_edges);
        // }
        // for (int edgeid : sources) {
        //     int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
        //     priorityVars.push_back(var);
        //     priorityVars.push_back(var + num_edges);
        // }
        // matrix = gaussianElimination(matrix, priorityVars);
        
        saveSparseMatrixToCSV(matrix, "../assets/" + path + "_reduced.csv");
    }
    
    vector<int> edges_from_node(int n) {
        vector<int> redges;
        bool flag = false;
        for (auto& [source, map] : edgelist) {
            for (auto& [target, list] : map) {
                if (target == n or source == n) {
                    for (int edgeid : list) {
                        flag = edges[edgeid].across_source or flag;
                        redges.push_back(edgeid);
                        if (flag) {break;}
                    }
                }
                if (flag) {break;}
            }
            if (flag) {break;}
        }
        if (!flag) {
            return redges;
        }
        vector<int> edges;
        return edges;
    }

    unordered_map<int, int> get_path(int from, int to) {
        unordered_set<int> visited; 
        unordered_map<int, int> parent;
        unordered_map<int, int> corresponding_edge;
        unordered_map<int, int> edge_values;
        stack<int> queue; 

        queue.push(from);
        visited.insert(from);
        parent[from] = -1;

        while (!queue.empty()) {
            int node = queue.top();
            queue.pop();
            if (node == to) break;
            if (normal_neighbors.find(node) == normal_neighbors.end()) continue;
            for (auto& [n, edgeid] : normal_neighbors[node]) {  
                if (visited.find(n) == visited.end()) {
                    visited.insert(n);
                    queue.push(n);
                    parent[n] = node;
                    corresponding_edge[n] = edgeid;
                }
            }
        }
        if (parent.find(to) == parent.end()) {
            return {};
        }
        int at = to;
        while (parent[at] != -1) {
            Edge edge = edges[corresponding_edge[at]];
            if (edge.target_node == at) {
                edge_values[edge.id] = 1;
            } else {
                edge_values[edge.id] = -1;
            }
            at = parent[at];
            
        }
        return edge_values;
    }

    private:
    int find_parent(int u) {
        int parent_u = parent[u];
        if (parent_u == u) {
            return u;
        }
        return find_parent(parent_u);
    }

    bool add_branch_to_normal(int u, int v) {
        int parent_u = find_parent(u);
        int parent_v = find_parent(v);
        if (parent_u == parent_v) { 
            return false; 
        }
        int rank_u = rank[parent_u];
        int rank_v = rank[parent_v];
        if (rank_u < rank_v) {
            rank[parent_v] = rank_v + 1;
            parent[parent_u] = parent_v; 
        } else {
            rank[parent_u] = rank_u + 1;
            parent[parent_v] = parent_u; 
        }
        return true;
    }

    

};

LinearGraph graph_from_json(string path) {
    json j;
    ifstream file(path);
    file >> j;
    file.close();
    vector<Edge> edges = j.get<vector<Edge>>();
    LinearGraph g(edges);
    return g;
}
