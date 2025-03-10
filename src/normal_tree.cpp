#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/SparseLU>
#include <unordered_set>
#include <stack>
#include <omp.h>
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
    vector<int> in_normal;
    bool normal_tree_built;
    unordered_map<int, unordered_map<int, vector<int>>> edgelist;
    vector<Edge> edges;
    unordered_map<int, unordered_map<int,int>> normal_neighbors;\
    unordered_map<int, int> elemental_rows;
    unordered_map<int, int> compatibility_rows; 
    unordered_map<int, int> continuity_rows;

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
        } else {
            elemental_rows[edge.id] = elemental_rows.size();
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
                    throw runtime_error("Across sources form a loop");
                };
                edge.in_normal = true;
                ++tree_branches;
                add_neighbors(edge.source_node, edge.target_node, edge.id);
                add_neighbors(edge.target_node, edge.source_node, edge.id);
                in_normal.push_back(edge.id);
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
                    in_normal.push_back(edge.id);
                    continuity_rows[edge.id] = continuity_rows.size()*2 + num_edges - num_sources_across - num_sources_through;
                } else {
                    compatibility_rows[edge.id] = compatibility_rows.size()*2 + num_edges + 2*num_nodes - 3*num_sources_across - num_sources_through - 2;
                }
            }
        }
        for (Edge& edge : edges) {
            // if (tree_branches == num_nodes - 1) {
            //     break;
            // }
            if (edge.type == "D") {
                if (add_branch_to_normal(edge.source_node, edge.target_node)) {
                    edge.in_normal = true;
                    ++tree_branches;
                    add_neighbors(edge.source_node, edge.target_node, edge.id);
                    add_neighbors(edge.target_node, edge.source_node, edge.id);
                    in_normal.push_back(edge.id);
                    continuity_rows[edge.id] = continuity_rows.size()*2 + num_edges - num_sources_across - num_sources_through;
                } else {
                    compatibility_rows[edge.id] = compatibility_rows.size()*2 + num_edges + 2*num_nodes - 3*num_sources_across - num_sources_through - 2;
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
                    in_normal.push_back(edge.id);
                    continuity_rows[edge.id] = continuity_rows.size()*2 + num_edges - num_sources_across - num_sources_through;
                } else {
                    independent.insert(edge.id);
                    compatibility_rows[edge.id] = compatibility_rows.size()*2 + num_edges + 2*num_nodes - 3*num_sources_across - num_sources_through - 2;
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

    void generate_state_eq(string path, bool verbose=false) {
        if (!normal_tree_built) {
            build_normal();
        }
        if (verbose) cout << "Normal tree built" << endl;
        
        // vector<string> primary;
        // vector<string> secondary;
        vector<Triplet<double>> triplets;  // Temporary storage
        vector<vector<Triplet<double>>> thread_triplets(omp_get_max_threads());
        if (verbose) cout << "Threads: " << omp_get_max_threads() << endl;
        

        // Elemental equations
        ElementalEq eqs; 
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < edges.size(); ++i) {
            Edge edge = edges[i];
            if (!edge.across_source and !edge.through_source) {
                int thread_id = omp_get_thread_num();
                int row = elemental_rows[edge.id];
                vector<Triplet<double>> &local_triplets = thread_triplets[thread_id];
                array<double,4> vals = eqs.equations[edge.constant_type] (edge.value);
                for (int i=0; i<4; ++i) {
                    local_triplets.emplace_back(row, edge.id + i*num_edges, vals[i]);
                }
            }
        }
        if (verbose) cout << "Elemental equations built" << endl;

        // Continuity equations
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < in_normal.size(); ++i) {
            int edgeid = in_normal[i];
            if (edges[edgeid].across_source) continue;
            int thread_id = omp_get_thread_num();
            int row = continuity_rows[edgeid];
            vector<Triplet<double>> &local_triplets = thread_triplets[thread_id];
            for (int n : children_nodes(edges[edgeid].source_node, edges[edgeid].target_node)) {
                for (int m : children_nodes(edges[edgeid].target_node, edges[edgeid].source_node)) {
                    if (edgelist.find(n) != edgelist.end() && edgelist[n].find(m) != edgelist[n].end()){
                        for (int edgeidx : edgelist[n][m]) {
                            local_triplets.emplace_back(row, edgeidx + 3*num_edges, 1);
                            local_triplets.emplace_back(row+1, edgeidx + 2*num_edges, 1);
                        }
                    }
                    if (edgelist.find(m) != edgelist.end() && edgelist[m].find(n) != edgelist[m].end()){
                        for (int edgeidx : edgelist[m][n]) {
                            local_triplets.emplace_back(row, edgeidx + 3*num_edges, -1);
                            local_triplets.emplace_back(row+1, edgeidx + 2*num_edges, -1);
                        } 
                    }
                    
                }
            }
        }
        if (verbose) cout << "Continuity equations built" << endl;


        // Compatibility
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < edges.size(); ++i) {
            const Edge edge = edges[i];
            if (edge.in_normal or edge.through_source) continue;
            int thread_id = omp_get_thread_num();
            int row = compatibility_rows[edge.id];
            vector<Triplet<double>> &local_triplets = thread_triplets[thread_id];
            local_triplets.emplace_back(row, edge.id + num_edges, -1);
            local_triplets.emplace_back(row+1, edge.id, -1);
            for (auto& [edgeid, value] : get_path(edge.source_node, edge.target_node)) {
                local_triplets.emplace_back(row, edgeid + num_edges, value);
                local_triplets.emplace_back(row+1, edgeid, value);
            }
        }
        for (const auto &local : thread_triplets) {
            triplets.insert(triplets.end(), local.begin(), local.end());
        }

        if (verbose) cout << "Compatibility equations built" << endl;

        Eigen::SparseMatrix<double, RowMajor> matrix(3*num_edges-3*(num_sources_across+num_sources_through), 4*num_edges);
        matrix.setFromTriplets(triplets.begin(), triplets.end());

        if (path != "") saveSparseMatrixToCSV(matrix, "../assets/" + path + "/eqs.csv");
        matrix.makeCompressed();
        if (verbose) cout << "Equations matrix built" << endl;
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
        if (verbose) cout << "Equations reduced" << endl;
        if (path != "") saveSparseMatrixToCSV(matrix, "../assets/" + path + "/reduced.csv");

        SparseMatrix<double> var_matrix(independent.size(), independent.size());
        SparseMatrix<double> source_matrix(independent.size(), sources.size()*2);
        unordered_map<int, int> var_cols;
        vector<string> x;
        for (int edgeid : independent) {
            int var = edgeid + 2 * num_edges * !edges[edgeid].in_normal;
            // var_cols[var] = var_cols.size();
            var_cols[var + num_edges] = var_cols.size();
            if (edges[edgeid].in_normal) {
                x.push_back(edges[edgeid].across); 
            } else {
                x.push_back(edges[edgeid].through); 
            }
        }
        unordered_map<int, int> source_cols;
        vector<string> s;
        for (int edgeid : sources) {
            int var = edgeid + 2 * num_edges * !edges[edgeid].in_normal;
            source_cols[var] = source_cols.size();
            source_cols[var + num_edges] = source_cols.size();
            if (edges[edgeid].in_normal) {
                s.push_back("d" + edges[edgeid].across); 
                s.push_back(edges[edgeid].across); 
                
            } else {
                s.push_back("d" + edges[edgeid].through);
                s.push_back(edges[edgeid].through);   
            }
        }
        saveVarIndex(x, "../assets/" + path + "/var_index.csv");
        saveVarIndex(s, "../assets/" + path + "/sources_index.csv");
        int new_row = 0;
        for (int edgeid : independent) {
            int var = edgeid + 2 * num_edges * !edges[edgeid].in_normal;
            for (int row = 0; row<matrix.rows(); row++) {
                if (abs(matrix.coeff(row, var)) > EPSILON) {
                    double factor = -matrix.coeff(row, var);
                    for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(matrix, row); it; ++it) {
                        if (abs(it.value()) < EPSILON) continue;
                        if (it.col() == var) continue;
                        
                        if (source_cols.find(it.col()) == source_cols.end()) {
                            if (var_cols.find(it.col()) == var_cols.end()) {
                                throw runtime_error("Could not reduce equations");
                            }
                            var_matrix.insert(new_row, var_cols[it.col()]) = it.value() / factor;
                        } else if (var_cols.find(it.col()) == var_cols.end()) {
                            if (source_cols.find(it.col()) == source_cols.end()) {
                                throw runtime_error("Could not reduce equations");
                            }
                            source_matrix.insert(new_row, source_cols[it.col()]) = it.value() / factor;
                        }
                        
                    }
                    break;
                }
            }
            
            new_row++;
        }


        var_matrix.makeCompressed();
        source_matrix.makeCompressed();
        
        if (path != "") saveSparseMatrixToCSV(var_matrix, "../assets/" + path + "/var.csv");
        if (path != "") saveSparseMatrixToCSV(source_matrix, "../assets/" + path + "/sources.csv");
        // return matrix;
    }
    
    unordered_set<int> children_nodes(int n, int avoid) {
        unordered_set<int> children = {n}; 
        stack<int> line;
        line.push(n);
        while (!line.empty()) {
            int node = line.top();
            line.pop(); 
            if (normal_neighbors.find(node) == normal_neighbors.end()) throw runtime_error("Graph not connected");
            for (auto& [m, list] : normal_neighbors[node]) {
                if (children.find(m) == children.end() and m != avoid) {
                    children.insert(m); 
                    line.push(m);
                }
            }
        }
        return children;
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
