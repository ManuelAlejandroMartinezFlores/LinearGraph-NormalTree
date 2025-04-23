#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <unordered_set>
#include <stack>
#include <omp.h>
#include "state_eq.cpp"
#include <memory>
#include <filesystem>

using namespace std;
using json = nlohmann::json;


filesystem::path join_paths(const std::vector<std::string>& parts) {
    filesystem::path result;
    for (const auto& part : parts) {
        result /= part;
    }
    return result;
}


class Edge {
    public:
    int id;
    int source_node;
    int target_node;
    string type;
    string constant_type;
    string constant;
    double value;
    bool across_source;
    bool through_source;
    string across;
    string through;
    bool in_normal;

    Edge() {
        id = -1;
        source_node = 0; 
        target_node = 0;
        type = "A";
        constant_type = "";
        value = 0;
        across_source = false;
        through_source = false;
        across = "";
        through = "";
        constant = "";
        in_normal = false;
    }

    Edge copy() {
        Edge ans;
        ans.id = id;
        ans.source_node = source_node;
        ans.target_node = target_node;
        ans.type = type;
        ans.constant_type = constant_type;
        ans.constant = constant;
        ans.value = value;
        ans.across_source = across_source;
        ans.through_source = through_source;
        ans.across = across;
        ans.through = through;
        ans.in_normal = in_normal;
        return ans;
    }

    friend std::ostream& operator<<(std::ostream& os, const Edge& e) {
        return os << e.source_node << " - " << e.target_node << " - " << e.type << " - normal: " << e.in_normal;
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Edge, source_node, target_node, type, constant_type, value, across_source, through_source, across, through, in_normal, constant);

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
    unordered_map<int, unordered_map<int,int>> normal_neighbors;
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


    void remove_edge(int edgeid) {
        Edge edge = edges[edgeid];
        edgelist[edge.source_node][edge.target_node].erase(
            find(
                edgelist[edge.source_node][edge.target_node].begin(),
                edgelist[edge.source_node][edge.target_node].end(),
                edgeid
            )
        );
        edges.erase(edges.begin() + edgeid);
        for (int i=edgeid; i<edges.size(); i++) {
            Edge e = edges[i];
            edgelist[e.source_node][e.target_node].erase(
                find(
                    edgelist[e.source_node][e.target_node].begin(),
                    edgelist[e.source_node][e.target_node].end(),
                    e.id
                )
            );
            edges[i].id = i;
            edgelist[e.source_node][e.target_node].emplace_back(i);
        }
        num_edges = edges.size();
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
        edges.push_back(edge);
        normal_tree_built = false;
        num_nodes = rank.size();
        num_edges = edges.size();
    }

    void simplify_edges() {
        for (Edge& edge : edges) {
            if (edge.type == "A") {
                // Find A edges in series 
                // Finde next edge
                if (edgelist.find(edge.target_node) != edgelist.end()) {
                    if (edgelist[edge.target_node].size() == 1) {
                        for (auto& [target, list] : edgelist[edge.target_node]) {
                            if (list.size() == 1) {
                                int edgeid = list[0];
                                if (edges[edgeid].type == "A") {
                                    Edge newedge = edges[edgeid].copy();
                                    newedge.across += edge.across;
                                    newedge.through += edge.through;
                                    newedge.source_node = edge.source_node;
                                    newedge.value = 1 / (1 / newedge.value + 1 / edge.value);
                                    newedge.constant = "1/(1/" + newedge.constant + "+1/" + edge.constant + ")";
                                    int otherid = edge.id;
                                    if (edgeid < otherid) {
                                        otherid = otherid - 1;
                                    }
                                    remove_edge(edgeid);
                                    remove_edge(otherid);
                                    add_edge(newedge);
                                    return simplify_edges();
                                }
                            }
                        }
                    }
                }
                // Find previous edge
                int sourceid = -1;
                int cnt = 0;
                for (auto& [source, list] : edgelist) {
                    if (list.find(edge.source_node) != list.end()) {
                        sourceid = source;
                        cnt++;
                    }
                }
                if (cnt == 1) {
                    if (edgelist[sourceid][edge.source_node].size() == 1) {
                        int edgeid = edgelist[sourceid][edge.source_node][0];
                        if (edges[edgeid].type == "A") {
                            Edge newedge = edges[edgeid].copy();
                            newedge.across += edge.across;
                            newedge.through += edge.through;
                            newedge.target_node = edge.target_node;
                            newedge.value = 1 / (1 / newedge.value + 1 / edge.value);
                            newedge.constant = "1/(1/" + newedge.constant + "+1/" + edge.constant + ")";
                            int otherid = edge.id;
                            if (edgeid < otherid) {
                                otherid = otherid - 1;
                            }
                            remove_edge(edgeid);
                            remove_edge(otherid);
                            add_edge(newedge) ;
                            return simplify_edges();
                        }
                    }
                }
            }
            if (edge.type == "T") {
                // Find parallel T edges
                if (edgelist.find(edge.source_node) != edgelist.end()) {
                    if (edgelist[edge.source_node].find(edge.target_node) != edgelist[edge.source_node].end()) {
                        for (int edgeid : edgelist[edge.source_node][edge.target_node]) {
                            if (edgeid != edge.id && edges[edgeid].type == "T") {
                                edges[edgeid].value = 1 / (1 / edges[edgeid].value + 1 / edge.value);
                                edges[edgeid].constant = "1/(1/" + edges[edgeid].constant + "+1/" + edge.constant + ")";
                                remove_edge(edge.id);
                                return simplify_edges();
                            }
                        }
                    }
                }
            }
        }
    }

    void save_to_json(string path) {
        /*
        saves to json file
        */
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
        /*
        Builds the normal tree:
        - Include all across sources
        - Include all possible "A" types (those included become independent components)
        - Include all possible "D" types
        - Include all possible "T" types (those not included become independent components)
        - Check that no through sources can be included
        Moreover, it assigns row to each variable in the state equation
        */
        simplify_edges();
        int tree_branches = 0;
        independent = {};
        sources = {};
        elemental_rows = {};
        num_sources_across = 0;
        num_sources_through = 0;
        num_edges = edges.size();
        num_nodes = 0;

        unordered_set<int> diff_nodes;
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
            if (edge.across_source or edge.through_source) {
                sources.insert(edge.id);
            } else {
                elemental_rows[edge.id] = elemental_rows.size();
            }
            if (edge.across_source) {
                ++num_sources_across;
            }
            if (edge.through_source) {
                ++num_sources_through;
            }
            diff_nodes.insert(edge.source_node);
            diff_nodes.insert(edge.target_node);
        }
        num_nodes = diff_nodes.size();
        // All posible A types
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
        // All posible D types
        for (Edge& edge : edges) {
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
        // All possible "T" types
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

        // if (tree_branches != num_nodes - 1) {
        //     throw runtime_error("Could not build normal tree");
        // }
        normal_tree_built = true;
    }

    void generate_state_eq(string path, bool verbose=false) {
        /*
        Generates state equations 
        - Writes B-S elemental equations for all non-source edges
        - Writes N-1-SA continuity equations by getting all the child nodes from both source and target nodes and assuring the flux in between is constant
        - Writes B-N+1-ST compatibility equations by inserting a link into the tree and assuring the path from source to target in the tree equals the link
        Then it solves the system using Gaussian elimination and selects the independent variables to get dx = Ax + Bs
        */
        if (!normal_tree_built) {
            build_normal();
        }
        if (verbose) cout << "Normal tree built" << endl;
        
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
            // Get all children nodes and assert equal flux
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
            if (edge.in_normal || edge.through_source) continue;
            int thread_id = omp_get_thread_num();
            int row = compatibility_rows[edge.id];
            vector<Triplet<double>> &local_triplets = thread_triplets[thread_id];
            local_triplets.emplace_back(row, edge.id + num_edges, -1);
            local_triplets.emplace_back(row+1, edge.id, -1);
            // Get path and assure equality
            for (auto& [edgeid, value] : get_path(edge.source_node, edge.target_node)) {
                local_triplets.emplace_back(row, edgeid + num_edges, value);
                local_triplets.emplace_back(row+1, edgeid, value);
            }
        }
        for (const auto &local : thread_triplets) {
            triplets.insert(triplets.end(), local.begin(), local.end());
        }

        if (verbose) cout << "Compatibility equations built" << endl;

        // Build matrix
        SparseMatrix<double> matrix(3*num_edges-3*(num_sources_across+num_sources_through), 4*num_edges, 0, 1);
        matrix.fromTriplets(triplets);

        if (path != "") matrix.toCSV(join_paths({path, "eqs.csv"}));
        if (verbose) cout << "Equations matrix built" << endl;

        // Define priority of variable elimination in the Gaussian method
        vector<int> priorityVars;
        

        // Non-sources and non-independent primary 
        for (Edge& edge : edges) {
            if ((independent.find(edge.id) == independent.end()) && (sources.find(edge.id) == sources.end())) {
                int var = edge.id + 2 * !edge.in_normal * num_edges;
                priorityVars.push_back(var);
                priorityVars.push_back(var + num_edges);
            }  
        }

        

        // Non-sources and non-independent secondary
        for (Edge& edge : edges) {
            if ((independent.find(edge.id) == independent.end()) || (sources.find(edge.id) == sources.end())) {
                int var = edge.id + 2 * edge.in_normal * num_edges;
                priorityVars.push_back(var);
                priorityVars.push_back(var + num_edges);
            }  
        }

        


        // Independent secondary
        for (int edgeid : independent) {
            int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            priorityVars.push_back(var + num_edges);
        }

        // Independent primary
        for (int edgeid : independent) {
            int var = edgeid + 2 * !edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            // priorityVars.push_back(var + num_edges);
        }

        

        // Sources
        for (int edgeid : sources) {
            int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            priorityVars.push_back(var + num_edges);
            var = edgeid + 2 * !edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            priorityVars.push_back(var + num_edges);
        }
        
        

        // Execute gaussian elimination
        matrix.gaussianElimination(priorityVars);

        if (verbose) cout << "Equations reduced" << endl;
        if (path != "") matrix.toCSV(join_paths({path, "reduced.csv"}));



        SparseMatrix<double> var_matrix(independent.size(), independent.size(), 0, 1);
        SparseMatrix<double> source_matrix(independent.size(), sources.size()*2, 0, 1);
        unordered_map<int, int> var_cols;

        // Identify independent variables
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

        // Only include independent variables in state equations
        saveVarIndex(x, join_paths({path, "var_index.csv"}));
        saveVarIndex(s, join_paths({path, "sources_index.csv"}));
        int new_row = 0;
        for (int edgeid : independent) {
            int var = edgeid + 2 * num_edges * !edges[edgeid].in_normal;
            for (int row = 0; row<matrix.rows(); row++) {
                if (abs(matrix(row, var)) > EPSILON) {
                    auto factor = - matrix(row, var);
                    for (auto& [col, value] : matrix.iterator(row)) {
                        if (abs(value) < EPSILON) continue;
                        if (col == var) continue;
                        if (source_cols.find(col) == source_cols.end()) {
                            if (var_cols.find(col) == var_cols.end()) {
                                // cout << "var " << var << endl;
                                // cout << row << ", " << col << endl;
                                throw runtime_error("Could not reduce equations");
                            }
                            var_matrix.ref(new_row, var_cols[col]) = value / factor;
                        } else if (var_cols.find(col) == var_cols.end()) {
                            if (source_cols.find(col) == source_cols.end()) {
                                throw runtime_error("Could not reduce equations");
                            }
                            source_matrix.ref(new_row, source_cols[col]) = value / factor;
                        }
                    }
                    break;
                }
            }
            
            new_row++;
        }


        // var_matrix.makeCompressed();
        // source_matrix.makeCompressed();
        
        if (path != "") var_matrix.toCSV(join_paths({path, "var.csv"}));
        if (path != "") source_matrix.toCSV(join_paths({path, "sources.csv"}));
    }
    
    void generate_state_eq_symbolic(string path, bool verbose=false) {
       /*
        Generates state equations 
        - Writes B-S elemental equations for all non-source edges
        - Writes N-1-SA continuity equations by getting all the child nodes from both source and target nodes and assuring the flux in between is constant
        - Writes B-N+1-ST compatibility equations by inserting a link into the tree and assuring the path from source to target in the tree equals the link
        Then it solves the system using Gaussian elimination and selects the independent variables to get dx = Ax + Bs
        */
        if (!normal_tree_built) {
            build_normal();
        }
        if (verbose) cout << "Normal tree built" << endl;
        
        vector<Triplet<shared_ptr<Expression>>> triplets;  // Temporary storage
        vector<vector<Triplet<shared_ptr<Expression>>>> thread_triplets(omp_get_max_threads());
        if (verbose) cout << "Threads: " << omp_get_max_threads() << endl;
        

        // Elemental equations
        ElementalEqSymbols eqs; 
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < edges.size(); ++i) {
            Edge edge = edges[i];
            if (!edge.across_source and !edge.through_source) {
                int thread_id = omp_get_thread_num();
                int row = elemental_rows[edge.id];
                vector<Triplet<shared_ptr<Expression>>> &local_triplets = thread_triplets[thread_id];
                array<shared_ptr<Expression>,4> vals = eqs.equations[edge.constant_type](edge.constant);
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
            vector<Triplet<shared_ptr<Expression>>> &local_triplets = thread_triplets[thread_id];
            // Get all children nodes and assert equal flux
            for (int n : children_nodes(edges[edgeid].source_node, edges[edgeid].target_node)) {
                for (int m : children_nodes(edges[edgeid].target_node, edges[edgeid].source_node)) {
                    if (edgelist.find(n) != edgelist.end() && edgelist[n].find(m) != edgelist[n].end()){
                        for (int edgeidx : edgelist[n][m]) {
                            local_triplets.emplace_back(row, edgeidx + 3*num_edges, make_shared<Scalar>(1));
                            local_triplets.emplace_back(row+1, edgeidx + 2*num_edges, make_shared<Scalar>(1));
                        }
                    }
                    if (edgelist.find(m) != edgelist.end() && edgelist[m].find(n) != edgelist[m].end()){
                        for (int edgeidx : edgelist[m][n]) {
                            local_triplets.emplace_back(row, edgeidx + 3*num_edges, make_shared<Scalar>(-1));
                            local_triplets.emplace_back(row+1, edgeidx + 2*num_edges, make_shared<Scalar>(-1));
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
            if (edge.in_normal || edge.through_source) continue;
            int thread_id = omp_get_thread_num();
            int row = compatibility_rows[edge.id];
            vector<Triplet<shared_ptr<Expression>>> &local_triplets = thread_triplets[thread_id];
            local_triplets.emplace_back(row, edge.id + num_edges, make_shared<Scalar>(-1));
            local_triplets.emplace_back(row+1, edge.id, make_shared<Scalar>(-1));
            // Get path and assure equality
            for (auto& [edgeid, value] : get_path(edge.source_node, edge.target_node)) {
                local_triplets.emplace_back(row, edgeid + num_edges, make_shared<Scalar>(value));
                local_triplets.emplace_back(row+1, edgeid, make_shared<Scalar>(value));
            }
        }
        for (const auto &local : thread_triplets) {
            triplets.insert(triplets.end(), local.begin(), local.end());
        }

        if (verbose) cout << "Compatibility equations built" << endl;

        // Build matrix
        SparseMatrix<shared_ptr<Expression>> matrix(3*num_edges-3*(num_sources_across+num_sources_through), 4*num_edges, make_shared<Scalar>(0), make_shared<Scalar>(1));
        matrix.fromTriplets(triplets);

        if (path != "") matrix.toCSV(join_paths({path, "eqs.csv"}));
        if (verbose) cout << "Equations matrix built" << endl;

        // Define priority of variable elimination in the Gaussian method
        vector<int> priorityVars;
        

        // Non-sources and non-independent primary 
        for (Edge& edge : edges) {
            if ((independent.find(edge.id) == independent.end()) && (sources.find(edge.id) == sources.end())) {
                int var = edge.id + 2 * !edge.in_normal * num_edges;
                priorityVars.push_back(var);
                priorityVars.push_back(var + num_edges);
            }  
        }

        

        // Non-sources and non-independent secondary
        for (Edge& edge : edges) {
            if ((independent.find(edge.id) == independent.end()) || (sources.find(edge.id) == sources.end())) {
                int var = edge.id + 2 * edge.in_normal * num_edges;
                priorityVars.push_back(var);
                priorityVars.push_back(var + num_edges);
            }  
        }

        


        // Independent secondary
        for (int edgeid : independent) {
            int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            priorityVars.push_back(var + num_edges);
        }

        // Independent primary
        for (int edgeid : independent) {
            int var = edgeid + 2 * !edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            // priorityVars.push_back(var + num_edges);
        }

        

        // Sources
        for (int edgeid : sources) {
            int var = edgeid + 2 * edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            priorityVars.push_back(var + num_edges);
            var = edgeid + 2 * !edges[edgeid].in_normal * num_edges;
            priorityVars.push_back(var);
            priorityVars.push_back(var + num_edges);
        }
        
        

        // Execute gaussian elimination
        matrix.gaussianElimination(priorityVars);

        if (verbose) cout << "Equations reduced" << endl;
        if (path != "") matrix.toCSV(join_paths({path, "reduced.csv"}));



        SparseMatrix<shared_ptr<Expression>> var_matrix(independent.size(), independent.size(), make_shared<Scalar>(0), make_shared<Scalar>(1));
        SparseMatrix<shared_ptr<Expression>> source_matrix(independent.size(), sources.size()*2, make_shared<Scalar>(0), make_shared<Scalar>(1));
        unordered_map<int, int> var_cols;

        // Identify independent variables
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

        // Only include independent variables in state equations
        saveVarIndex(x, join_paths({path, "var_index.csv"}));
        saveVarIndex(s, join_paths({path, "sources_index.csv"}));
        int new_row = 0;
        for (int edgeid : independent) {
            int var = edgeid + 2 * num_edges * !edges[edgeid].in_normal;
            for (int row = 0; row<matrix.rows(); row++) {
                if (abs(matrix(row, var)) > EPSILON) {
                    auto factor = - matrix(row, var);
                    for (auto& [col, value] : matrix.iterator(row)) {
                        if (abs(value) < EPSILON) continue;
                        if (col == var) continue;
                        if (source_cols.find(col) == source_cols.end()) {
                            if (var_cols.find(col) == var_cols.end()) {
                                // cout << "var " << var << endl;
                                // cout << row << ", " << col << endl;
                                throw runtime_error("Could not reduce equations");
                            }
                            var_matrix.ref(new_row, var_cols[col]) = value / factor;
                        } else if (var_cols.find(col) == var_cols.end()) {
                            if (source_cols.find(col) == source_cols.end()) {
                                throw runtime_error("Could not reduce equations");
                            }
                            source_matrix.ref(new_row, source_cols[col]) = value / factor;
                        }
                    }
                    break;
                }
            }
            
            new_row++;
        }


        // var_matrix.makeCompressed();
        // source_matrix.makeCompressed();
        
        if (path != "") var_matrix.toCSV(join_paths({path, "var.csv"}));
        if (path != "") source_matrix.toCSV(join_paths({path, "sources.csv"}));
    }

    unordered_set<int> children_nodes(int n, int avoid) {
        /*
        Get all children nodes in the tree avoiding one node
        */
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
        /*
        Gets the path in the normal tree using DFS
        */
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
};


void solve_system(string path, bool verbose=false) {
    LinearGraph g = graph_from_json(join_paths({path, "edges.json"}));
    if (verbose) cout << "Loaded from json" << endl;
    g.build_normal();
    if (verbose) cout << "Built normal" << endl;
    g.save_to_json(join_paths({path, "normal_tree.json"}));
    g.generate_state_eq(path, verbose);
    if (verbose) cout << "Generate equations" << endl;
    // g.generate_state_eq_symbolic(join_paths({path, "symb"}), false);
    
};