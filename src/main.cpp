#include <cwchar>
// #include "build_tests.cpp"
#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <omp.h>
#include "gen_circuits.cpp"
#include <filesystem>

using namespace std;
using json = nlohmann::json;


int main(int argc, char** argv) {
    // test(argv[1]); 
    for (int i=0; i<100; i++) {
        // filesystem::create_directory(join_paths({"..", "assets", "data", "squares", "sq" + to_string(i)}));
        cout << i << endl;
        // gen_square(join_paths({"..", "assets", "data", "squares", "sq" + to_string(i), "edges.json"}));
        solve_system(join_paths({"..", "assets", "data", "squares", "sq" + to_string(i)}));
        
    }
    // solve_system(join_paths({"..", "assets", "data", "squares", "sq4"}), true);
    return 0;
}




