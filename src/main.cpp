#include <cwchar>
#include "build_tests.cpp"
#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>

using namespace std;
using json = nlohmann::json;


int main(int argc, char** argv) {
    // vector<int> sources = {0, 1, 1};
    // vector<int> targets = {1, 0, 0};
    // save_graph(sources, targets, "test10");
    test(argv[1]); 
    return 0;
}

