#include <iostream>
#include "normal_tree.cpp"
#include <vector>
#include <ctime>
#include <random>
#include <filesystem>

using namespace std;

class Randomizer {
    private:
        uniform_real_distribution<> dis;
        mt19937 gen;
    
    public:
        Randomizer() {
            random_device rd;
            gen = mt19937(rd());
            dis = uniform_real_distribution<>(0, 1);
        }
        double uniform(double low=0, double high=1) {
            return dis(gen) * (high - low) + low;
        }

        int random_choice(vector<double> probs) {
            int N = probs.size();
            double cumprob = 0; 
            double prob = uniform();
            for (int i=0; i<N; i++) {
                cumprob += probs[i];
                if (prob <= cumprob) {
                    return i;
                }
            }
            throw runtime_error("Probability distribution not valid");
        }
};

class Components {
    public:
        vector<string> comp;
        vector<string> types;
        vector<double> dist;
        Components() {
            comp = {"R", "L", "C"};
            types = {"D", "T", "A"};
            double p = 1 / static_cast<double>(comp.size());
            for (auto c : comp) {
                dist.emplace_back(p);
            }
        }
};


void gen_square(string path) {
    vector<Edge> ans;
    Randomizer rd;
    Components comp;
    for (int i=0; i<4; i++) {
        Edge e;
        e.source_node = i;
        e.target_node = (i + 1) % 4;
        if (i == 1) {
            e.across_source = true; 
            e.type = "S";
            e.constant_type = "S";
            e.value = rd.uniform(0, 10);
            e.across = "Vs";
            e.through = "is";
        } else {
            int id = rd.random_choice(comp.dist);
            e.type = comp.types[id];
            e.constant_type = comp.comp[id];
            e.constant = e.constant_type + to_string(i);
            e.value = rd.uniform(0, 10);
            e.across = "V" + e.constant;
            e.through = "i" + e.constant;
        }
        ans.emplace_back(e);
    }
    LinearGraph g(ans);
    g.save_to_json(path);

}

