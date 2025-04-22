#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <unordered_set>
#include <fstream>
#include <omp.h>
#include "symbols.cpp"
#include <stdexcept>


using namespace std;
using json = nlohmann::json;



const double EPSILON = 1e-15;

class ElementalEq {
    public:
    unordered_map<string, function<array<double,4>(double)>> equations;

    ElementalEq() {
        /*
        Elemental equations in the form [dA, A, dT, T] = 0
        */
        equations["m"] = [](double value) -> array<double,4> {
            return {value, 0, 0, -1};
        };
        equations["K"] = [](double value) -> array<double,4> {
            return {0.0, value, -1, 0};
        };
        equations["B"] = [](double value) -> array<double,4> {
            return {0.0, value, 0, -1};
        };
        equations["C"] = [](double value) -> array<double,4> {
            return {value, 0, 0, -1};
        };
        equations["L"] = [](double value) -> array<double,4> {
            return {0, 1, -value, 0};
        };
        equations["R"] = [](double value) -> array<double,4> {
            return {0, 1, 0, -value};
        };
    }
};

class ElementalEqSymbols {

    public:
    unordered_map<string, function<array<shared_ptr<Expression>,4>(string)>> equations;

    ElementalEqSymbols() {
        /*
        Elemental equations in the form [dv, v, dF, F]
        */
        equations["m"] = [] (string var) -> array<shared_ptr<Expression>,4> {
            return {make_shared<Symbol>(var), make_shared<Scalar>(0), make_shared<Scalar>(0), make_shared<Scalar>(-1)};
        }; 
        equations["K"] = [] (string var) -> array<shared_ptr<Expression>,4> {
            return {make_shared<Scalar>(0), make_shared<Symbol>(var), make_shared<Scalar>(-1), make_shared<Scalar>(0)};
        }; 
        equations["B"] = [] (string var) -> array<shared_ptr<Expression>,4> {
            return {make_shared<Scalar>(0), make_shared<Symbol>(var), make_shared<Scalar>(0), make_shared<Scalar>(-1)};
        }; 
        equations["C"] = [] (string var) -> array<shared_ptr<Expression>,4> {
            return {make_shared<Symbol>(var), make_shared<Scalar>(0), make_shared<Scalar>(0), make_shared<Scalar>(-1)};
        }; 
        equations["L"] = [] (string var) -> array<shared_ptr<Expression>,4> {
            return {make_shared<Scalar>(-1), make_shared<Scalar>(0), make_shared<Symbol>(var), make_shared<Scalar>(0)};
        };
        equations["R"] = [] (string var) -> array<shared_ptr<Expression>,4> {
            return {make_shared<Scalar>(0), make_shared<Scalar>(-1), make_shared<Scalar>(0), make_shared<Symbol>(var)};
        }; 
    }


};

template <typename T>
class Triplet {
    public:
        int row;
        int col;
        T value;

    Triplet(int r, int c, T v) : row(r), col(c), value(v) {}
};


template <typename T>
class SparseMatrix {
    private:
        // Row first
        unordered_map<int, unordered_map<int, T>> data_;
        int rows_;
        int cols_;
        const T default_value_;
        const T unit_;

    public:
        SparseMatrix(int rows, int cols, const T default_value, const T unit) : rows_(rows), cols_(cols), default_value_(default_value), unit_(unit) {
            for (int i=0; i<rows_; i++) {
                data_.emplace(i, unordered_map<int, T>());
            }
        }

        void fromTriplets(vector<Triplet<T>> triplets) {
            for (Triplet<T> t : triplets) {
                data_[t.row][t.col] = t.value;
            }
        }

        // look
        const T& operator()(int i, int j) {
            if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
                throw out_of_range("Index out of range");
            } 
            if (data_[i].find(j) == data_[i].end()) {
                return default_value_;
            }
            return data_[i][j];
        }

        // edit
        T& ref(int i, int j) {
            if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
                throw out_of_range("Index out of range");
            } 
            if (data_[i].find(j) == data_[i].end()) {
                data_[i][j] = default_value_;
                return data_[i][j];
            }
            return data_[i][j];
        }

        int rows() const {return rows_;}
        int cols() const {return cols_;}

        unordered_map<int, T> iterator(int i){
            return data_[i];
        }

        void toCSV(const string& filename) {
            /*
            Saves matrix to CSV in (row, col, val) format
            */
            ofstream file;
            file.open(filename);
        
            if (!file.is_open()) {
                throw runtime_error("Error opening the file! " + filename);
            }
        
            // Iterate over all non-zero elements
            for (int i = 0; i < rows(); ++i) {
                for (auto& [j, val] : iterator(i)) {
                    // Write row, column, and value
                    if (abs(val) < EPSILON) continue;
                    file << i << "," << j << "," << to_string(val) << endl;
                }
            }
        
            file.close();
        }

        
        void gaussianElimination(const vector<int>& priorityCols) {
    
            // Convert to a format more suitable for row operations
            
            for (int row = rows()-1; row >= 0; --row) {
                // Find the pivot column
                int pivotCol = -1;
                for (int id = 0; id < priorityCols.size(); ++id) {
                    if (abs(operator()(row, priorityCols[id])) > EPSILON) {
                        pivotCol = priorityCols[id];
                        break;  // Take the first suitable pivot (added break)
                    }
                }

                if (pivotCol == -1) continue; // No pivot in this column

                T pivotValue = operator()(row, pivotCol);
                for (int c=0; c<cols(); c++) {
                    if (c == pivotCol) {
                        ref(row, c) = unit_;
                        continue;
                    }
                    T val = operator()(row, c);
                    if (abs(val) > EPSILON) {
                        ref(row, c) = val / pivotValue;
                    }
                }
                for (int r = rows()-1; r >= 0; --r) {
                    if (r != row) {
                        T factor = operator()(r, pivotCol);
                        if (abs(factor) > EPSILON) {
                            for (int c=0; c<cols(); c++) {
                                if (c == pivotCol) {
                                    ref(r, c) = default_value_;
                                    continue;
                                }
                                T val = operator()(r, c);
                                if (true) {
                                    ref(r, c) = val - operator()(row, c) * factor;
                                }
                            }
                        }
                    }
                }
                // toCSV("../assets/trial"+to_string(row) + "_" + to_string(pivotCol) +".csv");
            }
        }


};

void saveVarIndex(const vector<string>& vars, const string& filename) {
    ofstream file;
    file.open(filename);

    if (!file.is_open()) {
        throw runtime_error("Error opening the file! " + filename);
    }
    for (int i=0; i<vars.size(); i++) {
        file << i << "," << vars[i] << endl ;
    }
    file.close();
}


  