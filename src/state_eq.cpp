#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/SparseLU>
#include <unordered_set>
#include <fstream>

using namespace std;
using json = nlohmann::json;
using namespace Eigen;

const double EPSILON = 1e-8;

class ElementalEq {
    public:
    unordered_map<string, function<array<double,4>(double)>> equations;

    ElementalEq() {
        equations["m"] = [](double value) -> array<double,4> {
            return {value, 0, 0, -1};
        };
        equations["K"] = [](double value) -> array<double,4> {
            return {0.0, value, -1, 0};
        };
        equations["B"] = [](double value) -> array<double,4> {
            return {0.0, value, 0, -1};
        };
    }

    // void load_equations(string path) {
    //     ifstream file("../assets/" + path);
    //     json j;
    //     file >> j;
    //     file.close();
    //     equations = j.get<unordered_map<string, function<array<double,4>(double)>>>();

    // }
};
 
    void saveSparseMatrixToCSV(const SparseMatrix<double>& mat, const string& filename) {
        ofstream file;
        file.open(filename);
    
        if (!file.is_open()) {
            cerr << "Error opening the file!" << endl;
            return;
        }
    
        // Iterate over all non-zero elements
        for (int k = 0; k < mat.outerSize(); ++k) {
            for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
                // Write row, column, and value
                if (abs(it.value()) < EPSILON) continue;
                file << it.row() << "," << it.col() << "," << it.value() << endl;
            }
        }
    
        file.close();
    }

// Function to save a SparseMatrix to a CSV file (with each row of the matrix as a row in the CSV)
    void saveSparseMatrixRowWise(const SparseMatrix<double>& mat, const string& filename) {
        ofstream file;
        file.open(filename);
        if (!file.is_open()) {
            cerr << "Error opening the file!" << endl;
            return;
        }
        // Iterate over all rows of the matrix
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                // Write the value (zero if not explicitly set)
                double value = mat.coeff(i, j); // Returns 0 if the element is zero
                file << value;  // Write the value to the CSV file
    
                // Avoid trailing comma after the last element in the row
                if (j < mat.cols() - 1) {
                    file << ","; // Separate values with commas
                }
            }
            file << endl; // End of the row in the CSV
        }
        file.close();
    }
    



 // Function to perform Gaussian elimination on a sparse matrix
SparseMatrix<double> gaussianElimination(const SparseMatrix<double>& A, const vector<int>& priorityCols) {
    SparseMatrix<double> mat = A; // Copy the matrix
    int rows = mat.rows();
    int cols = mat.cols();



    for (int row = rows-1; row >= 0; --row) {
        // Find the pivot column

        int pivotCol = -1;
        for (int id=0; id<priorityCols.size(); ++id) {
            if (abs(mat.coeff(row, priorityCols[id])) > EPSILON) {
                pivotCol = priorityCols[id];
            }
        }

        if (pivotCol == -1) continue; // No pivot in this column

        // Normalize the pivot row
        double pivotValue = mat.coeff(row, pivotCol);
        for (int c = 0; c < cols; ++c) {
            mat.coeffRef(row, c) /= pivotValue;
        }

        // Eliminate the current column in other rows
        for (int r = 0; r < rows; ++r) {
            double factor = mat.coeff(r, pivotCol);
            if (r != row && abs(factor) > EPSILON) {
                for (int c = 0; c < cols; ++c) {
                    mat.coeffRef(r, c) -= factor * mat.coeff(row, c);
                }
            }
        }
    }

    // Return the reduced matrix
    return mat;
}

  