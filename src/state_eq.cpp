#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/SparseLU>
#include <unordered_set>
#include <fstream>
#include <omp.h>

using namespace std;
using json = nlohmann::json;
using namespace Eigen;

const double EPSILON = 1e-8;

class ElementalEq {
    public:
    unordered_map<string, function<array<double,4>(double)>> equations;

    ElementalEq() {
        /*
        Elemental equations in the form [dv, v, dF, F]
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
    }
};

 
    void saveSparseMatrixToCSV(const SparseMatrix<double>& mat, const string& filename) {
        /*
        Saves matrix to CSV in (row, col, val) format
        */
        ofstream file;
        file.open(filename);
    
        if (!file.is_open()) {
            throw runtime_error("Error opening the file! " + filename);
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

    void saveVarIndex(const vector<string>& vars, const string& filename) {
        /*
        Saves the index of variables in a (id, var) format
        */
        ofstream file(filename);
        if (!file.is_open()) {
            throw runtime_error("Could not open file: " + filename);
        }
        for (int i=0; i<vars.size(); ++i) {
            file << i << "," << vars[i] << endl;
        }
        file.close();
    }



SparseMatrix<double> gaussianElimination(const SparseMatrix<double>& A, const vector<int>& priorityCols) {
    SparseMatrix<double> mat = A; // Copy the matrix
    int rows = mat.rows();
    int cols = mat.cols();
    
    // Convert to a format more suitable for row operations
    mat.makeCompressed();
    
    for (int row = rows-1; row >= 0; --row) {
        // Find the pivot column
        int pivotCol = -1;
        for (int id = 0; id < priorityCols.size(); ++id) {
            if (abs(mat.coeff(row, priorityCols[id])) > EPSILON) {
                pivotCol = priorityCols[id];
                break;  // Take the first suitable pivot (added break)
            }
        }

        if (pivotCol == -1) continue; // No pivot in this column

        // Extract the current row as a dense vector for efficient access
        VectorXd currentRow = mat.row(row);
        double pivotValue = currentRow(pivotCol);
        
        // Normalize the pivot row
        currentRow /= pivotValue;
        
        // Update the matrix row (all at once)
        for (int c = 0; c < cols; ++c) {
            mat.coeffRef(row, c) = currentRow(c);
        }

        // Process other rows in parallel, but each thread works on a different row
        #pragma omp parallel for
        for (int r = 0; r < rows; ++r) {
            if (r != row) {
                double factor = mat.coeff(r, pivotCol);
                if (abs(factor) > EPSILON) {
                    // Extract row to dense vector for efficient operations
                    VectorXd rowToModify = mat.row(r);
                    // Subtract appropriately scaled pivot row
                    rowToModify -= factor * currentRow;
                    
                    // Update the matrix row (all at once)
                    for (int c = 0; c < cols; ++c) {
                        mat.coeffRef(r, c) = rowToModify(c);
                    }
                }
            }
        }
    }
    
    // Compress the matrix again before returning
    mat.makeCompressed();
    return mat;
}

  