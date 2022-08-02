#ifndef LA_SOLVERS_H_
#define LA_SOLVERS_H_

#include "Headers.h"
 
// ===========================================================================

//                       Class Declarations

// ===========================================================================

class LUdcmp;

// Class for LU Decomposition with partial pivoting

// ===========================================================================

class LUdcmp {
    unsigned int n;                        // Size of the matrix
    FullMatrix<double> LU;        // Stores the decomposition
    std::vector<int> indx;        // Stores the permutation
    double d;                     // Used in determinant calculation
public:
    LUdcmp();                     // Default Constructor
    LUdcmp(FullMatrix<double>);               // Constructor
    void initialize(FullMatrix<double>);      // Initialize with a given Matrix
    void solve(Vector<double>, Vector<double>&);  // Solve for sing RHS
    void solve(FullMatrix<double>, FullMatrix<double>&);  // Solve for multiple RHS
    void inverse(FullMatrix<double>&);        // Invert the matrix A
    double det();
};

// ===========================================================================

//                   Class Function Definitions




#endif /* LA_SOLVERS_H_ */
