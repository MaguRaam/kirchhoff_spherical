#ifndef CLS_H_
#define CLS_H_

#include "Headers.h"  

// ===========================================================================

// Class Declarations 

// ===========================================================================
 
class Constrained_LS; 

// Solve a constrained least squares problem: 
//     *Minimize r(x) = ||Ax - b||^2
//     *Subject to: Cx = d 
// Using QR factorization 

// ===========================================================================

// Class Definitions 

// ===========================================================================

class Constrained_LS {

    unsigned int m, n, p;
    bool solve_least_squares;  
    
    FullMatrix<double> Q1; 
    FullMatrix<double> Q2;
    FullMatrix<double> R; 
    FullMatrix<double> Qbar; 
    FullMatrix<double> Rbar; 
    
public:
    Constrained_LS();                // Default constructor 
    
    Constrained_LS(FullMatrix<double>, FullMatrix<double>);  // Takes two matrices (Constrained Least Squares)
    Constrained_LS(FullMatrix<double>); // Take one matrix (Least-squares)  
    
    void initialize(FullMatrix<double>, FullMatrix<double>); // Initialize the objetc - same as constructor 2
    void initialize(FullMatrix<double>); // Least squares initializer 
    
    void solve(Vector<double>, Vector<double>, Vector<double>&);
    void solve(Vector<double>, Vector<double>&);  
};

#endif /* CLS_H_ */
