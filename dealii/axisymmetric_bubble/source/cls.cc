#include "../include/CLS.h"
#include "../include/Headers.h"

// Constrained Least Squares class function definitions 

// ===========================================================================

// Function Definitions 

void print_matrix(FullMatrix<double> A) {
    unsigned int m = A.m(); unsigned int n = A.n(); 
    
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            std::cout << A(i, j) << "    \t"; 
        }
        std::cout << std::endl; 
    }
}

// ===========================================================================

// Default constructor 

// ===========================================================================

Constrained_LS::Constrained_LS() {}

// ===========================================================================

// Initialize the object (For Constrained least squares)

// ===========================================================================

void Constrained_LS::initialize(FullMatrix<double> A, FullMatrix<double> C) {

    solve_least_squares = false; 
    
    assert(A.n() == C.n());
    
    m = A.m(); n = A.n(); p = C.m(); 
    
    R.reinit(n, n); Q1.reinit(m, n); Q2.reinit(p, n); Rbar.reinit(p, p); Qbar.reinit(n, p);
    
    // Concatenate A and C int qr array 

    gsl_matrix * qr = gsl_matrix_alloc ((m+p), n);
    gsl_matrix * q = gsl_matrix_alloc ((m+p), (m+p));
    gsl_matrix * r = gsl_matrix_alloc ((m+p), n);
    
    
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            gsl_matrix_set(qr, i, j, A(i, j));  // Put A into qr 
        }
    }
    
    for (unsigned int i = m; i < (m+p); i++) {
        for (unsigned int j = 0; j < n; j++) {
            gsl_matrix_set(qr, i, j, C(i-m, j));
        }
    }
    
    // Find the QR decomposition of concatenated matrix AC and store the R matrix 
    
    int k = std::min((m+p), n); 
    
    gsl_vector * tau = gsl_vector_alloc (k);
    gsl_linalg_QR_decomp (qr, tau); 
    gsl_linalg_QR_unpack (qr, tau, q, r); 

    
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            R(i, j) = gsl_matrix_get(r, i, j);
        }
    }
    
    // Get Q1 and Q2 matrices  
    
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            Q1(i, j) =  gsl_matrix_get(q, i, j); 
        }
    }
    
    for (unsigned int i = m; i < (m+p); i++) {
        for (unsigned int j = 0; j < n; j++) {
            Q2(i-m, j) =  gsl_matrix_get(q, i, j); 
        }
    }
    
    
    
    // Get Qbar and Rbar 
    
    gsl_matrix * qr_bar = gsl_matrix_alloc (n, p);
    gsl_matrix * q_bar = gsl_matrix_alloc (n, n);
    gsl_matrix * r_bar = gsl_matrix_alloc (n, p);
    
    // Copy Q2' into qr2
    
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < p; j++) {
            gsl_matrix_set(qr_bar, i, j, Q2(j, i)); 
        }
    }
    
    k = std::min(n, p); 
    
    gsl_vector * tau_bar = gsl_vector_alloc (k);
    gsl_linalg_QR_decomp (qr_bar, tau_bar); 
    gsl_linalg_QR_unpack (qr_bar, tau_bar, q_bar, r_bar); 
    
    for (unsigned int i = 0; i < p; i++) {
        for (unsigned int j = 0; j < p; j++) {
            Rbar(i, j) = gsl_matrix_get(r_bar, i, j);
        }
    }
    
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < p; j++) {
                Qbar(i, j) = gsl_matrix_get(q_bar, i, j); 
        }
    }
    
    // Clean up all the memory 
    
    gsl_vector_free (tau); gsl_vector_free (tau_bar);
    gsl_matrix_free (qr); gsl_matrix_free (qr_bar);
    gsl_matrix_free (q); gsl_matrix_free (q_bar);
    gsl_matrix_free (r); gsl_matrix_free (r_bar);
    
}

// ===========================================================================

// Initialize the object (For Least squares)

// ===========================================================================

void Constrained_LS::initialize(FullMatrix<double> A) {

    solve_least_squares = true; 
    
    m = A.m(); n = A.n(); 
    p = 0; 
    
    R.reinit(n, n); Q1.reinit(m, n);
    
    // For qr array 

    gsl_matrix * qr = gsl_matrix_alloc (m, n);
    gsl_matrix * q = gsl_matrix_alloc (m, m);
    gsl_matrix * r = gsl_matrix_alloc (m, n);
    
    
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            gsl_matrix_set(qr, i, j, A(i, j));  // Put A into qr 
        }
    }
    
    // Find the QR decomposition of matrix A and store the R matrix 
    
    int k = std::min(m, n); 
    
    gsl_vector * tau = gsl_vector_alloc (k);
    gsl_linalg_QR_decomp (qr, tau); 
    gsl_linalg_QR_unpack (qr, tau, q, r); 

    
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            R(i, j) = gsl_matrix_get(r, i, j);
        }
    }
    
    // Get Q1 matrix  
    
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            Q1(i, j) =  gsl_matrix_get(q, i, j); 
        }
    }
    
    // Clean up all the memory 
    
    gsl_vector_free (tau); gsl_matrix_free (q); 
    gsl_matrix_free (qr); gsl_matrix_free (r);
    
}

// ===========================================================================

// Constructor taking two arguments  

// ===========================================================================

Constrained_LS::Constrained_LS(FullMatrix<double> A, FullMatrix<double> C) {
    
    initialize(A, C); 
}

// ===========================================================================

// Constructor taking one argument 

// ===========================================================================

Constrained_LS::Constrained_LS(FullMatrix<double> A) {
    
    initialize(A); 
}

// ===========================================================================

// Solve the actual constrained least squares problem 

// ===========================================================================

void Constrained_LS::solve(Vector<double> b, Vector<double> d, Vector<double>& x) {
   
    assert(b.size() == m && d.size() == p && x.size() == n);
    
    Vector<double> u(p); Vector<double> c(p); Vector<double> w(p); Vector<double> y(n); Vector<double> y_dash(n);
    
    FullMatrix<double> RbarT(p, p); 
    
    for (unsigned int i = 0; i < p; i++) {
        for (unsigned int j = 0; j < p; j++) {
            RbarT(i, j) = Rbar(j, i); 
        }
    }
    
    RbarT.forward(u, d);

    
    FullMatrix<double> Inter(p, m); 
    
    Qbar.TmTmult(Inter, Q1); 
    
    Inter.vmult(c, b); 
    
    c.sadd(1.0, -1.0, u); 
    
    Rbar.backward(w, c);
    
    Q1.Tvmult(y, b); 
    Q2.Tvmult(y_dash, w);
    
    
    y.sadd(1.0, -1.0, y_dash);

    x = 0.0; 

    R.backward(x, y); 
    
}

void Constrained_LS::solve(Vector<double> b, Vector<double>& x) {
   
    assert(b.size() == m && x.size() == n);
    
    Vector<double> y(n); 
    
    Q1.Tvmult(y, b); 
    R.backward(x, y); 
    
}

// ===========================================================================
  
