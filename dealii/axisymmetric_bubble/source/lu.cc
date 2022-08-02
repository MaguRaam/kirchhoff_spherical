#include "../include/LU.h"
#include "../include/Headers.h"

// ===========================================================================

// LUdcmp

// ===========================================================================

// Default Constructor

// ===========================================================================

LUdcmp::LUdcmp() {}

// ===========================================================================

// Initialize Function

// ===========================================================================


void LUdcmp::initialize(FullMatrix<double> A) {

    assert(A.m() == A.n());

    n = A.m();

    // Variables
    const double TINY = 1.0e-40;
    unsigned int i, imax, j, k;
    double big, temp;
    Vector<double> vv(n);
    d = 1.0;
    FullMatrix<double> Dummy(n, n);
    indx.resize(n);
    LU.copy_from(Dummy);
    LU = A;

    // Loop over rows to get the implicit scaling information

    for (i=0; i<n; i++) {

        big = 0.0;

        for (j=0; j<n; j++) {

            if ((temp=std::abs(LU(i,j))) > big) {
                big=temp;
            }
        }

        if (big == 0.0) {
            std::cerr << "Singular matrix in LUdcmp" << std::endl;
            std::exit(EXIT_FAILURE); // No nonzero largest element.
        }

        vv[i] = 1.0/big; // Save the scaling

    }

    for (k = 0; k < n; k++) {

        big = 0.0;  // Initialize for the search for largest pivot element

        for (i = k; i < n; i++) {

            temp=vv[i]*std::abs(LU(i,k));

            if (temp > big) {

                big = temp;
                imax = i;
            }

        }

        if (k != imax) {

            for (j = 0; j < n; j++) {

                temp = LU(imax, j);
                LU(imax, j) = LU(k,j);
                LU(k, j) = temp;


            }

            d = -d;
            vv[imax] = vv[k];
        }

        indx[k] = imax;
        if (LU(k,k) == 0.0) {
            LU(k, k)=TINY;
        }

        for (i=k+1;i<n;i++)  {

            temp=LU(i,k)/=LU(k,k); // Divide by pivot element

            for (j=k+1;j<n;j++) {

                LU(i, j) -= temp*LU(k ,j);

            }

        }

    }

}

// ===========================================================================

// Constructor taking one matrix as argument

// ===========================================================================

LUdcmp::LUdcmp(FullMatrix<double> A) {

    initialize(A); 
}

// ===========================================================================

// Solve for single RHS

// ===========================================================================

void LUdcmp::solve(Vector<double> b, Vector<double>& x) {

    if (b.size() != n || x.size() != n) {
        std::cerr << "Matrix and Vector sizes donot match" << std::endl;
        std::exit(EXIT_FAILURE); // No nonzero largest element.
    }

    unsigned int i, ii=0, ip, j;

    double sum;

    for (i=0; i<n; i++) {

        x[i] = b[i];
    }

    for (i=0; i<n; i++) {

        ip=indx[i];
        sum=x[ip];
        x[ip]=x[i];

        if (ii != 0) {
            for (j = ii-1; j<i; j++) {
                sum -= LU(i,j)*x[j];
            }
        }

        else if (sum != 0.0) {
            ii=i+1;
        }

        x[i] = sum;
    }

    for (int k=n-1;k>=0;k--) { // Back substitution

        sum=x[k];
        for (j=k+1;j<n;j++) {
            sum -= LU(k,j)*x[j];
        }
        x[k]=sum/LU(k,k);
    }

}

// ===========================================================================

// Solve for Multiple RHS

// ===========================================================================

void LUdcmp::solve(FullMatrix<double> B, FullMatrix<double>& X) {

    if (B.m() != n || X.n() != n || B.n() != X.n()) {
        std::cerr << "Matrix sizes donot match" << std::endl;
        std::exit(EXIT_FAILURE); // No nonzero largest element.
    }

    Vector<double> xx(n);

    unsigned int i, j, m = B.m();

    for (j = 0;j < m; j++) {
        for (i=0; i<n; i++) {
            xx[i] = B(i,j);
        }

        solve(xx,xx);

        for (i=0;i<n;i++) {
            X(i,j) = xx[i];
        }
    }
}

// ===========================================================================

// Invert the matrix A

// ===========================================================================

void LUdcmp::inverse(FullMatrix<double>& A_inv) {

    unsigned int i,j;

    if (A_inv.n() != n || A_inv.m() != n) {
        std::cerr << "Matrix sizes donot match" << std::endl;
        std::exit(EXIT_FAILURE); // No nonzero largest element.
    }


    for (i=0;i<n;i++) {

        for (j=0;j<n;j++) {
            A_inv(i, j) = 0.;
        }

        A_inv(i, i) = 1.;
    }

    solve(A_inv,A_inv);
}

// ===========================================================================

// Determinant the matrix A

// ===========================================================================

double LUdcmp::det() {

    //WARNING: For large systems, the determinant of the matrix may overflow or underflow

    double dd = d;

    for (unsigned int i=0; i < n; i++) {
        dd *= LU(i, i);
    }

    return dd;

} 
