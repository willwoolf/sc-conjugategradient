#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <ratio>

#include "mvector.h"
#include "mmatrix.h"

using namespace std;
using namespace chrono;

MVector conjugateGradientSolve(const MMatrix&, const MVector&, const MVector&, const string filename);
void computeLaplacian(int n, const string& filename);

int main()
{  
    MVector b(5, 1.0/36.0);

    MMatrix A(5, 5);
    makeTridiagonal(A, 2, -1);

    // starting guess is elementary basis vector
    MVector x_0 = {1,0,0,0,0};

    //regular 5x5 solve (3.35-3.36)
    MVector x = conjugateGradientSolve(A, b, x_0, "tridiagonal_convergence_0611.txt");
    cout << x << endl;

    MMatrix A1(5, 5), A2(5, 5), A3(5, 5);
    makeTridiagonal(A1, 2+1, -1);
    makeTridiagonal(A2, 2+4, -1);
    makeTridiagonal(A3, 2+10, -1);

    MVector b1(5, 2.5);

    // perform CG for different values of alpha (3.37)
    x = conjugateGradientSolve(A1, b1, x_0, "tridiagonal_a1_conv.txt");
    x = conjugateGradientSolve(A2, b1, x_0, "tridiagonal_a4_conv.txt");
    x = conjugateGradientSolve(A3, b1, x_0, "tridiagonal_a10_conv.txt");

    MMatrix P(16, 16);
    makeLaplacian(P, 4, -1);
    cout << P << endl;

    // Laplace's equation (3.39)
    computeLaplacian(5, "lap_conv_n5_0611.txt");
    computeLaplacian(6, "lap_conv_n6_0611.txt");
    computeLaplacian(7, "lap_conv_n7_0611.txt");
    computeLaplacian(8, "lap_conv_n8_0611.txt");
    computeLaplacian(9, "lap_conv_n9_0611.txt");
    computeLaplacian(10, "lap_conv_n10_0611.txt");
    computeLaplacian(20, "lap_conv_n20_0611.txt");
    computeLaplacian(50, "lap_conv_n50_0611.txt");


    return 0;
}

// conjugate gradient method implementation
MVector conjugateGradientSolve(const MMatrix& A, const MVector& b, const MVector& x_0, const string filename = "iteration.txt")
{
    // initialisations
    int maxIterations = 1000;
    double tolerance = 1e-6;

    if (A.Rows() != A.Cols())
    {
        throw invalid_argument("error: matrix is not square");
    }
    int dim = A.Rows();
    
    MVector rPrev(dim), pPrev(dim), xPrev(dim), r(dim), p(dim), x(dim);
    double alpha, beta;

    // solver writes a file, an optional filename argument is given
    ofstream of;
    int w = 12;
    of.open(filename);
    if (!of)
    {
        throw invalid_argument("error: failed to open file");
    }
    of << setw(w) << "Time t" << " | "
    << setw(w) <<  "L2 norm" << " | "
    << setw(w) << "LInf norm" << "\n";

    xPrev = x_0;
    rPrev = b - A*xPrev;
    pPrev = rPrev;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    of << setw(w) << 0 << " | "
    << setw(w) <<  rPrev.L2Norm() << " | "
    << setw(w) << rPrev.LInfNorm() << "\n";

    for (int iter=0; iter<maxIterations; iter++)
    {
        // ...calculate new values for x and r here...
        alpha = dot(rPrev, rPrev)/dot(pPrev, A*pPrev);
        x = xPrev + alpha*pPrev;
        r = rPrev - alpha*(A*pPrev);

        // write the iteration
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

        of << setw(10) << time_span.count() << " | "
        << setw(10) <<  r.L2Norm() << " | "
        << setw(10) << r.LInfNorm() << "\n";
        
        // check if solution is accurate enough
        if (r.L2Norm() < tolerance) break;
        
        // ...calculate new conjugate vector p here...
        beta = dot(r, r)/dot(rPrev, rPrev);
        p = r + beta*pPrev;

        // step the placeholder variables for the next iteration
        pPrev = p;
        rPrev = r;
        xPrev = x;
    }

    of.close();
    return x;
}

// procedure to compute a solution using the conjugate gradient method, given the value n and the file to write to
void computeLaplacian(int n, const string& filename)
{   
    MMatrix H(n*n, n*n);
    makeLaplacian(H, 4, -1);

    MVector a(n*n, 1.0/pow(n+1, 2));
    // initial guess is always the unit vector
    MVector x(n*n);
    x[0] = 1;

    conjugateGradientSolve(H, a, x, filename);
}