#include <vector>
#include <iostream>
#include <iomanip>
#include <ios>
#include <algorithm>
#include <cmath>

#include "mmatrix.h"
#include "mvector.h"

using namespace std;

MVector operator*(const MMatrix& A, const MVector& v)
{
    // check compatibility
    if (A.Cols() != v.size())
    {
        cout << "wrong dims" << endl;
        exception e;
        throw e;
    }
    // initialise output vector
    int vout_size = A.Rows();
    MVector vout = MVector(vout_size);

    for (int i = 0; i < A.Rows(); i++)
    {
        // sum across j
        double sum = 0;
        for (int j = 0; j < A.Cols(); j++)
        {
            sum += A(i, j)*v[j];
        }
        vout[i] = sum;        
    }
    return vout;
}

ostream& operator<<(ostream& os, const MMatrix& A)
{
    for (int i = 0; i < A.Rows(); i++)
    {
        os << "|  ";
        for (int j = 0; j < A.Cols(); j++)
        {
            os << setw(3) << A(i, j) << "  ";
        }
        os << "|\n";
    }
    return os;
}

void makeTridiagonal(MMatrix& A, double d, double bd)
{
    /*
    takes the arguments for the diagonal and band, and passes the matrix by reference.
    A must already be the size we want, but its entries will be entirely overwritten.
    */
    if (A.Rows() != A.Cols())
    {
        throw invalid_argument("error: matrix not square");
    }
    for (int i = 0; i < A.Rows(); i++)
    {
        for (int j = 0; j < A.Cols(); j++)
        {
            if (i == j)
            {
                A(i, j) = d;
            }
            else if (abs(i - j) == 1)
            {
                A(i, j) = bd;
            }
            else
            {
                A(i, j) = 0;
            }
        }
    }  
}

void makeLaplacian(MMatrix& A, double d, double bd)
{
    /*
    behaves similarly to makeTridiagonal in its arguments. for solutions to the Laplacian equation.
    */
    if (A.Rows() != A.Cols())
    {
        throw invalid_argument("error: matrix not square");
    }
    double sz = sqrt(A.Rows());
    if (floor(sz) != sz)
    {
        // A must be n^2 by n^2
        throw invalid_argument("error: matrix is not of dimension n^2.n^2");
    }
    int n = sz;
    for (int i = 0; i < A.Rows(); i++)
    {
        for (int j = 0; j < A.Cols(); j++)
        {
            if (i == j)
            {
                A(i, j) = d;
            }
            else if (abs(i - j) == n)
            {
                A(i, j) = bd;
            }
            else if ((abs(i - j) == 1) && ((i + j) % 2*n != 2*n - 1))
            {
                A(i, j) = bd;
            }
            else
            {
                A(i, j) = 0;
            }
        }
    }
}
