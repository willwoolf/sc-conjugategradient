#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include "mvector.h"

using namespace std;

// constructors in header

MVector operator*(const double& lhs, const MVector& rhs)
{
    MVector temp = rhs;
    for (int i = 0; i < temp.size(); i++)
    {
        temp[i] *= lhs;
    }
    return temp;
}

MVector operator*(const MVector& lhs, const double& rhs)
{
    return rhs*lhs;
}

MVector operator+(const MVector& lhs, const MVector& rhs)
{
    if (lhs.size() != rhs.size())
    {
        exception e;
        throw e;
    }
    MVector temp(lhs.size());
    for (int i = 0; i < temp.size(); i++)
    {
        temp[i] = lhs[i] + rhs[i];
    }
    return temp;
}

MVector operator-(const MVector& lhs, const MVector& rhs)
{
    return lhs+(-1*rhs);
}

// alternative overload such that given MVector v, -v is allowed
MVector MVector::operator-()
{
    MVector vec(v.size());
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] = -v[i];
    }
    return vec;
}

MVector operator/(const MVector& lhs, const double& rhs)
{
    return (1/rhs)*lhs;
}

ostream& operator<<(ostream& os, const MVector& v)
{
	int n = v.size();
	os << "(";
	for (int i = 0; i < n-1; i++)
    {
        os << setw(10) << v[i] << ", ";
    }
    os << v[n-1];
	os << ")";
	return os;
}

double MVector::LInfNorm() const
{
    double maxAbs = 0;
    size_t s = size();
    for (int i=0; i<s; i++)
    {
        maxAbs = max(abs(v[i]), maxAbs);
    }
    return maxAbs;
}

double MVector::L2Norm() const
{
    double sum = 0;
    for (int i = 0; i < v.size(); i++)
    {
        sum += v[i]*v[i];
    }
    return sqrt(sum);
}

double dot(const MVector& lhs, const MVector& rhs)
{
    if (lhs.size() != rhs.size())
    {
        throw invalid_argument("error: dot product not defined for vectors of non-equal dimension");
    }
    double sum = 0;
    for (int i = 0; i < lhs.size(); i++)
    {
        sum += lhs[i]*rhs[i];
    }
    return sum;
}