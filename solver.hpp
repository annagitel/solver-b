#pragma once
#include <complex>
using namespace std;
using std::complex;

namespace solver
{
	class RealVariable{
	public:
		double _a,_b,_c;
		RealVariable(){
            _a=0;
            _b=1;
            _c=0;
        }

        RealVariable(double a,double b,double c){
            _a=a;
            _b=b;
            _c=c;
        }//RealVariable
		friend const RealVariable operator+ (const RealVariable& c1, const RealVariable& c2);
		friend const RealVariable operator+ (const RealVariable& c1, const double& c2);
		friend const RealVariable operator+ (const double& c1, const RealVariable& c2);
		friend const RealVariable operator- (const RealVariable& c1, const RealVariable& c2);
		friend const RealVariable operator- (const RealVariable& c1, const double& c2);
		friend const RealVariable operator- (const double& c1, const RealVariable& c2);
		friend const RealVariable operator- (const RealVariable& c1);
		friend const RealVariable operator* (const RealVariable& c1, const RealVariable& c2);
		friend const RealVariable operator* (const RealVariable& c1, const double& c2);
		friend const RealVariable operator* (const double& c1, const RealVariable& c2);
		friend const RealVariable operator/ (const RealVariable& c1, const double& c2);
		friend const RealVariable operator^(const RealVariable& c1, const double& c2);
		friend const RealVariable operator==(const RealVariable& c1, const RealVariable& c2);
		friend const RealVariable operator==(const RealVariable& c1, const double& c2);
		friend const RealVariable operator==(const double& c1, const RealVariable& c2);
		friend ostream& operator<< (ostream& os, const RealVariable& c);

	};

	class ComplexVariable
	{
	public:
		complex<double> _a,_b,_comp;
        ComplexVariable(){
            _a=complex<double>(0,0);
            _b=complex<double>(1,0);
            _comp=complex<double>(0,0);
        }
         ComplexVariable(complex<double> a,complex<double> b,complex<double> comp){
            _a=a;
            _b=b;
            _comp=comp;
        }
		friend const ComplexVariable operator+ (const ComplexVariable& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator+ (const ComplexVariable& c1, const complex<double>& c2);
		friend const ComplexVariable operator+ (const complex<double>& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator+ (const ComplexVariable& c1, const double& c2);
		friend const ComplexVariable operator+ (const double& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator- (const ComplexVariable& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator- (const ComplexVariable& c1, const complex<double>& c2);
		friend const ComplexVariable operator- (const complex<double>& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator- (const ComplexVariable& c1, const double& c2);
		friend const ComplexVariable operator- (const double& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator- (const ComplexVariable& c1);
		friend const ComplexVariable operator* (const ComplexVariable& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator* (const ComplexVariable& c1, const complex<double>& c2);
		friend const ComplexVariable operator* (const complex<double>& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator* (const ComplexVariable& c1, const double& c2);
		friend const ComplexVariable operator* (const double& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator/ (const ComplexVariable& c1, const double& c2);
		friend const ComplexVariable operator/ (const double& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator^ (const ComplexVariable& c1, const double& c2);
		friend const ComplexVariable operator== (const ComplexVariable& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator== (const ComplexVariable& c1, const complex<double>& c2);
		friend const ComplexVariable operator== (const complex<double>& c1, const ComplexVariable& c2);
		friend const ComplexVariable operator== (const ComplexVariable& c1, const double& c2);
		friend const ComplexVariable operator== (const double& c1, const ComplexVariable& c2);

		friend ostream& operator<< (ostream& os, const ComplexVariable& c);
	};

	double solve(RealVariable equation);

	std::complex<double> solve(ComplexVariable equation);
	string toString(complex<double> comp);
}

namespace complexExtend{
	const complex<double> operator+(double r,const complex<double>& c2);
	const complex<double> operator+(const complex<double>&c1,double r);
	const complex<double> operator-(double r,const complex<double>& c2);
	const complex<double> operator-(const complex<double>&c1,double r);
	const complex<double> operator*(double r,const complex<double>& c2);
	const complex<double> operator*(const complex<double>&c1,double r);
	const complex<double> operator/(const complex<double>&c1,double r);
}