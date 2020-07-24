#include <complex>
#include "solver.hpp"
#include <math.h>
using std::complex;
namespace solver {
	//plus
	const RealVariable operator+(const RealVariable &c1, const RealVariable &c2) {
		return RealVariable(c1._a+c2._a,c1._b+c2._b,c1._c+c2._c);
	}

	const RealVariable operator+(const RealVariable &c1, const double &c2) {
		return RealVariable(c1._a,c1._b,c1._c+c2);
	}

	const RealVariable operator+(const double &c1, const RealVariable &c2) {
		return RealVariable(c2._a,c2._b,c2._c+c1);
	}
	//minus
	const RealVariable operator-(const RealVariable &c1, const RealVariable &c2) {
		return RealVariable(c1._a-c2._a,c1._b-c2._b,c1._c-c2._c);
	}

	const RealVariable operator-(const RealVariable &c1, const double &c2) {
		return RealVariable(c1._a,c1._b,c1._c-c2);
	}

	const RealVariable operator-(const double &c1, const RealVariable &c2) {
		return RealVariable(-c2._a,-c2._b,c1-c2._c);
	}

	const RealVariable operator-(const RealVariable &c1) {
		return RealVariable(-c1._a,-c1._b,-c1._c);
	}
	//multi
	const RealVariable operator*(const RealVariable &c1, const RealVariable &c2) {
		if((c1._a!=0 && c2._a!=0) || (c1._a!=0 && c1._b!=0) || (c2._a!=0 && c1._b!=0))
                throw std::out_of_range{"degree dos not match\n"};
            return RealVariable(c1._c*c2._a+c1._b*c2._b+c2._a*c1._c,c1._b*c2._c+c2._b*c1._c,c1._c*c2._c);
	}

	const RealVariable operator*(const RealVariable &c1, const double &c2) {
		return RealVariable(c2*c1._a,c2*c1._b,c2*c1._c);
	}

	const RealVariable operator*(const double &c1, const RealVariable &c2) {
		return RealVariable(c1*c2._a,c1*c2._b,c1*c2._c);
	}
	//davide
	const RealVariable operator/(const RealVariable &c1, const double &c2) {
		if(c2==0)
			throw string("division by zero\n");
		return RealVariable(c1._a/c2,c1._b/c2,c1._c/c2);
	}
	//power
	const RealVariable operator^(const RealVariable &c1, const double &c2) {
		if(c2==0)
			return RealVariable(0,0,1);
		else if(c1._a==0 && c1._c!=0 && c2==2)
			return RealVariable(pow(c1._a,2),2*c1._a*c1._b,pow(c1._c,2));
		else if(c1._a==0 && c1._c==0 && c2==2)
			return RealVariable(pow(c1._b,2),0,0);
		else
			throw std::out_of_range{"degree does not match\n"};
	}
	//equals
	const RealVariable operator==(const RealVariable &c1, const RealVariable &c2) {
		return RealVariable(c1._a-c2._a,c1._b-c2._b,c1._c-c2._c);
	}

	const RealVariable operator==(const RealVariable &c1, const double &c2) {
		return RealVariable(c1._a,c1._b,c1._c-c2);
	}

	const RealVariable operator==(const double &c1, const RealVariable &c2) {
		return RealVariable(c2._a,c2._b,c2._c-c1);
	}
	//stream
	ostream &operator<<(ostream &os, const RealVariable &c) {
		if(c._a!=0  && c._b!=0 && c._c!=0){
                if(c._b>0 && c._c>0)
                    return (os <<c._a<<"x^2"<<'+'<<c._b<<'x'<<'+'<<c._c);
                else if(c._b>0 && c._c<0)
                    return (os <<c._a<<"x^2"<<'+'<<c._b<<'x'<<c._c);
                else if(c._b<0 && c._c<0)
                    return (os <<c._a<<"x^2"<<c._b<<'x'<<c._c);
                else
                    return (os <<c._a<<"x^2"<<c._b<<'x'<<'+'<<c._c);
            }
            else if(c._a!=0  && c._b!=0 && c._c==0){
                if(c._b>0)
                    return (os<<c._a<<"x^2"<<'+'<<c._b<<'x');
                else
                    return (os<<c._a<<"x^2"<<c._b<<'x');
            }
            else if(c._a!=0  && c._b==0 && c._c!=0){
                if(c._c>0)
                    return (os<<c._a<<"x^2"<<'+'<<c._c);
                else
                    return (os<<c._a<<"x^2"<<c._c);
            }
            else if(c._a!=0  && c._b==0 && c._c==0)
            	return (os<<c._a<<"x^2");

            else if(c._a==0  && c._b!=0 && c._c!=0){
                if(c._c>0)
                    return (os<<c._b<<'x'<<'+'<<c._c);
                else
                    return (os<<c._b<<'x'<<c._c);
            }
            else if(c._a==0  && c._b==0 && c._c!=0)
                    return (os<<c._c);
            else if(c._a==0  && c._b!=0 && c._c==0)
                    return (os<<c._b<<'x');
            else
                throw string("ERR: 0x^2+0x+0\n");
	}

	//plus
	const ComplexVariable operator+(const ComplexVariable &c1, const ComplexVariable &c2) {
		return ComplexVariable(c1._a+c2._a,c1._b+c2._b,c1._comp+c2._comp);
	}

	const ComplexVariable operator+(const ComplexVariable &c1, const complex<double> &c2) {
		return ComplexVariable(c1._a,c1._b,c2+c1._comp);
	}

	const ComplexVariable operator+(const complex<double> &c1, const ComplexVariable &c2) {
		return ComplexVariable(c2._a,c2._b,c2._comp+c1);
	}

	const ComplexVariable operator+(const ComplexVariable &c1, const double &c2) {
		return ComplexVariable(c1._a,c1._b,c1._comp+c2);
	}

	const ComplexVariable operator+(const double &c1, const ComplexVariable &c2) {
		return ComplexVariable(c2._a,c2._b,c2._comp+c1);
	}
	//minus
	const ComplexVariable operator-(const ComplexVariable &c1, const ComplexVariable &c2) {
		return ComplexVariable(c1._a-c2._a,c1._b-c2._b,c1._comp-c2._comp);
	}

	const ComplexVariable operator-(const ComplexVariable &c1, const complex<double> &c2) {
		return ComplexVariable(c1._a,c1._b,c1._comp-c2);
	}

	const ComplexVariable operator-(const complex<double> &c1, const ComplexVariable &c2) {
		return ComplexVariable(-c2._a,-c2._b,c1-c2._comp);
	}

	const ComplexVariable operator-(const ComplexVariable &c1, const double &c2) {
		return ComplexVariable(c1._a,c1._b,c1._comp-c2);
	}

	const ComplexVariable operator-(const double &c1, const ComplexVariable &c2) {
		return ComplexVariable(-c2._a,-c2._b,c1-c2._comp);
	}

	const ComplexVariable operator-(const ComplexVariable &c1) {
		return ComplexVariable(-c1._a,-c1._b,-c1._comp);
	}
	//multi
	const ComplexVariable operator*(const ComplexVariable &c1, const ComplexVariable &c2) {
		if((c1._a!=complex<double>(0,0) && c2._a!=complex<double>(0,0)) || (c1._a!=complex<double>(0,0) && c1._b!=complex<double>(0,0)) || (c2._a!=complex<double>(0,0) && c1._b!=complex<double>(0,0)))
			throw std::out_of_range{"degree do not match\n"};
		return ComplexVariable(c1._b*c2._b,c1._b*c2._comp+c1._comp*c2._b,c1._comp*c2._comp);
	}

	const ComplexVariable operator*(const ComplexVariable &c1, const complex<double> &c2) {
		return ComplexVariable(c1._a*c2,c1._b*c2,c1._comp*c2);
	}

	const ComplexVariable operator*(const complex<double> &c1, const ComplexVariable &c2) {
		return ComplexVariable(c1*c2._a,c1*c2._b,c1*c2._comp);
	}

	const ComplexVariable operator*(const ComplexVariable &c1, const double &c2) {
		return ComplexVariable(c1._a*c2,c1._b*c2,c1._comp*c2);
	}

	const ComplexVariable operator*(const double &c1, const ComplexVariable &c2) {
		return ComplexVariable(c1*c2._a,c1*c2._b,c1*c2._comp);
	}
	//davide
	const ComplexVariable operator/(const ComplexVariable &c1, const double &c2) {
		if(c2==complex<double>(0,0))
			throw string("division by 0\n");
		return ComplexVariable(c1._a/c2,c1._b/c2,c1._comp/c2);
	}

	//const ComplexVariable operator/(const double &c1, const ComplexVariable &c2) { return ComplexVariable(); }
	//power
	const ComplexVariable operator^(const ComplexVariable &c1, const double &c2) {
		if(c2==complex<double>(0,0))
			return ComplexVariable(complex<double>(0,0),complex<double>(0,0),complex<double>(1,0));
		else if(c1._a==complex<double>(0,0) && c1._comp!=complex<double>(0,0) && c2==2)
			return ComplexVariable(pow(c1._a,2),(2.0*c1._a)*c1._b,pow(c1._comp,2));
		else if(c1._a==complex<double>(0,0) && c1._comp==complex<double>(0,0) && c2==2)
			return ComplexVariable(pow(c1._b,2),complex<double>(0,0),complex<double>(0,0));
		else
			throw std::out_of_range{"degree do not match\n"};
	}
	//equals
	const ComplexVariable operator==(const ComplexVariable &c1, const ComplexVariable &c2) {
		return ComplexVariable(c1._a-c2._a,c1._b-c2._b,c1._comp-c2._comp);
	}

	const ComplexVariable operator==(const ComplexVariable &c1, const complex<double> &c2) {
		return ComplexVariable(c1._a,c1._b,c1._comp-c2);
	}

	const ComplexVariable operator==(const complex<double> &c1, const ComplexVariable &c2) {
		return ComplexVariable(c2._a,c2._b,c2._comp-c1);
	}

	const ComplexVariable operator==(const ComplexVariable &c1, const double &c2) {
		return ComplexVariable(c1._a,c1._b,c1._comp-c2);
	}

	const ComplexVariable operator==(const double &c1, const ComplexVariable &c2) {
		return ComplexVariable(c2._a,c2._b,c2._comp-c1);
	}
	//stream
	ostream &operator<<(ostream &os, const ComplexVariable &c) {
		if(c._a!=complex<double>(0,0)  && c._b!=complex<double>(0,0) && c._comp!=complex<double>(0,0)){
                    if(c._b.real()>0 && c._comp.real()>0)
                        return (os <<toString(c._a)<<"z^2"<<toString(c._b)<<'z'<<toString(c._comp));
                    else if(c._b.real()>0 && c._comp.real()<0)
                        return (os <<toString(c._a)<<"z^2"<<toString(c._b)<<'z'<<toString(c._comp));
                    else if(c._b.real()<0 && c._comp.real()<0)
                        return (os <<toString(c._a)<<"z^2"<<toString(c._b)<<'z'<<toString(c._comp));
                    else    //(c._b.real()<0 && c._comp.real()>0)
                        return (os <<toString(c._a)<<"z^2"<<toString(c._b)<<'z'<<toString(c._comp));
                 }
                else if(c._a!=complex<double>(0,0)  && c._b!=complex<double>(0,0) && c._comp==complex<double>(0,0)){
                    if(c._b.real()>0)
                        return (os<<toString(c._a)<<"z^2"<<toString(c._b)<<'z');
                    else
                        return (os<<toString(c._a)<<"z^2"<<toString(c._b)<<'z');
                }
                else if(c._a!=complex<double>(0,0)  && c._b==complex<double>(0,0) && c._comp!=complex<double>(0,0)){
                    if(c._comp.real()>0)
                        return (os<<c._a<<"z^2"<<toString(c._comp));
                    else
                        return (os<<toString(c._a)<<"z^2"<<toString(c._comp));
                }
                else if(c._a!=complex<double>(0,0)  && c._b==complex<double>(0,0) && c._comp==complex<double>(0,0))
                    return (os<<toString(c._a)<<"z^2");

                else if(c._a==complex<double>(0,0)  && c._b!=complex<double>(0,0) && c._comp!=complex<double>(0,0)){
                    if(c._comp.real()>0)
                        return (os<<toString(c._b)<<'z'<<toString(c._comp));
                    else
                        return (os<<toString(c._b)<<'z'<<toString(c._comp));
                }
                else if(c._a==complex<double>(0,0)  && c._b==complex<double>(0,0) && c._comp!=complex<double>(0,0))
                        return (os<<toString(c._comp));
                else if(c._a==complex<double>(0,0)  && c._b!=complex<double>(0,0) && c._comp==complex<double>(0,0))
                        return (os<<toString(c._b)<<'z');
                else
                    throw string("ERR: 0z^2+0z+0");
	}

	double solve(RealVariable equation) {
		if(equation._a!=0 && equation._b!=0 && equation._c!=0 ){
			printf("Equation form: ax^2+bx+c=0\t");
			double discriminant =pow(equation._b,2)-4*equation._a*equation._c;
			if(discriminant <0)
				throw string("ERR:No solution in the realnumbers\n");
			double solutionA=((-1*equation._b)+sqrt(discriminant ))/2.0;
			double solutionB=((-1*equation._b)-sqrt(discriminant ))/2.0;
			return solutionA;
		}
		else if(equation._a!=0 && equation._b==0 && equation._c!=0){
			printf("Equation form: ax^2+c=0\t");
			double preSol=(-1*equation._c)/equation._a;
			if(preSol<0)
				throw std::out_of_range{"ERR:No solution in realnumbers\n"};

			double solutionA=sqrt(preSol);
			double solutionB=-1*sqrt(preSol);
			return solutionA;
		}
		else if(equation._a==0 && equation._b!=0 && equation._c!=0){
			printf("Equation form: bx+c=0\t");
			double solution=(-1*equation._c)/equation._b;
			return solution;
		}
		else if(equation._a==0 && equation._b!=0 && equation._c==0){
			printf("Equation form: bx=0\t");
			return 0;
		}
		else if(equation._a==0 && equation._b==0 && equation._c!=0){
			printf("Equation form: c=0\t");
			if(equation._c==0)
				throw string("TRUE for all realvariable x");
			else
				throw string("NO SOLUTION!");
		}
		else if(equation._a!=0 && equation._b!=0 && equation._c==0){
			double solutionA=(-1*equation._b)/equation._a;
			double solutionB=0;
			return solutionA;
		}
		else
			throw string("ERR: unexpected error for the equation ");
	}

	std::complex<double> solve(ComplexVariable equation) {
		if(equation._a==complex<double>(0,0) && equation._b==complex<double>(0,0) && equation._comp==complex<double>(0,0))
			throw string("ERR:0=0 True statment for all complex variable.");
		else if(equation._a==complex<double>(0,0) && equation._b==complex<double>(0,0) && equation._comp!=complex<double>(0,0))
			throw string("ERR:No solution.");
		else if(equation._a==complex<double>(0,0) && equation._b!=complex<double>(0,0) && equation._comp==complex<double>(0,0))
			return complex<double>(0,0);
		else if(equation._a==complex<double>(0,0) && equation._b!=complex<double>(0,0) && equation._comp!=complex<double>(0,0)){
			complex<double> solution=(-equation._comp)/equation._b;
			return solution;
		}
		else if(equation._a!=complex<double>(0,0) && equation._b==complex<double>(0,0) && equation._comp==complex<double>(0,0))
			return complex<double>(0,0);
		else if(equation._a!=complex<double>(0,0) && equation._b==complex<double>(0,0) && equation._comp!=complex<double>(0,0)){
			complex<double> solutionA= sqrt(-equation._comp/equation._a);
			complex<double> solutionB= -sqrt(-equation._comp/equation._a);
			return solutionA;
		}
		else if(equation._a!=complex<double>(0,0) && equation._b!=complex<double>(0,0) && equation._comp==complex<double>(0,0))
			return complex<double>(0,0);

		else if(equation._a!=complex<double>(0,0) && equation._b!=complex<double>(0,0) && equation._comp!=complex<double>(0,0))
			return complex<double>(0,0);
		return complex<double>(0,0);
	}
	string toString(complex<double> comp){
            if(comp.imag()>0){
                if(comp.real()>0)
                    return "+("+to_string(comp.real())+"+"+to_string(comp.imag())+"i)";
                else if(comp.real()==0)
                    return "+"+to_string(comp.imag())+"i";
                else
                    return "-("+to_string(comp.real())+"+"+to_string(comp.imag())+"i)";
            }
            else if(comp.imag()==0){
                if(comp.real()>0)
                    return to_string(comp.real());
                else if(comp.real()==0)
                     return "";
                else
                    return to_string(comp.real());
            }
            else{
                if(comp.real()>0)
                    return "+("+to_string(comp.real())+to_string(comp.imag())+"i)";
                else if(comp.real()==0)
                    return to_string(comp.imag())+"i";
                return "-("+to_string(comp.real())+to_string(comp.imag())+"i)";
            }
        }
}

namespace complexExtend {
	const complex<double> operator+(double r,const complex<double>& c2){
		return complex<double>(r+c2.real(),c2.imag());
	}

	const complex<double> operator+(const complex<double>&c1,double r){
		return complex<double>(r+c1.real(),c1.imag());
	}

	const complex<double> operator-(double r,const complex<double>& c2){
		return complex<double>(r-c2.real(),-c2.imag());
	}
	const complex<double> operator-(const complex<double>&c1,double r){
		return complex<double>(c1.real()-r,c1.imag());
	}

	const complex<double> operator*(double r,const complex<double>& c2){
		return complex<double>(r*c2.real(),r*c2.imag());
	}

	const complex<double> operator*(const complex<double>&c1,double r){
		return complex<double>(c1.real()*r,c1.imag()*r);
	}

	const complex<double> operator/(const complex<double>&c1,double r){
		return complex<double>(c1.real()/r,c1.imag()/r);
	}
}