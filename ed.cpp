#include <iostream>
#include <complex>
#include <string>

using namespace std;

// Read geometry
// Calculating Aij, Bi
// Решение системы
// Постобработка

class Vertex {
private:
	double x,y,z;
public:
	Vertex(){
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
	Vertex(double a, double b, double c){
		x = a;
		y = b;
		z = c;
	}
	Vertex operator+(const Vertex& right) const{
       		return Vertex(x+right.x,y+right.y,z+right.z);
	}
	Vertex operator-(const Vertex& right) const{
		return Vertex(x-right.x,y-right.y,z-right.z);
	}
	Vertex operator/(double scalar) const{
		return Vertex(x/scalar,y/scalar,z/scalar);
	}
	double scalart_product(const Vertex& right){
		double sum_ = x * right.x + y * right.y + z * right.z;
 		return sum_;
	}
	double calc_S(const Vertex& b, const Vertex& c){
		double S = 0.0;
 		S = S + (b.x - x) * (c.y - y);
 		S = S - (c.x - x) * (b.y - y);
 		S = 0.5 * S;
 		return S;
	}
};

struct Face {
    int v1, v2, v3;
};

extern "C"
{
	extern void zgesv_(int *N, int *Nrhs, complex<double> *A, int *ldA ,     int *Ipvt, complex<double> *B, int *ldB, int *info);
};

// func for D_i
void calc_Di(double &D_i, int i, Vertex &a, Vertex &b, Vertex &c, Vertex &d){
	bool cond1 = false;
	bool cond2 = false;
	//|-----a-------|
	//c     |	d
	//|---- b-------|
	if (cond1){
		D_i = - 2 * a.calc_S(b,d);
	}
	else if (cond2){
		D_i = 2 * a.calc_S(b,c);
	}
	else {
		D_i = 0;
	}
}

// func for calc e_i
void calc_ei(Vertex  &e_i, int i, Vertex &X, Vertex &a, Vertex &b, Vertex &c, Vertex &d){
	bool cond1 = false;
	bool cond2 = false;
	if (cond1){
		e_i = d - X;
		e_i = e_i/a.calc_S(b,d);
	}
	else if (cond2){
		e_i = X - c;
		e_i = e_i/a.calc_S(b,c);
	}
	else{
		e_i = Vertex();
	}	
}

// func for calc A_ij
void calc_Aij(){
	//double integral1 = k * k * scalar_product(ei,ej) - Di, Dj; // (k^2 * (ei,ej) - Di * Dj) F(x - y)dydt
	double integral2 = 0;
	double integral3 = 0;
	double integral4 = 0;
}

int main(){
	cout << "Geometry" << endl;
	Vertex a(0.1,0.1,0.1);
	Vertex b(0.2,0.3,0.5);
	Vertex c(0.0,0.7,0.9);
	double SSS = 0.0;
	SSS = a.calc_S(b,c);
	cout << "SSS: " << SSS << endl;
	return 0;
}
