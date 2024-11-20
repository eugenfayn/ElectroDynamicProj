#include <iostream>
#include <complex>
#include <string>

using namespace std;

// Read geometry
// Calculating Aij, Bi
// Решение системы
// Постобработка

struct Vertex {
    double x, y, z;
};

struct Face {
    int v1, v2, v3;
};

extern "C"
{
	extern void zgesv_(int *N, int *Nrhs, complex<double> *A, int *ldA ,     int *Ipvt, complex<double> *B, int *ldB, int *info);
};

double scalar_product(Vertex &a, Vertex &b){
	double sum_ = a.x * b.x + a.y * b.y + a.z * b.z;
	return sum_
}

// Calculation S by 3 points
double calc_S(Vertex a,Vertex b,Vertex c){
	double S = 0.0;
	S = S + (b.x - a.x) * (c.y - a.y);
	S = S - (c.x - a.x) * (b.y - a.y);
	S = 0.5 * S;
	return S;
}	

// func for D_i
void calc_Di(double &D_i, int i, Vertex &a, Vertex &b, Vertex &c, Vertex &d){
	bool cond1 = false;
	bool cond2 = false;
	//|-----a-------|
	//c     |	d
	//|---- b-------|
	if (cond1){
		D_i = - 2 * calc_S(a,b,d);
	}
	else if (cond2){
		D_i = 2 * calc_S(a,b,c);
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
		e_i.x  = (d.x - X.x) / calc_S(a,b,d);
		e_i.y  = (d.y - X.y) / calc_S(a,b,d);
		e_i.z  = (d.z - X.z) / calc_S(a,b,d);
	}
	else if (cond2){
		e_i.x = (X.x - c.x) / calc_S(a,b,c);
		e_i = (X.y - c.y) / calc_S(a,b,c);
		e_i = (X.z - c.z) / calc_S(a,b,c);
	}
	else{
		e_i.x = 0;
		e_i.y = 0;
		e_i.z = 0;
	}	
}

// func for calc A_ij
void calc_Aij(){
	double integral1 = k * k * scalar_product(ei,ej) - Di, Dj; // (k^2 * (ei,ej) - Di * Dj) F(x - y)dydt
	double integral2 = 0;
	double integral3 = 0;
	double integral4 = 0;
}

int main(){
	cout << "Geometry" << endl;
	double d = 0.0;
	calc_Di(d,3);
	return 0;
}
