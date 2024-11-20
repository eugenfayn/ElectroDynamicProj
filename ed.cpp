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
		D_i = 2 * calc_S(a,b,c);
	}
	else if (cond2){
		D_i = 2 * calc_S(a,b,c);
	}
	else {
		D_i = 0;
	}
}

// func for calc A_ij
void Aij_func(){
	cout << endl;
}

int main(){
	cout << "Geometry" << endl;
	double d = 0.0;
	calc_Di(d,3);
	return 0;
}
