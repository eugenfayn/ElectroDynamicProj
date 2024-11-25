#include "geometry/geometry.h"

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