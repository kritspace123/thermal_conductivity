#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double u(double x,  double t){
	return pow(x+0.1, 2) - (sin(2*M_PI*t))/2+x - 3.5*t;
}

double f(double x, double t){
	return -M_PI*cos(2*M_PI*t) - 3.164+0.56*x;
}
double fi(double x){
	return u(x, 0);//pow(x+0.1, 2) + x
}

double g1(double t){
	return u(0, t); // - (sin(2*M_PI*t))/2 -3.5t
}
double Courant_Сondition(double tao){
	return 0;
}

int main(){
	int N, J;
	ofstream file("value.txt");
	ofstream fout("value_function.txt");


	double tao = 0.01;
	double h;

	cout << "Введите h>0 удовлетворяющее условию кронекера:"<< endl;
	cin >> h;


	double a = 0.032;

	N = 1 / h;
	J = 1 / tao;


	double A = (a*tao) / h;
	double B  = 1 - A;

	double* array = new double[N+1]; 

	double old, now1, fun, s;
	double maxim = 0;
	
	for(int j = 0; j <= N; j++){
		if(j == 0){
			for(int k = 0; k <= N; k++){ 
				array[k] = fi(k*h);
				file << array[k] << " ";
				fout << u(h*k, j*tao)<< " ";
			}
			file << endl;
			fout << endl;
			continue;
		}
		for(int k = 0; k <= N; k++){
			if(k == 0){
				old = array[k];
				array[k] = g1(j*tao);
				file << array[k] << " ";
				continue;
			}
			now1 = B * array[k] + A * old + f(h*k, tao*j);
			old = array[k];
			array[k] = now1;
			fun = u(h*k, j*tao);

			fout << fun << " ";
			file << now1 << " ";

			s = abs(fun - now1);
			if (s > maxim){
				maxim = s;
			}
		}
		file << endl;
		fout << endl;
	}
}