#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

void PrintMassiv(double* array, int N){
	for(int i = 0; i <= N; i++){
		cout << array[i] << " ";
	}
	cout << endl;
}

double u(double x,  double t){
	return -0.5*pow(x, 4) + pow(x, 2) - x + t*x + 2 * pow(t,2) - t * exp(x);
}

double du_dt(double x, double t){
	return -6*pow(x, 2) + 2 - t * exp(x);
}

double mu(double x){
	return -0.5*pow(x, 4) + pow(x, 2) - x;
}

double mu1(double t){
	return 2*pow(t, 2) - t;
}

double mu2(double t){
	return -0.5 + t + 2*pow(t, 2) - t * exp(1);
}

double f(double x, double t){
	return x + 4*t-exp(x)+0.192*pow(x, 2) - 0.064 + 0.032*t*exp(x);
}

void sweep_method(int N, double A, double B, double* vector_B, double* vector_X){

	//straight stroke
	double* vector_V = new double[N];
	double* vector_U = new double[N];

	vector_V[1] = - (B) / (A);
	vector_U[1] = vector_B[1] / A;

	for (int i = 2; i < N-1; i++){
		vector_V[i] = (B) / (-A - (B * vector_V[i-1]));
		vector_U[i] = (B*vector_U[i-1] - vector_B[i]) / (-A - B* vector_V[i-1]);
	}
	
	vector_V[N-1] = 0;
	vector_U[N-1] = (B*vector_U[N-2] - vector_B[N-1]) / (- A - B* vector_V[N-2]);

	//reverse stroke

	vector_X[N-1] = vector_U[N-1];
	for(int i = N-2; i>= 1; i--){
		vector_X[i] = vector_V[i] * vector_X[i+1] + vector_U[i];
	}
}

int main(){
	int N, J;
	ofstream file("value.txt");
	ofstream fout("value_function.txt");


	double tao = 0.01;
	double h = 0.01;
	double a = 0.032;

	N = 1 / h;
	J = 1 / tao;


	double A = (J + 2* a*pow(N, 2));
	double B = (-a * pow(N, 2));

	
	double* vector_V = new double[N];
	double* vector_U = new double[N];
	double* vector_b = new double[N+1];
	double maxim = 0;

	for(int j = 0; j <= N; j++){
		if(j == 0){
			for(int k = 0; k <= N; k++){ 
				array[k] = mu(k*h);
				file << array[k] << " ";
				fout << u(h*k, j*tao)<< " ";
			}
			file << endl;
			fout << endl;
			continue;
		}	
	array[0] = mu1(j*tao);
	array[N] = mu2(j*tao);

	//вычисляем вектор b: Ax=b
	vector_b[1] = f(h, tao*j) + J*array[1] + a*pow(N, 2)*array[0];
	for(int k = 2; k < N-1; k++){
		vector_b[k] = f(h*k,tao*j)+ J*array[k];
	}
	vector_b[N-1] = f(h*N-1, tao * j)+J*array[N-1] + a*pow(N, 2)*array[N];
	sweep_method(N, A, B, vector_b, array);



	// //straight stroke
	// vector_V[1] = - B / A;
	// vector_U[1] =  (f(h,tao*j) + J* array[1] + a*pow(J, 2)*array[0])/ A;

	// for (int k = 2; k < N-1; k++){
	// 	vector_V[k] = (B) / (-A - (B * vector_V[k-1]));
	// 	vector_U[k] = (B*vector_U[k-1] - (f(h*k, tao*j)+J*array[k])) / (-A - B* vector_V[k-1]);
	// }
	
	// vector_V[N-1] = 0;
	// vector_U[N-1] = (B*vector_U[N-2] - (f(N-1*h, tao*j)+J*array[N-1] +a*pow(J, 2)*array[N])) / (- A - B* vector_V[N-2]);

	// //reverse stroke
	// array[N-1] = vector_U[N-1];
	// for(int i = N-2; i>= 0; i--){
	// 	array[i] = vector_V[i] * array[i+1] + vector_U[i];
	// }
	// array[0] = mu1(j*tao);	

	//записываем значения функции в узлах сетки на j-ом слое;
	
	double a;
	double t;
	double s;
	for(int k = 0; k <= N; k++){
		a = u(h*k, j*tao);
		t = array[k];
		fout << a<< " "; 
		file << t << " ";
		s = abs(a - t);
		if(s > maxim){
			maxim = s;
		}
	}
	file << endl;
	fout << endl;
	
	}


	file.close();
	fout.close();
	cout << maxim<< endl;
	cout << fixed << setprecision(10) << tao + pow(h, 2) << endl;
	return 0 ;
}