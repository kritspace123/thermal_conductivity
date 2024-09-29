#include <iostream>
#include <fstream>
#include <cmath>

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

int main(){
	int N, J;
	ofstream file("value.txt");
	ofstream fout("value_function.txt");

	cout << "Enter value for N (N*h = 1) and J ( J*tao = 1):" << endl;
	cin >> N >> J;

	double tao = double(1)/J;
	double h =double(1)/N;
	double a = 0.032;

	double A = (J + 2* a*pow(N, 2));
	double B = (-a * pow(N, 2));

	cout << A << " " << B << endl;
	double* array = new double[N+1]; 
	double* vector_V = new double[N];
	double* vector_U = new double[N];


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
	//straight stroke
	vector_V[1] = - B / A;
	vector_U[1] =  (f(h,tao*j) + J* array[1] + a*pow(J, 2)*array[0])/ A;

	for (int k = 2; k < N-1; k++){
		vector_V[k] = (B) / (-A - (B * vector_V[k-1]));
		vector_U[k] = (B*vector_U[k-1] - (f(h*k, tao*j)+J*array[k])) / (-A - B* vector_V[k-1]);
	}
	
	vector_V[N-1] = 0;
	vector_U[N-1] = (B*vector_U[N-2] - (f(N-1*h, tao*j)+J*array[N-1] +a*pow(J, 2)*array[N])) / (- A - B* vector_V[N-2]);

	//reverse stroke
	array[N-1] = vector_U[N-1];
	for(int i = N-2; i>= 0; i--){
		array[i] = vector_V[i] * array[i+1] + vector_U[i];
	}
	array[0] = mu1(j*tao);	

	//записываем значения функции в узлах сетки на j-ом слое;
	for(int k = 0; k <= N; k++){
		fout << u(h*k, j*tao)<< " "; 
		file << array[k] << " ";
	}
	file << endl;
	fout << endl;
	}

	file.close();
	fout.close();
	return 0 ;
}