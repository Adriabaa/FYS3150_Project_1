#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

//Writing to file
ofstream ofile;

// Analytical solution
double closed_solution(double x){
    return 1.0 - (1.0 - exp(-10))*x - exp(-10*x);
}

double source_term(double x){
    return 100.0*exp(-10.0*x);
}

void problemB(int n, string filename){
  


// Diagonal vectors 
    double *a_vector = new double[n + 1];
    double *b_vector = new double[n + 1];
    double *c_vector = new double[n + 1];

// Numeric solution
    double *v = new double[n + 1];

// Solution matrix RHS
    double *b_tilde = new double[n + 1];

// Analytical solution
    double *u = new double[n + 1];

// Step size
    double h = 1.0 / (n + 1.0);

// Steps
    double *x = new double[n+1];


    for (int i = 0; i < n+1; i++){
        b_tilde[i] = h*h*source_term(x[i]);
        u[i] = closed_solution(x[i]);
    }

// Initializing the vectors
    for (int i = 0; i <n+1; i++){
        x[i] = i*h;
        a_vector[i] = c_vector[i] = -1 ;
        b_vector[i] = 2;
    }

// General algorithm
    for (int i = 1; i <n+1; i++){
        b_vector[i] = b_vector[i] - a_vector[i] * c_vector[i-1]/b_vector[i-1];
    }
    cout << b_vector[7];

}




int main(){
    problemB(2, "test.txt");

   return 0;
}

