#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>


using namespace std;
using namespace arma;

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
    double *u = new double[n + 1];

// Solution matrix RHS
    double *f_tilde = new double[n + 1];
    double *f = new double[n + 1];

// Analytical solution
    double *v = new double[n + 1];



// Step size
    double h = 1.0 / (n);
    double hh = h*h;

// Steps
    double *x = new double[n+1];

// Initializing the vectors
    for (int i = 0; i <n+1; i++){
        x[i] = i*h;
        a_vector[i] = c_vector[i] = -1 ;
        b_vector[i] = 2;
    }
    
    for (int i = 0; i < n+1; i++){
        f[i] = hh*source_term(x[i]);
        v[i] = closed_solution(x[i]);
    }

//Time start
    clock_t start, finish;
    start = clock();

// General algorithm: Forward substitution
    f_tilde[0] = f[0];
    for (int i = 1; i < n+1; i++){
        b_vector[i] = b_vector[i] - a_vector[i] * c_vector[i-1]/b_vector[i-1];
        f_tilde[i] = f[i] - a_vector[i] * f_tilde[i-1]/b_vector[i-1];
    }

// Backward substitution
    u[n] = f_tilde[n]/b_vector[n];
    for (int i = n-1; i >= 1; i--){
        u[i] =(f_tilde[i] - c_vector[i]*u[i+1]) /b_vector[i];
    }
//Time end
    finish = clock();
    double timeused = (finish - start)/((double) CLOCKS_PER_SEC);
    cout << "Time used for computation=" << timeused << endl;

//File writing
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "     x        u       v" << endl;
    for (int i = 1; i < n+1; i++){        
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << u[i];
        ofile << setw(15) << setprecision(8) << v[i] << endl;
    }

    ofile.close();

    delete [] a_vector;
    delete [] b_vector;
    delete [] c_vector;
    delete [] v;
    delete [] f_tilde;
    delete [] f;
    delete [] u;
    delete [] x;

}

void problemC(int n, string filename){
    /* The specific algorithm is very similar
     so most code is copied from problemB
     what should have been done is to
     have created templates or use inheretance
     but poor planning and inexperience 
     made that difficult.
    */ 
    
// Diagonal vector (c and a are constant at 1)
    double *b_vector = new double[n + 1];

// Numeric solution
    double *u = new double[n + 1];

// Solution matrix RHS
    double *f_tilde = new double[n + 1];
    double *f = new double[n + 1];

// Analytical solution
    double *v = new double[n + 1];

// Step size
    double h = 1.0 / (n);
    double hh = h*h;

// Steps
    double *x = new double[n+1];

// Initializing the vectors
    for (int i = 0; i <n+1; i++){
        x[i] = i*h;
        b_vector[i] = 2;
    }
    
    for (int i = 0; i < n+1; i++){
        f[i] = hh*source_term(x[i]);
        v[i] = closed_solution(x[i]);
    }

//Time start
    clock_t start, finish;
    start = clock();    

// General algorithm: Forward substitution
    f_tilde[0] = f[0];
    for (int i = 1; i < n+1; i++){
        b_vector[i] -= 1/b_vector[i-1];
        f_tilde[i] = f[i] - f_tilde[i-1]/b_vector[i-1];
    }

// Backward substitution
    u[n] = f_tilde[n]/b_vector[n];
    for (int i = n-1; i >= 1; i--){
        u[i] =(f_tilde[i] - u[i+1]) /b_vector[i];
    }

    finish = clock();
    double timeused = (finish - start)/((double) CLOCKS_PER_SEC);
    cout << "Time used for computation=" << timeused << endl;

//Solving for the relativ error
    double error[n+1];
    

    for (int i = 0; i < n; i++){
        error[i] = log10(abs((v[i]-u[i])/u[i]));
    }
    double max = error[0];
    for (int i = 0; i < n; i++){
        if(abs(error[i]) > abs(max)){
            max = error[i];
        }
    }

    cout << "Max value of relative error " << max << endl;

//File writing
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "     x        u       v" << endl;
    for (int i = 1; i < n+1; i++){        
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << u[i];
        ofile << setw(15) << setprecision(8) << v[i] << endl;
    }

    ofile.close();

    delete [] b_vector;
    delete [] v;
    delete [] f_tilde;
    delete [] f;
    delete [] u;
    delete [] x;
}

void problemE(int n){

// Diagonal vectors
    vec a(n+1);
    vec b(n+1);
    vec c(n+1);
    vec x(n+1);
    vec f(n+1);

    

    double h = 1.0 / (n);
    double hh = h*h;

 

// Matrix
    mat A = zeros<mat>(n+1,n+1);


    // Initializing matrix
    for(int i = 0; i <= n; i++){
        for(int j = 0; j <= n; j++){
            if (i = j){
                A(i,j) = b(i);
            } else if(i = j - 1){
                A(i, j) = a(i);
            } else if(i = j + 1){
                A(i, j) = c(i);
            }
        }
    }

    for (int i=0; i < n+1; i++){
        x(i)= i*h;
        f(i)= hh*source_term(x(i));
    }


}




int main(){
    cout << "Problem B:" << endl;
    problemB(10, "Bn1.txt");
    problemB(100, "Bn2.txt");
    problemB(1000, "Bn3.txt");
    problemB(10000, "Bn4.txt");
    problemB(100000, "Bn5.txt");
    problemB(1000000, "Bn6.txt");


    cout << "Problem C:" << endl;
    problemC(10, "n1.txt");
    problemC(100, "n2.txt");
    problemC(1000, "n3.txt");
    problemC(10000, "n4.txt");
    problemC(100000, "n5.txt");
    problemC(100000, "n6.txt");
    cout << "problem E:" << endl;
    problemE(100);
   return 0;
}

