#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <time.h>


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
    double *f_tilde = new double[n + 1];
    double *f = new double[n + 1];

// Analytical solution
    double *u = new double[n + 1];



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
        u[i] = closed_solution(x[i]);
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
    v[n] = f_tilde[n]/b_vector[n];
    for (int i = n-1; i >= 1; i--){
        v[i] =(f_tilde[i] - c_vector[i]*v[i+1]) /b_vector[i];
    }
//Time end
    finish = clock();
    double timeused = (finish - start)/((double) CLOCKS_PER_SEC);
    cout << "Time used for computation=" << timeused << endl;

//File writing
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "     x        v       u" << endl;
    for (int i = 1; i < n+1; i++){        
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << v[i];
        ofile << setw(15) << setprecision(8) << u[i] << endl;
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
    double *v = new double[n + 1];

// Solution matrix RHS
    double *f_tilde = new double[n + 1];
    double *f = new double[n + 1];

// Analytical solution
    double *u = new double[n + 1];

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
        u[i] = closed_solution(x[i]);
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
    v[n] = f_tilde[n]/b_vector[n];
    for (int i = n-1; i >= 1; i--){
        v[i] =(f_tilde[i] - v[i+1]) /b_vector[i];
    }

    finish = clock();
    double timeused = (finish - start)/((double) CLOCKS_PER_SEC);
    cout << "Time used for computation=" << timeused << endl;

//Solving for the relativ error


//File writing
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "     x        v       u" << endl;
    for (int i = 1; i < n+1; i++){        
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << v[i];
        ofile << setw(15) << setprecision(8) << u[i] << endl;
    }

    ofile.close();

    delete [] b_vector;
    delete [] v;
    delete [] f_tilde;
    delete [] f;
    delete [] u;
    delete [] x;
}




int main(){
    problemB(10, "big.txt");
    problemB(100, "bigger.txt");
    problemB(1000, "biggest.txt");

   return 0;
}

