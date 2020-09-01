#include <iostream>
#include <cmath>

using namespace std;

//Closed-form solution
double closed_solution(double x){
    return 1.0 - (1.0 - exp(-10))*x - exp(-10x);
}

double source_term(double x){
    return 100.0*exp(-10.0*x)
}




int main(){
    
    /* All vectors are at same length
     * n+1 due to last and first being boundry
     * 
     */

// Diagonal vectors 
    double *a_vector = new double[n + 1];
    double *b_vector = new double[n + 1];
    double *c_vector = new double[n + 1];

    for(int i = 0; i < n; i++){
        a_vector[i] = c_vector[i] = -1;
        b_vector[i] = 2;
    }

// Step size
    h = 1 / n
    




// General algorithm gauss




   return 0
}

