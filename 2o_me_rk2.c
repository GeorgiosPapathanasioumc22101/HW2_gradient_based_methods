#include <stdlib.h>
#include <math.h>
#include <stdio.h>
double v0(double x, double a, double b, double c)
{
    const double pi = 3.141592653589793;
        if (x < 0.2 || x > 0.8)
        return 0.0;
    double term = a*x*x + b*x + c*sin(pi*x/6.0);
    return -(term * term);
}   
double f(double x, double delta, double U, double nu, double a, double b, double c)
{
    const double pi = 3.141592653589793;
    double vo = v0(x, a, b, c);

    double num = pi/2.0 - (vo/nu)*delta;
    double den = ((2.0/pi) - 0.5) * (U/nu) * delta;

    return num / den;
}
    

double calculate_F(double a, double b , double c , int N)
    
{

    double U  = 5.0 + 1.0/12.0;  
    double mi = 1.7e-5;
    double rho= 1.15;
    double nu = mi/rho;

    double L = 1.0;
    double dx = L/(N-1);

    const double pi = 3.141592653589793;

    double *x_arr     = calloc(N, sizeof(double));
    double *delta_arr = calloc(N, sizeof(double));
    double *tau_arr   = calloc(N, sizeof(double));
    double *v0_arr    = calloc(N, sizeof(double));

    double x = 0.0;
    double delta = 1e-6;

    x_arr[0] = x;
    delta_arr[0] = delta;
    v0_arr[0] = v0(x, a,b,c);

    // RK4 in x
    for(int i=0; i < N-1; i++)
    {
        double k1 = dx * f(x,          delta,            U, nu, a,b,c);
        double k2 = dx * f(x+0.5*dx,   delta+0.5*k1,     U, nu, a,b,c);
        
        delta += k2;
        x += dx;

        x_arr[i+1] = x;
        delta_arr[i+1] = delta;
        v0_arr[i+1] = v0(x, a,b,c);
    }
    //Υπολογισμός της F
        // tau_w(x)
    for(int i=0; i < N; i++)
        tau_arr[i] = ((pi/2.0) * mi * U) / delta_arr[i];

    // Total drag F
    double Ftot = 0.0;
    for(int i=0; i < N-1; i++)
        Ftot += 0.5*(tau_arr[i] + tau_arr[i+1]) * dx;

    //printf("F = %.15e\n", Ftot);
    free(x_arr);
    free(delta_arr);
    free(tau_arr);
    free(v0_arr);
    return Ftot;

}

int main () {
    double L = 1.0;
    int N = 100001;
    double dx = L/((double)N-1.0);
    double a, b, c;
    a=b=1;
    c = -5.494219125;
    double eps = 1e-6;
    calculate_F(a,b,c,N);
    calculate_F(a+eps, b, c, N);
    calculate_F(a-eps,b,c,N);
    calculate_F(a, b+eps, c, N);
    calculate_F(a,b-eps,c,N);
    calculate_F(a, b, c+eps, N);
    calculate_F(a,b,c-eps,N);


    double F_a = (calculate_F(a+eps, b, c, N) - calculate_F(a-eps, b, c, N))/(2*eps);
    double F_b = (calculate_F(a, b+eps, c, N) - calculate_F(a, b-eps, c, N))/(2*eps);
    double F_c = (calculate_F(a, b, c+eps, N) - calculate_F(a, b, c-eps, N))/(2*eps);

    printf("Derivative of F with respect to a is %20.12e\n" , F_a);
    printf("Derivative of F with respect to b is %20.12e\n" , F_b);
    printf("Derivative of F with respect to c is %20.12e\n" , F_c);
    printf("%d\n" , N);
    printf("%.5e\n" , dx);



    return 0;
}
