#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void writeField1d(const char *filename, double *field, int N)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) { printf("Error opening %s\n", filename); return; }

    for (int i = 0; i < N; i++)
        fprintf(fp, "%.15e\n", field[i]);

    fclose(fp);
}
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
double Calculate_F(int N, double *tau_arr, double pi, double mi, double U, double *delta_arr, double dx){
        for(int i=0; i < N; i++)
        tau_arr[i] = ((pi/2.0) * mi * U) / delta_arr[i];

    
    double Ftot = 0.0;
    for(int i=0; i < N-1; i++)
        Ftot += 0.5*(tau_arr[i] + tau_arr[i+1]) * dx;

    printf("F = %.15e\n", Ftot);
}

int main()
{

    int N = 100001;                 
    double a=1.0, b=1.0, c=-5.494219125;

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

    // RK2 in x
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



    // u(y) at x=L/2
    //διακριτοποίηση του y
    int i_mid = (N-1)/2;
    double delta_mid = delta_arr[i_mid];

    int Ny = 100001;
    double *y_arr = calloc(Ny, sizeof(double));
    double *u_arr = calloc(Ny, sizeof(double));

    double dy = delta_mid/(Ny-1);
    for(int j=0; j < Ny; j++)
    {
        y_arr[j] = j*dy;
        u_arr[j] = U * sin( (pi/2.0) *( y_arr[j] / delta_mid ) );
    }

    
    Calculate_F( N, tau_arr, pi, mi, U, delta_arr, dx);
    
    // write files
    writeField1d("x_100001_rk2.txt", x_arr, N);
    writeField1d("delta_100001_rk2.txt", delta_arr, N);
    writeField1d("tau_100001._rk2.txt", tau_arr, N);
    writeField1d("v0_100001_rk2.txt", v0_arr, N);

    writeField1d("y_100001_rk2.txt", y_arr, Ny);
    writeField1d("u_xmid_100001_rk2.txt", u_arr, Ny);

    free(x_arr);
    free(delta_arr);
    free(tau_arr);
    free(v0_arr);
    free(y_arr);
    free(u_arr);

    return 0;
}
