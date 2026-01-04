#include <stdlib.h>
#include <math.h>
#include <stdio.h>

const double PI_CONST = 3.141592653589793;


void writeField1d(const char *filename, double *field, int N)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) { printf("Error opening %s\n", filename); return; }

    for (int i = 0; i < N; i++)
        fprintf(fp, "%.15e\n", field[i]);

    fclose(fp);
}
// Συνάρτηση ταχύτητας έγχυσης v0(x)
double v0(double x, double a, double b, double c) {
    if (x < 0.2 || x > 0.8) return 0.0;
    double term = a*x*x + b*x + c*sin(PI_CONST*x/6.0);
    return -(term * term);
}


double f_primal(double x, double delta, double U, double nu, double a, double b, double c) {
    double vo_val = v0(x, a, b, c);
    double num = PI_CONST/2.0 - (vo_val/nu)*delta;
    double den = ((2.0/PI_CONST) - 0.5) * (U/nu) * delta;
    return num / den;
}
//  Υπολογισμός δ 
void compute_delta_field(double a, double b, double c, int N, double *delta_arr) {
    double U   = 5.0 + 1.0/12.0;  
    double mi  = 1.7e-5;
    double rho = 1.15;
    double nu  = mi/rho;
    double L   = 1.0;
    double dx  = L / (N - 1);

    double x = 0.0;
    double delta = 1e-6;
    
    delta_arr[0] = delta;

    for(int i=0; i < N-1; i++) {
        double k1 = dx * f_primal(x,          delta,          U, nu, a, b, c);
        double k2 = dx * f_primal(x + 0.5*dx, delta + 0.5*k1, U, nu, a, b, c);
        delta += k2;
        x += dx;
        delta_arr[i+1] = delta;
    }
}
double omega_p_min(double omega_p, double omega_p_max){
    double omega = 0.0;
    if (omega_p < omega_p_max)
    {
         omega = omega_p;
    }
    else omega = omega_p_max;
    return omega;
}
// Εξίσωση Adjoint (FAE)
double FAE(double x, double psi, double delta, double U, double nu, double a, double b, double c, double mi) {
    double k = (2.0/PI_CONST - 0.5) * U;

    double source_term_v0 = 0.0;
    if (x >= 0.2 && x <= 0.8) {
         double P = a*x*x + b*x + c*sin(PI_CONST*x/6.0);
         source_term_v0 = -(P*P);
    }
    
    double num = (psi * source_term_v0) - ((PI_CONST * mi * U)/(2.0 * delta * delta)); 
    double den = k * delta;
    
    return num / den; 
}

int main() {
    int N = 100001;                 
    double L = 1.0;
    double dx = L/(N-1);
    double htta = 0.2;
    double U   = 5.0 + 1.0/12.0;  
    double mi  = 1.7e-5;
    double rho = 1.15;
    double nu  = mi/rho;
    double a,b,c;
    double *delta_arr = calloc(N, sizeof(double));
    double *psi_mtr   = calloc(N, sizeof(double));
    double dPhi_da ,dPhi_db, dPhi_dc; 
    double lamda; 
    double dcda ,dcdb ,dcdc; 
    double omega_p = 0.25;
    double omega_p_max = 10; 
    double c_periorismos;
    double gamma = 1.01;
    //---------------------------------------//
    // 1.Αρχικοποίηση μεταβλητών σχεδιασμού  //
    //---------------------------------------//
    a=1;
    b=1;
    c=-5.494219125;
    lamda = -0.5;
    int Max_iterations = 10000;
    printf("--Value of F--\t\t--Constraint--\t\t\t  a\t\t  b\t\t  c\n");
    FILE *F = fopen("results_F.txt", "w");
    FILE *Iterations = fopen("Iterations.txt", "w");
    //-----------------------------------------------------//
    // Επαναληπτική διαδικαδία αλγορίθμου βελτιστοποίησης  //
    //-----------------------------------------------------//
    for (int i = 0; i < Max_iterations; i++)
    {
    double x = 0.0;
    
    //----------------------------------------------//
    //   2. Λύση του Primal και αποθήκευση του δ    //
    //----------------------------------------------//
    
    compute_delta_field(a, b, c, N, delta_arr);

    //-------------------//
    //   Εκτύπωση της F  //
    //-------------------//
    
    double Ftot = 0.0;
    for(int i=0; i < N-1; i++) {
        double tau1 = ((PI_CONST/2.0) * mi * U) / delta_arr[i];
        double tau2 = ((PI_CONST/2.0) * mi * U) / delta_arr[i+1];
        Ftot += 0.5 * (tau1 + tau2) * dx;
    }
    
    if (i % 100 == 0) {
            printf("%-20.10e \t",Ftot);
            fprintf(F, "%.15e\n", Ftot);
            fprintf(Iterations, "%d\n", i);
        }
    //---------------------------------------------------//
    //   3. Λύση και αποθήκευση του συζυγούς πεδίου Ψ    //
    //---------------------------------------------------//
    //Reset του χ στο άκρο L για να λυθεί από 'αριστερά' προς τα 'δεξιά' η FAE
    x = L;
    double psi = 0.0;
    psi_mtr[N-1] = psi; 
    
    double neg_dx = -dx; 

    for(int i=N-2; i >=0; i--) {
        double d_curr = delta_arr[i+1];
        double d_mid  = 0.5 * (delta_arr[i+1] + delta_arr[i]);
        
        double k1d = neg_dx * FAE(x,          psi,          d_curr, U, nu, a, b, c, mi);
        double k2d = neg_dx * FAE(x + 0.5*neg_dx, psi + 0.5*k1d, d_mid,  U, nu, a, b, c, mi);
        
        psi += k2d;
        psi_mtr[i] = psi;
        x += neg_dx;
    }
    
    
    // B. Adjoint Integration
    double dFda_Adj = 0.0, dFdb_Adj = 0.0, dFdc_Adj = 0.0;
    x = 0.0;
    
    for(int i=0; i < N-1; i++) {
        double src_a_curr=0, src_b_curr=0, src_c_curr=0;
        double src_a_next=0, src_b_next=0, src_c_next=0;
        double x_next = x + dx;
        
        // Σημείο i
        if(x >= 0.2 && x <= 0.8) {
             double P = a*x*x + b*x + c*sin(PI_CONST*x/6.0);
             src_a_curr = P * x * x;     
             src_b_curr = P * x;
             src_c_curr = P * sin(PI_CONST*x/6.0);
        }
        
        // Σημείο i+1
        if(x_next >= 0.2 && x_next <= 0.8) {
             double P = a*x_next*x_next + b*x_next + c*sin(PI_CONST*x_next/6.0);
             src_a_next = P * x_next * x_next;
             src_b_next = P * x_next;
             src_c_next = P * sin(PI_CONST*x_next/6.0);
        }

        
        dFda_Adj += 0.5 * (-2.0*delta_arr[i]*psi_mtr[i]*src_a_curr + -2.0*delta_arr[i+1]*psi_mtr[i+1]*src_a_next) * dx;
        dFdb_Adj += 0.5 * (-2.0*delta_arr[i]*psi_mtr[i]*src_b_curr + -2.0*delta_arr[i+1]*psi_mtr[i+1]*src_b_next) * dx;
        dFdc_Adj += 0.5 * (-2.0*delta_arr[i]*psi_mtr[i]*src_c_curr + -2.0*delta_arr[i+1]*psi_mtr[i+1]*src_c_next) * dx;
        
        x += dx;
    }

    //------------------------------------------------------------//
    // 4.   Διορθωση μεταβλητών σχεδιασμού (Steepest Descent)     //
    //------------------------------------------------------------//

    // Περιορισμός με τη μέθοδο ALM//
    c_periorismos = 0.065472*a*a + 0.204*a*b + 0.168*b*b + 0.104738*a*c + 0.172816*b*c + 0.044445*c*c - 0.254166;
    dcda = 0.130944*a + 0.204*b + 0.104738*c;
    dcdb = 0.204*a + 0.336*b + 0.172816*c;
    dcdc = 0.194738*a + 0.172816*b + 0.08889*c;
    dPhi_da = dFda_Adj - lamda * dcda + 2 * omega_p * c_periorismos * dcda;
    dPhi_db = dFdb_Adj - lamda * dcdb + 2 * omega_p * c_periorismos * dcdb;
    dPhi_dc = dFdc_Adj - lamda * dcdc + 2 * omega_p * c_periorismos * dcdc;
    
    a = a - htta * dPhi_da;
    b = b - htta * dPhi_db;
    c = c - htta * dPhi_dc;

    c_periorismos = 0.065472*a*a + 0.204*a*b + 0.168*b*b + 0.104738*a*c + 0.172816*b*c + 0.044445*c*c - 0.254166;
    lamda = lamda - 2 * omega_p * c_periorismos;
    omega_p = omega_p * gamma;
    omega_p = omega_p_min(omega_p , omega_p_max);
    
     if (i % 100 == 0) {
           printf("%10.e\t iteration %d\t%lf\t%lf\t%lf\n" , c_periorismos,i , a,b,c);
        }
    }

    
    free(psi_mtr);
    free(delta_arr); 

return 0;
}