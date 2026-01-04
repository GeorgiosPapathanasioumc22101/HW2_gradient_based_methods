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


// Υπολογίζει μονο την αντίσταση F (Χρήση: Finite Differences)
double calculate_drag_only(double a, double b, double c, int N) {
    double U   = 5.0 + 1.0/12.0;  
    double mi  = 1.7e-5;
    double rho = 1.15;
    double nu  = mi/rho;
    double L   = 1.0;
    double dx  = L / (N - 1);

    double x = 0.0;
    double delta = 1e-6; //αρχική συνθήκη του δ

    // Προσωρινός πίνακας για τον υπολογισμό του F 
    double *temp_delta = malloc(N * sizeof(double));
    temp_delta[0] = delta;

    // RK2 
    for(int i=0; i < N-1; i++) {
        double k1 = dx * f_primal(x,          delta,          U, nu, a, b, c);
        double k2 = dx * f_primal(x + 0.5*dx, delta + 0.5*k1, U, nu, a, b, c);
        delta += k2;
        x += dx;
        temp_delta[i+1] = delta;
    }

    // Υπολογισμός F
    double Ftot = 0.0;
    for(int i=0; i < N-1; i++) {
        double tau1 = ((PI_CONST/2.0) * mi * U) / temp_delta[i];
        double tau2 = ((PI_CONST/2.0) * mi * U) / temp_delta[i+1];
        Ftot += 0.5 * (tau1 + tau2) * dx;
    }

    free(temp_delta);
    return Ftot;
}


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

// =========================================================
// 3. ΣΥΝΑΡΤΗΣΕΙΣ RHS (DD & ADJOINT)
// =========================================================

double dd_rhs_a(double x, double lamda, double delta, double U, double nu, double a, double b, double c)
{
    double k = (2.0/PI_CONST - 0.5) * (U/nu);
    double vo_val = v0(x, a, b, c);
    double ddelta_dx = (PI_CONST/2.0 - (vo_val * delta)/nu) / (k * delta);
    
    double source_term = 0.0;
    if (x >= 0.2 && x <= 0.8) {
        double P = a*x*x + b*x + c*sin(PI_CONST*x/6.0);
        source_term = 2.0 * ((P * x * x * delta) / nu); 
    }
    double num = - k * lamda * ddelta_dx + source_term - ((vo_val * lamda)/(nu)); 
    return num / (k * delta);
}

double dd_rhs_b(double x, double lamda, double delta, double U, double nu, double a, double b, double c)
{
    double k = (2.0/PI_CONST - 0.5) * (U/nu);
    double vo_val = v0(x, a, b, c);
    double ddelta_dx = (PI_CONST/2.0 - (vo_val * delta)/nu) / (k * delta);

    double source_term = 0.0;
    if (x >= 0.2 && x <= 0.8) {
        double P = a*x*x + b*x + c*sin(PI_CONST*x/6.0);
        source_term = 2.0 * ((P * x * delta) / nu);
    }
    double num = - k * lamda * ddelta_dx + source_term - ((vo_val * lamda)/(nu)); 
    return num / (k * delta);
}

double dd_rhs_c(double x, double lamda, double delta, double U, double nu, double a, double b, double c)
{
    double k = (2.0/PI_CONST - 0.5) * (U/nu);
    double vo_val = v0(x, a, b, c);
    double ddelta_dx = (PI_CONST/2.0 - (vo_val * delta)/nu) / (k * delta);

    double source_term = 0.0;
    if (x >= 0.2 && x <= 0.8) {
        double P = a*x*x + b*x + c*sin(PI_CONST*x/6.0);
        source_term = 2.0 * ((P * sin(PI_CONST*x/6.0) * delta) / nu);
    }
    double num = - k * lamda * ddelta_dx + source_term - ((vo_val * lamda)/(nu)); 
    return num / (k * delta);
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
    double a=1.0, b=1.0, c=-5.494219125;
    double L = 1.0;
    double dx = L/(N-1);
    
    double U   = 5.0 + 1.0/12.0;  
    double mi  = 1.7e-5;
    double rho = 1.15;
    double nu  = mi/rho;

    printf("Nodes: %d, dx: %.5e\n", N, dx);

    // -----------------------------------------------------
    // 1. FINITE DIFFERENCES (FD)
    // -----------------------------------------------------
    printf("Calculation FD gradients...\n");
    double h = 1e-4;
    double F_a = (calculate_drag_only(a+h, b, c, N) - calculate_drag_only(a-h, b, c, N))/(2*h);
    double F_b = (calculate_drag_only(a, b+h, c, N) - calculate_drag_only(a, b-h, c, N))/(2*h);
    double F_c = (calculate_drag_only(a, b, c+h, N) - calculate_drag_only(a, b, c-h, N))/(2*h);

    // -----------------------------------------------------
    // 2. ΥΠΟΛΟΓΙΣΜΟΣ ΠΕΔΙΟΥ DELTA (Για DD & Adjoint)
    // -----------------------------------------------------
    printf("Solving Primal field (storing delta)...\n");
    double *delta_arr = calloc(N, sizeof(double));
    compute_delta_field(a, b, c, N, delta_arr);

    // -----------------------------------------------------
    // 3. DIRECT DIFFERENTIATION (DD)
    // -----------------------------------------------------
    printf("Calculation DD gradients...\n");
    double *lam_a = calloc(N, sizeof(double));
    double *lam_b = calloc(N, sizeof(double));
    double *lam_c = calloc(N, sizeof(double));
    //Reset του χ και συνοριακή συνθήκη του λ
    double x = 0.0;
    double la = 0.0, lb = 0.0, lc = 0.0;
    //RK-2 για τον υπολογισμό τον λ1,λ2,λ3
    for(int i=0; i < N-1; i++) {
        double d_curr = delta_arr[i];
        double d_mid  = 0.5 * (delta_arr[i] + delta_arr[i+1]);

        // Lambda A
        double k1a = dx * dd_rhs_a(x, la, d_curr, U, nu, a, b, c);
        double k2a = dx * dd_rhs_a(x + 0.5*dx, la + 0.5*k1a, d_mid, U, nu, a, b, c);
        la += k2a; lam_a[i+1] = la;

        // Lambda B
        double k1b = dx * dd_rhs_b(x, lb, d_curr, U, nu, a, b, c);
        double k2b = dx * dd_rhs_b(x + 0.5*dx, lb + 0.5*k1b, d_mid, U, nu, a, b, c);
        lb += k2b; lam_b[i+1] = lb;

        // Lambda C
        double k1c = dx * dd_rhs_c(x, lc, d_curr, U, nu, a, b, c);
        double k2c = dx * dd_rhs_c(x + 0.5*dx, lc + 0.5*k1c, d_mid, U, nu, a, b, c);
        lc += k2c; lam_c[i+1] = lc;

        x += dx;
    }

    // -----------------------------------------------------
    // 4. ADJOINT METHOD
    // -----------------------------------------------------
    printf("Calculation Adjoint gradients...\n");
    double *psi_mtr = calloc(N, sizeof(double));
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

    // -----------------------------------------------------
    // 5. Οοκλήρωση για τον υπολογισμό των παραγώγων 
    // -----------------------------------------------------
    double common_const = - (PI_CONST * mi * U) / 2.0;

    // A. DD Integration
    double dFda_DD = 0.0, dFdb_DD = 0.0, dFdc_DD = 0.0;
    for(int i=0; i < N-1; i++) {
        dFda_DD += 0.5 * (lam_a[i]/(pow(delta_arr[i],2)) + lam_a[i+1]/(pow(delta_arr[i+1],2))) * dx;
        dFdb_DD += 0.5 * (lam_b[i]/(pow(delta_arr[i],2)) + lam_b[i+1]/(pow(delta_arr[i+1],2))) * dx;
        dFdc_DD += 0.5 * (lam_c[i]/(pow(delta_arr[i],2)) + lam_c[i+1]/(pow(delta_arr[i+1],2))) * dx;
    }
    dFda_DD *= common_const; dFdb_DD *= common_const; dFdc_DD *= common_const;

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
    writeField1d("Adjoint_field_psi_100001.txt" , psi_mtr , N);
    // -----------------------------------------------------
    // 6. ΕΚΤΥΠΩΣΗ ΑΠΟΤΕΛΕΣΜΑΤΩΝ
    // -----------------------------------------------------
    printf("\n--- FINAL RESULTS ---\n");
    printf("Param |    Finite Diff (FD)    | Direct Diff (DD)       | Adjoint Method       | Error (FD-Adj)\n");
    printf("--------------------------------------------------------------------------------------------\n");
    printf("  a   | %20.12e | %20.12e | %20.12e | %.2e\n", F_a, dFda_DD, dFda_Adj, fabs(F_a - dFda_Adj));
    printf("  b   | %20.12e | %20.12e | %20.12e | %.2e\n", F_b, dFdb_DD, dFdb_Adj, fabs(F_b - dFdb_Adj));
    printf("  c   | %20.12e | %20.12e | %20.12e | %.2e\n", F_c, dFdc_DD, dFdc_Adj, fabs(F_c - dFdc_Adj));

    // Cleanup
    free(delta_arr);
    free(lam_a); free(lam_b); free(lam_c);
    free(psi_mtr);

    return 0;
}