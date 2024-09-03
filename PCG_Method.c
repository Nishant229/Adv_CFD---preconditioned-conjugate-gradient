#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double matrix_Vector_multiplication(int N, int Nj, double A[N+1][5], double B[N+1], double C[N+1]);
double vector_vector_multiplication(int N, double G[N+1], double H[N+1]);

double Initializing_A (int N, int Ni, int Nj, double A[N+1][5], double b[N+1], double T_left, double T_right, double T_top, double T_bottom );

double conjugate_gradient_method(int N, int Nj, double A[N+1][5], double T_initial[N+1], double T[N+1], double b[N+1], double epsilon, double residuals[1000] );
double pre_conditioner_conjugate_gradient_method(int N, int Nj, double A[N+1][5], double T_initial[N+1], double T[N+1], double b[N+1], double epsilon, char Pre_Conditioner_type, double residuals[1000] );

double Jacobi_pre_conditioner( int N, int Nj, double A[N+1][5], double Q_0[N+1], double F_0[N+1] );
double general_LU_pre_conditioner( int N, int Nj, double L[N+1][3], double U[N+1][3], double Q[N+1], double F[N+1] );
double ILU_matrix( int N, int Nj, double A[N+1][5], double L[N+1][3], double U[N+1][3] );
double SIP_matrix( int N, int Nj, double alpha, double A[N+1][5], double L[N+1][3], double U[N+1][3] );
void print_1D_Array(int N, double V[N+1]);

int main() {
	int Ni, Nj;
	/*printf("Please enter the number of points in x-direction --> ");
	scnaf("%d",%Ni);
	printf("Please enter the number of points in y-direction --> ");
	scnaf("%d",%Nj);*/
	Ni=128;
	Nj=128;
	int  N = Ni*Nj;
	
	double dx = 1.0/Ni;
	double dy = 1.0/Nj;
	
	// Boundary Conditions
	double T_left = 0;
	double T_right = 0;
	double T_top = 1.0;
	double T_bottom = 0;
	
	double epsilon = pow(10,-6);
	
	double T_initial[N+1], T_actual[N+1];
	double T_CG[N+1], T_PCG_Jacobi[N+1], T_PCG_ILU[N+1],T_PCG_SIP[N+1];
	double A[N+1][5]; // A[n][0] = A_W[n], A[n][1] = A_S[n], A[n][2] = A_P[n], A[n][3] = A_N[n], A[n][4] = A_E[n]
	double b[N+1];
	
	for(int i=0;i<=N;i++) {
		T_initial[i] = 0;  // Starting with zero initial guess
	}
	
	// Initializing A Matrix according to Boundary Conditions
	Initializing_A ( N, Ni, Nj, A, b, T_left, T_right, T_top, T_bottom );
	
	// Actual Solution
	double x,y;
	for(int i=1;i<=Ni;i++) {
		for(int j=1;j<=Nj;j++) {
			int n = Nj*(i-1) + j;
			x = (i-1)*dx + dx/2.0;
			y = (j-1)*dy + dy/2.0;
			T_actual[n] = 0.0;
			for(int k=1;k<=50;k++) {
				T_actual[n] = T_actual[n] + (( pow(-1,k+1) + 1 )/k)*sin(k*M_PI*x)*sinh(k*M_PI*y)/sinh(k*M_PI);
			}
			T_actual[n] = T_actual[n]*2/M_PI;
			// printf("%lf \t",T_actual[n]); // Testing
		}
		// printf("\n"); // Testing
	}
	
// Testing 
/*for(int i=1;i<=N;i++) {
	printf("n=%d, A_P = %lf, A_E = %lf, A_W = %lf, A_N = %lf, A-S = %lf, Source = %lf \n",i,A_P[i],A_E[i],A_W[i],A_N[i],A_S[i], b[i]);
}*/
	// Executing Functions
	
	double CG_residuals[1000];
	double PCG_Jacobi_residuals[1000];
	double PCG_ILU_residuals[1000];
	double PCG_SIP_residuals[1000];
	
	// Calling Conjugate Gradient method
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("------------ Conjugate Gradient method--------------\n");
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	int CG_iterations = conjugate_gradient_method( N, Nj, A, T_initial, T_CG, b, epsilon, CG_residuals );
	printf("\n \t The number of Iterations for Conjugate Gradient method are %d \n\n",CG_iterations);
	
	
	// Calling Jacobi Pre_conditioner Conjugate Gradient method
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("------------Jacobi Pre_conditioner Conjugate Gradient method--------------\n");
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	int PCG_Jacobi_iterations = pre_conditioner_conjugate_gradient_method( N, Nj, A, T_initial, T_PCG_Jacobi, b, epsilon,  'J' , PCG_Jacobi_residuals);
	printf("\n \t The number of Iterations for Pre Conditioner Jacobi Conjugate Gradient method are %d \n\n",PCG_Jacobi_iterations);
	
	// Calling ILU Pre_conditioner Conjugate Gradient method
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("------------ILU Pre_conditioner Conjugate Gradient method--------------\n");
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	int PCG_ILU_iterations = pre_conditioner_conjugate_gradient_method( N, Nj, A, T_initial, T_PCG_ILU, b, epsilon,  'I' , PCG_ILU_residuals);
	printf("\n \t The number of Iterations for Pre Conditioner ILU Conjugate Gradient method are %d \n\n",PCG_ILU_iterations);
	
	// Calling SIP Pre_conditioner Conjugate Gradient method
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	printf("------------SIP Pre_conditioner Conjugate Gradient method--------------\n");
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	int PCG_SIP_iterations = pre_conditioner_conjugate_gradient_method( N, Nj, A, T_initial, T_PCG_SIP, b, epsilon,  'S' , PCG_SIP_residuals);
	printf("\n \t The number of Iterations for Pre Conditioner SIP Conjugate Gradient method are %d \n\n",PCG_SIP_iterations);
	
	// Printing Results
	FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6;
	fp1=fopen("Iterations_vs_residual.dat","w");
	fp2=fopen("Temperature_Contours.dat","w");
	fp3=fopen("Temp_centerline_x_vs_y.dat","w");
	fp4=fopen("Temp_centerline_y_vs_x.dat","w");
	fp5=fopen("L_and_U_Diagonals_ILU.dat","w");
	fp6=fopen("L_and_U_Diagonals_SIP.dat","w");
	
	// Iteration vs Residual
	fprintf(fp1,"Iteration \t Conj_Grad \t Jacobi_Pre_Conj_Grad \t ILU_Pre_Conj_Grad \t SIP_Pre_Conj_Grad \n\n");	
	for(int i=1;i<=CG_iterations;i++) {
		fprintf(fp1,"%d \t",i);
		if (i<=CG_iterations) { fprintf(fp1,"%lf \t",CG_residuals[i]); }
		else { fprintf(fp1,"0 \t"); }
		if (i<=PCG_Jacobi_iterations) { fprintf(fp1,"%lf \t",PCG_Jacobi_residuals[i]); }
		else { fprintf(fp1,"0 \t"); }
		if (i<=PCG_ILU_iterations) { fprintf(fp1,"%lf \t",PCG_ILU_residuals[i]); }
		else { fprintf(fp1,"0 \t"); }
		if (i<=PCG_SIP_iterations) { fprintf(fp1,"%lf \t",PCG_SIP_residuals[i]); }
		else { fprintf(fp1,"0 \t"); }
	}
	
	// Contour Plots
	fprintf(fp2,"ZONE I=%d, J=%d \n",Ni,Nj);
	for(int j=1;j<=Nj;j++) {
		y=(j-1)*(dy) + dy/2.0;
		for(int i=1;i<=Ni;i++) {
			x=(i-1)*dx + dx/2.0;
			int n = Nj*(i-1) + j;
	fprintf(fp2,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",x,y,T_actual[n],T_CG[n],T_PCG_Jacobi[n],T_PCG_ILU[n],T_PCG_SIP[n]);
		}
	}
	
	// Temperature variation with y along mid-vertical plane x = 0.5
	fprintf(fp3,"Point(n) \t y \t T_actual \t T_CG \t T_PCG_Jacobi \t T_PCG_ILU \t T_PCG_SIP \n\n");
	for(int j=1;j<=Nj;j++) {
		int n = Nj*((Ni/2)-1) + j;
		y = (j-1)*dy + dy/2.0;
	fprintf(fp3,"%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf\n",n,y,T_actual[n],T_CG[n],T_PCG_Jacobi[n],T_PCG_ILU[n],T_PCG_SIP[n]);
	}
	
	// Temperature variation with x along mid-vertical plane y = 0.5
	fprintf(fp4,"Point(n) \t x \t T_actual \t T_CG \t T_PCG_Jacobi \t T_PCG_ILU \t T_PCG_SIP \n\n");
	for(int i=1;i<=Ni;i++) {
		int n = Nj*(i-1) + Nj/2;
		x=(i-1)*dx + dx/2.0;
	fprintf(fp4,"%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf\n",n,x,T_actual[n],T_CG[n],T_PCG_Jacobi[n],T_PCG_ILU[n],T_PCG_SIP[n]);
	}
	
	// PLease input Ni and Nj for which you wish to get the L and U Diagonals from ILU & SIP factorization
	Ni=4, Nj=4; 
	N = Ni*Nj;
	double W[N+1][5];
	double v[N+1];
	double L_Sample[N+1][3], U_Sample[N+1][3];
	Initializing_A ( N, Ni, Nj,  W, v, T_left, T_right, T_top, T_bottom );
	
	// L & U Diagonals from ILU factorization
	ILU_matrix( N, Nj, W, L_Sample, U_Sample );
	fprintf(fp5,"Ni = %d , Nj = %d , N = %d \n\n",Ni,Nj,N);
	fprintf(fp5,"Point(n)  L_W \t\t L_S \t\t L_P \t\t U_N \t\t U_E \n\n");
	for(int n=1;n<=N;n++) {
		fprintf(fp5,"%d \t %lf \t %lf \t %lf \t %lf \t %lf \n",n,L_Sample[n][0],L_Sample[n][1],L_Sample[n][2],U_Sample[n][1],U_Sample[n][2]);
	}
	
	// L & U Diagonals from SIP factorization
	SIP_matrix( N, Nj, 0.5, W, L_Sample, U_Sample ); // alpha = 0.5
	fprintf(fp6,"Ni = %d , Nj = %d , N = %d \n\n",Ni,Nj,N);
	fprintf(fp6,"Point(n)  L_W \t\t L_S \t\t L_P \t\t U_N \t\t U_E \n\n");
	for(int n=1;n<=N;n++) {
		fprintf(fp6,"%d \t %lf \t %lf \t %lf \t %lf \t %lf \n",n,L_Sample[n][0],L_Sample[n][1],L_Sample[n][2],U_Sample[n][1],U_Sample[n][2]);
	}
	
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);
	
	return 1;
}

double Initializing_A ( int N, int Ni, int Nj, double A[N+1][5], double b[N+1], double T_left, double T_right, double T_top, double T_bottom ) {
	
	double A_P[N+1], A_E[N+1],A_W[N+1], A_N[N+1], A_S[N+1];
	double dx = 1.0/Ni;
	double dy = 1.0/Nj;
	
	// Initializing Arrays with  zero values
	for(int i=0;i<=N;i++) {
		A_P[i] = 0; // Point
		A_E[i] = 0; // East Neighbour
		A_W[i] = 0; // West Neighbour
		A_N[i] = 0; // North Neighbour
		A_S[i] = 0; // South Neighbour
	}
	
	
	// Initializing Arrays with Coefficients
	for(int i=1;i<=Ni;i++) {
		for(int j=1;j<=Nj;j++) {
			int n = Nj*(i-1) + j;
			double Sp=0, Su=0;
			if (j==1) { // Bottom Boundary
				A_S[n] = 0;
				Sp = Sp - 2*dx/(dy);
				Su = Su + 2*T_bottom*dx/(dy);
			}
			else {
				A_S[n] = dx/(dy);
			}
			if (j==Nj) { // Top Boundary
				A_N[n] = 0;
				Sp = Sp - 2*dx/(dy);
				Su = Su + 2*T_top*dx/(dy);
				//printf("%lf \t",Su);
			}
			else {
				A_N[n] = dx/(dy);
			}
			if (i==1) { // Left Boundary
				A_W[n] = 0;
				Sp = Sp - 2*dy/(dx);
				Su = Su + 2*T_left*dy/(dx);
			}
			else {
				A_W[n] = dy/(dx);
			}
			if (i==Ni) { // Right Boundary
				A_E[n] = 0;
				Sp = Sp - 2*dy/(dx);
				Su = Su + 2*T_right*dy/(dx);
			}
			else {
				A_E[n] = dy/(dx);
			}
			A_P[n] = -1.0 * ( A_E[n] + A_W[n] + A_N[n] + A_S[n] - Sp );
			b[n] = -1.0*Su;
		}
	}
	
	// Storing 5 columns vectors of A matrix in A 2D matrix for transferring it to the function
	for(int i=1;i<=N;i++) {
			A[i][0] = A_W[i];
			A[i][1] = A_S[i];
			A[i][2] = A_P[i];
			A[i][3] = A_N[i];
			A[i][4] = A_E[i];
	}
	
	return 1;
}

double conjugate_gradient_method(int N, int Nj, double A[N+1][5], double T_initial[N+1], double T[N+1], double b[N+1], double epsilon, double residuals[1000] ) {
	for(int i=1;i<=N;i++) {
		T[i] = T_initial[i]; // Initializing temperature at all points
	}
	double C[N+1];
	double d[N+1];
	double r[N+1], r_old[N+1];
	double alpha, beta;
	double L2_norm_residual;
	double error = 1;
	int iteration_count=0;
	// step 1
	matrix_Vector_multiplication( N, Nj, A, T, C);
	for(int i=1;i<=N;i++) {
		r[i] = b[i] - C[i];
		r_old[i] = b[i] - C[i];
		d[i] = b[i] - C[i];
	}
	
	while /*(iteration_count<=0) { */ (error>epsilon) {
		// step 2
		matrix_Vector_multiplication( N, Nj, A, d, C);
		alpha = vector_vector_multiplication(N,r_old,r_old) / vector_vector_multiplication(N,d,C); 
		// printf("Alpha = %lf\n",alpha); // Testing
		
		// step 3
		for(int i=1;i<=N;i++) {
			T[i] = T[i] + alpha*d[i];
		}
		
		// step 4
		if (iteration_count%50 == 0) { // To avoid floating point Iteration
			matrix_Vector_multiplication( N, Nj, A, T, C);
			for(int i=1;i<=N;i++) {
				r[i] = b[i] - C[i];
			}
		}
		else {
			matrix_Vector_multiplication( N, Nj, A, d, C);
			for(int i=1;i<=N;i++) {
				r[i] = r_old[i] - alpha*C[i];
			}
		}
		
		// step 5
		beta = vector_vector_multiplication(N,r,r) / vector_vector_multiplication(N,r_old,r_old);
		
		// step 6
		for(int i=1;i<=N;i++) {
			d[i] = r[i] + beta*d[i];
		}
		
		// Error Calculation
		L2_norm_residual = 0.0;
		for(int i=1;i<=N;i++) {
			L2_norm_residual = L2_norm_residual + pow(r[i],2);
		}
		L2_norm_residual = pow(L2_norm_residual,0.5);
		
		error = L2_norm_residual;
		
		iteration_count++;
		residuals[iteration_count] = L2_norm_residual;
		
		// Updating old residue
		for(int i=1;i<=N;i++) {
			r_old[i] = r[i]	;
		}
		
		// Testing
		if (iteration_count%20==0) {
			printf("Iteration: %d -> Error/Residual - L2_norm = %lf \n",iteration_count,L2_norm_residual);
		}
	}
	
	return iteration_count;
}

double pre_conditioner_conjugate_gradient_method(int N, int Nj, double A[N+1][5], double T_initial[N+1], double T[N+1], double b[N+1], double epsilon, char Pre_Conditioner_type, double residuals[1000] ) {
	for(int i=1;i<=N;i++) {
		T[i] = T_initial[i]; // Initializing temperature at all points
	}
	double C[N+1], G[N+1], H[N+1]; // Temporary Vectors for extracting vector reults from functions
	double d[N+1];
	double r[N+1], r_old[N+1];
	double alpha, beta;
	double L2_norm_residual;
	double error = 1;
	int iteration_count=0;
	
	double L[N+1][3], U[N+1][3]; 
	if (Pre_Conditioner_type == 'I') { 
		ILU_matrix( N, Nj, A, L, U ); // L & U are Lower and Upper triangular matrix of pre_conditioner_matrix(P)
	}
	if (Pre_Conditioner_type == 'S') { 
		SIP_matrix( N, Nj, 0.5, A, L, U ); // L & U are Lower and Upper triangular matrix of pre_conditioner_matrix(P)
	}
	
	// step 1
	matrix_Vector_multiplication( N, Nj, A, T, C);
	for(int i=1;i<=N;i++) {
		r[i] = b[i] - C[i];
		r_old[i] = b[i] - C[i];
	}
	// print_1D_Array( N, r_old); // Testing
	
	// Caling Pre_Conditioner
	if (Pre_Conditioner_type == 'J') { 
		Jacobi_pre_conditioner( N, Nj, A, r_old, d ); // d = inv(P)*r_old   &   P=inv(A_P)
	}
	
	if (Pre_Conditioner_type == 'I'  || Pre_Conditioner_type == 'S') { 
		general_LU_pre_conditioner( N, Nj, L, U, r_old, d); // d = inv(P)*r_old   &   P=L*U
	}
	// print_1D_Array( N, d); // Testing
	
	while  /*(iteration_count<=0) {  */ (error>epsilon) {
	
		if (Pre_Conditioner_type == 'J') {
			Jacobi_pre_conditioner( N, Nj, A, r_old, G ); // G = inv(P)*r_old   &   P=inv(A_P)
		}
		if (Pre_Conditioner_type == 'I' || Pre_Conditioner_type == 'S') { 
			general_LU_pre_conditioner( N, Nj, L, U, r_old, G ); // G = inv(P)*r_old & P=L*U
		}
		// print_1D_Array( N, G); // Testing		
		
		// step 2
		matrix_Vector_multiplication( N, Nj, A, d, C);
		alpha = vector_vector_multiplication(N,r_old,G) / vector_vector_multiplication(N,d,C); 
		// Testing
		// printf("Numerator = %lf\n",vector_vector_multiplication(N,r_old,G));
		// printf("Denominator = %lf\n",vector_vector_multiplication(N,d,C));
		// printf("Alpha = %lf\n",alpha); 
		
		// step 3
		for(int i=1;i<=N;i++) {
			T[i] = T[i] + alpha*d[i];
		}
		//print_1D_Array( N, T); // Testing
		
		// step 4
		//printf("step 4 \n");
		if (iteration_count%50 == 0) { // To avoid floating point Iteration error due to approximation
			matrix_Vector_multiplication( N, Nj, A, T, C);
			for(int i=1;i<=N;i++) {
				r[i] = b[i] - C[i];
			}
		}
		else {
			matrix_Vector_multiplication( N, Nj, A, d, C);
			for(int i=1;i<=N;i++) {
				r[i] = r_old[i] - alpha*C[i];
			}
		}
		// print_1D_Array( N, r); // Testing
		
		if (Pre_Conditioner_type == 'J') {
			Jacobi_pre_conditioner( N, Nj, A, r, H ); // H = inv(P)*r   &   P=inv(A_P)
		}
		if (Pre_Conditioner_type == 'I' || Pre_Conditioner_type == 'S') { 
			general_LU_pre_conditioner( N, Nj, L, U, r, H ); // H = inv(P)*r   &   P=L*U
		}
		// print_1D_Array( N, H); // Testing
		
		// step 5
		//printf("step 5 \n");
		beta = vector_vector_multiplication(N,r,H) / vector_vector_multiplication(N,r_old,G);
		
		// step 6
		for(int i=1;i<=N;i++) {
			d[i] = H[i] + beta*d[i];
		}
		
		// Error Calculation
		L2_norm_residual = 0.0;
		for(int i=1;i<=N;i++) {
			L2_norm_residual = L2_norm_residual + pow(r[i],2);
		}
		L2_norm_residual = pow(L2_norm_residual,0.5);
		
		error = L2_norm_residual;
		
		iteration_count++;
		residuals[iteration_count] = L2_norm_residual;
		
		// Updating old residue
		for(int i=1;i<=N;i++) {
			r_old[i] = r[i]	;
		}
		
		// Testing
		if (iteration_count%10 == 0) {
			printf("Iteration: %d -> Error/Residual - L2_norm = %lf \n",iteration_count,L2_norm_residual);
		}
	}
	
	return iteration_count;
}

double Jacobi_pre_conditioner( int N, int Nj, double A[N+1][5], double Q_0[N+1], double F_0[N+1] ) {
	for(int n=1;n<=N;n++) {
		F_0[n] = (1.0/A[n][2])*Q_0[n];
	}
	return 1;
}

double general_LU_pre_conditioner( int N, int Nj, double L[N+1][3], double U[N+1][3], double Q[N+1], double F[N+1] ) {
	// p_inv * r = x >>> P*x = r
	// P = M = L*U >>> L*U*x = r
	// U*x = y --> L*y = r >>> Get y >>> U*x = y --> Get x 
	// Q = r , F = x 
	double x[N+1], y[N+1];
	
		// Forward Substitution --> Finding y
		y[1] = Q[1]/L[1][2];
		for(int n=2;n<=N;n++) {
			y[n] = Q[n];
			if (n>Nj) { y[n] =  y[n]  -  L[n][0]*y[n-Nj]; }
			if (n>1) { y[n] =  y[n]  -  L[n][1]*y[n-1]; }
			y[n] = y[n] / L[n][2];
		}
				
		// Backward Substitution --> Finding x
		x[N] = y[N]/U[N][0];
		for(int n=N-1;n>=1;n--) {
			x[n] = y[n];
			if (n<N) { x[n] = x[n] - U[n][1]*x[n+1]; }
			if (n<=N-Nj) { x[n] = x[n] - U[n][2]*x[n+Nj]; }
			x[n] = x[n] / U[n][0];
		}
		for(int i=1;i<=N;i++) {
			F[i] = x[i];
		}
			
	return 1;
}

double ILU_matrix( int N, int Nj, double A[N+1][5], double L[N+1][3], double U[N+1][3] ) {
	double L_P[N+1],L_S[N+1], L_W[N+1], U_N[N+1], U_E[N+1];
	// Initializing all Columns vectors to zero
	for(int i=0;i<=N;i++) {
		L_P[i] = 0; // Point
		L_W[i] = 0; // West Neighbour
		L_S[i] = 0; // South Neighbour
		U_N[i] = 0; // North Neighbour
		U_E[i] = 0; // East Neighbour
	}
	for(int n=1;n<=N;n++) {
		// L_W
		if (n>Nj) {
			L_W[n] = A[n][0]; // A[n][0] = A_W[n]
		}
		// L_S
		if (n>1) {
			L_S[n] = A[n][1]; // A[n][1] = A_S[n]
		}
		// L_P
		if( n == 1 ) {
			L_P[n] = A[n][2]; // A[n][2] = A_P[n]
		}
		else if( n>1 && n<=Nj ) {
			L_P[n] = A[n][2] - L_S[n]*U_N[n-1]; 
		}
		else {
			L_P[n] = A[n][2] - L_S[n]*U_N[n-1] - L_W[n]*U_E[n-Nj]; 
		}
		// U_N
		if (n>=1 && n<=N-1) {
			U_N[n] = A[n][3] / L_P[n]; // A[n][3] = A_N[n]
		}
		// U_E
		if (n>=1 && n<=N-Nj) {
			U_E[n] = A[n][4] / L_P[n]; // A[n][4] = A_E[n]
		}
		
	}
// For Understanding Decomposition of M as L & U
// double M[N+1][5]; // M[n][0] = L_W[n], M[n][1] = L_S[n], M[n][2] = L_P[n], M[n][3] = U_N[n], M[n][4] = U_E[n]
// double L[N+1][3], U[N+1][3];
// L[n][0] = L_W[n], L[n][1] = L_S[n], L[n][2] = L_P[n], 
// U[n][0] = 1, U[n][1] = U_N[n], U[n][2] = U_E[n]

	for(int n=1;n<=N;n++) {
		L[n][0] = L_W[n];
		L[n][1] = L_S[n];
		L[n][2] = L_P[n];
		U[n][0] = 1.0;
		U[n][1] = U_N[n];
		U[n][2] = U_E[n];
		// Testing
		// printf("n=%d --> L_W = %lf , L_S = %lf , L_P = %lf , U_P = %lf , U_N = %lf , U_E = %lf  \n",n,L[n][0],L[n][1],L[n][2],U[n][0],U[n][1],U[n][2]);
		
	}
	return 1;
}

double SIP_matrix( int N, int Nj, double alpha, double A[N+1][5], double L[N+1][3], double U[N+1][3] ) {
	double L_P[N+1],L_S[N+1], L_W[N+1], U_N[N+1], U_E[N+1];
	// Initializing all Columns vectors to zero
	for(int i=0;i<=N;i++) {
		L_P[i] = 0; // Point
		L_W[i] = 0; // West Neighbour
		L_S[i] = 0; // South Neighbour
		U_N[i] = 0; // North Neighbour
		U_E[i] = 0; // East Neighbour
	}
	for(int n=1;n<=N;n++) {
		// L_W
		if (n>Nj) {
			L_W[n] = A[n][0] / ( 1.0 + alpha*U_N[n-Nj] ) ; // A[n][0] = A_W[n]
		}
		// L_S
		if (n>1) {
			L_S[n] = A[n][1] / ( 1.0 + alpha*U_E[n-1] ); // A[n][1] = A_S[n]
		}
		// L_P
		if( n == 1 ) {
			L_P[n] = A[n][2]; // A[n][2] = A_P[n]
		}
		else if( n>1 && n<=Nj ) {
			L_P[n] = A[n][2] - L_S[n]*( U_N[n-1] - alpha*U_E[n-1] ); 
		}
		else {
			L_P[n] = A[n][2] - L_S[n]*( U_N[n-1] - alpha*U_E[n-1])  - L_W[n]*( U_E[n-Nj]  -  alpha*U_N[n-Nj] ); 
		}
		// U_N
		if (n>=1 && n<=Nj) {
			U_N[n] = ( A[n][3] ) / L_P[n]; 
		}
		else if (n>=Nj && n<=N-1) {
			U_N[n] = ( A[n][3] - alpha*L_W[n]*U_N[n-Nj] ) / L_P[n]; // A[n][3] = A_N[n]
		}
		// U_E
		if (n>=1 && n<=N-Nj) {
			U_E[n] = ( A[n][4] - alpha*L_S[n]*U_E[n-1] ) / L_P[n]; // A[n][4] = A_E[n]
		}
	}
// double M[N+1][5]; // M[n][0] = L_W[n], M[n][1] = L_S[n], M[n][2] = L_P[n], M[n][3] = U_N[n], M[n][4] = U_E[n]
// double L[N+1][3], U[N+1][3];
// L[n][0] = L_W[n], L[n][1] = L_S[n], L[n][2] = L_P[n], 
// U[n][0] = 1, U[n][1] = U_N[n], U[n][2] = U_E[n]
	for(int n=1;n<=N;n++) {
		L[n][0] = L_W[n];
		L[n][1] = L_S[n];
		L[n][2] = L_P[n];
		U[n][0] = 1.0;
		U[n][1] = U_N[n];
		U[n][2] = U_E[n];
		// Testing
		//printf("n=%d --> L_W = %lf , L_S = %lf , L_P = %lf , U_P = %lf , U_N = %lf , U_E = %lf  \n",n,L[n][0],L[n][1],L[n][2],U[n][0],U[n][1],U[n][2]);
		
	}
	return 1;
}

double matrix_Vector_multiplication(int N, int Nj, double A[N+1][5], double P[N+1], double C[N+1]) {
	for(int n=1;n<=N;n++) {
		C[n]=A[n][2]*P[n];
		if(n>Nj)    { C[n] = C[n] + A[n][0]*P[n-Nj]; }
		if(n<=N-Nj) { C[n] = C[n] + A[n][4]*P[n+Nj]; }  
		if(n<=N-1)  { C[n] = C[n] + A[n][3]*P[n+1]; }  
		if(n>1)     { C[n] = C[n] + A[n][1]*P[n-1]; } 
	}
	return 1;
}

double vector_vector_multiplication(int N, double V[N+1], double W[N+1]) {
	double result=0.0;
	for(int i=1;i<=N;i++) {
		result = result + V[i]*W[i];
	}
	return result;
}

void print_1D_Array(int N, double k[N+1]) {
	for(int i=1;i<=N;i++) {
		printf("k[%d] =  %lf \t",i,k[i]);
	}
	printf("\n");
	printf("-----------------------------------\n");
}
