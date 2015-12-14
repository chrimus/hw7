#include <iostream>
#include <cmath>

using namespace std;

void ODE( double* In, double Out[][4], int ki, double mu);
double error (double* RK4, double* RK5);
void calc_pos(double* output, double* sol, const double* A, double k[][4], int N, double dt);

int main(){
	
//------------ The Butcher Table for Dormand-Prince --------------------------------------------------------------------------
	// the c-Values determining the evaluated time position t_eva = tn + c*dt NOT NECESSARY	
	// const double c[7] 	= {0.0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1.0, 1.0};
	
	// the aij-values determining the evaluated space position y_eva= yn + a[i] ki + a[i+1] k(i+1) + ...
	const double a0[1]	= {0};
	const double a1[1]	= {1.0/5};
	const double a2[2]	= {3.0/40, 		9.0/40};
	const double a3[3]	= {44.0/45, 	-56.0/15, 		32.0/9};
	const double a4[4]	= {19372.0/6561,-25360.0/2187, 	64448.0/6561, 	-212.0/729};
	const double a5[5]	= {9017.0/3168, -355.0/33, 		46732.0/5247, 	49.0/176, 	-5103.0/18656};
	const double a6[6]	= {35.0/384, 	0.0 , 			500.0/1113, 	125.0/192, 	-2187.0/6784, 	11.0/84};
	
	// the weights for the iterationstep once for RK4 and RK5
	const double bRK4[7] 	= {35.0/384, 	0.0 , 			500.0/1113, 	125.0/192, -2187.0/6784, 	11.0/84, 	0.0};
	const double bRK5[7]	= {5179.0/57600, 0.0, 			7571.0/16695, 	393.0/640, -92097.0/339200, 187.0/2100, 1.0/40};
	//-----------------------------------------------------------------------------------------------------------------------------
	
//-------------- USER INTERFACE -----------------------------------------------------
	const double mu = 0.012277471;
	const double tmin = 0.0;
	const double tmax = 22.0;
	const double p = 0.7;			//safety factor
	const double tol = pow(10,-5);
	//-------------------------------------------------------------------------------
	
	// 7 slopes need to be evaluated, each slope being a 4dim "vector"
	// ki [j] = k[i][j] and there are k0, k1, k2, ..., k6
	double k[7][4];
	
	// initial conditions for solution-vector {x, x', y, y'}
	double sol[4] = {0.99400000000, 0, 0, -2.00158510637908};
	double sol_RK5[4]= {0.99400000000, 0, 0, -2.00158510637908};
	
	double dt = 10e-6;
	double tn = tmin;
	double pos[4];
	
	//int N = 100;
	//int i = 0;
	while (tn < tmax){
		// print:	t	dt	x(t)	y(t)
		cout << tn << "\t" << dt << "\t" << sol[0] << "\t" << sol[2] <<"\n";
		
		// current time tn
		tn += dt;
		
		// Calculate new k vectors
		/* k0= */ calc_pos(pos, sol, a0, k, 0, dt); ODE( pos, k, 0, mu);
		/* k1= */ calc_pos(pos, sol, a1, k, 1, dt); ODE( pos, k, 1, mu);
		/* k2= */ calc_pos(pos, sol, a2, k, 2, dt); ODE( pos, k, 2, mu);
		/* k3= */ calc_pos(pos, sol, a3, k, 3, dt); ODE( pos, k, 3, mu);
		/* k4= */ calc_pos(pos, sol, a4, k, 4, dt); ODE( pos, k, 4, mu);
		/* k5= */ calc_pos(pos, sol, a5, k, 5, dt); ODE( pos, k, 5, mu);
		/* k6= */ calc_pos(pos, sol, a6, k, 6, dt); ODE( pos, k, 6, mu);
		
		// calculate new solution vectors sol and sol_RK5
		calc_pos(sol_RK5, 	sol, bRK5, k, 7, dt);
		calc_pos(sol, 		sol, bRK4, k, 7, dt);
		
		//Adapt Step width
		double err = error(sol, sol_RK5);		
		dt = p*dt* pow(tol/err, 0.2);
	}
	
	return 0;
}

double error (double* RK4, double* RK5){
	double max_err = 0;
	double err = 0;
	for ( int i=0; i<4; i++){
		err = abs(RK4[i]-RK5[i]);
		if ( err > max_err ) {max_err = err;}
	}
	return max_err;
}
//calculate sol + A[0]*k0*dt + A[1]*k1*dt + ... + A[N-1]*kN-1*dt, where sol is 4dim vector
//k[i][] contains the i.th slope (also 4dim vector)
void calc_pos(double* output, double* sol, const double* A, double k[][4], int N, double dt){
	
	for (int i=0; i<4; i++){
		
		output[i] = sol[i];
		for (int j=0; j<N; j++){
			output[i] += A[j]*k[j][i]*dt;
		}
	}
}

//In:  vector {x, x', y, y'} determining phase space position at which ODE should be evaluted
//Out: vector {x', x'', y', y''} calculated via ODE (hardcoded). 
//     Note that k[i][j] is the j.th component of ki
//ki: determines which k-vector will be overwritten
void ODE( double* In, double Out[][4], int ki, double mu){
	
	double r = sqrt( (In[0]+mu)*(In[0]+mu) + In[2]*In[2] );
	double s = sqrt( (In[0]-1+mu)*(In[0]-1+mu) + In[2]*In[2] );
	
	Out[ki][0] = In[1];
	Out[ki][1] = In[0] + 2*In[3] - (1.0-mu)*(In[0]+mu)/(r*r*r) - mu*(In[0]-1.0+mu)/(s*s*s);
	Out[ki][2] = In[3];
	Out[ki][3] = In[2] - 2*In[1] - (1.0-mu)*In[2]/(r*r*r) - mu*In[2]/(s*s*s);
}

