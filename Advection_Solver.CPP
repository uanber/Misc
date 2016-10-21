#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// initial conditions
void pde_advect_IC(vector<double>& x, vector<double>& u){
    
   size_t N = x.size(); //size of vector x
    for (int i=0; i<=N; i++) {
        u[i]=x[i];
    }
}

// boundary conditions
double pde_advect_BC(double t){
    
    double phi;
    
    return phi=sin(10*t*M_PI);    
    
}


// declaration and integration in time
int main()
{
    double a=0;
    double b=1;
    int    nx=40;
    double dx=(b-a)/nx;
//    double xx[nx+1]; //array xx with intervals 
	vector<double> xx (nx+1);
    double sigma=0.8;        // CFL condition: sigma<1
    double c=4;
    double dt=sigma*dx/c;
    double time=0;
    int    Maxtime=50;
    
// open file and set the format:
    std::ofstream write_file("advection.dat");
     // Write numbers as +x.<13digits>e+00 (width 20)
     write_file.setf(std::ios::scientific);
     write_file.setf(std::ios::showpos);
     write_file.precision(13);
 //
    
    // allocate memory for vectors of solutions u0 and u1
   // double* u0 = new double [nx+1];
   // double* u1 = new double [nx+1];
    vector<double> u0 (nx+1);
	vector<double> u1 (nx+1);
    
    //fill in array x
    for (int i=0; i<=nx;i++){
        xx[i]=a+i*dx;
    }
    
	for (int i=0; i<=nx; i++) {
         write_file <<xx[i]<<"  ";
              }
   write_file<<endl;
 
    pde_advect_IC(xx,u0); // u0=x (initial conditions)
   // copy(u0, u0+(nx+1), u1);// u1=u0;
    u1=u0;
    
	 //////////////////////////////////////
	// now perform integration over time// 
   //////////////////////////////////////

    for (int n=1; n<=Maxtime; n++) {
        time=time+dt;
        
       // fill in the x at j 
        for (int j=nx+1; j>=1; j--) {
            u1[j]= (1-sigma)*u0[j]+sigma*u0[j-1];
        }
        u1[0]= pde_advect_BC(time); //impose boundary condition at first element
    
    // now update u0 with u1:
    //    for (int i=0; i<=nx; i++) {
    //        u0[i]=u1[i];
    //    }
       u0=u1;
	 // output results every 10 time steps
        if(n % 10==0){
            for (int i=0; i<=nx; i++) {
                write_file <<u1[i]<<"  ";
            }
           write_file <<endl; 
        } 
    }

write_file.close();
    
    //de-allocate memory 
    u0 = vector<double>();
	u1 = vector<double>();

    //delete[] u0;
    //delete[] u1;
    
    return 0;
}


//old stuff with arrays instead of vectors
/*
void pde_advect_IC(double* x, double* u){  
  // note * x is a pointer, with no actual size (8 or 4), to get
  // the size of the array stored in it, refer to it by &   
    int N = sizeof(& x) / sizeof(& x[0]); //size of vector u
    for (int i=0; i<=N; i++) {
         u[i]=x[i];
     }
 }
*/
