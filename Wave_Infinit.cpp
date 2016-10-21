//
//  main.cpp
//  PDE_WAVE
//
//  Created by Usama Anber on 1/27/16.
//  Copyright (c) 2016 Usama Anber. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

//IC 1: u0(0,x)=sin(2*pi*x)+0.25*sin(10*pi*x);
vector<double>  PDE_wave_infstring_U0(vector<double>& x){
    
    size_t    N= x.size(); //size of vector x
    
    vector<double> u0 (N+1);
    //  u0.resize(N+1);
    //   double u0;
    for (int i=0; i<=N; i++) {
        u0[i]=sin(2*M_PI*x[i])+0.25*sin(10*M_PI*x[i]);
    }
    return u0;
    
}


//IC 2: u1(0,x)=0;

vector<double> PDE_wave_infstring_U1(vector<double>& x){
    
size_t    N= x.size();
    
    vector<double> u1(N+1);
    
    for (int i=0; i<=N; i++) {
        u1[i]=0;
    }
    return u1;
}



//BCs are periodic (assumes same as U) and U1 at first and last element) 

int main()
{
// salt and papper
    double a=0;
    double b=1;
    double c=2; //wave speed
    int nx=50; // number of grid points
    double time=0;
    int Maxtime=50; //time integration
    double dx= (b-a)/nx;
    double dt=0.0100;
    double sigma=c*dt/dx;
	double sigma2=sigma*sigma;
    double coeff= 2*(1-sigma2);
    
    vector<double> xx (nx+1);
// fill in vector xx:
    for (int i=0; i<=nx; i++) {
        xx[i]=a+i*dx;
    }
// Now apply initial conditions:
    vector<double>  U0(nx+1);
    vector<double>  U1(nx+1);
    
    U0=PDE_wave_infstring_U0(xx);
    U1=PDE_wave_infstring_U1(xx);
    
    
//    U1=PDE_wave_infstring_U1(xx);
//    U0=PDE_wave_infstring_U0(xx);    

    
    // open file and set the format:
    std::ofstream write_file("wave_infinite.dat");
    // Write numbers as +x.<13digits>e+00 (width 20)
    write_file.setf(std::ios::scientific);
    write_file.setf(std::ios::showpos);
    write_file.precision(13);
    //
// first row contains xx
    for (int i=0; i<=nx; i++){
        write_file <<xx[i]<<" ";
    }
        write_file <<endl;

  ////////////////////////////////////////////////    
 // Now we integrate in time to get solution U2//
////////////////////////////////////////////////    

    vector<double>  U2(nx+1);
    
    for (int n=2; n<=2*Maxtime; n++) {
        time=time+dt;
        
        //solve for interior of the domain:
        
        for (int j=1; j<nx; j++) {
            U2[j]=-U0[j]+coeff*U1[j]+sigma*(U1[j-1]+U1[j+1]);
        }
        // impose boundary condtions at one boundary only, then use periodicity:
        
        U2[0]=U0[0]+dt*U1[0];
        U2[nx]=U2[0];   // periodic BCs!
        
    
             // output results every 10 time steps
                if(n % 10==0){
                    for (int i=0; i<=nx; i++) {
                        write_file <<U2[i]<<" ";
                   }
                      write_file <<endl;
                }
        
        // update solutions:

        U0=U1;
        U1=U2;
    
    }
    
    write_file.close();
    
    //de-allocate memory
    U0 = vector<double>();
    U1 = vector<double>();
    U2 = vector<double>();
    
    
    return 0;
}

//void PDE_wave_infstring_U0(vector<double>& x, vector<double>& u0){
//
//    size_t    N= x.size(); //size of vector x
//
//    for (int i=0; i<=N; i++) {
//        u0[i]=sin(2*M_PI*x[i])+0.25*sin(10*M_PI*x[i]);
//    }
//
//
//}
//
//
////IC 2: u1(0,x)=0;
//
//void PDE_wave_infstring_U1(vector<double>& x, vector<double>& u1){
//
//    size_t    N= x.size();
//
//    for (int i=0; i<=N; i++) {
//        u1[i]=0;
//    }
//   }

