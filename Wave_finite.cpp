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

// Exact solution:
vector<double> PDE_wave_fstring_exact(vector<double>& x, double t, double c,\
                                      vector<double>& amplitude, vector<double>& wavenumber){
    
    size_t    N= x.size(); //size of vector x
    
    vector<double> y (N+1,0); // vector of length N+1 of zeros
    
    for (int k=0; k<wavenumber.size(); k++) {
        //fill in vector y
        for (int i=0; i<=N; i++) {
            y[i]+=amplitude[k]*cos(wavenumber[k]*M_PI*c*t)*sin(wavenumber[k]*M_PI*x[i]);
        }
    }
    
    return y;   
}

//IC 1: u0(0,x)=sin(1*pi*x)+0.25*sin(10*pi*x);
vector<double>  PDE_wave_fstring_in(vector<double>& x, vector<double>& amplitude, vector<double>& wavenumber){
    
    size_t    N= x.size(); //size of vector x
    
    vector<double> y (N+1,0); // vector of length N+1 of zeros
   
    for (int k=0; k<wavenumber.size(); k++) {
        //fill in vector y
        for (int i=0; i<=N; i++) {
            y[i]+=amplitude[k]*sin(wavenumber[k]*M_PI*x[i]);
        }
    }
    
    return y;    
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
    int Maxtime=125; //time integration
    double dx= (b-a)/nx;
    double dt= 0.0080;
    double sigma=c*dt/dx;
    double sigma2=sigma*sigma;
    double coeff= 2*(1-sigma2);
    
    vector<double> xx (nx);
    // fill in vector xx:
    for (int i=0; i<=nx; i++) {
        xx[i]=a+double(i)*dx;
    }

    //declear amplitude and wavenumber vectors:
//    vector<double> amp [2]= {1.,0.25}; //amplitude
//    vector<double> k [2]= {1.,10.};  //wavenumber
   
    double aa[]={1,0.25};
    vector<double> amp (&aa[0], &aa[0]+2);
    
    double kk[]={1,10};
    vector<double> wn (&kk[0], &kk[0]+2);

 
    // Now apply initial conditions:
    vector<double>  U0(nx);
    vector<double>  U1(nx);
    
    U0=PDE_wave_fstring_in(xx,amp,wn); //IC1
    U1=U0; //IC2
    
    
    
    // open file and set the format:
    std::ofstream write_file("wave_finite_string.dat");
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
    
     /////////////////////////////////////////////////
    // Now we integrate in time to get solution U2 //
    ////////////////////////////////////////////////
    
    vector<double>  U2(nx);
    
    U2=U1;

    
    for (int n=2; n<=Maxtime; n++) {
        time=time+dt;
        
        //solve for interior of the domain:
        
        for (int j=1; j<nx; j++) {
            U2[j]=-U0[j]+coeff*U1[j]+sigma2*(U1[j-1]+U1[j+1]);
        }

        // output results every nt time steps
        if(n % Maxtime==0){
            for (int i=0; i<=nx; i++) {
                write_file <<U2[i]<<" ";
                }
                write_file <<endl;
            }
    
        // update solutions:
        
        U0=U1;
        U1=U2;
        
    }
   
	// Exact solution for comparison:
	vector<double>  UE=PDE_wave_fstring_exact(xx,time,c,amp,wn);
 
	//write exact solution on last row
	for (int i=0; i<=nx; i++){
		write_file <<UE[i]<<" ";
	}
	write_file <<endl;

    write_file.close();
     
    //de-allocate memory
    U0 = vector<double>();
    U1 = vector<double>();
    U2 = vector<double>();
 	UE = vector<double>();   
    
    return 0;
}
