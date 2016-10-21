#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

// Exact solution:
vector<double> PDE_heat_exact(vector<double>& x, double t,\
                                      vector<double>& amplitude, vector<double>& wavenumber){
    
    size_t    N= x.size(); //size of vector x
    
    vector<double> y (N+1,0); // vector of length N+1 of zeros
    
    for (int k=0; k<wavenumber.size(); k++) {
        //fill in vector y
        for (int i=0; i<=N; i++) {
            y[i]+=amplitude[k]*exp(-pow(wavenumber[k]*M_PI,2)*t)*sin(wavenumber[k]*M_PI*x[i]);
        }
    }
    
    return y;
    
}

//IC for infinite case (at the entire domain including the boundary)
vector<double>  PDE_heate_u0(vector<double>& x, vector<double>& amplitude, vector<double>& wavenumber){
    
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


int main()
{
    
    double a=0;
    double b=1;
    int nx=50; // number of grid points
    double dx= (b-a)/double(nx);
    
    double sigma=0.5;
    double dt=sigma*dx*dx;
    double coeff= 1-2*sigma;
    double time=0;
    int Maxtime=1200; //time integration

    
    vector<double> xx (nx+1); //note n+1 is the length of the vector xx: xx[0],xx[1],..xx[nx]
    // fill in vector xx:
    for (int i=0; i<=nx; i++) {
        xx[i]=a+double(i)*dx;
        //     cout<<xx[i]<<endl;
    }

    
    //initialize two vectors for amplitude and wavenumber this way.
    double aa[]={1,0.25};
    vector<double> amp (&aa[0], &aa[0]+2);
    
    double kk[]={1,10};
    vector<double> wn (&kk[0], &kk[0]+2);
    
    
    // Now apply initial conditions depending on the case:
    // if key==0, then finite length rod where:
    //u[0]=1 and zeros everywhere.
    //and if key==1, then specific BCs defined by
    //function PDE_heate_u0
    
    vector<double>  U0(nx+1);
    
/////// key for finite =0 or infinit =1 case
    int key=1;
    
// if (key==0)        
//        for (int i=1; i<=nx; i++) {
//            U0[i]=0;           
//                    }
//        U0[0]=1;
        
        // open file and set the format:
//  	        std::ofstream write_file("heat_finite.dat");
        // Write numbers as +x.<13digits>e+00 (width 20)
// 	        write_file.setf(std::ios::scientific);
//            write_file.setf(std::ios::showpos);
//            write_file.precision(13);
        //first row contains xx
//           for (int i=0; i<=nx; i++){
//               write_file <<xx[i]<<" ";
//           }
//           write_file <<endl;
    
    
//    if (key==1) {
        U0=PDE_heate_u0(xx,amp,wn); //IC
        
        // open file and set the format:
        std::ofstream write_file("heat_infinite.dat");
        // Write numbers as +x.<13digits>e+00 (width 20)
        write_file.setf(std::ios::scientific);
        write_file.setf(std::ios::showpos);
        write_file.precision(13);
        //first row contains xx
        for (int i=0; i<=nx; i++){
            write_file <<xx[i]<<" ";
        }
        write_file <<endl;
//    }

    
    /////////////////////////////////////////////////
    // Now we integrate in time to get solution U1 //
    ////////////////////////////////////////////////
    
    vector<double>  U1(nx+1);
    
    U1=U0; //
    
    for (int n=1; n<=Maxtime; n++) {
        time=time+dt;
        
        //solve for interior of the domain:
        
        for (int j=1; j<nx; j++) {
            U1[j]=coeff*U0[j]+sigma*(U0[j-1]+U0[j+1]);
        }
        
        //periodic boundary conditions for the case of infinte domain
        if (key==1) {
            U1[nx]=U0[nx];
            U1[0]=U1[nx];
        }
        
        // output results every 300 time steps
                if(n % 300==0){
                    for (int i=0; i<=nx; i++) {
                        write_file <<U1[i]<<" ";
                        }
                        write_file <<endl;
                    }
        
        
        // update solutions:
        U0=U1;
        
    }
    
        write_file.close();
   
	//exact solution
	vector<double> UE=PDE_heat_exact(xx,time,amp,wn);

   //de-allocate memory
    U0 = vector<double>();
    U1 = vector<double>();
    UE = vector<double>();


    return 0;
}


