//
//  BloodVessel.cpp
//  Sim3D-v3
//
//  Created by Edoardo Milotti on 25/7/2017.
//
//

#include "sim.h"
#include "BloodVessel.h"
using namespace vbl;
// costruttori
// costruttore di default (nessun attributo specificato)
BloodVessel::BloodVessel()
{
    //a=.assign(3,0.); // std::array is zero initialized by default
    //b.assign(3,0.);
    //va.assign(3,0.);
    //vb.assign(3,0.);
    R=0.;
    vR=0.;
    O2start=0.;
    O2end=0.;
    CO2start=0.;
    CO2end=0.;
    G=0.;
    A=0.;
    AcL=0.;

}

// costruttore con valori specificati
BloodVessel::BloodVessel( const std::array<double,3> ca, const std::array<double,3> cb, const std::array<double,3> cva, const std::array<double,3> cvb, const double cR, const double cvR, const double cO2start, const double cO2end, const double cCO2start, const double cCO2end, const double cG, const double cA, const double cAcL )
{
    a = ca;
    b = cb;
    va = cva;
    vb = cvb;
    R = cR;
    vR = cvR;
    O2start = cO2start;
    O2end = cO2end;
    CO2start = cCO2start;
    CO2end = cCO2end;
    G = cG;
    A = cA;
    AcL = cAcL;

}

// costruttore copia
BloodVessel::BloodVessel(const BloodVessel& bv)
{
    a = bv.a;
    b = bv.b;
    va = bv.va;
    vb = bv.vb;
    R = bv.R;
    vR = bv.vR;
    O2start = bv.O2start;
    O2end = bv.O2end;
    CO2start = bv.CO2start;
    CO2end = bv.CO2end;
    G = bv.G;
    A = bv.A;
    AcL = bv.AcL;
}

// funzione distanza e posizione punto piu' vicino
double BloodVessel::DistanceFromVessel( const std::array<double,3> x1, double* x0 )
{
    double tp = 0;
    double nrm = 0;
    double bvd = 0;
    
    for(int i=0; i<3; i++) // calcolo di t0
    {
        tp += (x1[i]-a[i])*(b[i]-a[i]);
        nrm += (b[i]-a[i])*(b[i]-a[i]);
    }
    tp /= nrm;
    
    // cout << "nrm: " << nrm << endl;
    // cout << "tp: " << tp << endl;
    
    
    if(tp >= 0 && tp <= 1) // calcolo della distanza dalla parte cilindrica
        for(int k=0; k<3; k++)
        {
            x0[k] = tp*(b[k]-a[k]) + a[k];
            bvd += (x1[k] - x0[k])*(x1[k] - x0[k]);
        }
    else if(tp<0)   // calcolo della distanza dalla calotta sferica all'inizio
        for(int k=0; k<3; k++)
        {
            x0[k] = a[k];
            bvd += ( x1[k]-a[k])*( x1[k]-a[k]);
        }
    else if(tp>1)   // calcolo della distanza dalla calotta sferica alla fine
        for(int k=0; k<3; k++)
        {
            x0[k] = b[k];
            bvd += ( x1[k]-b[k] )*( x1[k]-b[k] );
        }
    
    bvd = sqrt(bvd);
    
    return(bvd);
    
}
double BloodVessel::DistanceFromVessel( const ANNpoint &x1)
{
    
    double tp = 0;
    double nrm = 0;
    double bvd = 0;
    //std::array<double,3> x0 = {0.0,0.0,0.0};
    double x0 [3] = {0.0,0.0,0.0};
    
    for(int i=0; i<3; i++) // calcolo di t0
    {
        tp += (x1[i]-a[i])*(b[i]-a[i]);
        nrm += (b[i]-a[i])*(b[i]-a[i]);
    }
    tp /= nrm;
    
    // cout << "nrm: " << nrm << endl;
    // cout << "tp: " << tp << endl;
    
    
    if(tp >= 0 && tp <= 1) // Calculating the distance from the cylindrical part
        for(int k=0; k<3; k++)
        {
            x0[k] = tp*(b[k]-a[k]) + a[k];
            bvd += (x1[k] - x0[k])*(x1[k] - x0[k]);
        }
    else if(tp<0)   // Calculating the distance from the spherical shell at the beginning
        for(int k=0; k<3; k++)
        {
            x0[k] = a[k];
            bvd += ( x1[k]-a[k])*( x1[k]-a[k]);
        }
    else if(tp>1)   // Calculating the distance from the spherical shell to the end
        for(int k=0; k<3; k++)
        {
            x0[k] = b[k];
            bvd += ( x1[k]-b[k] )*( x1[k]-b[k] );
        }
    
    bvd = sqrt(bvd);
    
    return(bvd);
    
}
// namespace vbl{
// BloodVessel& operator=(const BloodVessel& bv)
// {
//   std::cout<< "here error" <<std::endl;
// }
// }
// BloodVessel& operator=(const BloodVessel& bv)
// {
//   std::cout<< "here error" <<std::endl;
// }

