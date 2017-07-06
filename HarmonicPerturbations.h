// This is start of the header guard.  DRAGFUNCTIONS_H can be any unique name.  By convention, we use the name of the header file.
#ifndef HARMONICPERTURBATIONS_H
#define HARMONICPERTURBATIONS_H

// This is the content of the .h file, which is where the declarations go
void EarthHarmonics(double a[3], double Rx, double Ry, double Rz, double Rmag, double Lat, double Long);
double Legendre(int n, int m, double x);
void PotentialPartials(double *dUdr, double *dUdphi, double *dUdlambda, double Req, double Rmag, double phi, double lambda);
void EarthHarmonicsDetailed(double a[3], double Rx, double Ry, double Rz, double Rmag, double Lat, double Long);

// This is the end of the header guard
#endif