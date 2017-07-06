// This is start of the header guard.  ELEMENTCONVERTERS_H can be any unique name.  By convention, we use the name of the header file.
#ifndef LATLONGBETA_H
#define LATLONGBETA_H

// This is the content of the .h file, which is where the declarations go
double Range_2pi(double Angle);
double JulianDate(int year, int month, double day);
double Geo_Alt(double Latitude);
void Longitude_Latitude(int year, int month, double day, double x[6], double *Longitude, double *Latitude);
void SunVectorCalc(double JD);
double BetaCalcs(double RAAN, double I, double e, double JD);

// This is the end of the header guard
#endif