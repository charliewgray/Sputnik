// This is start of the header guard.  ELEMENTCONVERTERS_H can be any unique name.  By convention, we use the name of the header file.
#ifndef DENSITY_H
#define DENSITY_H

// This is the content of the .h file, which is where the declarations go
void findSolarFlux(int year, int month, int Atmosphere_Select);
double get_density(double day, double alt, double g_lat, double g_long, double f107A, double f107, double ap);

// This is the end of the header guard
#endif