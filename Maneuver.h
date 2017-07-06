// This is start of the header guard.  ELEMENTCONVERTERS_H can be any unique name.  By convention, we use the name of the header file.
#ifndef MANEUVER_H
#define MANEUVER_H

// This is the content of the .h file, which is where the declarations go
void ConvertAndApplydV(double dt, double A, double e, double I, double Ra, double Wp, double TA, double x[6], double dV);

// This is the end of the header guard
#endif