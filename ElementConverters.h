// This is start of the header guard.  ELEMENTCONVERTERS_H can be any unique name.  By convention, we use the name of the header file.
#ifndef ELEMENTCONVERTERS_H
#define ELEMENTCONVERTERS_H

// This is the content of the .h file, which is where the declarations go
double DegToRad(double Var);
double RadToDeg(double Var);
double PiBound(double value);
void Elements_To_Vector(double A, double e, double I, double Ra, double Wp, double TA, double x[6]);
void State_Vector_To_Elements(double x[6], double *A, double *e, double *I, double *Ra, double *Wp, double *TA);
double *eci_to_ecf(double Rx, double Ry, double Rz);
double *ecf_to_eci(double Rx, double Ry, double Rz);
double Calculate_Invariant_Elements(double x[6]);

// This is the end of the header guard
#endif