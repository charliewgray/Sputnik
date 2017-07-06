/* -------------------------------------------------------------------- */
/* ----------      CHARLIE GRAY DECAY CALCULATIONS 2015      ---------- */
/* -------------------------------------------------------------------- */

/* This file is part of the Decay Calc  C source code package
*
* The orbital decay model was developed by Charles Gray. The MSISE-00
* density model used was develoed by Mike Picone, Alan Hedin, and Doug Drob.
*/

/* ------------------------------------------------------------------- */
/* ------------------------------ INCLUDES --------------------------- */
/* ------------------------------------------------------------------- */
#include <math.h>          /* maths functions */
#include "LatLongBeta.h"

/* ------------------------------------------------------------------- */
/* ------------------------- Shared Variables ------------------------ */
/* ------------------------------------------------------------------- */
//I want to use these two global variables within this c file, so i intialize them here
double Pi; 
double Mu;
double JD;

/* ------------------------------------------------------------------- */
/* ------------------- Element Conversion Routines ------------------- */
/* ------------------------------------------------------------------- */
double DegToRad(double Var){
	return Var * Pi / 180.;
}

double RadToDeg(double Var){
	return Var * 180. / Pi;
}

double PiBound(double value){
	//add 2*Pi to values below zero
	while (value < 0){ value += 2 * Pi; }
	//subtract 2*Pi from values above 2*Pi
	while (value > 2 * Pi){ value -= 2 * Pi; }

	return value;
}

void Elements_To_Vector(double A, double e, double I, double Ra, double Wp, double TA, double x[6]){
	//Convert the input fields that are in degrees to radians
	I = DegToRad(I);
	Ra = DegToRad(Ra);
	Wp = DegToRad(Wp);
	TA = DegToRad(TA);

	/*Converts Keplerian Orbital elements to Cartesian Coordinates (x)*/
	double p = A * (1 - pow(e,2));
	double R = A * (1 - pow(e,2)) / (1 + e * cos(TA));

	//double Rpi = R * cos(TA);
	//double Rpj = R * sin(TA);
	//double Vpi = sqrt(Mu/p) * sin(TA);
	//double Vpj = sqrt(Mu/p) * (e + cos(TA));

	double T11 = cos(Ra) * cos(Wp) - sin(Ra) * cos(I) * sin(Wp);
	double T12 = -cos(Ra) * sin(Wp) - sin(Ra) * cos(I) * cos(Wp);
	//double T13 = sin(Ra) * sin(I);
	double T21 = sin(Ra) * cos(Wp) + cos(Ra) * cos(I) * sin(Wp);
	double T22 = -sin(Ra) * sin(Wp) + cos(Ra) * cos(I) * cos(Wp);
	//double T23 = -cos(Ra) * sin(I);
	double T31 = sin(I) * sin(Wp);
	double T32 = sin(I) * cos(Wp);
	//double T33 = cos(I);

	x[0] = T11 * R * cos(TA) + T12 * R * sin(TA);
	x[1] = T21 * R * cos(TA) + T22 * R * sin(TA);
	x[2] = T31 * R * cos(TA) + T32 * R * sin(TA);

	x[3] = T11 * -sqrt(Mu/p) * sin(TA) + T12 * sqrt(Mu/p) * (e + cos(TA));
	x[4] = T21 * -sqrt(Mu/p) * sin(TA) + T22 * sqrt(Mu/p) * (e + cos(TA));
	x[5] = T31 * -sqrt(Mu/p) * sin(TA) + T32 * sqrt(Mu/p) * (e + cos(TA));
}

void State_Vector_To_Elements(double x[6], double *A, double *e, double *I, double *Ra, double *Wp, double *TA){
	/*Convert State Vector x to Orbital elements a, e, I , RA, Wp, and TA
	  Define the State Vector in more easily readable terms R and V*/
	double Ri = x[0];
	double Rj = x[1];
	double Rk = x[2];
	double Vi = x[3];
	double Vj = x[4];
	double Vk = x[5];

	//Define the magnitudes of the R and V vectors
	double RdotV = Ri * Vi + Rj * Vj + Rk * Vk;
	double Rmag = sqrt(pow(Ri,2) + pow(Rj,2) + pow(Rk,2));
	double Vmag = sqrt(pow(Vi,2) + pow(Vj,2) + pow(Vk,2));

	//Determine the semimajor axis (a)
	*A = 1 / (2 / Rmag - pow(Vmag, 2) / Mu);

	//Define the angular momentum vector h, perpendicular to both r and v:
	double hi = Rj * Vk - Rk * Vj;
	double hj = Rk * Vi - Ri * Vk;
	double hk = Ri * Vj - Rj * Vi;
	double hmag = sqrt(pow(hi, 2) + pow(hj, 2) + pow(hk, 2));

	//To find eccentricity (e), first define the constant vector c
	//c = v x h - mu * (R/Rmag)
	double ci = (Vj * hk - Vk * hj) - Mu * Ri / Rmag;
	double cj = (Vk * hi - Vi * hk) - Mu * Rj / Rmag;
	double ck = (Vi * hj - Vj * hi) - Mu * Rk / Rmag;
	double cmag = sqrt(pow(ci,2) + pow(cj,2) + pow(ck,2));

	//Calculate eccentricity
	*e = cmag / Mu;

	//Define the constants for each quadrant of the direction cosine matrix
	double C11 = ci / Mu * *e;
	double C12 = cj / Mu * *e;
	double C13 = ck / Mu * *e;

	double C31 = hi / hmag;
	double C32 = hj / hmag;
	double C33 = hk / hmag;

	//double C21 = C32 * C13 - C33 * C12;
	//double C22 = C33 * C11 - C31 * C13;
	double C23 = C31 * C12 - C32 * C11;

	//Determine the Euler elements
	*Ra = atan2(C31, -C32); //acos(-hj / sqrt(pow(hi, 2) + pow(hj, 2)))
	*I = acos(C33);
	*Wp = atan2(C13, C23); 

	//Quadrant check
	//if (hi < 0) {*Ra = 2 * Pi - *Ra;}
	//if (ck < 0) {*Wp = 2 * Pi - *Wp;}

	//Defining the true anomaly (TA) requires the eccentricity vector
	double num = (ci * Ri + cj * Rj + ck * Rk) / Mu;
	double den = *e * Rmag;
	*TA = acos(num / den);

	if (RdotV < 0) {*TA = 2 * Pi - *TA;}

	*I = RadToDeg( PiBound(*I) );
	*Ra = RadToDeg( PiBound(*Ra) );
	*Wp = RadToDeg( PiBound(*Wp) );
	*TA = RadToDeg( PiBound(*TA) );
}

double *eci_to_ecf(double Rx, double Ry, double Rz){

	//Rotation matrix:
	// Inverse of ECI to ECEF:
	// [X]     [C  S  0][X]
	// [Y]  =  [-S C  0][Y]
	// [Z]ecf  [0  0  1][Z]eci
	double ECEF_r[3];

	//to convert from ECI to ECEF, i need to know the angle between greenwich meridian and the First point of Aries?
	//The MJD has a starting point of midnight on November 17, 1858
	double MJD = JD - 2400000.5;
	double D = (MJD - 51544.5) / 36525.;

	double alpha_m = 100.46061838 + D * (0.7700537 + 0.000388 * D) + 360 * (D * 100 - floor(D * 100));
	alpha_m = Range_2pi(alpha_m);

	//gmst from degrees to radians
	double gmst = alpha_m * Pi / 180.;

	ECEF_r[0] = Rx * cos(gmst) + Ry * sin(gmst);
	ECEF_r[1] = -Rx * sin(gmst) + Ry * cos(gmst);
	ECEF_r[2] = Rz;

	return ECEF_r;
}

double *ecf_to_eci(double Rx, double Ry, double Rz){

	//Rotation matrix:
	// [X]     [C -S  0][X]
	// [Y]  =  [S  C  0][Y]
	// [Z]eci  [0  0  1][Z]ecf
	double ECI_r[3];

	//to convert from ECI to ECEF, i need to know the angle between greenwich meridian and the First point of Aries?
	//The MJD has a starting point of midnight on November 17, 1858
	double MJD = JD - 2400000.5;
	double D = (MJD - 51544.5) / 36525.;

	double alpha_m = 100.46061838 + D * (0.7700537 + 0.000388 * D) + 360 * (D * 100 - floor(D * 100));
	alpha_m = Range_2pi(alpha_m);

	//gmst from degrees to radians
	double gmst = alpha_m * Pi / 180.;

	ECI_r[0] = Rx * cos(gmst) - Ry * sin(gmst);
	ECI_r[1] = Rx * sin(gmst) + Ry * cos(gmst);
	ECI_r[2] = Rz;

	return ECI_r;
}

double Calculate_Invariant_Elements(double x[6]) {
	//Calculate invariant elements as outputs
	double Ri = x[0];
	double Rj = x[1];
	double Rk = x[2];
	double Vi = x[3];
	double Vj = x[4];
	double Vk = x[5];

	//Define the magnitudes of the R and V vectors
	double RdotV = Ri * Vi + Rj * Vj + Rk * Vk;
	double Rmag = sqrt(pow(Ri, 2) + pow(Rj, 2) + pow(Rk, 2));
	double Vmag = sqrt(pow(Vi, 2) + pow(Vj, 2) + pow(Vk, 2));

	//Define the angular momentum vector h, perpendicular to both r and v:
	double hi = Rj * Vk - Rk * Vj;
	double hj = Rk * Vi - Ri * Vk;
	double hk = Ri * Vj - Rj * Vi;
	double hmag = sqrt(pow(hi, 2) + pow(hj, 2) + pow(hk, 2));
	
	double EarthRadius = 6370; // 6378.140;
	double C20 = -0.001082626;

	double rslr = pow(hmag, 2) / Mu; //Semi-latus rectum
	double i = acos(hk / hmag); //inclination
	double theta = asin((Rk / Rmag) / sin(i)); //argument of latitude
	//double theta = acos((hi*rj - hj*ri) / hmag*rmag*sin(i));

	//Gottlieb Formulae:
	double lambda = (3 / 2)*C20*pow(EarthRadius / rslr, 2)*(Mu / pow(rslr, 2));
	double rho = lambda * pow(sin(i), 2);
	double v = lambda - rho * (3 / 2 - cos(2 * theta));
	double nyay = hmag / pow(rslr, 2);
	double A = (Rmag - rslr) - (1 / pow(nyay,2))*(v - rho*cos(2 * theta) / 6);
	double B = (1 / nyay)*(RdotV / Rmag) - (rho*sin(2 * theta) / 3 * pow(nyay, 2));

	//double Inv_Elements[7];
	double Inv_SMA = rslr + v / pow(nyay, 2); //Invariant Semimajor Axis
	double Inv_Hsma = Inv_SMA - EarthRadius;

	return Inv_Hsma;
	}