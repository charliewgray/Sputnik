/* -------------------------------------------------------------------- */
/* ----------      CHARLIE GRAY DECAY CALCULATIONS 2015      ---------- */
/* -------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* ------------------------------ INCLUDES --------------------------- */
/* ------------------------------------------------------------------- */
#include <math.h>          /* maths functions */
double Pi;
double Mu;

/* ------------------------------------------------------------------- */
/* ------------------------- Shared Variables ------------------------ */
/* ------------------------------------------------------------------- */
void ConvertAndApplydV(double dt, double A, double e, double I, double Ra, double Wp, double TA, double x[6], double dV){
	//Routine applies the input dV from a local frame to the earth body coordinate frame x,y,z
	//Input dV(1) = Vx, dV(2) = Vy, dV(3) = Vz

	//Convert the input fields that are in degrees to radians
	I = I * Pi / 180.;
	Wp = Wp * Pi / 180.;
	Ra = Ra * Pi / 180.;
	TA = TA * Pi / 180.;

	//Calculate magnitude of dV
	double VMag = sqrt(pow(x[3], 2) + pow(x[4], 2) + pow(x[5], 2));
	
	//Create the identity matrix
	double Vx = x[3] / VMag;
	double Vy = x[4] / VMag;
	double Vz = x[5] / VMag;
	
	//add dV to the velocity magnitude
	VMag = VMag + dV / 1000.;

	//add magnitude back into the identity matrix
	x[3] = Vx * VMag;
	x[4] = Vy * VMag;
	x[5] = Vz * VMag;
}