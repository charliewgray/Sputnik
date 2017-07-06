/* -------------------------------------------------------------------- */
/* ----------      CHARLIE GRAY DECAY CALCULATIONS 2015      ---------- */
/* ------------------------------------------------------------------- */
/* ------------------------------ INCLUDES --------------------------- */
/* ------------------------------------------------------------------- */
#include <math.h>          /* maths functions */
/* ------------------------------------------------------------------- */
/* ------------------------- Shared Variables ------------------------ */
/* ------------------------------------------------------------------- */
//I want to use these two global variables within this c file, so i intialize them here
double Mu, MuS;                                // Universal Gravitational Constants (MG) = kg^3/s^2 
double We;                                     // Earth Angular velocity radians/s 
double J[7], C[5][5], S[5][5];
double a_Earth, b_Earth, REarth;              // Earth radii
double Pi;

//This version is the short version of earth harmonics, does not include tesseral or sectoral harmonic variations
void EarthHarmonics(double aH[3], double Rx, double Ry, double Rz, double Rmag, double Lat, double Long){
	//This subroutine is specifically to calculate Earth harmonics
	//The current version only covers Zonal harmonics...
	//**************** J2 - J6 Effects ********************************************************
	//double Req = a_Earth;

	//double RxRmag = Rx / Rmag;
	//double RyRmag = Ry / Rmag;
	//double RzRmag = Rz / Rmag;

	//**************** J2 ********************************************************************
	//double aJ2x = -1.5 * J[2] * Mu * pow(a_Earth, 2) * pow(Rmag, -4) * (1 - 5 * pow(Rz / Rmag, 2)) * Rx / Rmag;
	//double aJ2y = -1.5 * J[2] * Mu * pow(a_Earth, 2) * pow(Rmag, -4) * (1 - 5 * pow(Rz / Rmag, 2)) * Ry / Rmag;
	//double aJ2z = -1.5 * J[2] * Mu * pow(a_Earth, 2) * pow(Rmag, -4) * (3 - 5 * pow(Rz / Rmag, 2)) * Rz / Rmag;

	//**************** J3 ********************************************************************
	//double aJ3x = -0.5 * J[3] * Mu * pow(a_Earth, 3) * pow(Rmag, -5) * (15 * Rz / Rmag - 35 * pow(Rz / Rmag, 3)) * Rx / Rmag;
	//double aJ3y = -0.5 * J[3] * Mu * pow(a_Earth, 3) * pow(Rmag, -5) * (15 * Rz / Rmag - 35 * pow(Rz / Rmag, 3)) * Ry / Rmag;
	//double aJ3z = -0.5 * J[3] * Mu * pow(a_Earth, 3) * pow(Rmag, -5) * (30 * pow(Rz / Rmag, 2) - 35 * pow(Rz / Rmag, 4) - 3);

	//**************** J4 ********************************************************************
	//double aJ4x = 0.125 * J[4] * Mu * pow(a_Earth, 4) * pow(Rmag, -6) * (15 - 210 * pow(Rz / Rmag, 2) + 315 * pow(Rz / Rmag, 4)) * Rx / Rmag;
	//double aJ4y = 0.125 * J[4] * Mu * pow(a_Earth, 4) * pow(Rmag, -6) * (15 - 210 * pow(Rz / Rmag, 2) + 315 * pow(Rz / Rmag, 4)) * Ry / Rmag;
	//double aJ4z = 0.125 * J[4] * Mu * pow(a_Earth, 4) * pow(Rmag, -6) * (75 - 350 * pow(Rz / Rmag, 2) + 315 * pow(Rz / Rmag, 4)) * Rz / Rmag;

	//**************** J5 ********************************************************************
	//double aJ5x = 0.375 * J[5] * Mu * pow(a_Earth, 5) * pow(Rmag, -7) * (35 * RzRmag - 210 * pow(RzRmag, 3) + 210 * pow(RzRmag, 5)) * RxRmag;
	//double aJ5y = 0.375 * J[5] * Mu * pow(a_Earth, 5) * pow(Rmag, -7) * (35 * RzRmag - 210 * pow(RzRmag, 3) + 210 * pow(RzRmag, 5)) * RyRmag;
	//double aJ5z = 0.375 * J[5] * Mu * pow(a_Earth, 5) * pow(Rmag, -7) * pow(RzRmag, 2) * (105 - 315 * pow(RzRmag, 2) + 231 * pow(RzRmag, 4)) - 1.875 * J[5] * Mu * pow(Req, 5) * pow(Rmag, -7);

	//**************** J6 ********************************************************************
	//double aJ6x = -0.0625 * J[6] * Mu * pow(a_Earth, 6) * pow(Rmag, -8) * (35 - 945 * pow(RzRmag, 2) + 3465 * pow(RzRmag, 4) - 3003 * pow(RzRmag, 6)) * RxRmag;
	//double aJ6y = -0.0625 * J[6] * Mu * pow(a_Earth, 6) * pow(Rmag, -8) * (35 - 945 * pow(RzRmag, 2) + 3465 * pow(RzRmag, 4) - 3003 * pow(RzRmag, 6)) * RyRmag;
	//double aJ6z = -0.0625 * J[6] * Mu * pow(a_Earth, 6) * pow(Rmag, -8) * (245 - 2205 * pow(RzRmag, 2) + 4851 * pow(RzRmag, 4) - 3003 * pow(RzRmag, 6)) * RzRmag;

	//Add these effects to the accelleration vector
	//aH[0] = aJ2x + aJ3x + aJ4x;// +aJ5x + aJ6x;
	//aH[1] = aJ2y + aJ3y + aJ4y;// +aJ5y + aJ6y;
	//aH[2] = aJ2z + aJ3z + aJ4z;// +aJ5z + aJ6z;

	aH[0] = -1.5 * J[2] * Mu * pow(a_Earth, 2) * pow(Rmag, -4) * (1 - 5 * pow(Rz / Rmag, 2)) * Rx / Rmag + -0.5 * J[3] * Mu * pow(a_Earth, 3) * pow(Rmag, -5) * (15 * Rz / Rmag - 35 * pow(Rz / Rmag, 3)) * Rx / Rmag + 0.125 * J[4] * Mu * pow(a_Earth, 4) * pow(Rmag, -6) * (15 - 210 * pow(Rz / Rmag, 2) + 315 * pow(Rz / Rmag, 4)) * Rx / Rmag;// +aJ5x + aJ6x;
	aH[1] = -1.5 * J[2] * Mu * pow(a_Earth, 2) * pow(Rmag, -4) * (1 - 5 * pow(Rz / Rmag, 2)) * Ry / Rmag + -0.5 * J[3] * Mu * pow(a_Earth, 3) * pow(Rmag, -5) * (15 * Rz / Rmag - 35 * pow(Rz / Rmag, 3)) * Ry / Rmag + 0.125 * J[4] * Mu * pow(a_Earth, 4) * pow(Rmag, -6) * (15 - 210 * pow(Rz / Rmag, 2) + 315 * pow(Rz / Rmag, 4)) * Ry / Rmag;// +aJ5y + aJ6y;
	aH[2] = -1.5 * J[2] * Mu * pow(a_Earth, 2) * pow(Rmag, -4) * (3 - 5 * pow(Rz / Rmag, 2)) * Rz / Rmag + -0.5 * J[3] * Mu * pow(a_Earth, 3) * pow(Rmag, -5) * (30 * pow(Rz / Rmag, 2) - 35 * pow(Rz / Rmag, 4) - 3) + 0.125 * J[4] * Mu * pow(a_Earth, 4) * pow(Rmag, -6) * (75 - 350 * pow(Rz / Rmag, 2) + 315 * pow(Rz / Rmag, 4)) * Rz / Rmag;// +aJ5z + aJ6z;
}

double Legendre(int l, int m, double phi){
	//Legendre polynomials
	double P;

	switch (l){
	case 0:
		P = 1;
		break;
	case 1:
		switch (m){
		case 0:
			//P(1,0)
			P = sin(phi);
			break;
		case 1:
			//P(1,1)
			P = cos(phi);
			break;
		case 2:
			//P(1,2)
			P = 0;
			break;
		}
		break;
	case 2:
		switch (m){
		case 0:
			//P(2,0)
			P = 0.5 * (3 * pow(sin(phi), 2) - 1);
			break;
		case 1:
			//P(2,1)
			P = 3 * sin(phi) * cos(phi);
			break;
		case 2:
			//P(2,2)
			P = 3 * pow(cos(phi), 2);
			break;
		case 3:
			//P(2,3)
			P = 0;
			break;
		}
		break;
	case 3:
		switch (m){
		case 0:
			//P(3,0)
			P = 0.5 * (5 * pow(sin(phi), 3) - 3 * sin(phi));
			break;
		case 1:
			//P(3,1)
			P = 0.5 * cos(phi) * (15 * pow(sin(phi), 2) - 3);
			break;
		case 2:
			//P(3,2)
			P = 15 * pow(cos(phi), 2) * sin(phi);
			break;
		case 3:
			//P(3,3)
			P = 15 * pow(cos(phi), 3);
			break;
		case 4:
			//P(3,4)
			P = 0;
			break;
		}
		break;
	case 4:
		switch (m){
		case 0:
			//P(4,0)
			P = 0.125 * (35 * pow(sin(phi), 4) - 30 * pow(sin(phi), 2) + 3);
			break;
		case 1:
			//P(4,1)
			P = 2.5 * cos(phi) * (7 * pow(sin(phi), 3) - 3 * sin(phi));
			break;
		case 2:
			//P(4,2)
			P = 7.5 * pow(cos(phi), 2) * (7 * pow(sin(phi), 2) - 1);
			break;
		case 3:
			//P(4,3)
			P = 105 * pow(cos(phi), 3) * sin(phi);
			break;
		case 4:
			//P(4,4)
			P = 105 * pow(cos(phi), 4);
			break;
		case 5:
			//P(4,5)
			P = 0;
			break;
		}
		break;
	}
	return P;
}

void PotentialPartials(double *dUdr, double *dUdphi, double *dUdlambda, double Req, double Rmag, double phi, double lambda){
	//Sum the nonspherical partial derivatives
	int LandM[15][2] = { { 2, 0 }, { 2, 1 }, { 2, 2 }, { 3, 0 }, { 3, 1 }, { 3, 2 }, { 3, 3 }, { 4, 0 }, { 4, 1 }, { 4, 2 }, { 4, 3 }, { 4, 4 } };
	//int LandM[9][2] = {{ 2, 1 }, { 2, 2 }, { 3, 1 }, { 3, 2 }, { 3, 3 }, { 4, 1 }, { 4, 2 }, { 4, 3 }, { 4, 4 } };
	int m, l;

	//Initialize the values of the partial derivatives
	*dUdr = 0;
	*dUdphi = 0;
	*dUdlambda = 0;

	for (int count = 0; count < 15; count++){
		l = LandM[count][0];
		m = LandM[count][1];

		double P = Legendre(l, m, phi);
		double P1 = Legendre(l, m + 1, phi);

		*dUdr += pow(Req / Rmag, l) * (l + 1) * P * (C[l][m] * cos(m * lambda) + S[l][m] * sin(m * lambda));
		*dUdphi += pow(Req / Rmag, l) * (P1 - m * tan(phi) * P) * (C[l][m] * cos(m * lambda) + S[l][m] * sin(m * lambda));
		*dUdlambda += pow(Req / Rmag, l) * m * P * (S[l][m] * cos(m * lambda) - C[l][m] * sin(m * lambda));
	}

	//Add the mu over r terms
	*dUdr *= -Mu / pow(Rmag, 2);
	*dUdphi *= Mu / Rmag;
	*dUdlambda *= Mu / Rmag;
}

//This version is the new version of earth harmonics, much more detailed
void EarthHarmonicsDetailed(double aH[3], double Rx, double Ry, double Rz, double Rmag, double Lat, double Long){
	double Req = a_Earth;

	//*********** To describe the J2-J6 effects, i need the partial potential function values first**************
	double dUdphi, dUdlambda, dUdr;
	double phi = Lat * Pi / 180;
	double lambda = Long * Pi / 180;

	PotentialPartials(&dUdr, &dUdphi, &dUdlambda, Req, Rmag, phi, lambda);

	//Determine the perturbations to the cartesian accelerations
	double ai = ((1 / Rmag) * dUdr - dUdphi * Rz / (pow(Req, 2) * sqrt(pow(Rx, 2) + pow(Ry, 2)))) * Rx - (dUdlambda / sqrt(pow(Rx, 2) + pow(Ry, 2))) * Ry;
	double aj = ((1 / Rmag) * dUdr - dUdphi * Rz / (pow(Req, 2) * sqrt(pow(Rx, 2) + pow(Ry, 2)))) * Ry - (dUdlambda / sqrt(pow(Rx, 2) + pow(Ry, 2))) * Rx;
	double ak = (1 / Rmag) * dUdr * Rz + dUdphi * sqrt(pow(Rx, 2) + pow(Ry, 2)) / pow(Rmag, 2);

	//now that i have the potential functions, i can define the components of nonspherical accelleration
	aH[0] = ai;
	aH[1] = aj;
	aH[2] = ak;
}
