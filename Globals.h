// This is start of the header guard.  ELEMENTCONVERTERS_H can be any unique name.  By convention, we use the name of the header file.
#ifndef GLOBALS_H
#define GLOBALS_H

/*define constant values*/
extern double Pi = 3.141592653589793238463;

extern double Mu = 398600.4418;       // Earth Standard Gravitational Parameter (MG) = km^3/s^2 
extern double MuS = 132712440018;     // Sun Standard Gravitational Parameter (MG) = km^3/s^2 
extern double MuMo = 4902.794;        // Moon Standard Gravitational Parameter (MG) = km^3/s^2 

extern double We = 0.000072921156001; // Earth Angular velocity radians/s 
extern double a_Earth = 6378.137;     // Earth equatorial radius 
extern double b_Earth = 6356.7523;    // Earth polar radius 
extern double REarth = 6370;                 // Earth radius at an instantaneous moment, to be calculated in the code

// Earth Ablateness Constants J2 - J6
extern double J[7] = { 0, 0, 0.001082635, -0.000002531, -0.00000160, -0.00000015, 0.00000057 };
extern double C[5][5] = { { 0, 0, 0, 0, 0 },
						{ 0, 0, 0, 0, 0 },
						{ 0, -2.414E-10, 1.57446037456E-6, 0, 0 },
						{ 0, 2.19263852717E-6, 3.08989206881E-7, 1.00548778064E-7, 0 },
						{ 0, -5.08799360404E-7, 7.84175859844E-8, 5.92099402629E-8, -3.98407411766E-9 } };

extern double S[5][5] = { { 0, 0, 0, 0, 0 },
						{ 0, 0, 0, 0, 0 },
						{ 0, 1.5431E-9, -9.03803806639E-7, 0, 0 },
						{ 0, 2.68424890297E-7, -2.11437612437E-7, 1.97222559006E-7, 0 },
						{ 0, -4.49144872839E-7, 1.48177868296E-7, -1.20077667634E-8, 6.52571425370E-9 } };

extern double JD;					  //Julian Date (needs to be accessible to a large number of calculations)
extern int Month_Array[13];           //Array to keep the number of days in each month (changes as a function of year)
extern char SFUmatrix[1500][8][10];    //Array to keep the input SFU and Geomagnetic Index predictions (change as a function of month / year)
extern float f107A;                   // 81 day average of F10.7 flux (centered on doy) 
extern float f107;                    // daily F10.7 flux for previous day 
extern float ap;                      // magnetic index(daily) 

extern char Month_Names[13][4] = { "\0", "JAN\0", "FEB\0", "MAR\0", "APR\0", "MAY\0", "JUN\0", "JUL\0", "AUG\0", "SEP\0", "OCT\0", "NOV\0", "DEC\0" };
//extern int Atmosphere_Select;         //A way to select the atmosphere to use: 1 = M95, 2 = M50, 3 = M5
extern char inputMatrix[19][40];
extern char ReboostsStringMatrix[50][40];
extern float Reboosts[50][8];
extern char *atmoStr;
extern char Configmatrix[500][13][20];
extern int ConfigDates[500][4];
extern int ConfigBetas[10];
extern int currentSFU_row;
int rows;

extern double SunVector[3];         //Describes the position of the sun as a function of julian date
extern double RSun;                 //Magnitude of the sun vector
extern double ECEF_r[3];

// This is the end of the header guard
#endif

/*extern double J2 = 0.001082635;       // Earth Ablateness Constant J2
extern double J3 = -0.000002531;      // Earth Ablateness Constant J3
extern double J4 = -0.00000160;       // Earth Ablateness Constant J4
extern double J5 = -0.00000015;       // Earth Ablateness Constant J5
extern double J6 = 0.00000057;        // Earth Ablateness Constant J6 */