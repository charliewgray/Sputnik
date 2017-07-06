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

/* ------------------------------------------------------------------- */
/* ------------------------- Shared Variables ------------------------ */
/* ------------------------------------------------------------------- */
//I want to use these two global variables within this c file, so i intialize them here
double Pi;
double Mu;
double We;
double a_Earth;     // Earth equatorial radius 
double b_Earth;    // Earth polar radius 
double SunVector[3], RSun;
//double JD;

/* ------------------------------------------------------------------- */
/* ------------ Latitude / Longitude Calculation Routine ------------- */
/* ------------------------------------------------------------------- */

double Range_2pi(double Angle) {
	//Force the RADIAN angle to fall between 0 and 2*Pi (0 to 360 degrees)
	//Create an integer value of the angle over 2Pi
	int IntegerVar = Angle / 2. / Pi;
	return Angle - 2. * Pi * IntegerVar;
}

double JulianDate(int year, int month, double day){
	//The Julian date (JD) of any instant is the Julian day number for the 
	//preceding noon in Greenwich Mean Time plus the fraction of the day since that instant.

	//http://en.wikipedia.org/wiki/Julian_day
	//Compute the Julian Day first
	double a = floor((14 - month) / 12.);
	double y = year + 4800 - a;
	double m = month + 12 * a - 3;

	double JDN = day + 365 * y + floor(y/4) - floor(y/100) + floor(y/400) - 32045;
	return JDN;
}

double Geo_Alt(double Latitude){
	//Calculate the Geodetic altitude 
	double Geodetic_Lat = Latitude * Pi / 180.;

	//it's not the angle from the center of the earth, so we can't compute the radius with this value alone
	double altRatio = pow(b_Earth, 2) / pow(a_Earth, 2);

	//Calculate the geocentric latitude, which is the angle from the center of the earth to the R vector
	double Geocentric_Lat = atan(altRatio * tan(Geodetic_Lat));

	//Solve for R using ellipse equation [p^2/a^2 + z^2/b^2 = 1] & p = r*cos(geocentricLat) & z = r*sin(geocentricLat)
	//double Eq1 = pow(cos(Geocentric_Lat) / a_Earth, 2);
	//double Eq2 = pow(sin(Geocentric_Lat) / b_Earth, 2);

	return sqrt(1 / (pow(cos(Geocentric_Lat) / a_Earth, 2) + pow(sin(Geocentric_Lat) / b_Earth, 2)));
}

void Longitude_Latitude(int year, int month, double day, double x[6], double *Longitude, double *Latitude){
	//add these variables to make the rest easier to comprehend
	double Ri = x[0], Rj = x[1], Rk = x[2];

	//Latitude is the easy part...
	double R_mag = sqrt(pow(Ri, 2) + pow(Rj, 2) + pow(Rk, 2));
	*Latitude = asin(Rk / R_mag) * 180. / Pi;

	//Convert earth rotation to degrees and julian date
	double we_degree = We * 180. / Pi;
	double JD = JulianDate(year, month, day);

	//*****************Modified Julian Date***********************
	//The MJD has a starting point of midnight on November 17, 1858
	double MJD = JD - 2400000.5;
	double D = (MJD - 51544.5) / 36525.;

	//Determine the angle since midnight of the current date
	double alpha_m = 100.46061838 + D * (0.7700537 + 0.000388 * D) + 360 * (D * 100 - floor(D * 100));
	alpha_m = Range_2pi(alpha_m);

	//Determine how much time has passed since midnight of the current date (and how much rotation has occurred)
	int intDay = day;
	double t_0 = (day - intDay) * 3600 * 24;
	double prog_today = t_0 * we_degree;

	//Alpha = Alpha @ midnight + Alpha rotation since midnight
	double ALPHA = alpha_m + prog_today;
	ALPHA = Range_2pi(ALPHA);

	//Alpha describes the current Greenwich 0 longitude line, still need the
	//ascending node to calculate the exact location of ISS w.r.t. Greenwich 0

	//Calculate Right Ascension and Declination
	double Right_asc = acos(Ri / sqrt(pow(Ri,2) + pow(Rj,2))) * 180. / Pi;
	int Declination = asin(Rj / sqrt(pow(Ri,2) + pow(Rj,2))) * 180. / Pi;
	if (Declination <= 0){Right_asc = -Right_asc;}

	//Find the difference between ISS longitude location and Greenwich 0
	*Longitude = Right_asc - ALPHA;

	//Make sure the numbers fall within the 0 - 360
	if (*Longitude < -180){*Longitude = *Longitude + 360;}
	else if (*Longitude > 180){*Longitude = *Longitude - 360;}
}

void SunVectorCalc(double JD){
	//******************************************************************************************
	//*****************Angles associated with the position of the sun **************************
	//******************************************************************************************
	
	double n = JD - 2451545.0;                          //Compute: time since Jan 1, 2000: n

	double RAANsun = 2.1429 - 0.0010394594 * n;
	RAANsun = Range_2pi(RAANsun);

	//mean longitude of sun at JD 2451545.0 (deg) + rate of mean longitude of sun (deg/day)
	double Lsun = 4.895062994 + 0.017202791698 * n;     //Compute: Mean longitude of sun
	Lsun = Range_2pi(Lsun);

	//mean anomaly of sun at JD 2451545.0 (deg) + rate of mean anomaly of sun (deg/day)
	double Msun = 6.2400600 + 0.0172019699 * n;         //Compute: Mean anomaly of sun
	Msun = Range_2pi(Msun);

	//ecliptic longitude of the sun
	double Lsun_Ecl = Lsun + 0.03341607 * sin(Msun) + 0.00034894 * sin(2 * Msun) - 0.0001134 - 0.0000206 * sin(RAANsun);  //Compute: Ecliptic longitude of sun
	Lsun_Ecl = Range_2pi(Lsun_Ecl);

	//obliquity of the ecliptic at JD 2451545.0 (deg) + rate of obliquity (deg/day)
	double Obl_Ecl = 0.4090928 - 6.2140 * pow(10, -9) * n + 0.0000396 * cos(RAANsun);  //Compute: Obliquity of ecliptic
	Obl_Ecl = Range_2pi(Obl_Ecl);

	//Create the sun (unit) vector in ECI coordinates
	SunVector[0] = cos(Lsun_Ecl);
	SunVector[1] = cos(Obl_Ecl) * sin(Lsun_Ecl);
	SunVector[2] = sin(Obl_Ecl) * sin(Lsun_Ecl);
	//******************************************************************************************

	//Now to find the radius magnitude (distance from earth to sun)
	double T = n / 36525.;
	double eEarth = 0.016708617 - 0.000042037 * T - 0.0000001236 * pow(T, 2);  //Earth orbit eccentricity
	double TAEarth = Msun + Lsun_Ecl - Lsun;                                   //true anomaly of earth
	double aEarth = 149.6 * pow(10, 6);                                        //Earth semimajor axis
	
	//R = a * (1/e^2) / (1 + e*cos(TA))
	RSun = aEarth * (1 - pow(eEarth, 2)) / (1 + eEarth * cos(TAEarth));
}

double BetaCalcs(double RAAN, double I, double e, double JD){
	// Analytic Calculation of Solar Beta - Taken from STRAT / TPS
	//     References:  "Astronomical Almanac 1995"
	//                  "BGFCON-SEDA User's Guide"

	//Beta angle is the angle between the satellite orbit plane and the sun vector,
	//so we need to define those two vectors first

	//Convert the Ra, I, and e to radians
	RAAN = RAAN * Pi / 180.;
	I = I * Pi / 180.;
	e = e * Pi / 180.;

	//use the sunvector subroutine to calculate the sun vector
	SunVectorCalc(JD);

	//h vector = r x v, perpendicular to satellite orbit plane
	double h[3] = { sin(RAAN)*sin(I), -cos(RAAN)*sin(I), cos(I) };

	//The dot product between s and v vectors since: sin(beta) = s_dot_v
	double sdotv = h[0] * SunVector[0] + h[1] * SunVector[1] + h[2] * SunVector[2];
	double SolarBeta = asin(sdotv) * 180. / Pi;  

	return SolarBeta;
}