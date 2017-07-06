/* -------------------------------------------------------------------- */
/* ------------------  MAIN ATMOSPHERIC DECAY MODEL  ------------------ */
/* -------------------------------------------------------------------- */

/*---------------------------------------------------------------------*/
/*-----This file is part of the Decay Calc  C source code package------*/
/*-------The orbital decay model was developed by Charles Gray.--------*/
/*---------------------------------------------------------------------*/

/* ------------------------------------------------------------------- */
/* ------------------------------ INCLUDES --------------------------- */
/* ------------------------------------------------------------------- */
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Globals.h"
#include "InputReaders.h"
#include "TimeFunctions.h"
#include "ElementConverters.h"
#include "LatLongBeta.h"
#include "Density.h"
#include "Propagator.h"
#include "Maneuver.h"

/* ------------------------------------------------------------------- */
/* -----------------Global Variables to Include ---------------------- */
/* ------------------------------------------------------------------- */
int reboostFlag;
int boostCount;
float f107A, f107, ap;  
float Reboosts[50][8];
char inputMatrix[19][40];
char *atmoStr;
double SunVector[3], RSun;
int currentSFU_row;

/* ------------------------------------------------------------------- */
/* -------------------------------- MAIN ----------------------------- */
/* ------------------------------------------------------------------- */

void DecayProgram(int runlengthinseconds, int dt, int outputincrement, int Atmosphere_Select, int caseNumber){
	//------------------------------INPUTS------------------------------------------------------
	//Define inputs to the entire code (user inputs eventually)
	char *objName;
	objName = strtok(inputMatrix[0], "\t");

	//initialize time
	int year = atoi(inputMatrix[1]);    
	double day = strtod(inputMatrix[2], NULL);
	int month;

	double A = strtod(inputMatrix[3], NULL);    //semimajor axis in m
	double e = strtod(inputMatrix[4], NULL);    //eccentricity
	double I = strtod(inputMatrix[5], NULL);    //inclination
	double Wp = strtod(inputMatrix[6], NULL);   //argument of perigee
	double Ra = strtod(inputMatrix[7], NULL);   //longitude of ascending node
	double TA = strtod(inputMatrix[8], NULL);   //true anomaly

	//Placeholder values for ECI coordinates and geodetic coordinate variables
	double x[6] = { 0, 0, 0, 0, 0, 0 };
	double Long, Lat, rho;

	//Calculate Ballistic coefficient (beta) and Ballistic Number (BN) 
	double BN = strtod(inputMatrix[9], NULL);   //Ballistic Number (most familiar to me)

	//Calculate Apogee and perigee
	double Ha = A * (1 + e);
	double Hp = A * (1 - e);

	//Define the number of days in each month of this year (leap year calculation)
	defineMonthArray(year);

	//Calculate GMT
	int GMT = day; //- only used if user inputs month and day instead of GMT

	//Determine the month
	find_month_from_GMT(&month, GMT);

	//Initialize solar flux on this date
	findSolarFlux(year, month, Atmosphere_Select);

	//Convert the input variables into cartesian coordinates
	Elements_To_Vector(A, e, I, Ra, Wp, TA, x);

	//Initialize the latitude and longitude of the vehicle at this time
	Longitude_Latitude(year, month, day, x, &Long, &Lat);

	//initialize the geodetic altitude
	double Rmag = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
	
	//Calculate the earth radius at this latitude
	double Geodetic_Alt = Rmag - REarth; // REarth;

	//Calculate Julian Date and solar beta angle
	double JD = JulianDate(year, month, day);
	double SolarBeta = BetaCalcs(Ra, I, e, JD);

	//Run the density function
	rho = get_density(day, Geodetic_Alt, Lat, Long, f107A, f107, ap);
	rho *= pow(1000, 3);       //Convert from meters to km
	
	//open the output file
	char outputName[30]; 
	snprintf(outputName, sizeof outputName, "output_%s.txt", objName);
	FILE* f = fopen(outputName, "w");
	if (f == NULL) {perror("Error opening output file!\n");}

	//Print the header in the output file
	fprintf(f, "year,day,A,e,I,Wp,Ra,TA,Density,BN,Hsma,Ha,Hp\n");

	//create a progress bar!!! - This one is a comparison bar to show how much "complete" will be
	printf("Case #%d: %s\n|___________________________________________________________|\n|", caseNumber, objName);

	//variables for the progress bar counter!!!
	double progressincrement = runlengthinseconds / 60;
	double incStep = progressincrement;

	//Initialize the variable for the output counter (will change as the run progresses)
	double outputincStep = outputincrement;
	double elapsedseconds = 0;     //Essentially just a loop counter for the upcoming while statement
	double Inv_Hsma;               // = Calculate_Invariant_Elements(x);

	//Loop through the time frames until i've hit the run length or the apogee falls below ISS
	while (elapsedseconds < runlengthinseconds && Ha > REarth + 150) {

		//--------------------------------------Print the results -----------------------------------------
		if (elapsedseconds == outputincStep || elapsedseconds == 0) {
			//Convert the x vector back to orbital elements so i can calculate apogee and perigee
			State_Vector_To_Elements(x, &A, &e, &I, &Ra, &Wp, &TA);
			
			//calculate apogee and perigee
			Ha = A * (1 + e);
			Hp = A * (1 - e);
			
			//Calcualate Invariant SMA
			Inv_Hsma = Calculate_Invariant_Elements(x);
			
			//Print the data from the analysis
			fprintf(f, "%d, %f", year, day);
			fprintf(f, ",%f,%f,%f,%f,%f,%f", A, e, I, Wp, Ra, TA);
			fprintf(f, ",%1.3e", rho);
			fprintf(f, ",%f", BN);
			fprintf(f, ",%f", Inv_Hsma);
			fprintf(f, ",%f,%f", Ha-REarth, Hp-REarth);
			fprintf(f, "\n");

			//Increment the output step so we know which line to output next time
			if (elapsedseconds != 0) { outputincStep += outputincrement; }
		}

		//There are 40 underscores to fill to show completion of the program. 
		//every time i hit a 1/40 increment, add another underscore
		if (elapsedseconds == incStep) {
			printf("_");
			incStep += progressincrement;
		}
		//------------------------------------------------------------------------------------------------------

		//Increment the time
		elapsedseconds = elapsedseconds + dt;

		//Redefine year, GMT, month, day, hour, minute, & second variables
		TimeStep(dt, &GMT, &year, &month, &day, Atmosphere_Select);

		//Recalculate Julian Date 
		JD = JulianDate(year, month, day);

		//Run the density function 
		rho = get_density(GMT, Geodetic_Alt, Lat, Long, f107A, f107, ap);
		rho *= pow(1000, 3); //Convert from meters to km

		//Apply Perturbations on the vehicle and its effect over the chosen time step
		RK4(dt, BN, rho, Lat, Long, x);

		//Determine the new latitude and longitude of the vehicle at this time
		Longitude_Latitude(year, month, day, x, &Long, &Lat);

		//Determine the geodetic altitude
		Rmag = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));

		//Calculate the radius of earth at this latitude
		Geodetic_Alt = Rmag - REarth;

		//Calculate Beta angle
		//SolarBeta = BetaCalcs(Ra, I, e, JD);
	}

	//close the output file
	fprintf(f, "# of Mvrs: %1d", boostCount-1);
	fprintf(f, "\n The object has deorbited within %f days or ~%f months or ~%f years \n", elapsedseconds / 3600 / 24, elapsedseconds / 3600 / 24 / 30.5, elapsedseconds / 3600 / 24 / 365);
	fclose(f);
	
	//If the object has completedly deorbited, output the number of days in orbit, otherwise notify user that the program is complete
	if (runlengthinseconds != elapsedseconds) {
		int elapseddays = elapsedseconds / (3600 * 24);
		printf("\n The object has crossed the set lower Alt limit within %d days \n", elapseddays);
	}
	else{ 
		printf("|\n Program Complete. Final Invariant SMA: %f km \n", Inv_Hsma);
	}

}

void main() {
	//Read in the input file & record the data into a public variable matrix
	int caseNumber = 0;
	while (caseNumber < 5000) {
		openInputFile(caseNumber);

		//Read in the atmosphere file & record the data into a public variable SFUMatrix
		int Atmosphere_Select = atoi(inputMatrix[10]);
		openAtmosphereFile(Atmosphere_Select);

		//Assign a string to define what atmosphere is being used
		if (Atmosphere_Select == 0) { atmoStr = "M95"; }
		if (Atmosphere_Select == 1) { atmoStr = "M50"; }
		if (Atmosphere_Select == 2) { atmoStr = "M05"; }
		if (Atmosphere_Select == 3) { atmoStr = "M75"; }

		//Define the time increments for the analysis
		int dt = atoi(inputMatrix[11]);                                //time step in seconds
		int runlengthinseconds = atoi(inputMatrix[12]) * 3600 * 24;    //How long with the program run (seconds)?
		int outputincrement = atoi(inputMatrix[13]);			       //how frequent should i output data (seconds)
																	   //Reset the SFU row
		currentSFU_row = 2;

		//Check if the proram should continue running
		if (dt == 0) { 
			goto endProgram; 
		}

		//Run the program
		DecayProgram(runlengthinseconds, dt, outputincrement, Atmosphere_Select, caseNumber);

		//move onto the next case in the input file
		caseNumber++;
	}
endProgram:;
}