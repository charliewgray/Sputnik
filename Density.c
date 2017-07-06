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
#define _CRT_SECURE_NO_DEPRECATE
#include <math.h>          /* maths functions */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrlmsise-00.h"

float f107A;                 // 81 day average of F10.7 flux (centered on doy) 
float f107;                  // daily F10.7 flux for previous day 
float ap;                    // magnetic index(daily) 

char SFUmatrix[1500][8][10];   // matrix with SFU predictions
char Month_Names[13][4];
int currentSFU_row;
int rows;

//int Atmosphere_Select;

/* ------------------------------------------------------------------- */
/* ------------------------------ RUN_GTD7 DENSITY-------------------- */
/* ------------------------------------------------------------------- */
void findSolarFlux(int year, int month, int Atmosphere_Select) {
	//Searches the public SFUmatrix for the current year and month and updates the f107A, f107, and ap values

	//Read in the data from the file that's just been opened
	//Function atoi converts a string to an integer
	//function strtod(string to conver, variable for the new double value) converts a string to a double
	char temp[40];
	char currentMonth[4];
	char MatrixMonth[4];

	strcpy(currentMonth, Month_Names[month]);

	//int currentSFU_row = 2;

	//assign the year string to the temp variable, then convert it to an actual numerical value
	strcpy(temp, SFUmatrix[currentSFU_row][0]);
	int MatrixYear = atoi(temp);

	//Define the number of rows and columns in the SFUmatrix
	if (rows == 0) {
		rows = sizeof SFUmatrix / sizeof SFUmatrix[0];
	}

	//Loop currentSFU_row variable until i find the row that starts on the correct year
	while (MatrixYear != year) {
		currentSFU_row++;
		strcpy(temp, SFUmatrix[currentSFU_row][0]);
		MatrixYear = atoi(temp);

		//break out of this loop if i've run out of SFU data
		if (currentSFU_row >= rows) break;
	}

	//Find the right month
	if (currentSFU_row < rows) {
		int returnVal = 1;
		while (returnVal != 0) {
			strcpy(MatrixMonth, SFUmatrix[currentSFU_row][1]);
			//do the two strings have the same first 3 characters?
			if (strcmp(MatrixMonth, currentMonth) == 0) {
				returnVal = 0; //get out of this loop
			}
			else {
				currentSFU_row++; //increment to the next row
			}
		}
	}

	//Record the f107 value to the global variable
	strcpy(temp, SFUmatrix[currentSFU_row][Atmosphere_Select + 2]);
	f107A = strtod(temp, NULL);      // 81 day average of F10.7 flux (centered on doy) 
	f107 = f107A;                    // daily F10.7 flux for previous day 
			
	//Record the ap value to the global variable
	strcpy(temp, SFUmatrix[currentSFU_row][Atmosphere_Select + 5]);
	ap = strtod(temp, NULL);         // magnetic index(daily) 

}

double get_density(double day, double alt, double g_lat, double g_long, double f107A, double f107, double ap) {
	struct nrlmsise_output output;
	struct nrlmsise_input input;
	struct nrlmsise_flags flags;
	struct ap_array aph;

	/* input values */
	for (int i = 0; i < 7; i++){
		aph.a[i] = ap;
	}

	int GMT = day;

	//Set the flags
	for (int i = 0; i<24; i++)
		flags.switches[i] = 1;

	flags.switches[0] = 1; /*1 makes the output in kilograms and meters, - is cm and grams*/

	input.doy = GMT;                                 /* day of year */
	input.year = 0;                                  /* year, currently ignored */
	input.sec = (day - GMT)*24*3600;                 /* seconds in day (UT) */
	input.alt = alt;                                 /* altitude in kilometers */
	input.g_lat = g_lat;                             /* geodetic latitude */
	input.g_long = g_long;                           /* geodetic longitude */
	input.lst = input.sec / 3600. + g_long / 15.;    /* local apparent solar time (hours), see note below */
	input.f107A = f107A;                             /* 81 day average of F10.7 flux (centered on doy) */
	input.f107 = f107;                               /* daily F10.7 flux for previous day */
	input.ap = ap;                                   /* magnetic index(daily) */

	/* evaluate the atmosphere using the MSISE-00 model */
	gtd7(&input, &flags, &output);
	//gts7(&input, &flags, &output);

	/*in this case, i'm only concerned with the final atmospheric density*/
	return output.d[5];
}