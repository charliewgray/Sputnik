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
#include "LatLongBeta.h"

char SFUmatrix[1500][8][10];   // matrix with SFU predictions
char inputMatrix[19][40];
char ReboostsStringMatrix[50][40];
float Reboosts[50][8];
char Configmatrix[500][13][20];
int ConfigDates[500][4];
int ConfigBetas[10];
//int Atmosphere_Select;
double JD;
/* ------------------------------------------------------------------- */
/* -----------------------OPEN THE ATMOSPHERE FILE-------------------- */
/* ------------------------------------------------------------------- */
void ReboostConverter(){
	//Loops through the rest of the input file and creates reboost data based on the maneuver section
	int indiceVar = 0;
	int counter;

	char *pch;
	char line[80];

	//Start at the beginning of the manuever rows (they start at 21)
	strcpy(line, ReboostsStringMatrix[indiceVar]);

	//Loop through all the manuever rows and parse out the date and dV values, storing them in Reboosts
	while (line[0] != 0){
		counter = 0;

		//parse each reboost line, initializing pch as the first value separated by
		pch = strtok(line, "/");

		//Check that the value being pulled isn't null, since i need to operate on it more
		while (pch != NULL){
			Reboosts[indiceVar][counter] = atof(pch);

			//move on to the next value separated by /
			pch = strtok(NULL, "/");
			counter += 1;
		}

		//Record the Julian Date in the reboost matrix
		//JulianDate(Reboosts[indiceVar][0], Reboosts[indiceVar][1], Reboosts[indiceVar][2], Reboosts[indiceVar][3], Reboosts[indiceVar][4], Reboosts[indiceVar][5]);
		//Reboosts[indiceVar][7] = JulianDate(Reboosts[indiceVar][0], Reboosts[indiceVar][1], Reboosts[indiceVar][2], Reboosts[indiceVar][3], Reboosts[indiceVar][4], Reboosts[indiceVar][5]);

		indiceVar += 1;
		strcpy(line, ReboostsStringMatrix[indiceVar]);
	}
}

void openInputFile(int caseNumber){
	
	//Once i get to doing this multiple times, i can use the colVar variable as my loop counter!
	
	//open the output file
	FILE* f = fopen("InputFiles/input.txt", "r");
	if (f == NULL) { perror("Error opening input file!\n"); }

	char *pch;
	char line[800000];

	//Set indice variables
	int indiceVar = 0;
	int colVar = 0;

	//Check the input file line by line!
	while (fgets(line, 800000, f) != NULL){

		//parse each line, initializing pch as the first value separated by \t
		pch = strtok(line, "\t");

		//Check that the value being pulled isn't null, since i need to operate on it more
		while (pch != NULL){
			//Record into the input matrix if the previous value was "="
			if (colVar == caseNumber){
				strcpy(inputMatrix[indiceVar], pch);
				//if (indiceVar == 0) {
				//	objName = pch;
				//}
			}

			//move on to the next value separated by \t
			pch = strtok(NULL, "\t");
			colVar += 1;
		}
		//Once the end of the line is reached, move to the next row and reset the column values
		indiceVar += 1;
		colVar = 0;
	}

	//close the input file
	fclose(f);

	//Read through and store the Reboost data in the reboost matrix
	//ReboostConverter();
}

void openAtmosphereFile(int Atmosphere_Select){
	//This function reads in an atmosphere .DAS file and inputs the data into a matrix i can reference later
	FILE *fp;
	char *file_name = "InputFiles/";

	if (Atmosphere_Select == 3){
		file_name = "InputFiles/MSFC75.DAT";
		Atmosphere_Select -= 3;
	}else{
		file_name = "InputFiles/MSFC.DAT";
	}

	char *pch;
	char *monthName;
	char line[80];

	//Creete a matrix to store the solar flux data
	int currentRow = 0;
	int currentCol = 0;

	//open the chosen file in read only mode
	fp = fopen(file_name, "r"); // read mode

	//Error check in case the file can't be opened
	if (fp == NULL){
		perror("Error while opening the atmosphere file.\n");
		//exit();
	}

	//Check the input file line by line!
	while (fgets(line, 80, fp) != NULL){

		//print the entire line for reference first
		//printf("\n%s", line);

		//parse each line, initializing pch as the first value separated by \t
		int monthNumber = 0;
		pch = strtok(line, "\t");

		//Check that the value being pulled isn't null, since i need to operate on it more
		while (pch != NULL){
			//Record into the structure
			strcpy(SFUmatrix[currentRow][currentCol], pch);

			//move on to the next value separated by \t
			pch = strtok(NULL, "\t");
			currentCol += 1;
		}
		//Once the end of the line is reached, move to the next row and reset the column values
		currentRow += 1;
		currentCol = 0;
	}

	//Close the file
	fclose(fp);
}

void openConfigFile(){
	//This function reads in an Config file and inputs the data into a matrix i can reference later
	FILE *fp;
	char *file_name = "InputFiles/ConfigInput.txt";

	char *pch;
	char line[160];

	//Create a matrix to store the solar flux data
	int currentRow = 0;
	int currentCol = 0;
	int dateCol = 0;
	int betaCol = 0;

	//open the chosen file in read only mode
	fp = fopen(file_name, "r"); // read mode

	//Error check in case the file can't be opened
	if (fp == NULL){
		perror("Error while opening the Config file.\n");
		//exit();
	}

	//Check the input file line by line!
	while (fgets(line, 160, fp) != NULL){

		//parse each line, initializing pch as the first value separated by \t
		pch = strtok(line, "\t");

		//Check that the value being pulled isn't null, since i need to operate on it more
		while (pch != NULL){
			//Record into the structure
			strcpy(Configmatrix[currentRow][currentCol], pch);

			if(currentRow == 0){
				//Record the beta ranges on the first row
				if (betaCol > 0) {
					ConfigBetas[betaCol-1] = atoi(pch);
				}
				//move on to the next value separated by \t
				betaCol += 1;
			}

			//move on to the next value separated by \t
			pch = strtok(NULL, "\t");
			currentCol += 1;
		}
		//Once the end of the line is reached, move to the next row and reset the column values
		currentRow += 1;
		currentCol = 0;
	}

	//Split up the dates and record them into a new matrix (actual numbers instead of strings)
	int thisRow = 1;
	while (thisRow < currentRow) {
		//Split up the date into corresponding values
		dateCol = 0;
		char *test = Configmatrix[thisRow][0];
		pch = strtok(test, "/");
		
		while (pch != NULL){
			ConfigDates[thisRow][dateCol] = atoi(pch);

			//move on to the next value separated by \t
			pch = strtok(NULL, "/");
			dateCol += 1;
		}

		//Record the Julian date for this event
		//Only record the date if there is a data point...
		if (ConfigDates[thisRow][2] != 0) {
			ConfigDates[thisRow][3] = JulianDate(ConfigDates[thisRow][2], ConfigDates[thisRow][0], ConfigDates[thisRow][1], 0, 0, 0);
		}

		thisRow += 1;
	}

	//Close the file
	fclose(fp);
}