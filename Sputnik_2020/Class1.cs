using OrbitalElement_Tools; //THIS PULLS FROM THE SPUTNIK ACTIVE FOLDER, SO THE FILE MUST BE COPIED OVER
using System;
using System.Collections.Generic;
using System.Diagnostics.Eventing.Reader;
using System.IO;
using System.Linq;
using System.Management.Instrumentation;
using System.Reflection.Emit;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.ComTypes;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Schema;

namespace Sputnik_2020
{
    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.AutoDual)]
    public class Global
    {
        public static string find_directory()
        {
            string this_dir = Directory.GetCurrentDirectory();
            string[] dir_split = this_dir.Split(Path.DirectorySeparatorChar);

            //Initialize the new directory
            string new_dir = dir_split[0] + Path.DirectorySeparatorChar;

            int i = 1;
            while (i < 3)
            {
                new_dir = Path.Combine(new_dir, dir_split[i]);
                i++;
            }

            return new_dir;
        }

        public static string this_dir = find_directory();

        public static string repodir = Path.Combine(this_dir, "Documents", "source", "repos");
        public static string SpXStatusDir = Path.Combine(repodir, "SpaceXStatus");
        public static string Sputnikdir = Path.Combine(repodir, "Sputnik_2020");
        public static string inputdir = Path.Combine(Sputnikdir, "InputFiles");
        public static string outputdir = Path.Combine(repodir, "Output");
        public static string TLEdir = Path.Combine(repodir, "SatelliteTLEs");
    }
    public class Sputnik
    {
        public static void write_to_file(StringBuilder output, string savePath, string fileName)
        {
            //Save to a file so I can check
            using (var writer = new StreamWriter(savePath + fileName.Split('.')[0] + ".txt"))
            {
                writer.WriteLine(output);
            }
        }

        static double Sqrt(double input) { return Math.Pow(input, 0.5); }
        static double Pow(double x, double y) { return Math.Pow(x, y); }
        static double Range_2pi(double Angle)
        {
            //Force the RADIAN angle to fall between 0 and 2*Pi (0 to 360 degrees)
            //Create an integer value of the angle over 2Pi
            int IntegerVar = Convert.ToInt32(Angle / (2 * Math.PI));
            return Angle - 2 * Math.PI * IntegerVar;
        }
        static double find_target_TA(double burndV, double thrust, double mass, double preset_dt, double burnTAtarget, double period_m, double prev_TA, double[] x, double BN, double Current_SFU, double Current_AP, double JDN, DateTime thisDate)
        {
            //tried several smart solutions that predict TA in the future that just don't work, so it's brute fucking force
            //We're going to propagate forward until we find the TA we want, then record the TA target that's burntime/2 before that
            OE_Tools OETools = new OE_Tools();

            //find how many time steps the burn will take
            double burn_time = (Math.Abs(burndV) * mass) / (thrust * 9.82); //burn time in seconds
            double burn_time_steps = burn_time / Math.Abs(preset_dt); //number of time steps burn will take
            int steps_from_TA = Convert.ToInt32(burn_time_steps / 2);                   //divide by 2 so that we center the burn

            //This is the dumb solution - it's OK but assumes TA changes linearly which isn't true and ends up misplacing the burns sometimes
            double solution = burnTAtarget - (360.0 / period_m) * (burn_time / 60.0);
            //double solution = burnTAtarget - TA_delta * Math.Abs(burndV) * 2; //the burndV needs to be time increments
            while (solution > 360) { solution -= 360; }
            while (solution < 0) { solution += 360; }

            //propagate forward
            bool TAFound = false;
            double[] Osc_Elements = OETools.ECI_Vector_to_Elements(x);
            double[] TA_Array = new double[50000];
            TA_Array[0] = Osc_Elements[6];
            int iterate_Count = 1;
            double[] SunVector = new double[4] { 0, 0, 0, 0 };

            while (TAFound == false & iterate_Count < 5000)
            {
                //Find latitude and longitude so i can calculate density
                thisDate = thisDate.Add(TimeSpan.FromSeconds(preset_dt));
                double thisDay_Fraction = thisDate.Day + thisDate.Hour / 24.0 + thisDate.Minute / (60.0 * 24) + thisDate.Second / (24.0 * 3600);

                DateTime yearStart = new DateTime(Convert.ToInt32(thisDate.Year), 1, 1);
                double currentGMT = (thisDate - yearStart).TotalDays + 1;

                //Calculate Julian Date
                JDN = OETools.JulianDate(thisDate.Year, thisDate.Month, thisDay_Fraction);
                double[] LongLat = OETools.Long_Lat(JDN, thisDay_Fraction, x);

                //Find Densities
                double this_Density = Pow(1000, 3) * OETools.VBA_Density_Finder(currentGMT, LongLat[2], LongLat[3], LongLat[0], Current_SFU, Current_AP);

                x = RK4_5(preset_dt, BN, this_Density, SunVector, x);

                //Calculate new TA
                Osc_Elements = OETools.ECI_Vector_to_Elements(x);
                TA_Array[iterate_Count] = Osc_Elements[6];

                //Check whether we've crossed the target
                if (CrossCheck(Osc_Elements[6], prev_TA, burnTAtarget, preset_dt) && iterate_Count > steps_from_TA)
                {
                    TAFound = true;
                    solution = TA_Array[iterate_Count - steps_from_TA];
                }

                prev_TA = TA_Array[iterate_Count];
                iterate_Count++;
            }

            return solution;
        }
        static DateTime find_Burn_End_date(double deorbitBurnTime, DateTime thisDate, double dt)
        {
            DateTime burnEndDate = new DateTime(4000, 1, 1);
            TimeSpan duration = new TimeSpan(0, 0, 0, Convert.ToInt32(deorbitBurnTime));
            if (dt > 0)
            {
                burnEndDate = thisDate.Add(duration);
            }
            else
            {
                burnEndDate = thisDate.Subtract(duration);
            }
            return burnEndDate;
        }

        //The following DLL call only works if i dynamically change the directory before the density function is called in the code
        //[DllImport(@"C:\Users\cwgray1\Documents\source\repos\densityDLL\x64\Debug\densityDLL.dll")]
        //[DllImport("densityDLL.dll")]
        //public static extern double Density_Finder(double day, double Geodetic_Alt, double Lat, double Long, double f107A, double f107, double ap);
        static public double[] Accelerations(double BN, double rho, double[] SunVector, double[] x)
        {
            double Rx = x[0];
            double Ry = x[1];
            double Rz = x[2];
            double Vx = x[3];
            double Vy = x[4];
            double Vz = x[5];

            double We = 0.000072921156001; // Earth Rotation Rate in radians
            double Mu = 398600.5;          // Earth Standard Gravitational Parameter (MG) = km^3/s^2 
            //double MuS = 132712440018;     // Sun Standard Gravitational Parameter (MG) = km^3/s^2 
            double a_Earth = 6378.137;     // Earth equatorial radius 
            double[] J = new double[] { 0, 0, 0.001082635, -0.000002531, -0.00000160, -0.00000015, 0.00000057 };

            double Rmag = Sqrt(Pow(Rx, 2) + Pow(Ry, 2) + Pow(Rz, 2)); //magnitude of the position vector

            //**************Earth Gravitational Potential**********************************************
            double Eaccel = -Mu / Pow(Rmag, 3);
            double aEx = Eaccel * Rx;
            double aEy = Eaccel * Ry;
            double aEz = Eaccel * Rz;

            // Earth Harmonic effects: **************************************************************
            double[] aH = new double[3]; // , aHdetailed[3];
            //aH[0] = -1.5 * J[2] * Mu * Pow(a_Earth, 2) * Pow(Rmag, -4) * (1 - 5 * Pow(Rz / Rmag, 2)) * Rx / Rmag + -0.5 * J[3] * Mu * Pow(a_Earth, 3) * Pow(Rmag, -5) * (15 * Rz / Rmag - 35 * Pow(Rz / Rmag, 3)) * Rx / Rmag + 0.125 * J[4] * Mu * Pow(a_Earth, 4) * Pow(Rmag, -6) * (15 - 210 * Pow(Rz / Rmag, 2) + 315 * Pow(Rz / Rmag, 4)) * Rx / Rmag;// +aJ5x + aJ6x;
            //aH[1] = -1.5 * J[2] * Mu * Pow(a_Earth, 2) * Pow(Rmag, -4) * (1 - 5 * Pow(Rz / Rmag, 2)) * Ry / Rmag + -0.5 * J[3] * Mu * Pow(a_Earth, 3) * Pow(Rmag, -5) * (15 * Rz / Rmag - 35 * Pow(Rz / Rmag, 3)) * Ry / Rmag + 0.125 * J[4] * Mu * Pow(a_Earth, 4) * Pow(Rmag, -6) * (15 - 210 * Pow(Rz / Rmag, 2) + 315 * Pow(Rz / Rmag, 4)) * Ry / Rmag;// +aJ5y + aJ6y;
            //aH[2] = -1.5 * J[2] * Mu * Pow(a_Earth, 2) * Pow(Rmag, -4) * (3 - 5 * Pow(Rz / Rmag, 2)) * Rz / Rmag + -0.5 * J[3] * Mu * Pow(a_Earth, 3) * Pow(Rmag, -5) * (30 * Pow(Rz / Rmag, 2) - 35 * Pow(Rz / Rmag, 4) - 3) + 0.125 * J[4] * Mu * Pow(a_Earth, 4) * Pow(Rmag, -6) * (75 - 350 * Pow(Rz / Rmag, 2) + 315 * Pow(Rz / Rmag, 4)) * Rz / Rmag;// +aJ5z + aJ6z;

            double RxRmag = Rx / Rmag;
            double RyRmag = Ry / Rmag;
            double RzRmag = Rz / Rmag;

            //**************** J2 ********************************************************************
            double aJ2x = -1.5 * J[2] * Mu * Pow(a_Earth, 2) * Pow(Rmag, -4) * (1 - 5 * Pow(Rz / Rmag, 2)) * Rx / Rmag;
            double aJ2y = -1.5 * J[2] * Mu * Pow(a_Earth, 2) * Pow(Rmag, -4) * (1 - 5 * Pow(Rz / Rmag, 2)) * Ry / Rmag;
            double aJ2z = -1.5 * J[2] * Mu * Pow(a_Earth, 2) * Pow(Rmag, -4) * (3 - 5 * Pow(Rz / Rmag, 2)) * Rz / Rmag;

            //****************J3 * *******************************************************************
            double aJ3x = -0.5 * J[3] * Mu * Pow(a_Earth, 3) * Pow(Rmag, -5) * (15 * Rz / Rmag - 35 * Pow(Rz / Rmag, 3)) * Rx / Rmag;
            double aJ3y = -0.5 * J[3] * Mu * Pow(a_Earth, 3) * Pow(Rmag, -5) * (15 * Rz / Rmag - 35 * Pow(Rz / Rmag, 3)) * Ry / Rmag;
            double aJ3z = -0.5 * J[3] * Mu * Pow(a_Earth, 3) * Pow(Rmag, -5) * (30 * Pow(Rz / Rmag, 2) - 35 * Pow(Rz / Rmag, 4) - 3);

            //****************J4 * *******************************************************************
            double aJ4x = 0.125 * J[4] * Mu * Pow(a_Earth, 4) * Pow(Rmag, -6) * (15 - 210 * Pow(Rz / Rmag, 2) + 315 * Pow(Rz / Rmag, 4)) * Rx / Rmag;
            double aJ4y = 0.125 * J[4] * Mu * Pow(a_Earth, 4) * Pow(Rmag, -6) * (15 - 210 * Pow(Rz / Rmag, 2) + 315 * Pow(Rz / Rmag, 4)) * Ry / Rmag;
            double aJ4z = 0.125 * J[4] * Mu * Pow(a_Earth, 4) * Pow(Rmag, -6) * (75 - 350 * Pow(Rz / Rmag, 2) + 315 * Pow(Rz / Rmag, 4)) * Rz / Rmag;

            //****************J5 * *******************************************************************
            double aJ5x = 0.375 * J[5] * Mu * Pow(a_Earth, 5) * Pow(Rmag, -7) * (35 * RzRmag - 210 * Pow(RzRmag, 3) + 210 * Pow(RzRmag, 5)) * RxRmag;
            double aJ5y = 0.375 * J[5] * Mu * Pow(a_Earth, 5) * Pow(Rmag, -7) * (35 * RzRmag - 210 * Pow(RzRmag, 3) + 210 * Pow(RzRmag, 5)) * RyRmag;
            double aJ5z = 0.375 * J[5] * Mu * Pow(a_Earth, 5) * Pow(Rmag, -7) * Pow(RzRmag, 2) * (105 - 315 * Pow(RzRmag, 2) + 231 * Pow(RzRmag, 4)) - 1.875 * J[5] * Mu * Pow(a_Earth, 5) * Pow(Rmag, -7);

            //****************J6 * *******************************************************************
            double aJ6x = -0.0625 * J[6] * Mu * Pow(a_Earth, 6) * Pow(Rmag, -8) * (35 - 945 * Pow(RzRmag, 2) + 3465 * Pow(RzRmag, 4) - 3003 * Pow(RzRmag, 6)) * RxRmag;
            double aJ6y = -0.0625 * J[6] * Mu * Pow(a_Earth, 6) * Pow(Rmag, -8) * (35 - 945 * Pow(RzRmag, 2) + 3465 * Pow(RzRmag, 4) - 3003 * Pow(RzRmag, 6)) * RyRmag;
            double aJ6z = -0.0625 * J[6] * Mu * Pow(a_Earth, 6) * Pow(Rmag, -8) * (245 - 2205 * Pow(RzRmag, 2) + 4851 * Pow(RzRmag, 4) - 3003 * Pow(RzRmag, 6)) * RzRmag;

            //Combine these effects into one Earth Harmonics accelleration vector
            aH[0] = aJ2x + aJ3x + aJ4x + aJ5x + aJ6x;
            aH[1] = aJ2y + aJ3y + aJ4y + aJ5y + aJ6y;
            aH[2] = aJ2z + aJ3z + aJ4z + aJ5z + aJ6z;

            //****************Acceleration due to Drag***********************************************************
            double coef = (-1 * rho) / (2.0 * BN * Pow(1000, 2));
            double VxRel = Vx + We * Ry;
            double VyRel = Vy + We * Rx;
            double VzRel = Vz;
            double Vmag = Sqrt(Pow(VxRel, 2) + Pow(VyRel, 2) + Pow(VzRel, 2));

            double aDx = coef * Pow(Vmag, 2) * (VxRel / Vmag);
            double aDy = coef * Pow(Vmag, 2) * (VyRel / Vmag);
            double aDz = coef * Pow(Vmag, 2) * (VzRel / Vmag);

            ////**************Sun Gravitational Potential************************************************
            ////the current sun vector is w.r.t.the earth, so change that to the satellite
            //double RsunX = SunVector[0] * SunVector[3] - Rx;
            //double RsunY = SunVector[1] * SunVector[3] - Ry;
            //double RsunZ = SunVector[2] * SunVector[3] - Rz;
            //double rSun = Sqrt(Pow(RsunX, 2) + Pow(RsunY, 2) + Pow(RsunZ, 2));

            ////Find the ratio of the distance from the central body (earth) to the third body (sun) over the distance from the satellite to the third body (sun)
            //double Bk = SunVector[3] / rSun - 1;
            //double bk = 3 * Bk + 3 * Pow(Bk, 2) + Pow(Bk, 3);

            //double Saccel = -MuS / Pow(SunVector[3], 3);
            //double aSx = Saccel * (SunVector[0] - bk * RsunX / rSun);
            //double aSy = Saccel * (SunVector[1] - bk * RsunY / rSun);
            //double aSz = Saccel * (SunVector[2] - bk * RsunZ / rSun);

            //Add each accelleration effect to the accelleration vector -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
            double[] a = new double[3];
            a[0] = aEx + aH[0] + aDx;// + aSx;
            a[1] = aEy + aH[1] + aDy;// + aSy;
            a[2] = aEz + aH[2] + aDz;// + aSz;

            //////return a;
            double[] xdot = new double[6];
            xdot[0] = Vx;
            xdot[1] = Vy;
            xdot[2] = Vz;
            xdot[3] = a[0];
            xdot[4] = a[1];
            xdot[5] = a[2];
            return xdot;
        }
        static public Boolean CrossCheck(double thisValue, double previousValue, double targetValue, double dt)
        {
            //if (Math.Abs(thisValue - previousValue) > 60) {
            //    if (dt < 0){
            //        previousValue += 360;
            //    } else {
            //        thisValue += 360;
            //    }
            //}

            if (Math.Abs(previousValue - thisValue) > 60)
            {
                if (thisValue < previousValue)
                {
                    thisValue += 360;
                }
                else
                {
                    previousValue += 360;
                }
            }

            if (thisValue >= targetValue && previousValue < targetValue)
            {
                return true;
            }
            else if (previousValue >= targetValue && thisValue < targetValue)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        static public double[] RK4_5(double h, double BN, double rho, double[] SunVector, double[] x)
        {
            //4(5) Fehlberg Runge Kunge - Kutta Method to solve fundamental equation
            //Define k1 * ***********************************************************
            double[] y = new double[6] { 0, 0, 0, 0, 0, 0 };
            for (int j = 0; j < 6; j++) { y[j] = x[j]; }
            double[] k1 = Accelerations(BN, rho, SunVector, y);
            for (int j = 0; j < 6; j++) { k1[j] = h * k1[j]; }

            //Define k2 ************************************************************
            for (int j = 0; j < 6; j++) { y[j] = x[j] + k1[j] * 1.0 / 4.0; }
            double[] k2 = Accelerations(BN, rho, SunVector, y);
            for (int j = 0; j < 6; j++) { k2[j] = h * k2[j]; }

            //Define k3 ************************************************************
            for (int j = 0; j < 6; j++) { y[j] = x[j] + k1[j] * 3.0 / 32.0 + k2[j] * 9.0 / 32.0; }
            double[] k3 = Accelerations(BN, rho, SunVector, y);
            for (int j = 0; j < 6; j++) { k3[j] = h * k3[j]; }

            //Define k4 ************************************************************
            for (int j = 0; j < 6; j++) { y[j] = x[j] + k1[j] * 1932.0 / 2197.0 - k2[j] * 7200.0 / 2197.0 + k3[j] * 7296.0 / 2197.0; }
            double[] k4 = Accelerations(BN, rho, SunVector, y);
            for (int j = 0; j < 6; j++) { k4[j] = h * k4[j]; }

            //Define k5 ************************************************************
            for (int j = 0; j < 6; j++) { y[j] = x[j] + k1[j] * 439.0 / 216.0 - k2[j] * 8 + k3[j] * 3680.0 / 513.0 - k4[j] * 845.0 / 4104.0; }
            double[] k5 = Accelerations(BN, rho, SunVector, y);
            for (int j = 0; j < 6; j++) { k5[j] = h * k5[j]; }

            //Define k6 ************************************************************
            for (int j = 0; j < 6; j++) { y[j] = x[j] - k1[j] * 8.0 / 27.0 + k2[j] * 2 - k3[j] * 3544.0 / 2565.0 + k4[j] * 1859.0 / 4104.0 - k5[j] * 11.0 / 40.0; }
            double[] k6 = Accelerations(BN, rho, SunVector, y);
            for (int j = 0; j < 6; j++) { k6[j] = h * k6[j]; }

            //Define the new position and velocity components according to the k factors
            double[] new_x = new double[6];
            for (int j = 0; j < 6; j++)
            {
                new_x[j] = x[j] + 16.0 / 135.0 * k1[j] + 6656.0 / 12825.0 * k3[j] + 28561.0 / 56430.0 * k4[j] - 9.0 / 50.0 * k5[j] + 2.0 / 55.0 * k6[j];
            }

            return new_x;
        }
        static public double[] SunVectorCalc(double JD)
        {
            //******************************************************************************************
            //*****************Angles associated with the position of the sun **************************
            //******************************************************************************************
            double[] SunVector = new double[4] { 0, 0, 0, 0 };
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
            double Lsun_Ecl = Lsun + 0.03341607 * Math.Sin(Msun) + 0.00034894 * Math.Sin(2 * Msun) - 0.0001134 - 0.0000206 * Math.Sin(RAANsun);  //Compute: Ecliptic longitude of sun
            Lsun_Ecl = Range_2pi(Lsun_Ecl);

            //obliquity of the ecliptic at JD 2451545.0 (deg) + rate of obliquity (deg/day)
            double Obl_Ecl = 0.4090928 - 6.2140 * Pow(10, -9) * n + 0.0000396 * Math.Cos(RAANsun);  //Compute: Obliquity of ecliptic
            Obl_Ecl = Range_2pi(Obl_Ecl);

            //Create the sun (unit) vector in ECI coordinates
            SunVector[0] = Math.Cos(Lsun_Ecl);
            SunVector[1] = Math.Cos(Obl_Ecl) * Math.Sin(Lsun_Ecl);
            SunVector[2] = Math.Sin(Obl_Ecl) * Math.Sin(Lsun_Ecl);

            //Now to find the radius magnitude (distance from earth to sun)
            double T = n / 36525.0;
            double eEarth = 0.016708617 - 0.000042037 * T - 0.0000001236 * Pow(T, 2);  //Earth orbit eccentricity
            double TAEarth = Msun + Lsun_Ecl - Lsun;                                   //true anomaly of earth
            double aEarth = 149.6 * Pow(10, 6);                                        //Earth semimajor axis

            //R = a * (1/e^2) / (1 + e*cos(TA))
            double RSun = aEarth * (1 - Pow(eEarth, 2)) / (1 + eEarth * Math.Cos(TAEarth));
            SunVector[3] = RSun;

            return SunVector;
        }
        static public double[] Reboost(double dt, double dV_ratio, double[] x)
        {
            //Apply reboost dV to the velocity vector
            double dV = dV_ratio * Math.Abs(dt);    //find the dV total for this specific time increment

            //Calculate magnitude of dV
            double VMag = Math.Sqrt(Math.Pow(x[3], 2) + Math.Pow(x[4], 2) + Math.Pow(x[5], 2));

            //Create the identity matrix
            double Vx = x[3] / VMag;
            double Vy = x[4] / VMag;
            double Vz = x[5] / VMag;

            //add dV to the velocity magnitude
            VMag = VMag + dV / 1000;

            //add magnitude back into the identity matrix
            x[3] = Vx * VMag;
            x[4] = Vy * VMag;
            x[5] = Vz * VMag;
            return x;
        }
        static public double[,] Read_Config(DateTime[] DateArray, string this_dir)
        {
            //string config_dir = Path.Combine(this_dir, "Documents\\source\\repos\\Sputnik_2020\\InputFiles\\config.txt");
            string config_dir = Path.Combine(this_dir,"config.txt");

            string[] configArray = System.IO.File.ReadAllLines(config_dir);
            //string[] configArray = System.IO.File.ReadAllLines(@"C:\Users\cwgray1\Documents\source\repos\Sputnik_2020\InputFiles\config.txt");

            int configCount = configArray.Length;
            double[,] Config_Array = new double[configCount, 5];

            int config_index = 0;
            while (config_index < configCount)
            {
                string this_Date = configArray[config_index];
                int config_year = Convert.ToInt32(configArray[config_index].Split(',')[0].Split('/')[2]);
                int config_day = Convert.ToInt32(configArray[config_index].Split(',')[0].Split('/')[1]);
                int config_month = Convert.ToInt32(configArray[config_index].Split(',')[0].Split('/')[0]);//.Split('\"')[1]);

                DateTime ConfigDate = new DateTime(config_year, config_month, config_day);
                //Add time in fractions of a day
                ConfigDate = ConfigDate.Add(TimeSpan.FromDays(Convert.ToDouble(configArray[config_index].Split(',')[0].Split('/')[3])));

                double ConfigBN = Convert.ToDouble(configArray[config_index].Split(',')[1]);

                string reboost_Info = configArray[config_index].Split(',')[2].Split('\"')[0];
                double ConfigdV = Convert.ToDouble(reboost_Info.Split('/')[0]);
                double ConfigReboost_Time = Convert.ToDouble(reboost_Info.Split('/')[1]);
                double ConfigThrust = Convert.ToDouble(reboost_Info.Split('/')[2]);

                DateArray[config_index] = ConfigDate;
                Config_Array[config_index, 0] = ConfigBN;
                Config_Array[config_index, 1] = ConfigdV;
                Config_Array[config_index, 2] = ConfigReboost_Time;
                Config_Array[config_index, 3] = ConfigThrust;
                config_index++;
            }
            return Config_Array;
        }
        static public double Smart_Average(double[] inputArray)
        {
            //This function is a simple average, but it doesn't take 0's into account
            double totalSum = 0;
            double Count = 0;

            for (int i = 0; i < inputArray.Length; i++)
            {
                if (inputArray[i] != 0)
                {
                    totalSum += inputArray[i];
                    Count += 1;
                }
            }

            return totalSum / Count;
        }

        static public DateTime makeStartDate(string year, double EPOCH)
        {
            //Converts year and EPOCH (or just a year input that's a date) into an actual dateTime date I can work with
            DateTime startDate;
            if (year.IndexOf("/") != -1)
            {
                //it's a date
                startDate = Convert.ToDateTime(year);
            }
            else
            {
                double yearDub = Convert.ToDouble(year);
                if (yearDub > 10000)
                {
                    startDate = new DateTime(Convert.ToInt32(1900), 1, 1);
                    startDate = startDate.Add(TimeSpan.FromDays(yearDub - 2));
                }
                else
                {
                    startDate = new DateTime(Convert.ToInt32(year), 1, 1);
                    startDate = startDate.Add(TimeSpan.FromHours(EPOCH * 24 - 24));
                }
            }
            return startDate;
        }

        public static void Propagator(string[] inputArray)
        {
            OE_Tools OETools = new OE_Tools();
            //string this_dir = Regex.Split(Directory.GetCurrentDirectory(), @"Documents")[0];

            double XKMPER = 6378.135; //km per AE
            double GM = 398600.435507; //Earth Gravitational Paramater
            double We = 0.000072921156001; //Earth Rotation in radians

            //Grab Satellite Name and Write update to console
            string SatID = inputArray[0];

            //Grab Epoch
            string year = inputArray[1];
            double EPOCH = Convert.ToDouble(inputArray[2]);

            //Get Oscillating Elements (should already be converted to Osc from Mean TLE elements during the SpacePull process
            double A = Convert.ToDouble(inputArray[3]);
            double e = Convert.ToDouble(inputArray[4]);                 //eccentricity
            double I = Convert.ToDouble(inputArray[5]);                 //Inclination
            double Wp = Convert.ToDouble(inputArray[6]);                //Argument of Perigee
            double Ra = Convert.ToDouble(inputArray[7]);                //Right Ascension
            double MA = Convert.ToDouble(inputArray[8]);                //Mean Anomaly
            //double Mean_Motion = Convert.ToDouble(inputArray[9].Split('\t')[i]);     //Mean Motion (in orbits per day)
            //double BSTAR = Convert.ToDouble(inputArray[10].Split('\t')[i]);          //B*
            //double XNDT2O = Convert.ToDouble(inputArray[11].Split('\t')[i]);         //first derivative of mean motion
            //double XNDD6O = Convert.ToDouble(inputArray[12].Split('\t')[i]);         //second derivative of mean motion
            double BN = Convert.ToDouble(inputArray[13]);
            int Atmosphere_Select = Convert.ToInt32(inputArray[14]);
            double preset_dt = Convert.ToDouble(inputArray[15]);              //time increment in seconds
            double runLength = Convert.ToDouble(inputArray[16]) * 86400;      //total run Length in seconds
            double Output_dt = Convert.ToDouble(inputArray[17]);             //output Increment in seconds

            //Work to get the date + time manageable
            DateTime startDate = makeStartDate(year, EPOCH);

            //Check if the burn trigger is an altitude or date
            int inputBurnCount = 0;
            string deboostDateInput = inputArray[18];
            string modifiedDeboostDateInput = deboostDateInput.Split(';')[0];
            double burnAlt = 0;
            DateTime specifiedBurnDate = new DateTime(4000, 1, 1);

            //Determine if there are multiple inputs
            if (modifiedDeboostDateInput.IndexOf("/") != -1 && deboostDateInput.IndexOf(";") != -1)
            {
                //if (deboostDateInput.IndexOf("/") != -1 && deboostDateInput.IndexOf(";") != -1) {
                while (Convert.ToDateTime(modifiedDeboostDateInput) < startDate)
                {
                    inputBurnCount++;
                    modifiedDeboostDateInput = deboostDateInput.Split(';')[inputBurnCount];
                }
                modifiedDeboostDateInput = deboostDateInput.Split(';')[inputBurnCount];
            }

            //Check if the dV is a value or "M"
            double burndV = 0;
            int inputBurnMax = 0;
            string burndVinput = inputArray[19];
            if (burndVinput != "M")
            {
                //it's a value or series of values
                if (burndVinput.IndexOf(";") != -1)
                {
                    inputBurnMax = burndVinput.Split(';').Length;
                    burndV = Convert.ToDouble(burndVinput.Split(';')[inputBurnCount]);
                    inputBurnCount++;
                }
                else
                {
                    burndV = Convert.ToDouble(burndVinput);
                }
            }

            //Determine whether the input is a date or altitude
            if (modifiedDeboostDateInput.IndexOf("/") != -1)
            {
                //it's a date
                specifiedBurnDate = Convert.ToDateTime(modifiedDeboostDateInput);
            }
            else
            {
                //it's an altitude
                burnAlt = Convert.ToDouble(modifiedDeboostDateInput);
            }

            double burnTAtarget = Convert.ToDouble(inputArray[20]);

            //Check if the breakup altitude is an altitude or actually a date to stop reboosting
            string breakupAltInput = inputArray[21];
            double breakupAlt = 0;
            //DateTime specifiedSTOPDate = new DateTime(4000, 1, 1);
            DateTime BurnSTOPDate = new DateTime(4000, 1, 1);

            if (breakupAltInput.IndexOf("/") != -1)
            {
                //it's a date
                BurnSTOPDate = Convert.ToDateTime(breakupAltInput);
            }
            else
            {
                //it's an altitude
                breakupAlt = Convert.ToDouble(breakupAltInput);
            }

            double breakupBN = Convert.ToDouble(inputArray[22]);
            double mass = Convert.ToDouble(inputArray[23]);
            double thrust = Convert.ToDouble(inputArray[24]);
            if (mass == 0) { mass = 480000; }
            if (thrust == 0) { thrust = 814.6; }

            int config_Check = Convert.ToInt32(inputArray[25]);

            //Make some parameters such that BN can be updated
            double Cd = 2.07;
            double massAreaRatio = BN * Cd;

            //Create a boolean to check whether i want to output every LAN, or only  on specified output increments. 
            //If the output_dt has a decimal (86400.1, for example), the code will also output every time Lat crosses 0
            bool LAN_Output;
            if (Output_dt - Math.Round(Output_dt) > 0)
            {
                LAN_Output = true;
                Output_dt = Math.Round(Output_dt);
            }
            else
            {
                LAN_Output = false;
            }

            //dt can change if config file requires it, but start off with the preset value from the input file
            double dt = preset_dt;

            //initialize the config variables
            int config_index = 0;

            DateTime[] DateArray = new DateTime[1000];
            double[,] config_Array = new double[1000, 4];
            DateTime configDate = DateArray[0];

            //Convert Mean motions to radians per second - don't really use these anymore
            //double n = Mean_Motion * 2 * Math.PI / 86400; //1440;
            //XNDT2O = XNDT2O * 2 * Math.PI / Pow(86400,2); // (1440 * 1440);
            //XNDD6O = XNDD6O * 2 * Math.PI / 86400; // 1440;

            //Pull in the atmosphere prediction file so i can assign SFU and AP values for density lookup
            string solarPredictString = "MSFC.DAT";
            if (Atmosphere_Select == 3)
            {
                solarPredictString = "MSFC75.DAT";
            }

            // Pull in the SFU data as convert to a List for formatting purposes
            string[] SFUFile = System.IO.File.ReadAllLines(@"C:\strap\data\" + solarPredictString);
            string[] oldSFUFile = System.IO.File.ReadAllLines(@"C:\strap\data\MSFC_Pre1998.txt");
            var SFUList = new List<string>(SFUFile);
            var oldSFUList = new List<string>(oldSFUFile);
            SFUList.RemoveAt(0); //Delete the first two lines
            SFUList.RemoveAt(0); //Delete the first two lines

            oldSFUList.AddRange(SFUList);

            //Convert back to array for speed
            string[] SFUArray = oldSFUList.ToArray();
            SFUList = null; //Clear the List to save memory
            oldSFUList = null;//Clear the List to save memory
            string[] SFU_Start_Row = SFUArray[2].Split('\t');
            int SFU_Start_Year = Int16.Parse(SFU_Start_Row[0]);

            //If this is a case with BN changes and reboosts, pull in the config file
            if (config_Check == 1)
            {
                //Read the config input file (this also updates the DateArray)
                config_Array = Read_Config(DateArray, Global.inputdir);

                //Initialize the configuration index (find the starting point)
                if (dt > 0)
                {
                    //if we're going forward
                    while (configDate <= startDate)
                    {
                        config_index++;
                        configDate = DateArray[config_index];
                    }
                }
                else
                {
                    //if we're going backward, 
                    configDate = DateArray[config_index];
                }
            }
            else
            {
                configDate = new DateTime(4000, 1, 1);
            }

            //DRIVE the program with the correct dt values to ensure i get the correct output increments
            double current_dt = 0;
            double Output_Check = 0;
            double[] decay_Array = new double[10];

            double[] x = OETools.Elements_To_Vector_ECI(A, e, I, Ra, Wp, MA);
            double R = Sqrt(Pow(x[0], 2) + Pow(x[1], 2) + Pow(x[2], 2));
            double V = Sqrt(Pow(x[3], 2) + Pow(x[4], 2) + Pow(x[5], 2));

            //Set up a status bar
            Console.WriteLine("Running Propagation for " + SatID + ", Atmosphere: " + Atmosphere_Select + ", BN:" + BN);
            Console.WriteLine("|___________________________________________________________|");
            Console.Write("|");
            double statusBarLen = runLength / 60;
            double currentStatus = statusBarLen;

            //Initialize the output string (header row)
            string FileName = "A" + Atmosphere_Select + "-Out_" + SatID + ".txt";

            //Resets the clean string
            StringBuilder cleanResponse = new StringBuilder("SatID,year,day,A,e,I,Wp,Ra,MA,TA,SFU,AP,Density,BN,Hsma,Ha,Hp,Long,Lat,Beta,Date,SpatialDen,reboostdV,Daily_Orbit,LAN,Time Since Start (days), Distance Travelled (km), Avg Population, Nation ID,AirPop,ShipPop,,,");

            //Add altitude ranges to the header row for spatial density calculations
            for (int j = 0; j < 1; j++)
            {
                double altBin = 0;
                int altStepSize = 10;
                while (altBin < 650)
                {
                    cleanResponse.Append(altBin + ",");
                    altBin += altStepSize;
                }
                cleanResponse.Append(",");
                cleanResponse.Append("x1,x2,x3,v1,v2,v3,,Rmag,Vmag");
            }
            cleanResponse.Append("\r\n");

            string writeDir = Path.Combine(Global.outputdir, FileName);

            System.IO.File.WriteAllText(writeDir, cleanResponse.ToString());
            cleanResponse = new StringBuilder();

            double Hp = A * (1 - e) - XKMPER;
            double Ha = A * (1 + e) - XKMPER;
            double Inv_SMA = A;
            double Inv_hSMA = Inv_SMA - XKMPER;
            bool special_output = false;
            double reported_dV;
            double dV_ratio = 0;
            int decay_count = 0;
            double prev_hSMA = Inv_hSMA;
            double targetTA = 500;
            double prev_TA = 500;
            double prev_dt = 0;
            double previous_Lat = 0;
            double previous_Lon = 0;
            double LAN = 0;
            DateTime burnEndDate = new DateTime(4000, 1, 1);
            if (dt < 0) { burnEndDate = new DateTime(1000, 1, 1); }
            int burnFlag = 0;
            bool TA_Check = false;

            //initialize Time and Lat/Long values
            int DailyOrbit = 0;
            DateTime thisDate = startDate;
            double thisDay_Fraction = thisDate.Day + thisDate.Hour / 24.0 + thisDate.Minute / (60.0 * 24) + thisDate.Second / (24.0 * 3600) + thisDate.Millisecond / (24.0 * 3600000);
            double JDN = OETools.JulianDate(thisDate.Year, thisDate.Month, thisDay_Fraction);
            double[] LongLat = OETools.Long_Lat(JDN, thisDay_Fraction, x);
            int startBurn = 0;

            //Here's the loop
            while (Math.Abs(current_dt) < runLength)
            {
                //Recalculate values based on new x
                double[] Osc_Elements = OETools.ECI_Vector_to_Elements(x);

                //Check if the object has deorbited, break from this loop if so
                //if (R - 6378.134 < 50 || Inv_hSMA < 0) {
                if (R - XKMPER < 50 || Osc_Elements[0] < 0)
                {
                    break;
                }
                else if (dt < 0)
                {
                    if (Inv_hSMA >= 416)
                    {
                        break;
                    }
                }

                Hp = Osc_Elements[0] * (1 - Osc_Elements[1]) - XKMPER;
                Ha = Osc_Elements[0] * (1 + Osc_Elements[1]) - XKMPER;

                R = Sqrt(Pow(x[0], 2) + Pow(x[1], 2) + Pow(x[2], 2));
                V = Sqrt(Pow(x[3], 2) + Pow(x[4], 2) + Pow(x[5], 2));
                //double Geo_Alt = R - XKMPER;

                thisDate = startDate.Add(TimeSpan.FromSeconds(current_dt));
                thisDay_Fraction = thisDate.Day + thisDate.Hour / 24.0 + thisDate.Minute / (60.0 * 24) + thisDate.Second / (24.0 * 3600);

                DateTime yearStart = new DateTime(Convert.ToInt32(thisDate.Year), 1, 1);
                double currentGMT = (thisDate - yearStart).TotalDays + 1;

                //Calculate Julian Date
                JDN = OETools.JulianDate(thisDate.Year, thisDate.Month, thisDay_Fraction);

                //When i need to find SFU for a specific year/month combination, this is the algorithm
                int SFU_Row = 12 * (thisDate.Year - SFU_Start_Year) + thisDate.Month - 1;

                //check if i'm past the end of the array
                while (SFU_Row >= SFUArray.Length)
                {
                    //If we're past the end of existing data, just go back ~10 years and redo the previous solar cycle
                    SFU_Row = SFU_Row - 129;
                }
                string[] SFU_Current_Row = SFUArray[SFU_Row].Split('\t');// SFUList[SFU_Row].Split('\t');
                double Current_SFU = 0;
                double Current_AP = 0;

                //In some cases i just need to set the SFU - I can just input an SFU value instead of an index for M95/M75/etc.
                if (Atmosphere_Select > 3)
                {
                    Current_SFU = Atmosphere_Select;
                    Current_AP = Atmosphere_Select / 10;
                }
                else
                {
                    if (Atmosphere_Select == 3)
                    {
                        Current_SFU = Convert.ToDouble(SFU_Current_Row[2]);
                        Current_AP = Convert.ToDouble(SFU_Current_Row[5]);
                    }
                    else
                    {
                        Current_SFU = Convert.ToDouble(SFU_Current_Row[Atmosphere_Select + 2]);
                        Current_AP = Convert.ToDouble(SFU_Current_Row[Atmosphere_Select + 5]);
                    }
                }

                //Apply any dV as necessary
                reported_dV = 0;
                if (dV_ratio != 0)
                {
                    special_output = true; //output this time increment since it's a reboost
                    x = Reboost(dt, dV_ratio, x);
                    reported_dV = dV_ratio * Math.Abs(dt);
                }

                //Calculate Invariant Elements
                double[] Inv_Elements = OETools.Calculate_Invariant_Elements_ECI(x);

                Inv_SMA = Inv_Elements[0];
                Inv_hSMA = Inv_SMA - XKMPER;
                double Inv_e = Inv_Elements[1];
                //double Inv_i = Inv_Elements[2];
                //double Inv_Wp = Inv_Elements[3];
                double Inv_Ha = Inv_SMA * (1 + Inv_e) - XKMPER;
                double Inv_Hp = Inv_SMA * (1 - Inv_e) - XKMPER;
                double period = 2 * Math.PI * Math.Sqrt(Math.Pow(Inv_SMA, 3) / (GM));
                double period_m = period / 60.0;

                //Calculate an estimated date to stop burning to hit target Stop Date
                double decay_Rate = Smart_Average(decay_Array);

                //Create a reboost or deorbit burn from the input file (not from the config file)
                if (burnFlag == 0 && thisDate.Date < BurnSTOPDate.Date)
                { //Just to make sure i'm not already performing a burn
                  //Single Burns only when dV is specified
                    if (burndV != 0)
                    {
                        //Single burns based on a target altitude 
                        if (Inv_Hp <= burnAlt)
                        {
                            startBurn = 1;
                        }
                        //Single burn based on target date
                        //else if (dt > 0 && thisDate.date >= specifiedBurnDate.date || dt < 0 && thisDate.date <= specifiedBurnDate.date)
                        else if (dt > 0 && thisDate >= specifiedBurnDate || dt < 0 && thisDate <= specifiedBurnDate)
                        {
                            startBurn = 1;
                        }
                    }
                    else if (burndVinput == "M" && current_dt >= 86400 * 12)
                    {
                        //Average the last X day's worth of decay to estimate how high to reboost
                        //If we're within 30 days of the stop date, only reboost enough to hit the target altitude ON the stop date
                        double target_duration = (BurnSTOPDate - thisDate).TotalDays;
                        if (target_duration > 30)
                        {
                            target_duration = 30;
                        }
                        double target_Alt_Gain = decay_Rate * target_duration * 86400; //total estimated decay over 30 days

                        //Maintenance burns based on a target altitude 
                        if (dt > 0 && Inv_hSMA <= burnAlt || dt < 0 && Inv_hSMA >= (burnAlt - target_Alt_Gain))
                        {
                            startBurn = 1;
                            //Average the last X day's worth of decay to estimate how high to reboost
                            if (target_Alt_Gain != 0)
                            {
                                burndV = target_Alt_Gain / 1.852;
                            }
                        }
                        //Maintenance burns based on target date
                        else if (thisDate.Date == specifiedBurnDate.Date)
                        {
                            //Set up the next burn
                            specifiedBurnDate = thisDate.AddDays(30);
                            startBurn = 1;

                            //Average the last X day's worth of decay to estimate how high to reboost
                            burndV = target_Alt_Gain / 1.852;
                        }
                    }
                }

                //If flagged above, turn on the burn
                if (startBurn == 1)
                {
                    //Only go through the trouble of checking for the target TA once,
                    //after that we're just here to see if we've actually found the target TA
                    if (targetTA == 500)
                    {
                        targetTA = find_target_TA(burndV, thrust, mass, preset_dt, burnTAtarget, period_m, prev_TA, x, BN, Current_SFU, Current_AP, JDN, thisDate);
                    }

                    TA_Check = CrossCheck(Osc_Elements[6], prev_TA, targetTA, dt);

                    //We now know we're at the right date/altitude to burn, but have to check TA too
                    if (TA_Check && burndV != 0)
                    {
                        //turn on the burn and set timer to end the burn
                        double deorbitBurnTime = (Math.Abs(burndV) * mass) / (thrust * 9.82); // Math.Abs(60.0 * burndV);
                        dV_ratio = burndV / deorbitBurnTime;

                        // If the burn is > 250 m/s, it's treated as an ALTITUDE instead
                        if (burndV > 100)
                        {
                            burndV = (burndV - Inv_hSMA) / 1.852;

                            deorbitBurnTime = (Math.Abs(burndV) * mass) / (thrust * 9.82); // Math.Abs(60.0 * burndV);
                            dV_ratio = burndV / deorbitBurnTime;
                        }

                        burnEndDate = find_Burn_End_date(deorbitBurnTime, thisDate, dt);
                        configDate = burnEndDate; //this should make another output on exactly the burn stop time

                        burndV = 0;
                        burnFlag = 1;
                        startBurn = 0;
                        targetTA = 500;
                    }
                }

                //Turn off the burn if we're past the burn end date
                if (dt > 0 && thisDate >= burnEndDate || dt < 0 && thisDate <= burnEndDate)
                {
                    dV_ratio = 0;
                    burnFlag = 0;

                    //Incremenet to the next burn
                    if (inputBurnCount < inputBurnMax)
                    {
                        burndV = Convert.ToDouble(burndVinput.Split(';')[inputBurnCount]);
                        if (deboostDateInput.IndexOf("/") != -1)
                        {
                            specifiedBurnDate = Convert.ToDateTime(deboostDateInput.Split(';')[inputBurnCount]);
                        }
                        else
                        {
                            burnAlt = Convert.ToDouble(deboostDateInput.Split(';')[inputBurnCount]);
                        }

                        inputBurnCount++;

                        //Also re-initialize burn end date
                        burnEndDate = new DateTime(4000, 1, 1);
                        if (dt < 0) { burnEndDate = new DateTime(1000, 1, 1); }
                    }
                }

                //Check if Config needs updating or new reboost dV needs applying
                if (thisDate == configDate)
                {
                    double ConfigdV = 0; // config_Array[config_index, 1];
                    double ConfigReboost_Time = 0; // config_Array[config_index, 2];
                    double ConfigThrust = 0; //config_Array[config_index, 3];

                    //Apply the config info
                    if (config_index < config_Array.GetLength(0))
                    {
                        if (config_Array[config_index, 0] != 0)
                        {
                            BN = config_Array[config_index, 0];
                            massAreaRatio = BN * 2.07;
                        }

                        ConfigdV = config_Array[config_index, 1];
                        ConfigReboost_Time = config_Array[config_index, 2];
                        ConfigThrust = config_Array[config_index, 3];

                        if (ConfigdV > 0)
                        {
                            burnEndDate = find_Burn_End_date(ConfigReboost_Time, thisDate, dt);
                        }
                    }

                    //Turn off dV if this is the end of the burn
                    if (ConfigThrust == 0)
                    {
                        dV_ratio = 0;
                    }
                    else if (ConfigThrust == 1)
                    {
                        //in this case the dV is defined specifically
                        dV_ratio = ConfigdV / ConfigReboost_Time;
                    }
                    else if (ConfigThrust == 2)
                    {
                        //in this case the target altitude is defined rather than dV, so i need to find target dV
                        //dV_ratio = ((ConfigdV - (Inv_SMA - 6378.134)) / 1.852) / dt; //this method reboosts directly to the target altitude immediately

                        //in this case, we're actually bouncing off the target altitude - i.e. i need to find the 
                        //altitude of the next target reboost date and the decay rate, then reboost so that i'm 
                        //at the target altitude on the target date
                        int k = config_index + 2;
                        while (config_Array[k, 1] == 0) { k++; }
                        DateTime nextReboost = DateArray[k];
                        //if the next reboost is an altitude target, this works great. If it's a dV target just have the target altitude be current altitude
                        double targetAlt = Inv_SMA - XKMPER;
                        if (config_Array[k, 3] == 2)
                        {
                            targetAlt = config_Array[k, 1];
                        }

                        double time_to_Reboost = nextReboost.Subtract(thisDate).TotalSeconds;
                        //Average the last X day's worth of decay to estimate how high to reboost
                        double target_Alt_Gain = targetAlt - (Inv_SMA - XKMPER) + (decay_Rate * time_to_Reboost);
                        dV_ratio = (target_Alt_Gain / 1.852) / Math.Abs(dt);
                    }

                    //update the config index and config values
                    config_index++;
                    configDate = DateArray[config_index];
                }

                //Check if I need to adjust BN due to low altitudes. Should be a more complex bridging function, but for now
                //I'm just assuming Cd is linear from 2.07 above 200 km down to 0.9197 at 100 km
                Cd = 2.07;
                //if (dt > 0){
                //    if (R - 6378.134 < 100){
                //        Cd = 0.9197;
                //    }else if (R - 6378.134 < 200){
                //        Cd = 0.9197 + ((2.07 - 0.9197) / 100) * (R - 6378.134 - 100);
                //    }
                //}

                //Flag if the object is breaking up and modify BN accordingly
                if (R - XKMPER <= breakupAlt && breakupBN > 0)
                {
                    massAreaRatio = breakupBN * 2.07;
                }

                //Regularly recalculate BN in case the Cd has changed
                BN = massAreaRatio / Cd;

                //Find latitude and longitude so i can calculate density
                double lastLong = LongLat[0];
                double lastLat = LongLat[3];
                LongLat = OETools.Long_Lat(JDN, thisDay_Fraction, x);

                //Find Densities
                double this_Density = Pow(1000, 3) * OETools.VBA_Density_Finder(currentGMT, LongLat[2], LongLat[3], LongLat[0], Current_SFU, Current_AP);

                //Calculate Beta angle 
                double Beta = OETools.BetaCalcs(Osc_Elements[3], Osc_Elements[2], Osc_Elements[1], JDN);

                bool LatCheck; //this is a boolean to check if we're on the ascending node
                if (previous_Lat < 0 && LongLat[3] >= 0)
                {
                    if (DailyOrbit != 0)
                    {
                        DailyOrbit++;
                    }
                    LatCheck = true;

                    //Interpolate Longitude of Ascending Node
                    LAN = lastLong - lastLat * (LongLat[0] - lastLong) / (LongLat[3] - lastLat);

                    //Find this Daily Orbit
                    double orbit_Longitude_movement = period * We * 180.0 / Math.PI;
                    double DO1Long = 20;
                    double testLAN = LAN;

                    if (testLAN > DO1Long)
                    {
                        testLAN -= 360;
                    }

                    DailyOrbit = 1;
                    double checkLon = DO1Long - orbit_Longitude_movement;
                    while (testLAN < checkLon)
                    {
                        checkLon -= orbit_Longitude_movement;
                        DailyOrbit++;
                    }

                    //check daily orbit
                    if (previous_Lon > DO1Long && LongLat[0] <= DO1Long)
                    {
                        //Define daily orbit parameters
                        DailyOrbit = 1;
                    }
                    previous_Lon = LongLat[0];
                }
                else
                {
                    LAN = 0;
                    LatCheck = false;
                }

                //Check whether to output this time step
                if (Math.Abs(current_dt) >= Output_Check || LatCheck && LAN_Output == true || Inv_Hp < 150 && Inv_Ha < 500 || special_output)
                {
                    //Append the data to my output string
                    cleanResponse.Append(SatID + ",");
                    cleanResponse.Append(thisDate.Year + ",");
                    cleanResponse.Append(currentGMT + ",");
                    cleanResponse.Append(Osc_Elements[0] + ","); //SMA
                    cleanResponse.Append(Osc_Elements[1] + ","); //e
                    cleanResponse.Append(Osc_Elements[2] + ","); //i
                    cleanResponse.Append(Osc_Elements[4] + ","); //Wp
                    cleanResponse.Append(Osc_Elements[3] + ","); //Ra
                    cleanResponse.Append(Osc_Elements[5] + ","); //Ma
                    cleanResponse.Append(Osc_Elements[6] + ","); //TA
                    cleanResponse.Append(Current_SFU + ",");
                    cleanResponse.Append(Current_AP + ",");
                    cleanResponse.Append(this_Density + ",");
                    cleanResponse.Append(BN + ",");
                    cleanResponse.Append(Inv_hSMA + ",");
                    cleanResponse.Append(Inv_Ha + ",");
                    cleanResponse.Append(Inv_Hp + ",");
                    cleanResponse.Append(LongLat[0] + ",");
                    cleanResponse.Append(LongLat[3] + ",");
                    //cleanResponse.Append(LongLat[1] + ",");
                    cleanResponse.Append(Beta + ",");
                    cleanResponse.Append(thisDate + ",");
                    cleanResponse.Append(",");
                    //cleanResponse.Append(dV_ratio * dt + ",");
                    cleanResponse.Append(reported_dV + ",");
                    cleanResponse.Append(DailyOrbit + ",");
                    if (LAN != 0)
                    {
                        cleanResponse.Append(LAN + ",");
                    }
                    else
                    {
                        cleanResponse.Append(",");
                    }
                    double thisDuration = thisDate.Subtract(startDate).TotalDays;
                    cleanResponse.Append(thisDuration + ",");
                    cleanResponse.Append(",");// thisDistance + ",");
                    cleanResponse.Append(",,,,,,");

                    //Add Spatial Density 
                    double altBin = 0;
                    int altStepSize = 10;
                    while (altBin < 650)
                    {
                        double SpatialDensity = OETools.spatialDensity(Hp + XKMPER, Ha + XKMPER, altBin + XKMPER, altBin + XKMPER + altStepSize);
                        cleanResponse.Append(SpatialDensity + ",");
                        altBin += altStepSize;
                    }

                    double inst_alt = R - XKMPER;
                    cleanResponse.Append("," + x[0] + "," + x[1] + "," + x[2] + "," + x[3] + "," + x[4] + "," + x[5] + ",," + inst_alt + "," + V);
                    cleanResponse.Append("\r\n");

                    //Record this line in the data
                    System.IO.File.AppendAllText(writeDir, cleanResponse.ToString());

                    //Clear the string
                    cleanResponse = new StringBuilder();

                    //Increment to the next output check
                    if (special_output == false) { Output_Check += Output_dt; }
                }

                previous_Lat = LongLat[3];

                //Increment the TAindex or wrap back around to 0
                prev_TA = Osc_Elements[6];

                //Calculate the Sun Vector
                double[] SunVector = new double[4] { 0, 0, 0, 0 };
                //SunVector = SunVectorCalc(JDN);

                //Apply Perturbations on the vehicle and its effect over the chosen time step
                //x = RK4(dt, BN, this_Density, SunVector, x);
                x = RK4_5(dt, BN, this_Density, SunVector, x);

                //Increment Time for next loop
                double delta = configDate.Subtract(thisDate).TotalSeconds;
                //if (delta < Math.Abs(preset_dt) && delta > 0){
                if (Math.Abs(delta) < Math.Abs(preset_dt))
                {
                    special_output = true;
                    if (delta > 0)
                    {
                        dt = delta;
                    }
                    else
                    {
                        dt = -delta;
                    }
                }
                else
                {
                    dt = preset_dt;
                    special_output = false;
                }
                current_dt += dt;

                //Need to find altitude decay over time for reboost values - this creates an array with the last X day's worth of decay rates
                double delta_time = Math.Abs(prev_dt - current_dt);
                if (delta_time >= 86400 && LatCheck)
                {
                    double altDelta = prev_hSMA - Inv_hSMA;
                    if (dt > 0 && altDelta > 0 || dt < 0 && altDelta < 0)
                    {
                        //Loop over the same 10 array indices to prevent saving too much useless data
                        if (decay_count > decay_Array.Length - 1) { decay_count = 0; }

                        decay_Array[decay_count] = altDelta / 86400;// Output_dt;
                        decay_count++;
                    }
                    prev_dt = current_dt;
                    prev_hSMA = Inv_hSMA;
                }

                //Check if i need a status bar update
                if (Math.Abs(current_dt) >= currentStatus)
                {
                    Console.Write("_");
                    currentStatus += statusBarLen;
                }
            }

            ////Record the output - already tried it this way, it's not faster
            //write_to_file(cleanResponse, this_dir + "Documents\\source\\repos\\Output\\", FileName);
            //cleanResponse = new StringBuilder();
        }
    

        //This program does orbital propagation using inputs that have been converted from TLE mean elements to 
        //oscillating keplerian elements
        static void Main() {

            OE_Tools OETools = new OE_Tools();

            //Read in the input file
            string inputfile = Path.Combine(Global.inputdir, "input.txt"); 
            string[] inputArray = System.IO.File.ReadAllLines(inputfile);
            int caseCount = inputArray[0].Split(',').Length - 1;
            
            int i = 0;
            //Satellite loop
            while (i <= caseCount)
            {
                string[] thisInput = new string[26];
                for (int j = 0; j < inputArray.Length; j++)
                {
                    thisInput[j] = inputArray[j].Split(',')[i];
                }
                Console.WriteLine("\nLoop " + i + " of " + caseCount);
                Propagator(thisInput);

                //Loop to next satellite input
                i += 1;
            }
        }
    }
}
