/*
 * General Code Structure (GCS) for one-dimensional simulation
 * 
 * A. Code sructure is re-organised using standard C++ to fit into my project.
 * B. The friendly style makes the project easier understood.
 * C. The version is more extendable whatever to further 1D,2D or another single cell model. 
 * 
 * Under Intellectual Property Protection.
 * 
 * 
 * Author      : Shugang Zhang <zsg@ouc.edu.cn>
 * Last update : 14-09-2022
 */

// #include "SingleCell/NeonatalRatAtria.cc"
// #include "SingleCell/RatSAN.cc"
#include "SingleCell/TP06.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "omp.h"

using namespace std;

// ONLY ONE LINE CAN BE KEPT HERE
// #define SAN
// #define ATRIA
#define VENT
// #define HETEROGENEITY 


// #define INIT // no initiazation file, give successive S1 simulation to generate initialization file.
#define S1S2 // for normal S1S2 2D simulation with S2 for spiral wave detection, supposing 2D initialization file exists.



int main(int args, char *argv[])
{
	// 3mm * 10mm 2D sheet, i.e.,30(transmural) * 100(latitude) = 3000 cells
	// spatial resolution: 0.1mm
	//
	// Diagram
	// *********** EPI  CELLS ***********
	// *               T    			*
	// *               R    			*
	// *           	   A			  	*
	// *               S    			*
	// *               M    			*
	// *********** ENDO CELLS ***********

	double dx = 0.15; // might not be exactly the same in different cell types. remain fully understand
	double dt = 0.02; // ms
	int BCL = 1000;
	double startTime = 0.0;
	double time = 0.0;
	int step = 0;
	double stimDuration = 3; // unit: s, i.e.,1ms
	double stimStrength = -52.0; // unit: nA or 10pA/pF
	double coeffx = 0.154; // mm2/ms from Henggui's paper 
	double coeffy = 0.154; // mm2/ms from Henggui's paper 

	// default.
	#ifndef S1S2
	int numS1 = 2;
	double stopTime = numS1 * BCL;
	#endif

	#ifdef S1S2
	int numS1 = 1; // 待确认
	double s2stimStrength = -104;
	// s222222222
	double s2startTime = 352; //        ( 370.3 CON    291 SQT1    352 QUIN)             
	double s2stimDuration = 3;
	double stopTime = numS1 * BCL;

	cout << "s2startTime = " << s2startTime << endl;
	#endif

	// create Volt file or recreate file and delete previous existing version
	FILE *FF;
	char filename[200];
	sprintf(filename, "Outputs/Transmural2D.dat");
	FF = fopen(filename,"w+");
	// fclose(FF);



	// 2D sheet
	int D1 = 100; // 30 cells in transmural orientation
	int D2 = 500; // 100 cells in latitude orientation
	int endoCellNum = 37; // homogenous tissue!
	int mCellNum = 26;
	int epiCellNum = 37;
	TP06 ***sheet = new TP06**[D1]; // 2D sheet
	double oldV[D1][D2];

	FILE *initfile;

	for(int i = 0; i < D1; i++)
	{
		sheet[i] = new TP06*[D2];
		for(int j = 0; j < D2; j++)
		{
			#ifdef CON
			if(i < endoCellNum)
			{
				sheet[i][j] = new TP06(ENDO);
				initfile = fopen("SingleCell/TP06InitialValues_CON_ENDO.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			else if(i < endoCellNum + mCellNum)
			{
				sheet[i][j] = new TP06(MCELL);
				initfile = fopen("SingleCell/TP06InitialValues_CON_MCELL.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			else
			{
				sheet[i][j] = new TP06(EPI);
				initfile = fopen("SingleCell/TP06InitialValues_CON_EPI.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			#endif


			#ifdef SQT1
			if(i < endoCellNum)
			{
				sheet[i][j] = new TP06(ENDO);
				initfile = fopen("SingleCell/TP06InitialValues_SQT1_ENDO.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			else if(i < endoCellNum + mCellNum)
			{
				sheet[i][j] = new TP06(MCELL);
				initfile = fopen("SingleCell/TP06InitialValues_SQT1_MCELL.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			else
			{
				sheet[i][j] = new TP06(EPI);
				initfile = fopen("SingleCell/TP06InitialValues_SQT1_EPI.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			#endif

			#ifdef QUIN
			if(i < endoCellNum)
			{
				sheet[i][j] = new TP06(ENDO);
				initfile = fopen("SingleCell/TP06InitialValues_QUIN_ENDO.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			else if(i < endoCellNum + mCellNum)
			{
				sheet[i][j] = new TP06(MCELL);
				initfile = fopen("SingleCell/TP06InitialValues_QUIN_MCELL.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			else
			{
				sheet[i][j] = new TP06(EPI);
				initfile = fopen("SingleCell/TP06InitialValues_QUIN_EPI.dat","r");
				sheet[i][j]->readinAllStates(initfile);
				fclose(initfile);
			}
			#endif

			sheet[i][j]->setDt(dt);
		}
	}




	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		if(step%2000 == 0) // 2e4 * 5e-6 = 0.1s
			cout << "Progress = " << 100.0*time/stopTime << "\%." << endl;
			// cout << "s2startTime = " << 1000*s2startTime << "ms, Progress = " << 100.0*time/stopTime << "\%." << endl;
		// Open file and write time at line start
		// BEGIN! FOR! 1D!
		
		// #pragma omp parallel for
		// for(int i = 0; i < transCellNum; i++)
		// {
		// 	for(int j = 0; j < latiCellNum; j++)
		// 	{
		// 		ratLVSheet[i][j]->setIstim(0.0);
		// 	}
		// }


		for(int i = 0; i < D1; i++)
			for(int j = 0; j < D2; j++)
				oldV[i][j] = sheet[i][j]->getV();



		#pragma omp parallel for schedule(static)
		for(int i = 0; i < D1; i++)
		{			
			for(int j = 0; j < D2; j++)
			{
				sheet[i][j]->setIstim(0);
				// if( i >= 20 && i < 30 && j >= 20 && j < 30 )  // give stimulus to first 3 row cells
				if(i < 3 )  // give stimulus to first 3 row cells
				{// first 4 cells get stimulation
			    	if( // floor(time/BCL) <= numS1 &&
			    		time - floor(time/BCL)*BCL >= startTime && 
			   			time - floor(time/BCL)*BCL < startTime + stimDuration)
					{// cells get stimulation in certain duration
						// as ratLVSheet[i][j] is a pointer, A->B should be used instead of A.B
						// ref. http://blog.csdn.net/qq457163027/article/details/54237782
						sheet[i][j]->setIstim(stimStrength);
					}
					else
					{
						sheet[i][j]->setIstim(0.0);
					}
				}

				// s222222222222222222222222

				#ifdef S1S2  //   (314 CON  247 SQT6  377 QUIN)
				else if( (i >= 66 && i < 69) && (j >= 0 && j < 377) ) //   312 short   314 ok
				{// index from 10 to 15 cells get S2 stimulation
					if(time - (numS1-1)*BCL >= s2startTime && // s2 is an advanced stimulation, cthus s2startTime ranges between 0~(BCL)
			   			time - (numS1-1)*BCL < s2startTime + s2stimDuration) // further constraint s2startTime to 0~(BCL - stimulation)
					{// cells get stimulation in certain duration
						// as ratLVStrand[i] is a pointer, A->B should be used instead of A.B
						// ref. http://blog.csdn.net/qq457163027/article/details/54237782
						sheet[i][j]->setIstim(s2stimStrength);
					}

				}
				#endif 

				// calculate dVgapx
				double dVgapx = 0;
				if(i == 0) 
				{
					// use ratLVSheet[0][j] instead of "ratLVSheet[-1]"
					dVgapx = coeffx*(oldV[i+1][j] - 2.0*oldV[i][j] + oldV[i][j])/(dx*dx);
				}
				else if(i > 0 && i < D1 - 1) 
				{
					// normal case
					dVgapx = coeffx*(oldV[i+1][j] - 2.0*oldV[i][j] + oldV[i-1][j])/(dx*dx);
					// the coefficient was decreased 5-fold at the position 60
					// if(i == 59)
					// 	dVgapx = (-0.4*coeffx)*(oldV[i+1][j] - oldV[i-1][j])/(4.0*dx*dx) + 1.0*coeffx*(oldV[i+1][j] - 2.0*oldV[i][j] + oldV[i-1][j])/(dx*dx);
					if(i == 60)
						dVgapx = 0.2*coeffx*(oldV[i+1][j] - 2.0*oldV[i][j] + oldV[i-1][j])/(dx*dx);
					// if(i == 61)
					// 	dVgapx = ( 0.4*coeffx)*(oldV[i+1][j] - oldV[i-1][j])/(4.0*dx*dx) + 1.0*coeffx*(oldV[i+1][j] - 2.0*oldV[i][j] + oldV[i-1][j])/(dx*dx);
				}
				else if(i == D1 - 1) 
				{
					// use ratLVSheet[cellNum-1] instead of "ratLVSheet[cellNum]" as the latter is out of index
					dVgapx = coeffx*(oldV[i][j] - 2.0*oldV[i][j] + oldV[i-1][j])/(dx*dx);
				}

				



				// calculate dVgapy
				double dVgapy = 0;
				if(j == 0) 
				{
					// use ratLVSheet[0][j] instead of "ratLVSheet[-1]"
					dVgapy = coeffy*(oldV[i][j+1] - 2.0*oldV[i][j] + oldV[i][j])/(dx*dx);
				}
				else if(j > 0 && j < D2 - 1) 
				{
					// normal case
					dVgapy = coeffy*(oldV[i][j+1] - 2.0*oldV[i][j] + oldV[i][j-1])/(dx*dx);
				}
				else if(j == D2 - 1) 
				{
					// use ratLVSheet[cellNum-1] instead of "ratLVSheet[cellNum]" as the latter is out of index
					dVgapy = coeffy*(oldV[i][j] - 2.0*oldV[i][j] + oldV[i][j-1])/(dx*dx);
				}
				//if(step%1000 == 0) cout<<ratLVSheet[i]->getVgap()<<endl;
				sheet[i][j]->setDVgap_dt(dVgapx + dVgapy);

				// update V
				sheet[i][j]->update();

			} // end latitude loop
		} // end transmural loop






		// write volt of ratLVSheet[i]
		// if(step%200 == 0 && floor(time/BCL) >= numS1 - 1) // 1ms once
		int count = 0;
		if(step%250 == 0) // 0.02*250 = 5ms once
		{
			fprintf(FF,"%4.10f\t",time);
			for(int i = 0; i < D1; i++)
			{
				for(int j = 0; j < D2; j++)
				{
					count ++;
					// FF = fopen(filename,"a");
					if(i == 0) // write '\n' at line start
					{
						// fprintf(FF,"%4.10f\t",1000*time); // unit: ms
						// fprintf(FF,"%4.10f\t",ratLVSheet[i]->getStimu()); // unit: ms
					}				
					fprintf(FF,"%4.10f\t",sheet[i][j]->getV()); // unit: mV						
					if(i == D1 - 1 && j == D2 - 1) // write '\n' at line end
						fprintf(FF,"\n");
					// fclose(FF);
				}
			}
		}// end writing file
		// cout << "total V: " << count << endl;

		
	}// end of timeloop
	fclose(FF);
}