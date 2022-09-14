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


int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.
	
	// 1D parameters for all cell types
	double dx = 0.15; // might not be exactly the same in different cell types. remain fully understand
	double dt = 0.02; // ms
	int sanCellNum = 0;
	int atrialCellNum = 0;
	int epiCellNum = 0;
	int mCellNum = 0;
	int endoCellNum = 0;

	#if defined(VENT)
	endoCellNum = 0;
	mCellNum = 0;
	epiCellNum = 100;	
	#endif
	// following parameters are for 1D or 2D.3D
	double ventCoeff = 0.154; // 0.0106

	// for ventricle
	#ifdef VENT
	int numS1 = 10;
	double BCL = 800; // ms   325(2:1)   250(1:1)  1000(EAD) 500(normal) 750(normal)
	double stopTime = numS1*BCL; //ms
	double stimStrength = -80.0;//8.78; //8.78;//-8.78; // pA
	double stimDuration = 1.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	#endif



	// --------start simulation--------	
	// CV calculation stuff
	double cvStartTime = 0;
	double cvEndTime = 0;
	int cvStartFlag = 0; // 0 for not start yet
	int cvEndFlag = 0;
	double cv = 0;

	// parallel stuff
	int coreNum = 8;//omp_get_num_procs();
	omp_set_num_threads(2 * coreNum);

	// strand initilization, diffusion stuff
	int cellNum = sanCellNum + atrialCellNum + epiCellNum + mCellNum + endoCellNum; // number of cells in OneD strand
	typedef Cell* CellPointer;
	TP06* strand[cellNum]; // note that constructor contains the initializer
	double coeff[cellNum]; // diffusion parameters for each cell
	double dcoeff_dx[cellNum]; // first order derivative for each cell
	double oldV[cellNum];

	// assign coeff according to cell type
	for(int i = 0; i < cellNum; i++)
	{
		// Step 1: general case without yet considering transition in the conjunction of two heterogenous tissues
		#ifdef VENT // set coeff to 'ventCoeff' whatever it was if VENT defined.
		if (i == 60)
			coeff[i] = ventCoeff;
		else
			coeff[i] = ventCoeff;
		#endif
	}

	// for(int j = 0; j < 3; j++)
		// coeff[j] = 1*ventCoeff;
	// coeff[35] = 1*ventCoeff;
	// coeff[36] = 1*ventCoeff;
	// coeff[65] = 1*ventCoeff;
	// coeff[66] = 1*ventCoeff;

	// Calculate the dcoeff/dx(i.e. dcoeff_dx in the code) in the 1D strand
	for(int i = 0; i < cellNum; i++)
	{
		if (i == 0) 
			dcoeff_dx[i] = (coeff[i+1] - coeff[i])/dx;
		else if (i == cellNum-1) 
			dcoeff_dx[i] = (coeff[i] - coeff[i-1])/dx;
		else
			dcoeff_dx[i] = (coeff[i+1] - coeff[i-1])/(2.0*dx);
	}



	// Initialization
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < cellNum; i++)
	{

		#ifdef CON
		FILE *initfile;
		if(i >= 0 && i < endoCellNum)
		{
			strand[i] = new TP06(ENDO);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_CON_ENDO.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new TP06(MCELL);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_CON_MCELL.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else // i < total cellnum
		{
			strand[i] = new TP06(EPI);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_CON_EPI.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}	
		#endif

		#ifdef SQT1
		FILE *initfile;
		if(i >= 0 && i < endoCellNum)
		{
			strand[i] = new TP06(ENDO);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_SQT1_ENDO.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new TP06(MCELL);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_SQT1_MCELL.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else // i < total cellnum
		{
			strand[i] = new TP06(EPI);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_SQT1_EPI.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}	
		#endif

		#ifdef QUIN
		FILE *initfile;
		if(i >= 0 && i < endoCellNum)
		{
			strand[i] = new TP06(ENDO);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_QUIN_ENDO.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new TP06(MCELL);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_QUIN_MCELL.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else // i < total cellnum
		{
			strand[i] = new TP06(EPI);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_QUIN_EPI.dat","r");
			// strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}	
		#endif


		// apply user configuration about dt
		strand[i]->setDt(dt);
	}


	#ifdef VENT
	FILE *datafile = fopen("Outputs/VentOneDResults.dat","w+");
	#endif



	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		// 1. progress stats
		if(step%25000 == 0) // 25000 * dt ms = 0.5s 
			cout << "S1 Progress = " << 100.0*time/stopTime << "\%." << endl;

		// set oldV
		for(int i = 0; i < cellNum; i++)
		{
			oldV[i] = strand[i]->getV();
		}
		
		// 2. cell loop
		// #pragma omp parallel for schedule(static)
		for(int i = 0; i < cellNum; i++)
		{
			// ---stimulation applying or not, based on the time and cell location---
			// initialize all stim to 0.0
			strand[i]->setIstim(0.0);
			// Apply stimulus according to user configuration (about time)
			if(time - floor(time/BCL)*BCL >= stimStart && 
		   	   time - floor(time/BCL)*BCL < stimStart + stimDuration)
			{// first 5 cells get S1 stimulation
		    	if(i < 3)
				{// cells get stimulation in certain duration
					// as strand[i] is a pointer, A->B should be used instead of A.B
					// ref. http://blog.csdn.net/qq457163027/article/details/54237782
					// cout << "here" << endl;
					strand[i]->setIstim(stimStrength);
				}
			}
			// following is S2 protocal, temporally commented out
			/**
			else if(time - (numS1-1)*BCL >= s2startTime && // s2 is an advanced stimulation, thus s2startTime ranges between 0~(BCL)
		   			time - (numS1-1)*BCL < s2startTime + stimDuration) 
			{// index from 10 to 15 cells get S2 stimulation
				if(i >= 10 && i < 16)
				// further constraint s2startTime to 0~(BCL - stimulation)
				{// cells get stimulation in certain duration
					// as strand[i] is a pointer, A->B should be used instead of A.B
					// ref. http://blog.csdn.net/qq457163027/article/details/54237782
					strand[i]->setIstim(s2stimStrength);
				}
			}
			**/


			// ---------calculate diffusion, i.e. dVgap---------
		
			double dVgap_dt = 0;
			double first_order;
			double second_order;

			// Step 1: calculate first and second order of membrane potential
			if(i == 0) 
			{
				// use oldV[0] instead of "oldV[-1]"
				first_order = (oldV[i+1] - oldV[i])/(1.0*dx);
				second_order = (oldV[i+1] + oldV[i] - 2.0*oldV[i])/(dx*dx);
			}
			else if(i > 0 && i < cellNum - 1) 
			{
				// normal case
				first_order = (oldV[i+1] - oldV[i-1])/(2.0*dx);
				second_order = (oldV[i+1] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}
			else if(i == cellNum - 1)
			{
				// use oldV[cellNum-1] instead of "oldV[cellNum]" as the latter is out of index
				first_order = (oldV[i] - oldV[i-1])/(1.0*dx);
				second_order = (oldV[i] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}

			// Step 2: calculate dVgap according to equations
			dVgap_dt = dcoeff_dx[i]*first_order + coeff[i]*second_order;

			//if(step%1000 == 0) cout<<strand[i]->getVgap()<<endl;
			strand[i]->setDVgap_dt(dVgap_dt);
				
			// update
			strand[i]->update();
		}// end cell loop


		// 3. output file. Unfortunately, this part cannot run in paralell
		// if( floor(time/BCL) >= numS1 - 20) // output final cycle only
		if(step%10 == 0 && floor(time/BCL) == numS1-1) // 50*dt = 1 ms once
		{
			for(int j = 0; j < cellNum; j++)
			{
				// write time before each time step
				if(j == 0)
					fprintf(datafile,"%4.10f\t", time);
				// write volt for each cell
				fprintf(datafile,"%4.10f\t", strand[j]->getV()); // unit: mV
				// write enter after each time step
				if(j == cellNum - 1)
					fprintf(datafile,"\n");
			}
		}// end Membrane Potential recording


	}// end of s1 timeloop
	fclose(datafile);
	

	// save s1 results!!!
	FILE *s2_initfile = fopen("Outputs/S2_INITFILE_1D.dat","w+");
	for(int i = 0; i < cellNum; i++)
	{
		strand[i] -> outputAllStates(s2_initfile);
	}


	// loop s2
	double apd90;
	#ifdef CON
	FILE *rcfile = fopen("Outputs/CVRestCurve_CON_DI.dat","w+");
	apd90 = 308.94;// APD90 can be put here.
	#endif

	#ifdef SQT1
	FILE *rcfile = fopen("Outputs/CVRestCurve_SQT1_DI.dat","w+");
	apd90 = 235.16;// APD90 can be put here.
	#endif

	#ifdef QUIN
	FILE *rcfile = fopen("Outputs/CVRestCurve_QUIN_DI.dat","w+");
	apd90 = 284.74;// APD90 can be put here.
	#endif

	double DI = 5000;
	     
	while(DI > 50)
	{
		// init!!!!
		rewind(s2_initfile);
		for(int i = 0; i < cellNum; i++)
		{
			strand[i] -> readinAllStates(s2_initfile);
		}

		// reset.
		cvStartFlag = 0;
		cvEndFlag = 0;

		// calculate DI for s2 restitution.
		if(DI > 4000)
		{
			DI = DI - 1000;
		}
		else if(DI > 3000)
		{
		  	DI = DI - 1000;
		}
		else if(DI > 2000)
		{
		  	DI = DI - 1000;
		}
		else if(DI > 1000)
		{
		  	DI = DI - 1000;
		}
		else if(DI > 500)
		{	     
		  	DI = DI - 250;
		}
		else if(DI > 400)
		{
		  	DI = DI - 50;
		}
		else if(DI > 300)
		{
		  	DI = DI - 10;
		}
		else if(DI > 250)
		{
		  	DI = DI - 5;
		}
		else if(DI > 50)
		{
		  	DI = DI - 1;
		}

		// start s2
		for(time = 0.0, step = 0; time < apd90 + DI + BCL; time += dt, step++)
		{
			// set oldV
			for(int i = 0; i < cellNum; i++)
			{
				oldV[i] = strand[i]->getV();
			}
			// 2. cell loop
			// #pragma omp parallel for schedule(static)
			for(int i = 0; i < cellNum; i++)
			{
				// ---stimulation applying or not, based on the time and cell location---
				// initialize all stim to 0.0
				strand[i]->setIstim(0.0);
				// Apply stimulus according to user configuration (about time)
				if( (time >= 0 && time < stimDuration) ||  // s1
					(time >= stimDuration + apd90 + DI && time < apd90 + DI + 2*stimDuration))  // s2  // 加stimDuration是因为APD90是从stim结束后开始算起的
				{// first 5 cells get S1 stimulation
			    	if(i < 3 && i >= 0)
					{// cells get stimulation in certain duration
						// as strand[i] is a pointer, A->B should be used instead of A.B
						// ref. http://blog.csdn.net/qq457163027/article/details/54237782
						// cout << "here" << endl;
						strand[i]->setIstim(stimStrength);
					}
				}



				// ---------calculate diffusion, i.e. dVgap---------
			
				double dVgap_dt = 0;
				double first_order;
				double second_order;

				// Step 1: calculate first and second order of membrane potential
				if(i == 0) 
				{
					// use strand[0] instead of "strand[-1]"
					first_order = (oldV[i+1] - oldV[i])/(1.0*dx);
					second_order = (oldV[i+1] + oldV[i] - 2.0*oldV[i])/(dx*dx);
				}
				else if(i > 0 && i < cellNum - 1) 
				{
					// normal case
					first_order = (oldV[i+1] - oldV[i-1])/(2.0*dx);
					second_order = (oldV[i+1] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
				}
				else if(i == cellNum - 1)
				{
					// use oldV[cellNum-1] instead of "oldV[cellNum]" as the latter is out of index
					first_order = (oldV[i] - oldV[i-1])/(1.0*dx);
					second_order = (oldV[i] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
				}

				// Step 2: calculate dVgap according to equations
				dVgap_dt = dcoeff_dx[i]*first_order + coeff[i]*second_order;

				//if(step%1000 == 0) cout<<strand[i]->getVgap()<<endl;
				strand[i]->setDVgap_dt(dVgap_dt);
					
				// update
				strand[i]->update();
			}// end cell loop


			// 4. CV stats. Note that this can be only applied to non-selfpacing cells!!!
			// only record the last S1 data
			if (time > apd90 + DI)
			{ 
				// record propagation start and end time for CV calculation
				if(strand[10]->getV() > -30 && cvStartFlag == 0)
				{
					cvStartTime = time;
					cvStartFlag = 1;
				}
				if(strand[90]->getV() > -30 && cvEndFlag == 0)
				{
					cvEndTime = time;
					cvEndFlag = 1;
					cv = (dx * (90 - 10)) / (cvEndTime - cvStartTime);
				}
			} // end CV calculation

		} // end s2.
		if(cvStartFlag == 1 && cvEndFlag == 1)
		{
			cout << "DI = " << DI << " ms    CV = " << cv << " m/s." << endl;
			fprintf(rcfile,"%.5f\t",DI); // unit: ms
			fprintf(rcfile,"%.5f\n",100*cv); // unit: cm/s
		}
		else
		{
			cout << "DI = " << DI << " ms    Conduction failure!" << endl;
			fprintf(rcfile,"%.5f\t",DI); // unit: ms
			fprintf(rcfile,"inf\n"); // unit: ms
		}


	}// end DI loop



	fclose(rcfile);
	fclose(s2_initfile);

	printf("All done.\n");

	return 0;
}