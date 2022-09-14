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
// #include "SingleCell/TP06.h"
#include "SingleCell/TP06.h"
// #include "SingleCell/TPORd.cc"
// #include "SingleCell/ORd.cc"
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
	double dx = 0.15; // mm
	double dt = 0.02; // ms
	int sanCellNum = 0;
	int atrialCellNum = 0;
	int epiCellNum = 0;
	int mCellNum = 0;
	int endoCellNum = 0;
	#if defined(SAN) || defined(HETEROGENEITY)
	sanCellNum = 50; // cell numbers in OneD strand
	#endif
	#if defined(ATRIA) || defined(HETEROGENEITY) 
	atrialCellNum = 50;
	#endif
	#if defined(VENT)	
	endoCellNum = 37;  // modified according to ORd
	mCellNum = 26;
	epiCellNum = 37;
	#endif
	// following parameters are for 1D or 2D.3D
	double atrialCoeff = 1.0*0.0195; // coefficient parameter in OneD strand, but needed to be validated by fitting CV.
	double sanCoeff = 0.1*0.0195; // coefficient parameter in OneD strand, but needed to be validated by fitting CV.
	double ventCoeff = 0.154; // 0.0106      0.154 for TP06


	// for atria
	#ifdef ATRIA
	int numS1 = 5;
	double BCL = 1000; // ms
	double stopTime = numS1*BCL; //ms
	double stimStrength = -10.0;//8.78; //8.78;//-8.78; // pA
	double stimDuration = 5.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	#endif
	
	// for SAN
	#ifdef SAN
	int numS1 = 50; // actually its useless in SAN, just to keep a consistent equations form
	double BCL = 500; // ms, actually its useless in SAN, just to keep a consistent equations form
	double stopTime = 2000; //ms
	double stimStrength = -0.0; // pA
	double stimDuration = 0.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	#endif

	// for ventricle
	#ifdef VENT
	int numS1 = 2;
	double BCL = 800; // ms   325(2:1)   250(1:1)  1000(EAD) 500(normal) 750(normal)
	double stopTime = numS1*BCL; //ms
	double stimStrength = -52.0;//8.78; //8.78;//-8.78; // pA
	double stimDuration = 1.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	#endif


	// for heteregeneity 1D strand
	#ifdef HETEROGENEITY
	int numS1 = 2; // actually its useless in SAN+Atria, just to keep a consistent equations form
	double BCL = 1000; // ms, actually its useless in SAN+Atria, just to keep a consistent equations form
	double stopTime = 1500; //ms
	double stimStrength = -0.0; // pA
	double stimDuration = 0.0;	// ms
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
		#ifdef VENT // set coeff to 'ventCoeff' whatever it was if VENT defined.
		if(i == 63)
			coeff[i] = 0.2*ventCoeff;
		else
			coeff[i] = ventCoeff;
		#endif
	}


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



	




	int right = BCL;
	int left = 0;
	int s2startTime = 0;
	int erpflag = 0;
	double lastV = -1000;
	double time = 0;
	int step = 0;
	FILE *initfile;
	FILE *datafile;

	while(right - left > 1)
	{
		// 
		s2startTime = (right + left)/2; // 472 小 473 ok
		erpflag = 0;
		lastV = -1000;
		cout << endl << "s2 = " << s2startTime << "; left = " << left << "; right = " << right << endl;
		datafile = fopen("Outputs/VentOneDResults.dat","w+");



		// initialization
		for (int i = 0; i < cellNum; i++)
		{
			#ifdef CON
			if(i >= 0 && i < endoCellNum)
			{
				strand[i] = new TP06(ENDO);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_CON_ENDO.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			else if (i < endoCellNum + mCellNum)
			{
				strand[i] = new TP06(MCELL);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_CON_MCELL.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			else // i < total cellnum
			{
				strand[i] = new TP06(EPI);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_CON_EPI.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			#endif


			#ifdef SQT1
			if(i >= 0 && i < endoCellNum)
			{
				strand[i] = new TP06(ENDO);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_SQT1_ENDO.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			else if (i < endoCellNum + mCellNum)
			{
				strand[i] = new TP06(MCELL);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_SQT1_MCELL.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			else // i < total cellnum
			{
				strand[i] = new TP06(EPI);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_SQT1_EPI.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			#endif


			#ifdef QUIN
			if(i >= 0 && i < endoCellNum)
			{
				strand[i] = new TP06(ENDO);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_QUIN_ENDO.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			else if (i < endoCellNum + mCellNum)
			{
				strand[i] = new TP06(MCELL);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_QUIN_MCELL.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			else // i < total cellnum
			{
				strand[i] = new TP06(EPI);
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				initfile = fopen("SingleCell/TP06InitialValues_QUIN_EPI.dat","r");
				strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
				fclose(initfile);
			}
			#endif

			// apply user configuration about dt
			strand[i]->setDt(dt);
		}


		time = 0;
		step = 0;
		for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
		{
			// 1. progress stats
			if(step%25000 == 0) // 25000 * dt ms = 0.125s 
				cout << "Progress = " << 100.0*time/stopTime << "\%." << endl;

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
			    	if(i < 3 && i >= 0)
					{// cells get stimulation in certain duration
						// as strand[i] is a pointer, A->B should be used instead of A.B
						// ref. http://blog.csdn.net/qq457163027/article/details/54237782
						// cout << "here" << endl;
						strand[i]->setIstim(stimStrength);
					}
				}
				// following is S2 protocal			
				else if(time - (numS1-1)*BCL >= s2startTime && // s2 is an advanced stimulation, thus s2startTime ranges between 0~(BCL)
			   			time - (numS1-1)*BCL < s2startTime + stimDuration) 
				{// index from 10 to 15 cells get S2 stimulation
					if(i < 3 && i >= 0)
					// further constraint s2startTime to 0~(BCL - stimulation)
					{// cells get stimulation in certain duration
						// as strand[i] is a pointer, A->B should be used instead of A.B
						// ref. http://blog.csdn.net/qq457163027/article/details/54237782
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
					// use strand[cellNum-1] instead of "strand[cellNum]" as the latter is out of index
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
			if(step%10 == 0) // 50*dt = 1 ms once
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

			// 4. CV stats. Note that this can be only applied to non-selfpacing cells!!!
			// only record the last S1 data
			if (floor(time/BCL) == numS1-1) 
			{ 
				if(strand[99]->getV() >= -30 && lastV < -30)
				{
					erpflag++;
				}
				lastV = strand[99]->getV();
			} // end CV calculation

		}// end of timeloop
		fclose(datafile);

		if(erpflag == 2)
		{
			cout << "beyond ERP!" << endl;
			right = s2startTime;
		}
		else if(erpflag == 1)
		{
			cout << "within ERP!" << endl;
			left = s2startTime;
		}
		else
		{
			cout << "conduction failure!" << endl;
		}
	}

	cout << "Tissue ERP = " << left << " (ms)." << endl;


	printf("All done.\n");

	return 0;
}