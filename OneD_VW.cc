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


#include "SingleCell/TP06.h"
// #include "SingleCell/ORd.cc"
// #include "SingleCell/TPORd.cc"
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

	#ifdef VENT	
	endoCellNum = 37;
	mCellNum = 26;
	epiCellNum = 37;
	#endif

	double ventCoeff = 0.154;  // 0.154 for TNNP06      0.3 for ORd (0.653)



	// for ventricle
	#ifdef VENT
	int numS1 = 2;
	double BCL = 800; // ms
	double stopTime = numS1*BCL; //ms
	double stimStrength = -52.0;//8.78; //8.78;//-8.78; // pA
	double stimDuration = 3.0;	// ms
	double s2stimStrength = -104.0;
	double s2stimDuration = 3.0;
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
	Cell* strand[cellNum]; // note that constructor contains the initializer
	double coeff[cellNum]; // diffusion parameters for each cell
	double dcoeff_dx[cellNum]; // first order derivative for each cell
	double oldV[cellNum];

	// assign coeff according to cell type
	for(int i = 0; i < cellNum; i++)
	{
		// Step 1: general case without yet considering transition in the conjunction of two heterogenous tissues
		#ifdef VENT // set coeff to 'ventCoeff' whatever it was if VENT defined.
		if (i == 63)
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
			strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new TP06(MCELL);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_CON_MCELL.dat","r");
			strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else // i < total cellnum
		{
			strand[i] = new TP06(EPI);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_CON_EPI.dat","r");
			strand[i]->readinAllStates(initfile);
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
			strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new TP06(MCELL);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_SQT1_MCELL.dat","r");
			strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else // i < total cellnum
		{
			strand[i] = new TP06(EPI);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_SQT1_EPI.dat","r");
			strand[i]->readinAllStates(initfile);
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
			strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new TP06(MCELL);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_QUIN_MCELL.dat","r");
			strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}
		else // i < total cellnum
		{
			strand[i] = new TP06(EPI);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			initfile = fopen("SingleCell/TP06InitialValues_QUIN_EPI.dat","r");
			strand[i]->readinAllStates(initfile);
			fclose(initfile);
		}	
		#endif

		// apply user configuration about dt
		strand[i]->setDt(dt);
	}


	double time = 0;
	int step = 0;
	FILE* FF;

	double left = 269;  //351.1 ~ 356
	double right = 269;  //351.1 ~ 356       332.8 337.5
	// #define MODE1
	// #define MODE2
	#define MODE3
	int tempstart = 60;  // position

	double boundary = 10; //negative
	double s2startTime = 0.5*(left + right);

	// used for VW determine
	double maxEndo = -100;
	double maxEpi = -100;

	// flag
	int flag = 0;


	while(right - left >= 0.1 || right == left)
	{
		// NOTE THAT PREVIOUS FILENAMES ARE USELESS NOW.
		char VoltFileName[200];
		sprintf(VoltFileName, "Outputs/Transmural1D_s2@%.2f.dat", s2startTime);
		FF = fopen(VoltFileName,"w");
		fclose(FF);

		// re-initialization!!
		#pragma omp parallel for schedule(static)
		for(int i = 0; i < cellNum; i++)
		{
			strand[i]->init(strand[i]->getCellType());

			FILE *initfile;			
			if( strand[i]->getCellType() == ENDO)
			{

				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				// initfile = fopen("SingleCell/ORdInitialValues_SQT6_ENDO.dat","r");
				// strand[i]->readinAllStates(initfile);
				// fclose(initfile);
			}
			else if( strand[i]->getCellType() == MCELL)
			{

				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				// initfile = fopen("SingleCell/ORdInitialValues_SQT6_MCELL.dat","r");
				// strand[i]->readinAllStates(initfile);
				// fclose(initfile);
			}
			else
			{
				// read in initial values (this is because the original init values is not stable yet)
				// if the initfile is not available, run the initialization.cc first
				// initfile = fopen("SingleCell/ORdInitialValues_SQT6_EPI.dat","r");
				// strand[i]->readinAllStates(initfile);
				// fclose(initfile);
			}
		}


		for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
		{
			// cout << "????" << endl;
			if(step%20000 == 0) // 2e4 * 5e-6 = 0.1s 
			{
				cout << "s2startTime = " << s2startTime << "ms, Progress = " << 100.0*time/stopTime << "\%." << endl;
			}


			for(int i = 0; i < cellNum; i++)
			{
				oldV[i] = strand[i]->getV();
			}
			
			// 1. progress stats
			// if(step%25000 == 0) // 25000 * dt ms = 0.125s 
				// cout << "Progress = " << 100.0*time/stopTime << "\%." << endl;
			
			// 2. cell loop
			#pragma omp parallel for schedule(static)
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
				// following is S2 protocal			
				else if(time - (numS1-1)*BCL >= s2startTime && // s2 is an advanced stimulation, thus s2startTime ranges between 0~(BCL)
			   			time - (numS1-1)*BCL < s2startTime + s2stimDuration) 
				{// index from 10 to 15 cells get S2 stimulation
					if(i >= tempstart && i < tempstart+3)
					// further constraint s2startTime to 0~(BCL - stimulation)
					{// cells get stimulation in certain duration
						// as strand[i] is a pointer, A->B should be used instead of A.B
						// ref. http://blog.csdn.net/qq457163027/article/details/54237782
						strand[i]->setIstim(s2stimStrength);
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

			// 3. output file. Unfortunately, this part cannot run in paralell
			if( floor(time/BCL) == numS1 - 1) // output final cycle only
			// if(step%50 == 0) // 50*dt = 1 ms once
			{
				FF = fopen(VoltFileName,"a");
				for(int j = 0; j < cellNum; j++)
				{					
					if(step%50 == 0) // 1ms once
					{									
						fprintf(FF,"%4.10f\t",strand[j]->getV()); // unit: mV						
						if(j == cellNum-1) // write '\n' at line end
							fprintf(FF,"\n");
					}	
				}
				fclose(FF);
			}// end Membrane Potential recording



			// get maxEpi and maxEndo during s2
			if(time - (numS1-1)*BCL >= s2startTime ) 
			{
				if(strand[cellNum-1]->getV() >= maxEpi) maxEpi = strand[cellNum-1]->getV();
				if(strand[0]->getV() >= maxEndo) maxEndo = strand[0]->getV();
			}


		}// end of timeloop
		// fclose(datafile);
		

		cout << "maxEndo = " << maxEndo << endl;
		cout << "maxEpi = " << maxEpi << endl;
		cout << endl;


		// *************mode 1: divide and find*****************

		#ifdef MODE1
		if(maxEndo+boundary < 0 && maxEpi+boundary < 0)
		{
			left = s2startTime;
			s2startTime = 0.5*(left + right);
		}
		else if (maxEndo+boundary > 0 && maxEpi+boundary > 0)
		{
			right = s2startTime;
			s2startTime = 0.5*(left + right);
		}
		else if (  (maxEndo + boundary) * (maxEpi + boundary) <= 0 )
		{
			flag = 1;
		}

		if (flag == 1)	
			break;
		#endif


		// *************mode 2: find left boundary (left must equal to right and must left to VW)*****************
		#ifdef MODE2
		if (  (maxEndo + boundary) * (maxEpi + boundary) <= 0 )
		{
			flag = 1;
		}
		else
		{
			left += 0.1;//0.00001;
			right += 0.1;//0.00001;
			s2startTime = left;
		}

		if (flag == 1)	
			break;
		#endif


		// *************mode 3: find right boundary (left must equal to right and must right to VW)*****************
		#ifdef MODE3
		if (  (maxEndo + boundary) * (maxEpi + boundary) <= 0 )
		{
			flag = 1;
		}
		else
		{
			left -= 0.1;//0.00001;
			right -= 0.1;//0.00001;
			s2startTime = left;
		}

		if (flag == 1)	
			break;
		#endif

		// **************VW search finished here**************


		// reinitialization
		maxEndo = -100;
		maxEpi = -100;
	}

	if(flag == 1) cout << "VM found!" << endl;
	else cout << "No VM!" << endl;

	printf("All done.\n");
	return 0;
}