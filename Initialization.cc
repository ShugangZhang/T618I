/*
 * General Code Structure (GCS) for obtaining steady-state single cell initial values
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
// #include "SingleCell/GrandiCaMK.cc"
// #include "SingleCell/TP06.h"
#include "SingleCell/ToRORD.cc"

using namespace std;

// ONLY ONE LINE CAN BE KEPT HERE
// #define SAN
// #define ATRIA
#define VENT

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.

	#ifdef VENT
	double BCL = 1000; // ms
	double dt = 0.02; //ms
	double stopTime = 200*BCL; //ms
	double stimStrength = -53; // pA/pF
	double stimDuration = 1.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	typedef ToRORD CellType;
	typedef ToRORD* CellPointer;
	#endif


	// --------start simulation--------
	// note that constructor contains the initializer
	CellPointer cell = new CellType(EPI);
	#ifdef VENT
	FILE *initfile = fopen("SingleCell/ToRORDInitialValues_CON_EPI.dat","w+");
	#endif
	// apply user configuration about dt
	cell->setDt(dt);

	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		// 1. progress stats
		if(step%4000000 == 0) // 300000 * dt ms = 4.5s    0.005ms*4000000 = 20s
			cout << "Generating initial values, current progress = " << 100.0*time/stopTime << "\%." << endl;

		// Apply stimulus according to user configuration
		if(time - floor(time/BCL)*BCL >= stimStart && 
		   time - floor(time/BCL)*BCL < stimStart + stimDuration)	
		{
			cell->setIstim(stimStrength);
		}
		else
		{
			cell->setIstim(0.0);
		}

		// update all states and currents
		cell->update();

		// write file for each 1ms 
		if(time + dt > stopTime) // 50*dt = 1ms once
		{
			cell->outputAllStates(initfile);
			cout <<"Initial values of all states have been generated!" << endl;
		}
	}
	fclose(initfile);
	return 0;
}




