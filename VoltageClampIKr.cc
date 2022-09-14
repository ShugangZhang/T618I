/*
 * General Code Structure (GCS) for single cell simulation
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

using namespace std;

// ONLY ONE LINE CAN BE KEPT HERE
// #define SAN
// #define ATRIA
#define VENTRICLE

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.
	
	double BCL = 800; // ms
	double numS1 = 5;
	double dt = 0.02; //ms
	double stopTime = numS1*BCL; //ms
	double stimStrength = -0.0; // pA/pF   -6.0pA/pF(-0.6nA) for RatAtrial // -12.5pA/pF for CaMKII; -8.78pA for neonatalRatAtrial
	double stimDuration = 5;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	// typedef GrandiCaMKII CellType;
	// typedef GrandiCaMKII* CellPointer;
	// typedef NeonatalRatAtria CellType;
	// typedef NeonatalRatAtria* CellPointer;
	typedef TP06 CellType;
	typedef TP06* CellPointer;
	// typedef ORdHumanVentricle CellType;
	// typedef ORdHumanVentricle* CellPointer;

	// statistics for single cell
	double apd20;
	double apd25;
	double apd50;
	double apd80;
	double apd90;
	double dvdt_max;
	double rest_potential;
	double amplitude;
	double overshoot;


	// --------start simulation--------
	// note that constructor contains the initializer
	CellPointer cell = new CellType(MCELL);

	FILE *datafile = fopen("Outputs/VoltageClampResults.dat","w+");
	FILE *ivfile = fopen("Outputs/IVcurve_Drug.dat","w+");

	
	// uncomment following two lines to read in initial values (this is because the original init values is not stable yet)
	// if the initfile is not available, run the initialization.cc first
	// FILE *initfile = fopen("SingleCell/NeonatalRatAtriaInitialValues.dat","r");
	// cell->readinAllStates(initfile);

	// apply user configuration about dt
	cell->setDt(dt);

	double time = 0;
	int step = 0;
	double tail = 0;
	double v;
	for(double i = -49; i <= 62 ; i += 2)
	{
		if(i == 15) v = 15.0001;
		else v = i;


		cell->init(MCELL);

		time = 0;
		tail = 0;
		step = 0;
		cout << "Test potential = " << v << "; ";

		for(time = 0.0, step = 0; time <= 400+3500+3500; time += dt, step++)
		{
			if(time <= 400) cell->setVolt(-80);
			else if(time <= 400 + 3500) cell->setVolt(v);
			else cell->setVolt(-50);

			cell->update();
			if(step%(int(0.1/dt)) == 0) 
			{
				fprintf(datafile,"%4.10f\t", time);
				fprintf(datafile,"%4.10f\t", cell->getV()); // unit: mV
				fprintf(datafile,"%4.10f\t", cell->getIKr()); // unit: mV
				fprintf(datafile,"\n");
			}


			if(time <= 400+3500 && time+dt > 400+3500)
			{
				tail = cell->getIKr();
			}
		}	

		cout << "Tail = " << tail << "." << endl;

		fprintf(ivfile,"%4.10f\t", v);
		fprintf(ivfile,"%4.10f\t", tail);
		fprintf(ivfile,"\n");
	}

	
	fclose(datafile);
	fclose(ivfile);
	return 0;
}




