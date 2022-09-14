#
# General Code Structure (GCS) for Makefile
# 
# A. Code sructure is re-organised using standard C++ to fit into my project.
# B. The friendly style makes the project easier understood.
# C. The version is more extendable whatever to further 1D,2D or another single cell model. 
# 
# Under Intellectual Property Protection.
# 
# 
# Author      : Shugang Zhang <zsg@ouc.edu.cn>
# Last update : 14-09-2022
# 

all: SingleCell ERP VoltageClampIKr OenD RestCurveS1S2 OneD_VW CVRestCurve CVRestCurve_PCL TwoD Initialization

common = SingleCell/TP06.cc  SingleCell/Cell.cc

CC  	=	g++

CFLAGS	=	-w -O3 #-g:warning output to screen   -w:warning ignored

CFLAGS2	=	-fopenmp

CFLAGS3 =   -arch sm_61 -Xptxas -dlcm=cg


SingleCell: $(common) SingleCell.cc
	$(CC) $(CFLAGS) -o model_single_cell $(common) SingleCell.cc

Initialization: $(common) Initialization.cc
	$(CC) $(CFLAGS) -o model_initialization $(common) Initialization.cc

VoltageClampIKr: $(common) VoltageClampIKr.cc
	$(CC) $(CFLAGS) -o model_voltage_clamp_ikr $(common) VoltageClampIKr.cc

ERP: $(common) ERP.cc
	$(CC) $(CFLAGS) -o model_ERP $(common) ERP.cc

OneD: $(common) OneD.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_oned $(common) OneD.cc

OneD_VW: $(common) OneD_VW.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_oned_vw $(common) OneD_VW.cc

TissueERP: $(common) TissueERP.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_tissue_erp $(common) TissueERP.cc

RestCurveS1S2: $(common) RestCurveS1S2.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_rest_curve_s1s2 $(common) RestCurveS1S2.cc

RestCurve_PCL: $(common) RestCurve_PCL.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_rest_curve_pcl $(common) RestCurve_PCL.cc

CVRestCurve: $(common) CVRestCurve.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_cv_rest_curve $(common) CVRestCurve.cc

CVRestCurve_PCL: $(common) CVRestCurve_PCL.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_cv_rest_curve_pcl $(common) CVRestCurve_PCL.cc

TwoD: $(common) TwoD.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_twod $(common) TwoD.cc

TwoD_EPI: $(common) TwoD_EPI.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_twod_epi $(common) TwoD_EPI.cc

clean:
	rm model_*
