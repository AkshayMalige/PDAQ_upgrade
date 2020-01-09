#ifndef PDAQ_RES_H
#define PDAQ_RES_H

#include "TFile.h"
#include <TGraph.h>
#include <TH1F.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "SttEvent.h"
#include "SttHit.h"
#include "TF1.h"
#include "TRandom.h"
#include "TTree.h"
#include <TH2D.h>
#include <TLinearFitter.h>
#include <TNamed.h>
#include <TObject.h>
#include <cstdlib>
#include <math.h>

#include "TCanvas.h"
#include <TChain.h>
#include <TROOT.h>
#include "TObject.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TLine.h"

using namespace std;

//Bool_t PDAQ_Spl_Res(void);

Bool_t PDAQ_Res(char* intree, char* outtree, int maxEvents);

struct histograms{
    TH1F* hx = new TH1F ( "hx", "Spatial_Resolution_X", 500, -5, 5 );
    TH1F* hdx = new TH1F ( "hdx", "res_non_corrected", 500, -5, 5 );
    TH1F* h_dr = new TH1F ( "h_dr", "drift_radius", 350, 0, 350 );
    TH1F* hx_theeta = new TH1F ( "h_theeta_X", "h_theeta_X", 400, -2, 2 );    
    TH1F* hx_slope = new TH1F ( "hx_slope", "Slope_X", 1000, -0.25, 0.25 );
    TH1F* hx_const = new TH1F ( "hx_const", "Constant_X", 1000, 0, 50 );
    TH2F* h_XfvsZ = new TH2F ( "h_XfvsZ", "h_XfvsZ;X of track fit X_{f} [mm]; Z [mm]", 200, 150, 350,7,-10,600 );
    TH2F* Dt_vs_dr = new TH2F ("Dt_vs_dr", "Dt_vs_dr;Drift time [ns];dr [mm]",200,0,200,100,-1,1 );
    TH1F* hx_chi = new TH1F ("hx_chi","hx_chi",1000,-0.05,0.05);
    
    TH1F* h_dx[18];
    TH1F* h_res_plane[18];
    TH2F* h_Dr_vs_dr[18];
    TH1F* h_plane_str_no[18];
    TH2F* h_dt_vs_dr[18];
    TH2F* h_Dt_vs_dr_Lay[8];
    TH1F* h_ch_sq_od6;
    TH1F* h_cal_chi = new TH1F ( "h_cal_chi", "h_cal_chi", 50, 0, 50 );;
};



struct vec{
    std::vector<double>* vec_Drifttime = 0;
    std::vector<double>* vec_x = 0;
    std::vector<double>* vec_y = 0;
    std::vector<double>* vec_z = 0;
    std::vector<double>* vec_layer = 0;
    std::vector<double>* vec_straw = 0;
    std::vector<double>* vec_plane = 0;
    std::vector<double>* vec_tot = 0;
    
};


#endif
