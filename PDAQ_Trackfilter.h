#ifndef PDAQ_TRACK_FILTER_H
#define PDAQ_TRACK_FILTER_H

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include <TF1.h>
#include <TGraph.h>
#include <TLinearFitter.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <cctype>
#include <sstream>
#include <stdio.h>

#include "SciHit.h"
#include "SttEvent.h"
#include "SttHit.h"
#include "SttRawHit.h"
#include "SttTrackEvent.h"

#include "panda_stt_cal.h"
#include "panda_stt_track.h"
#include "panda_subsystem.h"
#include "panda_subsystem_sb.h"
#include "panda_subsystem_sci.h"
#include "panda_subsystem_stt.h"

#include "FTGeo.h"
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>

#include "TObject.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

struct histograms{
  
    TCanvas *HitMultiplicity = new TCanvas ( "HitMultiplicity","HitMultiplicity" );
    TCanvas *TimeOverThreshold = new TCanvas ( "TimeOverThreshold","TimeOverThreshold" );
	
    TH1F* h_hitmultiplicity0 = new TH1F ( "h_hitmultiplicity0", "h_hitmultiplicity0;Number of hits", 110, -10, 100 );
    TH1F* h_hitmultiplicity1 = new TH1F ( "h_hitmultiplicity1", "h_hitmultiplicity1;Number of hits", 110, -10, 100 );
    TH1F* h_hitmultiplicity2 = new TH1F ( "h_hitmultiplicity2", "h_hitmultiplicity2;Number of hits", 110, -10, 100 );
    TH1F* h_hitmultiplicity3 = new TH1F ( "h_hitmultiplicity3", "h_hitmultiplicity3;Number of hits", 110, -10, 100 );
    TH1F* h_hitmultiplicity4 = new TH1F ( "h_hitmultiplicity4", "h_hitmultiplicity4;Number of hits", 110, -10, 100 );
    
    TH1F* h_pLayer_eff = new TH1F ( "h_pLayer_eff", "h_pLayer_eff;Number of fired straws", 30, 0, 30 );
    TH1F* h_Layer_eff = new TH1F ( "h_Layer_eff", "h_Layer_eff;Number of fired straws", 30, 0, 30 );
    TH1F* h_Layer_eff2 = new TH1F ( "h_Layer_eff2", "h_Layer_eff2;Number of fired straws", 30, 0, 30 );
    TH1F* h_sq_ch = new TH1F ( "h_sq_ch", "h_sq_ch;channel", 300, 0, 300 );

    
    TH1F* h_cluster_size =  new TH1F ( "h_cluster_size", "h_cluster_size;Number of hits in the cluster", 110, -10, 100 );
    
    TH2F* h_LTvsLayer0 = new TH2F ( "h_LTvsLayer0", "h_LTvsLayer0;Layer;Lead Time", 8, 0, 8, 10000, -10000,0 );
    TH2F* h_LTvsLayer1 = new TH2F ( "h_LTvsLayer1", "h_LTvsLayer1;Layer;Lead Time", 8, 0, 8, 10000, -10000,0 );
    
    TH2F* h_drifttimevstot = new TH2F("h_drifttimevstot", "h_drifttimevstot;Drift Time [ns];TOT [ns]", 700, 0, 700,700,0,700);
    
    TH2F* h_p_drifttimevstot = new TH2F("h_p_drifttimevstot", "h_p_drifttimevstot;Drift Time [ns];TOT [ns]", 500, 0, 500,700,0,700);
    
    TH1F* h_tot0 = new TH1F("h_tot0", "h_tot0;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot1 = new TH1F("h_tot1", "h_tot1;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot2 = new TH1F("h_tot2", "h_tot2;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot3 = new TH1F("h_tot3", "h_tot3;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot4 = new TH1F("h_tot4", "h_tot4;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_double_hit = new TH1F("h_double_hit", "h_double_hit;Number of hits", 5, 0, 5);
    
    TH1F* h_p_drifttime = new TH1F("h_p_drifttime", "h_p_drifttime;Drift Time [ns]", 710, -10,700 );
    TH1F* h_drifttime = new TH1F("h_drifttime", "h_drifttime;Drift Time [ns]", 800, -100,700 );
  
    TH1F* h_p_tot = new TH1F("h_p_tot", "h_p_tot;Time Over Threshold [ns]", 710, -10,700 );
    TH1F* h_tot = new TH1F("h_tot", "h_tot;Time Over Threshold [ns]", 710, -10,700 );
    
    
    TH2F* h_drifttimevsLayer = new TH2F("h_drifttimevsLayer", "h_drifttimevsLayer;Drift Time [ns];Layer", 360, -10, 350,20,0,10);
    
    TH1F* h_scint_timediff = new TH1F("h_scint_timediff", "h_scint_timediff;Number of hits", 100000, -50000, 50000);
    TH1F* h_scint_timediffa = new TH1F("h_scint_timediffa", "h_scint_timediffa;Number of hits", 50000, 0, 50000);
    TH1F* h_scint_timediffb = new TH1F("h_scint_timediffb", "h_scint_timediffb;Number of hits", 50000, 0, 50000);

    TH1F* h_raw_leadtimes = new TH1F("h_raw_leadtimes", "h_raw_leadtimes;", 100000, -100000, 0);
    
    TH1F* h_drifttimeTRB1 = new TH1F("h_drifttimeTRB1", "h_drifttimeTRB1;Drift Time [ns]", 710, -10,700 );
    TH1F* h_drifttimeTRB2 = new TH1F("h_drifttimeTRB2", "h_drifttimeTRB2;Drift Time [ns]", 710, -10,700 );
    
    TH1F* h_scnt_plup_cnts =  new TH1F ( "h_scnt_plup_cnts", "h_scnt_plup_cnts;Number of hits ", 3, 0, 3 );
    TH1F* h_scint_mult_b = new TH1F("h_scint_mult_b", "h_scint_mult_b;Number of hits", 25, 0, 25);
    TH1F* h_scint_mult_a = new TH1F("h_scint_mult_a", "h_scint_mult_a;Number of hits", 25, 0, 25);

    
    TH1F* h_L5_S9LTdiff = new TH1F("h_L5_S9LTdiff", "h_L5_S9LTdiff;Time diff [ns]", 1100, -100, 1000);
    TH1F* h_pL5_S9LTdiff = new TH1F("h_pL5_S9LTdiff", "h_pL5_S9LTdiff;Time diff [ns]", 1100, -100, 1000);

    
    TH2F* h_pLayerMult = new TH2F("h_pLayerMult", "h_pLayerMult;Layer;Multiplicity", 8, 0, 8,75,0,75);
    TH2F* h_LayerMult = new TH2F("h_LayerMult", "h_LayerMult;Layer;Multiplicity", 8, 0, 8,75,0,75);
    TH1F* h_0LMultiplicity = new TH1F("h_0LMultiplicity", "h_0LMultiplicity;Multiplicity", 10, 0, 10);
    TH1F* h_pLMultiplicity = new TH1F("h_pLMultiplicity", "h_pLMultiplicity;Multiplicity", 10, 0, 10);
    TH1F* h_LMultiplicity = new TH1F("h_LMultiplicity", "h_LMultiplicity;Multiplicity", 10, 0, 10);
    
    TH1F* h_hitBlock = new TH1F("h_hitBlock", "h_hitBlock;Number of hits", 5, 0, 5);
    TH1F* h_ShitBlock = new TH1F("h_ShitBlock", "h_ShitBlock;Number of hits", 5, 0, 5);
    
    TH2F* h_ppL_dtvstot[8];
    
    TH1F* h_pL_layerDT[8];
    TH2F* h_pL_dtvstot[8];
    TH1F* h_pL_TOT[8];
    TH2F* h_pL_channel_mult[8];
    TH1F* h_pChDiff[8];
    TH1F* h_straw[8];
    
    
    TH1F* h_L_layerDT[8];    
    TH2F* h_L_dtvstot[8];
    TH1F* h_L_TOT[8];
    TH2F* h_L_channel_mult[8];
    TH1F* h_Ch_Dt[256];
    TH1F* h_pLT_Diff[256];
//     for ( int i = 0; i < 8; i++ ) {
//       
//        h_layerDT[i] = new TH1F ( Form ( "Layer_%d_DT", i ) , Form ( "Layer_%d_DT", i ), 600, -100, 500 );
//     }
    

    
};

typedef std::vector<SttHit> VecSttHit;
bool f_sttHitCompareLeadTime(SttHit a, SttHit b);

bool f_sttHitCompareCell(SttHit a, SttHit b);

VecSttHit GetPairs(VecSttHit vec_get_pairs);

std::vector<VecSttHit> clusterfinder(VecSttHit vec_flayer);


void add_tuples(const std::vector<std::vector<VecSttHit>>& vectors,
                std::size_t pos, std::vector<VecSttHit> prefix,
                std::vector<std::vector<VecSttHit>>& result);

std::vector<std::vector<VecSttHit>>
make_tuples(const std::vector<std::vector<VecSttHit>>& vectors);


bool PDAQ_Event_Finder(VecSttHit vec_stthits, int i,
                       TTree* PDAQ_tree, Stt_Track_Event* stt_event,
                       MFTGeomPar* ftGeomPar, PandaSubsystemSCI* SCI_CAL,  histograms* h);

#endif















