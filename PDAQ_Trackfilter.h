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
    
    TH1F* h_cluster_size =  new TH1F ( "h_cluster_size", "h_cluster_size;Number of hits in the cluster", 110, -10, 100 );
    
    TH2F* h_LTvsLayer0 = new TH2F ( "h_LTvsLayer0", "h_LTvsLayer0;Layer;Lead Time", 8, 0, 8, 10000, -10000,0 );
    TH2F* h_LTvsLayer1 = new TH2F ( "h_LTvsLayer1", "h_LTvsLayer1;Layer;Lead Time", 8, 0, 8, 10000, -10000,0 );
    
    TH2F* h_drifttimevstot = new TH2F("h_drifttimevstot", "h_drifttimevstot;Drift Time [ns];TOT [ns]", 360, -10, 350,350,0,350);
    
    TH1F* h_tot0 = new TH1F("h_tot0", "h_tot0;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot1 = new TH1F("h_tot1", "h_tot1;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot2 = new TH1F("h_tot2", "h_tot2;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot3 = new TH1F("h_tot3", "h_tot3;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_tot4 = new TH1F("h_tot4", "h_tot4;Time Over Threshold [ns]", 710, -10, 700);
    TH1F* h_double_hit = new TH1F("h_double_hit", "h_double_hit;Number of hits", 5, 0, 5);
    
    TH1F* h_drifttime = new TH1F("h_drifttime", "h_drifttime;Drift Time [ns]", 710, -10,700 );
    
    TH2F* h_drifttimevsplane = new TH2F("h_drifttimevsplane", "h_drifttimevsplane;Drift Time [ns];Layer", 360, -10, 350,20,0,10);
    
    TH1F* h_scint_mult = new TH1F("h_scint_mult", "h_scint_mult;Number of hits", 11, -1, 10);
    TH1F* h_scint_timediff = new TH1F("h_scint_timediff", "h_scint_timediff;Number of hits", 5010, -10, 5000);
    TH1F* h_raw_leadtimes = new TH1F("h_raw_leadtimes", "h_raw_leadtimes;", 100000, -100000, 0);
    TH1F* h_TRB_ref_diff = new TH1F ( "h_TRB_ref_diff", "h_TRB_ref_diff;Time diff [ns]", 1000, -5, 5 );



};

bool f_sttHitCompareLeadTime(SttHit* a, SttHit* b);

bool f_sttHitCompareCell(SttHit* a, SttHit* b);

std::vector<SttHit*> GetPairs(std::vector<SttHit*> vec_get_pairs);

std::vector<std::vector<SttHit*>> clusterfinder(std::vector<SttHit*> vec_flayer);


void add_tuples(const std::vector<std::vector<std::vector<SttHit*>>>& vectors,
                std::size_t pos, std::vector<std::vector<SttHit*>> prefix,
                std::vector<std::vector<std::vector<SttHit*>>>& result);

std::vector<std::vector<std::vector<SttHit*>>>
make_tuples(const std::vector<std::vector<std::vector<SttHit*>>>& vectors);


bool PDAQ_Event_Finder(std::vector<SttHit*> vec_stthits, int i,
                       TTree* PDAQ_tree, Stt_Track_Event* stt_event,
                       MFTGeomPar* ftGeomPar, PandaSubsystemSCI* SCI_CAL,  histograms* h);

#endif














