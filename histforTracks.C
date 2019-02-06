#include <fstream>
#include <TF1.h>
#include <TLinearFitter.h>
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <TGraph.h>
#include <math.h>
#include <cstdlib>
#include <TLinearFitter.h>


#include "SttRawHit.h"
#include "SttHit.h"
#include "SttTrackHit.h"
#include "SttEvent.h"
#include "SttTrackEvent.h"

                          
#include "panda_subsystem.h"
#include "panda_subsystem_stt.h"
#include "panda_subsystem_sb.h"
#include "panda_stt_cal.h"
#include "panda_stt_track.h"

using namespace std;


Bool_t histforTracks(void)
{
  cout<<"Open"<<endl;
  
  PandaSubsystemSTT* STT_RAW = 0;
  PandaSubsystemSTT* CAL = 0;
  PandaSttCal* STT_CAL = 0;
  PandaSttTrack* STT_TRACK = new PandaSttTrack();
  
  
    TFile * file = TFile::Open("c.root", "READ");

    TTree* tree = 0;
    file->GetObject("PDAQ_tree", tree);
    if (!tree) {
        std::cerr << "Tree doesn't exists" << std::endl;
        return 1;
    }
    tree->Print();
 
    std::vector<Double_t> vec_event;
    std::vector<Double_t> vec_test;
    std::vector<Double_t> vec_driftTime;
    std::vector<Double_t> vec_roundoff;
    std::vector<Double_t> vec_occurance;
    std::vector<Double_t> vec_pos_DT;
    std::vector<Double_t> vec_cumsum;
    std::vector<Double_t> vec_drift_radius;
    std::vector<SttHit*> vec_All_tracks;
    
    std::vector<double> vec_o_test;
    std::vector<double> vec_o_test1;

    std::vector<double> vec_x;
    std::vector<double> vec_y;
    std::vector<double> vec_z;
    std::vector<double> vec_layer;
    std::vector<double> vec_module;
    std::vector<double> vec_fee;
    std::vector<double> vec_fee_ch;
    std::vector<double> vec_tdc_ch;
    
    int counterofdt =0;
    int driftTimeCounter2 =0;

  
    tree->SetBranchAddress("STT_TRACKS", &STT_TRACK);
    TFile* ftree = new TFile("histforTracks.root", "RECREATE");
    // TTree* DR_Tree = new TTree("DR_Tree", "DR_Tree");
    
    // DR_Tree->Branch("vec_Driftradius", &vec_o_test);
    // DR_Tree->Branch("vec_x", &vec_x);
    // DR_Tree->Branch("vec_y", &vec_y);
    // DR_Tree->Branch("vec_z", &vec_z);
    // DR_Tree->Branch("vec_layer", &vec_layer);

    TH2F* h_track[100];

    
    Int_t iev = (Int_t)tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl
         << endl;
	 
	 cout<<STT_TRACK->stt_track_can.total_track_NTDCHits<<endl;
	
    for (Int_t i = 0; i < iev; i++) 
    {
        tree->GetEntry(i);
		if (i%100==0){
		  cout<<i<<endl;
		}
		
		for (int n = 0; n < STT_TRACK->stt_track_can.total_track_NTDCHits; n++)
		{
	  		h_track[i] = new TH2F (Form("Track No %d - Straw Vs Layer",i),Form ("Track %d Straw VsLayer; Straw;Layer",i ),18.0,-1.0,17.0,20,0.0,10.0);


		    std::vector<SttHit*> vec_track_can;
		    //SttHit* cal_hit  = (SttHit*)STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt(n);
		    SttTrackHit* track_hit  = (SttTrackHit*)STT_TRACK->stt_track_can.tdc_track_hits->ConstructedAt(n); 	    
		    vec_track_can = track_hit->vec_Track;
		    	    
		    for (int t =0; t< vec_track_can.size(); t++)
		    {
				// vec_o_test1.push_back(vec_track_can[t]->drifttime);
				// vec_o_test.push_back(vec_track_can[t]->drifttime);
				// vec_x.push_back(vec_track_can[t]->x);
				// vec_y.push_back(vec_track_can[t]->y);
				// vec_z.push_back(vec_track_can[t]->z);
				// vec_layer.push_back(vec_track_can[t]->layer);
		    	if (vec_track_can[t]->straw%2 ==0)
		    	{
		    		h_track[i]->Fill((floor(vec_track_can[t]->straw)/2)-1,vec_track_can[t]->layer+0.5);
		    	}
		    	else
		    		h_track[i]->Fill(vec_track_can[t]->straw/2,vec_track_can[t]->layer);

				counterofdt++;		
	      
		    }
		    	h_track[i]->Write();
		 	  //DR_Tree->Fill();
			  // vec_o_test.clear();
			  // vec_x.clear();
			  // vec_y.clear();
			  // vec_z.clear();

	    
		    	  
		} 

	  
      }

    return kTRUE;
    
 //////////////////// 
}