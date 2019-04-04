#include <stdio.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"

#include "base/EventProc.h"
#include "base/Event.h"
#include "base/ProcMgr.h"
#include "hadaq/TdcSubEvent.h"
#include "hadaq/definess.h"
#include "hadaq/HldProcessor.h"
#include "hadaq/TrbProcessor.h"
#include "hadaq/TdcProcessor.h"

#include <sstream>
#include <string>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include <fstream>

#include <algorithm>
#include <TError.h>
#include <TRint.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TView.h>
#include <TGeoManager.h>
#include <TGeoTrack.h>
#include <TVirtualGeoTrack.h>
#include <TLine.h>


#include <TFile.h>
#include <TMath.h>
#include <TH2F.h>
#include <TTree.h>
#include <TChain.h>

#include <TStopwatch.h>
#include <fstream>
#include <stdio.h>

#define NUMBER_OF_TDCS 8
#define CHANNELS_OFFSET 49

#define MAX_HITS 50
#define CHANNELS 300
#define LAYERS 8

using namespace std;

vector<string> &split(const string &s, char delim, vector<string> &elems) {

        stringstream ss(s);

        string item;

        while(getline(ss, item, delim)) {

                elems.push_back(item);

        }

        return elems;

}

 

vector<string> split(const string &s, char delim) {

        vector<string> elems;

        split(s, delim, elems);

        return elems;

}



struct detLoc {
  detLoc(){};
 
  detLoc(int p_station, int p_layer, int p_straw){
    layer = p_layer;
    station = p_station;
    straw = p_straw;
  }

  int layer;
  int station;
  int straw;
};

class Detector {
    public:
        map<int, detLoc> detMap;

        Detector();

        detLoc GetDetectorLocFromTDCChannel(int channel);

};

class LPetProcessor : public base::EventProc {

    private:
        int internalEventCtr;

    public:

        Detector* det;

        vector<hadaq::TdcSubEvent*> tdcs;
       
        TH1F* h_channelMult;
        TH2F* h_leadTimeVsChannel;
        TH1F* h_Tot;
        TH2F* h_TotVsChannel;
        TH1F* h_leadTime;
        TH1F* h_trailTime;
        TH2F* h_leadTimeVsTot;
        TH2F* h_channelVsHits;

        TH2D* h_Ft_Geo;
        TH2D* h_FT_Geo;
        TH2D* h_FT_Geo2;
        TH2F* hL_DT_Mult[LAYERS];

        TH2F* h_layerChannelVsLeadTime[LAYERS];
        TH2F* h_layerChannelVsTot[LAYERS];

        TH1F* h_refTimeTRB1;
        TH1F* h_refTimeTRB2;
        TH1F* h_refTimeTDC[6];
        TH1F* h_layerMultiplicity;
        TH1F* h_driftTime;
        TH1F* h_driftTimeTRB1;
        TH1F* h_driftTimeTRB2;
        TH2F* h_driftTimeVsChannel;
        TH2F* hL_driftTimeVsChannel[LAYERS];
        TH2F* hL_driftTimeVsTOT[LAYERS];
        TH2F* h_driftTimeVsTOT;

	TH2F* h_geo;
        TH1F* h_Tot_All;

        LPetProcessor() : base::EventProc("LPET")
        {



    	det = new Detector();

	std::map<UInt_t,int> stt_channel_offsets;
	stt_channel_offsets[6400] = 0 * 49;
	stt_channel_offsets[6410] = 1 * 49;
	stt_channel_offsets[6411] = 2 * 49;
	stt_channel_offsets[6420] = 3 * 49;
	stt_channel_offsets[6430] = 4 * 49;
	stt_channel_offsets[6431] = 5 * 49;
	stt_channel_offsets[6500] = 6 * 49;



    	fstream infile("/home/pandastraws/go4/new_release/trb3/beamtime201902/tdc_map_cosysetup.txt");

	string line;

	stringstream ss;

	//vector <vector <int> > container; 

	if (!infile) { cout<<"Error reading file"<<endl; exit(0); }

	while(getline(infile, line)) {

	    vector<string> elems = split(line, '\t');
	    vector<int> int_elems;
	    for (int a=0; a<elems.size(); a++)
		{
			int num = atoi(elems.at(a).c_str());
			int_elems.push_back(num);
		}
//tdc channel elems[0] + offset
	//detMap[int_elems.at(0) + offset  ] = detLoc(int_elems.at(0), 1, 1);
	

		det->detMap[stt_channel_offsets[int_elems.at(0)]+int_elems.at(1)] = detLoc(int_elems.at(2), int_elems.at(3), int_elems.at(4));

	    //container.push_back(int_elems);

	}


            internalEventCtr = 0;
           
            tdcs = vector<hadaq::TdcSubEvent*>();   
    //TH1F* h1_width = new TH1F("h1_width","Baseline_scan;Width [mV]",32,0,32);

            h_channelMult = new TH1F("Channel_Mult", "channelMult;Channel No", 401, 0, 400);
            h_Tot = new TH1F("TOT", "TOT;Time Over Threshold [ns]", 1000, 0, 1000);
            h_leadTimeVsChannel = new TH2F("LeadTime_vs_Channel", "LeadTime_vs_Channel;Lead time;Channel No", 10000, -10000, 10000, 300, 0, 300);
            h_TotVsChannel = new TH2F("TOT_vs_Channel", "TOT_vs_Channel;Time over threshold [ns];Channel No", 600, 0, 600, 300, 0, 300);
            h_leadTime = new TH1F("leadTime", "leadTime;Lead time", 10000, -10000, 10000);
            h_layerMultiplicity = new TH1F("Layer_Mult", "Layer;Hits", 10, 0, 10);
            h_trailTime = new TH1F("Trail_Time", "Trail_Time; Trail time", 10000, -10000, 10000);
            h_leadTimeVsTot = new TH2F("LeadTime_vs_TOT", "LeadTime_vs_TOT;Leadtime;Time over threshold [ns]", 10000, -10000, 1000, 300, 0, 300);
            h_channelVsHits = new TH2F("Channels_vs_Hits", "Channels_vs_Hits;channel no;No.hits/event", 300, 0, 300, 70, -35, 35);
            h_Tot_All = new TH1F("TOT_All", "TOT_All;Time Over Threshold All[ns]", 1000, 0, 1000);

            
            
            
            h_Ft_Geo = new TH2D("h_Ft_Geo", "h_Ft_Geo;Straw No;Sub Layer no", 18, -1, 17, 52, -1, 51);
            h_FT_Geo = new TH2D("h_FT_Geo", "h_FT_Geo;Straw No;Sub Layer no", 36, -1, 35, 20, 0, 10);
            h_FT_Geo2 = new TH2D("h_FT_Geo2", "h_FT_Geo2;Straw No;check Sub Layer no", 36, -1, 35, 20, 0, 10);
            
          
            
            h_driftTime = new TH1F("DriftTime", "DriftTime; Drift Time", 5000, -50000, 50000);
            h_driftTimeTRB1 = new TH1F("DriftTime_TRB1", "DriftTime_TRB1; Drift Time", 5000, -50000, 50000);
            h_driftTimeTRB2 = new TH1F("DriftTime_TRB2", "DriftTime_TRB2; Drift Time", 5000, -50000, 50000);
            h_driftTimeVsChannel = new TH2F("DriftTime_vs_Channel", "Drift_Time_vs_Channel;Drift Time;Channel", 500, 0, 500, 300, 0, 300);
            h_driftTimeVsTOT = new TH2F("DriftTime_vs_TOT", "DriftTime_vs_TOT;Drift Time;Time Over Threshold", 500, 0, 500, 300, 0, 300);

            for (int i = 0; i < LAYERS; i++) {
                h_layerChannelVsLeadTime[i] = new TH2F(Form("Layer%d_Channel_vs_LeadTime", i + 1), Form("Layer_%d_Channel_vs_LeadTime;Lead time;Channel No", i + 1), 10000, -10000, 10000, 32, 0, 32 );
                h_layerChannelVsTot[i] = new TH2F(Form("Layer%d_Channel_vs_TOT", i + 1), Form("Layer%d_Channel_vs_TOT;Time Over Threshold [ns];Channel No", i + 1), 600, 0, 400, 32, 0, 32 );
                hL_driftTimeVsChannel[i] = new TH2F(Form("Layer%d_DriftTime_vs_Channel",i+1),Form( "Layer%d_DriftTime_vs_Channel;Drift Time;Channel",i+1), 500, 0, 500, 32, 0, 32);
                hL_driftTimeVsTOT[i] = new TH2F(Form("Layer%d_DriftTime_vs_TOT",i+1),Form( "Layer%d_DriftTime_vs_TOT;Drift Time;Time Over Threshold",i+1), 500, 0, 500, 600, 0, 600);


		hL_DT_Mult[i] = new TH2F(Form("Layer%d_DriftTime_Mult",i+1),Form( "Layer%d__DriftTime_Mult;Channel;Multiplicity",i+1), 33, 0, 33, 10, 0, 10);

            }

            h_refTimeTRB1 = new TH1F("Ref_Time_TRB1", "Ref_Time_TRB1", 500, -25, 25);
            h_refTimeTRB2 = new TH1F("Ref_Time_TRB2", "Ref_Time_TRB2", 500, -25, 25);

            for (int i = 0; i < 6; i++) {
                h_refTimeTDC[i] = new TH1F(Form("Ref_Time_TDC%d", i) , Form("Ref_Time_TDC%d", i), 500, -2500, 2500);
            }

		h_geo = new TH2F("geo", "geo", 32, 0, 32, 9, 0, 9);

            RegisterObject(h_channelMult, "general");
            RegisterObject(h_Tot, "general");
            RegisterObject(h_TotVsChannel, "general");
            RegisterObject(h_leadTimeVsChannel, "general");
            RegisterObject(h_leadTime, "general");
            RegisterObject(h_trailTime, "general");
            RegisterObject(h_leadTimeVsTot, "general");
            RegisterObject(h_channelVsHits, "general");

RegisterObject(h_geo, "general");

            RegisterObject(h_Ft_Geo, "general");
            RegisterObject(h_FT_Geo, "general");
            RegisterObject(h_FT_Geo2, "general");
            RegisterObject(h_layerMultiplicity, "general");
            RegisterObject(h_driftTime, "general");
            RegisterObject(h_driftTimeTRB1, "general");
            RegisterObject(h_driftTimeTRB2, "general");
            RegisterObject(h_driftTimeVsChannel, "general");
            RegisterObject(h_driftTimeVsTOT, "general");
            RegisterObject(h_Tot_All, "general");


            for (int i = 0; i < LAYERS; i++) {
                RegisterObject(h_layerChannelVsLeadTime[i], "LayerInfo");
                RegisterObject(h_layerChannelVsTot[i], "LayerInfo");
                RegisterObject(hL_driftTimeVsChannel[i], "LayerInfo");
                RegisterObject(hL_driftTimeVsTOT[i], "LayerInfo");
	        RegisterObject(hL_DT_Mult[i], "LayerInfo");
                     
            }

            for (int i = 0; i < 6; i++) {
                RegisterObject(h_refTimeTDC[i], "tdcs");
            }

            RegisterObject(h_refTimeTRB1, "tdcs");
            RegisterObject(h_refTimeTRB2, "tdcs");
        }



        virtual bool Process(base::Event* ev)
        {

            internalEventCtr++;

//cout<<"Event: "<<internalEventCtr<<endl;
            
            if (internalEventCtr % 10000 == 0) cout<<"Event: "<<internalEventCtr<<endl;

            tdcs.clear();
            // order of TRBs in this vector determines the channels offsets
/*            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_E100"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_E101"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_E102"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_E200"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_E201"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_E202"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_E203"))) ;*/


            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_6400"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_6410"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_6411"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_6420"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_6430"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_6431"))) ;
            tdcs.push_back( dynamic_cast<hadaq::TdcSubEvent*> (ev->GetSubEvent("TDC_6500"))) ;


            // reference times
            double refTime[NUMBER_OF_TDCS];

            double leadTimes[CHANNELS][MAX_HITS];
            double trailTimes[CHANNELS][MAX_HITS];
            int leadTimesNum[CHANNELS];
            int trailTimesNum[CHANNELS];


            double Tot[CHANNELS][MAX_HITS];
            int TotNum[CHANNELS];


            detLoc loc;
            double hit_time =0;
            double scint_time =0;
            double scint_ref =0;
            double ref_diff =0;
            double drifttime =0;
            bool scint_hit =false;
            
            
            for (int i = 0; i < CHANNELS; i++) {
                leadTimesNum[i] = 0;
                trailTimesNum[i] = 0;
                TotNum[i] = 0;
            }

	        int hitMultOnLayer[LAYERS+1];
		int dtHitsMult[LAYERS + 1][CHANNELS+1];
            for (int l = 0; l <= LAYERS; l++) {
                hitMultOnLayer[l] = 0;
                for (int m =0; m<= CHANNELS; m++)
                {
                    dtHitsMult[l][m] = 0;
                }
            }

		//layerMult.clear();


            // gather hit times
            int tdc_ptr = 0;
            vector<hadaq::TdcSubEvent*>::iterator tdcs_it;
            for (tdcs_it = tdcs.begin(); tdcs_it != tdcs.end(); tdcs_it++) {
//cout<<"TDC: "<<tdc_ptr<<endl;
                tdc_ptr++;
                if ( (*tdcs_it) == 0 ) continue;
//cout<<"SS TDC: "<<tdc_ptr<<endl;                
                for (unsigned cnt = 0; cnt < (*tdcs_it)->Size(); cnt++) {
                    const hadaq::TdcMessageExt& ext = (*tdcs_it)->msg(cnt);

                    int local_ch = ext.msg().getHitChannel() + CHANNELS_OFFSET * (tdc_ptr - 1);
//cout<<"channel: "<<local_ch<<endl;
                    
                    /*if (ext.msg().getHitChannel() == 37 && tdc_ptr == 2) {
                            cout<<"RAW HIT ON 38"<<endl;
                    }*/

                    if (ext.msg().isHitRisingEdge() == true) {
                        leadTimes[local_ch][leadTimesNum[local_ch]] = ext.GetGlobalTime() * 1e9;
                        leadTimesNum[local_ch]++;
                   
                    }
                    else if (ext.msg().isHitRisingEdge() == false) {
                        trailTimes[local_ch][trailTimesNum[local_ch]] = ext.GetGlobalTime() * 1e9;
                        trailTimesNum[local_ch]++;
                    }



       
                    if (ext.msg().isHitRisingEdge() == true) {
                        if (ext.msg().getHitChannel() == 0) {
                            refTime[tdc_ptr - 1] = ext.GetGlobalTime() * 1e9;
                            
                            //if ((tdc_ptr - 1)==0)cout<<"ref zero"<<refTime[tdc_ptr - 1]<<endl;

                        }
                        h_channelMult->Fill(ext.msg().getHitChannel() + CHANNELS_OFFSET * (tdc_ptr - 1));
                       
                        h_leadTimeVsChannel->Fill((ext.GetGlobalTime() * 1e9), ext.msg().getHitChannel() + CHANNELS_OFFSET * (tdc_ptr - 1));
                   
                        loc = det->GetDetectorLocFromTDCChannel(ext.msg().getHitChannel() + CHANNELS_OFFSET * (tdc_ptr - 1) );
			int tdc_ad = ext.msg().getHitChannel() + CHANNELS_OFFSET * (tdc_ptr - 1);
			
//			printf("hit %d %d %d\n", loc.layer, loc.station, loc.straw);

                        h_layerChannelVsLeadTime[loc.layer - 1]->Fill((ext.GetGlobalTime() * 1e9),  loc.straw);
                     /*   hit_time = (ext.GetGlobalTime() * 1e9);			
                        //cout << "hit time  "<<hit_time<<"\t"<<tdc_ptr<<endl;
                        if (tdc_ptr == 7)
                        {
                            scint_hit = true;
                            scint_ref = refTime[6];
                            scint_time = (ext.GetGlobalTime() * 1e9);
                           // cout << "scint time  "<<scint_time<<"\t"<<tdc_ptr<<endl;

                        }  
                        
                        if (scint_hit == true)
                        {
                           ref_diff =  (hit_time - scint_ref); 
                           drifttime = (ext.GetGlobalTime() * 1e9) - (scint_time) - (ref_diff);
                          // cout << "ref_diff  "<<ref_diff<<"\t"<<drifttime<<endl;

                        }   */                     
                         //cout << drifttime <<endl;                     
                        
                    }

                }
            }
                       
                   
               
        h_refTimeTRB1->Fill(refTime[0] - refTime[1]);
        h_refTimeTRB1->Fill(refTime[0] - refTime[2]);
        h_refTimeTRB2->Fill(refTime[3] - refTime[4]);
        h_refTimeTRB2->Fill(refTime[3] - refTime[5]);
        h_refTimeTRB2->Fill(refTime[3] - refTime[6]);

        h_refTimeTDC[0]->Fill(refTime[0] - refTime[1]);
        h_refTimeTDC[1]->Fill(refTime[0] - refTime[2]);
        h_refTimeTDC[2]->Fill(refTime[0] - refTime[3]);
        
        h_refTimeTDC[3]->Fill(refTime[3] - refTime[4]);
        h_refTimeTDC[4]->Fill(refTime[3] - refTime[5]);
        h_refTimeTDC[5]->Fill(refTime[3] - refTime[6]);

        scint_ref = refTime[6];
        //cout<<" ref 0"<<refTime[0]<<endl;
        
        if (leadTimesNum[295] > 0) {
            scint_time = leadTimes[295][0];
	}

            
            
            
//cout<<"Processing begins"<<endl;
        int tdc = 0;
        for (int k = 0; k < CHANNELS; k++)
        {

            tdc = floor(k / 49);
    
            if (tdc < 3) {
                ref_diff = (scint_ref - refTime[tdc]) - (refTime[3] - refTime[0]);
                //cout<< "trb2\t"<<scint_ref<<"\t"<<refTime[tdc]<<"\t"<<ref_diff<<endl;
            }
            else {
                ref_diff = scint_ref - refTime[tdc];
                //cout<<"trb2-1\t"<<scint_ref<<"\t"<<refTime[tdc]<<"\t"<<ref_diff<<endl;
            }
            

            

            if (k % 49 != 0) {

                loc = det->GetDetectorLocFromTDCChannel(k);
                //cout<<ref_diff<<"\t"<<refTime[tdc]<<endl;

//printf("%d %d\n", loc.straw, loc.layer);
                

                for (int l = 0; l < leadTimesNum[k]; l++)
                {
                    h_leadTime->Fill(leadTimes[k][l]);
                    //cout<<"leadTimes : "<<leadTimes[k][l]<<endl;

                    h_trailTime->Fill(trailTimes[k][l]);



                    if (trailTimesNum[k] == leadTimesNum[k] && leadTimes[k][l] < trailTimes[k][l])  {

//    printf("tot: channel %d hit number %d value %f\n", k, TotNum[k], (trailTimes[k][l] - leadTimes[k][l]));

                        // calculate times
                        Tot[k][TotNum[k]] = (trailTimes[k][l] - leadTimes[k][l]) ;
                        drifttime = ( scint_time - leadTimes[k][l] ) - ref_diff;

                        h_Tot_All->Fill(Tot[k][TotNum[k]]);
                        
                       if ( (drifttime >= 0 and drifttime <= 500) )
                        {
                            h_Tot->Fill(Tot[k][TotNum[k]]);
                        
                        if ( (Tot[k][TotNum[k]] >= 0 and Tot[k][TotNum[k]] <= 600) ) {
                        

                        // fill histograms
                        h_TotVsChannel->Fill(Tot[k][TotNum[k]], k);
                        h_leadTimeVsTot->Fill(leadTimes[k][l], Tot[k][TotNum[k]]);

			hitMultOnLayer[loc.layer]++;
                                   //printf("HERE %d\n",hitMultOnLayer[0]);
              //          printf("tot: %d %d %d %d %d\n", loc.layer, loc.station, loc.straw,  loc.straw);
                                   
                        h_layerChannelVsTot[loc.layer -1]->Fill(Tot[k][TotNum[k]],  loc.straw);
                        
                        if (tdc < 3) {
                            h_driftTimeTRB1->Fill(drifttime);
                        
                            //printf("%f %f %f\n", ref_diff, (scint_time - leadTimes[k][l]), drifttime);
                        }
                        if (tdc > 3) {
                            h_driftTimeTRB2->Fill(drifttime);
                        
                            //printf("%f %f %f\n", ref_diff, (scint_time - leadTimes[k][l]), drifttime);
                        }  
                        h_driftTime->Fill(drifttime);
                        h_driftTimeVsChannel->Fill(drifttime,k);
                        hL_driftTimeVsChannel[loc.layer -1]->Fill(drifttime,loc.straw);
                        hL_driftTimeVsTOT[loc.layer -1]->Fill(drifttime,Tot[k][TotNum[k]]);
                        //cout << drifttime<<"\t"<<Tot[k][TotNum[k]] <<endl;
	//Float_t layer_gap[8] = { 3.0, 9.0, 15.0, ... };
        //h_Ft_Geo->Fill((floor(loc.straw)/2)-1,layer_gap[loc.layer]);
			//////////////////////////////////////////////////////////////////////////////////
                        if (loc.layer==1)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,3.0);
				}
				else 
					h_Ft_Geo->Fill(loc.straw/2,1.0);
			}
			else if (loc.layer==2)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,9.0);
				}
				else
					h_Ft_Geo->Fill(loc.straw/2,7.0);
			}
			else if (loc.layer==3)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,15.0);
				}
				else
					h_Ft_Geo->Fill(loc.straw/2,13.0);
			}
			else if (loc.layer==4)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,21.0);
				}
				else
					h_Ft_Geo->Fill(loc.straw/2,19.0);
			}
			else if (loc.layer==5)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,27.0);
				}
				else
					h_Ft_Geo->Fill(loc.straw/2,25.0);
			}
			else if (loc.layer==6)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,33.0);
				}
				else
					h_Ft_Geo->Fill(loc.straw/2,31.0);
			}
			else if (loc.layer==7)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,39.0);
				}
				else
					h_Ft_Geo->Fill(loc.straw/2,37.0);
			}
			else if (loc.layer==8)
			{ 			
				if (loc.straw%2 ==0)
				{
					h_Ft_Geo->Fill((floor(loc.straw)/2)-1,45.0);
				}
				else
					h_Ft_Geo->Fill(loc.straw/2,43.0);
			}
///////////////////////////////////////////////////////////////////////////////////////////
			
				if (loc.straw%2 ==0)
				{
					h_FT_Geo->Fill( (loc.straw ) , loc.layer+0.5);
				}
				else
					h_FT_Geo->Fill( loc.straw, loc.layer);
                                
				if (loc.straw%2 ==0)
				{
					h_FT_Geo2->Fill( (loc.straw ) , loc.layer+0.5);
				}
				else
					h_FT_Geo2->Fill( loc.straw +2, loc.layer);

                                h_geo->Fill(loc.straw, loc.layer);
		
///////////////////////////////////////////////////////////////////////////////////////////
                                if (drifttime <=500 && drifttime >0)
                                {
                                    dtHitsMult[loc.layer][loc.straw]++;
                                }
                                
                        TotNum[k]++;
                        }
                        }
                    }
                }

                h_channelVsHits->Fill(k,TotNum[k]);

		

            }



        }
    //}

		int mult = 0;
		for (int a=1; a<=LAYERS; a++){
			if (hitMultOnLayer[a]> 0 )
                        {mult++;}
			//printf("dt multiplicity%i \n",dtHitsMult[a]);
                        for (int b=0; b<CHANNELS; b++)
                        {
                            if (dtHitsMult[a][b] >0)
                            hL_DT_Mult[a - 1]->Fill(b,dtHitsMult[a][b]); //Fill(b,dtHitsMult[a][b]);
                        }
		}

		// if (mult > 0) {
			h_layerMultiplicity->Fill(mult);
		// }

// post analysis




//for (leadTimes[local_ch] = 0; leadTimes[local_ch] < 200; leadTimes[local_ch]++)
//{
//for (leadTimes[leadTimesNum[local_ch]] = 0; leadTimes[leadTimesNum[local_ch]] < 10; leadTimes[leadTimesNum[local_ch]]++)
//{
//            Tot[local_ch][leadTimesNum[local_ch]]  << "  ";
//        }
//        cout << endl;
//    }
// post analysis

            return true;
        }

                                                       
    virtual void UserPostLoop() {

    }
};


void second()
{
    gSystem->Load("libGeom");
    LPetProcessor* proc = new LPetProcessor();
}

detLoc Detector::GetDetectorLocFromTDCChannel(int channel) {
    return detMap[channel];
}


Detector::Detector() {}

/*

////////////////////// Layer 1 ////////////////////////////

// TDC 6400

//fee1
    detMap[1] = detLoc(1, 1, 1);
    detMap[2] = detLoc(1, 1, 2);
    detMap[3] = detLoc(1, 1, 3);
    detMap[4] = detLoc(1, 1, 4);

    detMap[5] = detLoc(1, 1, 5);
    detMap[6] = detLoc(1, 1, 6);
    detMap[7] = detLoc(1, 1, 7);
    detMap[8] = detLoc(1, 1, 8);

    detMap[9] = detLoc(1, 1, 9);
    detMap[10] = detLoc(1, 1, 10);
    detMap[11] = detLoc(1, 1, 11);
    detMap[12] = detLoc(1, 1, 12);


    detMap[13] = detLoc(1, 1, 13);
    detMap[14] = detLoc(1, 1, 14);
    detMap[15] = detLoc(1, 1, 15);
    detMap[16] = detLoc(1, 1, 16);



//fee2
    detMap[17] = detLoc(1, 1, 17);
    detMap[18] = detLoc(1, 1, 18);
    detMap[19] = detLoc(1, 1, 19);
    detMap[20] = detLoc(1, 1, 20);

    detMap[21] = detLoc(1, 1, 21);
    detMap[22] = detLoc(1, 1, 22);
    detMap[23] = detLoc(1, 1, 23);
    detMap[24] = detLoc(1, 1, 24);

    detMap[25] = detLoc(1, 1, 25);
    detMap[26] = detLoc(1, 1, 26);
    detMap[27] = detLoc(1, 1, 27);
    detMap[28] = detLoc(1, 1, 28);

    detMap[29] = detLoc(1, 1, 29);
    detMap[30] = detLoc(1, 1, 30);
    detMap[31] = detLoc(1, 1, 31);
    detMap[32] = detLoc(1, 1, 32);

////////////////////// Layer 2 ////////////////////////////

//fee3
    detMap[33] = detLoc(1, 2, 1);
    detMap[34] = detLoc(1, 2, 2);
    detMap[35] = detLoc(1, 2, 3);
    detMap[36] = detLoc(1, 2, 4);

    detMap[37] = detLoc(1, 2, 5);
    detMap[38] = detLoc(1, 2, 6);
    detMap[39] = detLoc(1, 2, 7);
    detMap[40] = detLoc(1, 2, 8);

    detMap[41] = detLoc(1, 2, 9);
    detMap[42] = detLoc(1, 2, 10);
    detMap[43] = detLoc(1, 2, 11);
    detMap[44] = detLoc(1, 2, 12);

    detMap[45] = detLoc(1, 2, 13);
    detMap[46] = detLoc(1, 2, 14);
    detMap[47] = detLoc(1, 2, 15);
    detMap[48] = detLoc(1, 2, 16);



// TDC 6410

//fee4
    detMap[50] = detLoc(1, 2, 17);
    detMap[51] = detLoc(1, 2, 18);
    detMap[52] = detLoc(1, 2, 19);
    detMap[53] = detLoc(1, 2, 20);

    detMap[54] = detLoc(1, 2, 21);
    detMap[55] = detLoc(1, 2, 22);
    detMap[56] = detLoc(1, 2, 23);
    detMap[57] = detLoc(1, 2, 24);

    detMap[58] = detLoc(1, 2, 25);
    detMap[59] = detLoc(1, 2, 26);
    detMap[60] = detLoc(1, 2, 27);
    detMap[61] = detLoc(1, 2, 28);

    detMap[62] = detLoc(1, 2, 29);
    detMap[63] = detLoc(1, 2, 30);
    detMap[64] = detLoc(1, 2, 31);
    detMap[65] = detLoc(1, 2, 32);


////////////////////// Layer 3 ////////////////////////////

//fee5
    detMap[66] = detLoc(1, 3, 1);
    detMap[67] = detLoc(1, 3, 2);
    detMap[68] = detLoc(1, 3, 3);
    detMap[69] = detLoc(1, 3, 4);

    detMap[70] = detLoc(1, 3, 5);
    detMap[71] = detLoc(1, 3, 6);
    detMap[72] = detLoc(1, 3, 7);
    detMap[73] = detLoc(1, 3, 8);

    detMap[74] = detLoc(1, 3, 9);
    detMap[75] = detLoc(1, 3, 10);
    detMap[76] = detLoc(1, 3, 11);
    detMap[77] = detLoc(1, 3, 12);

    detMap[78] = detLoc(1, 3, 13);
    detMap[79] = detLoc(1, 3, 14);
    detMap[80] = detLoc(1, 3, 15);
    detMap[81] = detLoc(1, 3, 16);



//fee6
    detMap[82] = detLoc(1, 3, 17);
    detMap[83] = detLoc(1, 3, 18);
    detMap[84] = detLoc(1, 3, 19);
    detMap[85] = detLoc(1, 3, 20);

    detMap[86] = detLoc(1, 3, 21);
    detMap[87] = detLoc(1, 3, 22);
    detMap[88] = detLoc(1, 3, 23);
    detMap[89] = detLoc(1, 3, 24);

    detMap[90] = detLoc(1, 3, 25);
    detMap[91] = detLoc(1, 3, 26);
    detMap[92] = detLoc(1, 3, 27);
    detMap[93] = detLoc(1, 3, 28);

    detMap[94] = detLoc(1, 3, 29);
    detMap[95] = detLoc(1, 3, 30);
    detMap[96] = detLoc(1, 3, 31);
    detMap[97] = detLoc(1, 3, 32);

////////////////////// Layer 4 ////////////////////////////

// TDC 6411

//fee7
    detMap[99] = detLoc(1, 4, 1);
    detMap[100] = detLoc(1, 4, 2);
    detMap[101] = detLoc(1, 4, 3);
    detMap[102] = detLoc(1, 4, 4);

    detMap[103] = detLoc(1, 4, 5);
    detMap[104] = detLoc(1, 4, 6);
    detMap[105] = detLoc(1, 4, 7);
    detMap[106] = detLoc(1, 4, 8);

    detMap[107] = detLoc(1, 4, 9);
    detMap[108] = detLoc(1, 4, 10);
    detMap[109] = detLoc(1, 4, 11);
    detMap[110] = detLoc(1, 4, 12);

    detMap[111] = detLoc(1, 4, 13);
    detMap[112] = detLoc(1, 4, 14);
    detMap[113] = detLoc(1, 4, 15);
    detMap[114] = detLoc(1, 4, 16);



//fee8

    detMap[115] = detLoc(1, 4, 17);
    detMap[116] = detLoc(1, 4, 18);
    detMap[117] = detLoc(1, 4, 19);
    detMap[118] = detLoc(1, 4, 20);

    detMap[119] = detLoc(1, 4, 21);
    detMap[120] = detLoc(1, 4, 22);
    detMap[121] = detLoc(1, 4, 23);
    detMap[122] = detLoc(1, 4, 24);

    detMap[123] = detLoc(1, 4, 25);
    detMap[124] = detLoc(1, 4, 26);
    detMap[125] = detLoc(1, 4, 27);
    detMap[126] = detLoc(1, 4, 28);

    detMap[127] = detLoc(1, 4, 29);
    detMap[128] = detLoc(1, 4, 30);
    detMap[129] = detLoc(1, 4, 31);
    detMap[130] = detLoc(1, 4, 32);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////// Layer 5 ////////////////////////////

// TDC 6420

//fee9
    detMap[148] = detLoc(1, 5, 1);
    detMap[149] = detLoc(1, 5, 2);
    detMap[150] = detLoc(1, 5, 3);
    detMap[151] = detLoc(1, 5, 4);

    detMap[152] = detLoc(1, 5, 5);
    detMap[153] = detLoc(1, 5, 6);
    detMap[154] = detLoc(1, 5, 7);
    detMap[155] = detLoc(1, 5, 8);

    detMap[156] = detLoc(1, 5, 9);
    detMap[157] = detLoc(1, 5, 10);
    detMap[158] = detLoc(1, 5, 11);
    detMap[159] = detLoc(1, 5, 12);

    detMap[160] = detLoc(1, 5, 13);
    detMap[161] = detLoc(1, 5, 14);
    detMap[162] = detLoc(1, 5, 15);
    detMap[163] = detLoc(1, 5, 16);



//fee10
    detMap[164] = detLoc(1, 5, 17);
    detMap[165] = detLoc(1, 5, 18);
    detMap[166] = detLoc(1, 5, 19);
    detMap[167] = detLoc(1, 5, 20);

    detMap[168] = detLoc(1, 5, 21);
    detMap[169] = detLoc(1, 5, 22);
    detMap[170] = detLoc(1, 5, 23);
    detMap[171] = detLoc(1, 5, 24);

    detMap[172] = detLoc(1, 5, 25);
    detMap[173] = detLoc(1, 5, 26);
    detMap[174] = detLoc(1, 5, 27);
    detMap[175] = detLoc(1, 5, 28);

    detMap[176] = detLoc(1, 5, 29);
    detMap[177] = detLoc(1, 5, 30);
    detMap[178] = detLoc(1, 5, 31);
    detMap[179] = detLoc(1, 5, 32);


////////////////////// Layer 6 ////////////////////////////

//fee11
    detMap[180] = detLoc(1, 6, 1);
    detMap[181] = detLoc(1, 6, 2);
    detMap[182] = detLoc(1, 6, 3);
    detMap[183] = detLoc(1, 6, 4);

    detMap[184] = detLoc(1, 6, 5);
    detMap[185] = detLoc(1, 6, 6);
    detMap[186] = detLoc(1, 6, 7);
    detMap[187] = detLoc(1, 6, 8);

    detMap[188] = detLoc(1, 6, 9);
    detMap[189] = detLoc(1, 6, 10);
    detMap[190] = detLoc(1, 6, 11);
    detMap[191] = detLoc(1, 6, 12);

    detMap[192] = detLoc(1, 6, 13);
    detMap[193] = detLoc(1, 6, 14);
    detMap[194] = detLoc(1, 6, 15);
    detMap[195] = detLoc(1, 6, 16);



// TDC E6430

//fee12


    detMap[197] = detLoc(1, 6, 17);
    detMap[198] = detLoc(1, 6, 18);
    detMap[199] = detLoc(1, 6, 19);
    detMap[200] = detLoc(1, 6, 20);

    detMap[201] = detLoc(1, 6, 21);
    detMap[202] = detLoc(1, 6, 22);
    detMap[203] = detLoc(1, 6, 23);
    detMap[204] = detLoc(1, 6, 24);

    detMap[205] = detLoc(1, 6, 25);
    detMap[206] = detLoc(1, 6, 26);
    detMap[207] = detLoc(1, 6, 27);
    detMap[208] = detLoc(1, 6, 28);

    detMap[209] = detLoc(1, 6, 29);
    detMap[210] = detLoc(1, 6, 30);
    detMap[211] = detLoc(1, 6, 31);
    detMap[212] = detLoc(1, 6, 32);


////////////////////// Layer 7 ////////////////////////////

//fee13

    detMap[213] = detLoc(1, 7, 1);
    detMap[214] = detLoc(1, 7, 2);
    detMap[215] = detLoc(1, 7, 3);
    detMap[216] = detLoc(1, 7, 4);

    detMap[217] = detLoc(1, 7, 5);
    detMap[218] = detLoc(1, 7, 6);
    detMap[219] = detLoc(1, 7, 7);
    detMap[220] = detLoc(1, 7, 8);

    detMap[221] = detLoc(1, 7, 9);
    detMap[222] = detLoc(1, 7, 10);
    detMap[223] = detLoc(1, 7, 11);
    detMap[224] = detLoc(1, 7, 12);

    detMap[225] = detLoc(1, 7, 13);
    detMap[226] = detLoc(1, 7, 14);
    detMap[227] = detLoc(1, 7, 15);
    detMap[228] = detLoc(1, 7, 16);


//fee14


    detMap[229] = detLoc(1, 7, 17);
    detMap[230] = detLoc(1, 7, 18);
    detMap[231] = detLoc(1, 7, 19);
    detMap[232] = detLoc(1, 7, 20);

    detMap[233] = detLoc(1, 7, 21);
    detMap[234] = detLoc(1, 7, 22);
    detMap[235] = detLoc(1, 7, 23);
    detMap[236] = detLoc(1, 7, 24);

    detMap[237] = detLoc(1, 7, 25);
    detMap[238] = detLoc(1, 7, 26);
    detMap[239] = detLoc(1, 7, 27);
    detMap[240] = detLoc(1, 7, 28);

    detMap[241] = detLoc(1, 7, 29);
    detMap[242] = detLoc(1, 7, 30);
    detMap[243] = detLoc(1, 7, 31);
    detMap[244] = detLoc(1, 7, 32);


////////////////////// Layer 8 ////////////////////////////

// TDC 6431

//fee15

    detMap[246] = detLoc(1, 8, 1);
    detMap[247] = detLoc(1, 8, 2);
    detMap[248] = detLoc(1, 8, 3);
    detMap[249] = detLoc(1, 8, 4);

    detMap[250] = detLoc(1, 8, 5);
    detMap[251] = detLoc(1, 8, 6);
    detMap[252] = detLoc(1, 8, 7);
    detMap[253] = detLoc(1, 8, 8);

    detMap[254] = detLoc(1, 8, 9);
    detMap[255] = detLoc(1, 8, 10);
    detMap[256] = detLoc(1, 8, 11);
    detMap[257] = detLoc(1, 8, 12);

    detMap[258] = detLoc(1, 8, 13);
    detMap[259] = detLoc(1, 8, 14);
    detMap[260] = detLoc(1, 8, 15);
    detMap[261] = detLoc(1, 8, 16);


//fee16


    detMap[262] = detLoc(1, 8, 17);
    detMap[263] = detLoc(1, 8, 18);
    detMap[264] = detLoc(1, 8, 19);
    detMap[265] = detLoc(1, 8, 20);

    detMap[266] = detLoc(1, 8, 21);
    detMap[267] = detLoc(1, 8, 22);
    detMap[268] = detLoc(1, 8, 23);
    detMap[269] = detLoc(1, 8, 24);

    detMap[270] = detLoc(1, 8, 25);
    detMap[271] = detLoc(1, 8, 26);
    detMap[272] = detLoc(1, 8, 27);
    detMap[273] = detLoc(1, 8, 28);

    detMap[274] = detLoc(1, 8, 29);
    detMap[275] = detLoc(1, 8, 30);
    detMap[276] = detLoc(1, 8, 31);
    detMap[277] = detLoc(1, 8, 32);

}*/
