#include "SttEvent.h"
#include "SttHit.h"
#include "SttRawHit.h"
#include "SttTrackEvent.h"
#include "SttTrackHit.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include <TCanvas.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TLinearFitter.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TObject.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>

#include "panda_stt_cal.h"
#include "panda_stt_track.h"
#include "panda_subsystem.h"
#include "panda_subsystem_sb.h"
#include "panda_subsystem_stt.h"

using namespace std;

Bool_t histforTracks(void)
{
    cout << "Open" << endl;

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;
    PandaSttCal* STT_CAL = 0;
    PandaSttTrack* STT_TRACK = new PandaSttTrack();

    Double_t w = 2000; // ToT Canvas width
    Double_t h = 1000; // ToT Canvas hight

    TCanvas* cToT = new TCanvas("ToT_Canvas", "ToT_Canvas", w, h);
    cToT->Divide(2, 2);

    TCanvas* cDMult = new TCanvas("DMult_Canvas", "DMult_Canvas", w, h);
    cDMult->Divide(2, 2);

    TFile* file = TFile::Open("N1800clus.root", "READ");

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

    int counterofdt = 0;
    int driftTimeCounter2 = 0;

    tree->SetBranchAddress("STT_TRACKS", &STT_TRACK);
    TFile* ftree = new TFile("c.root", "RECREATE");


    TH1F* h_x;
    TH1F* h_STT_Hit_Diff =
        new TH1F("h_STT_Hit_Diff", "h_STT_Hit_Diff", 800, -100, 700);
    TH1I* L_STT_Hit_Diff[8];
    for (int i = 0; i < 8; i++)
    {
        L_STT_Hit_Diff[i] =
            new TH1I(Form("Layer%dSTT_Hit_Diff", i + 1),
                     Form("Layer%dSTT_Hit_Diff", i + 1), 800, -100, 700);
    }

	TH1I* L_STT_Hit_Mult[8];
    for (int i = 0; i < 8; i++)
    {
        L_STT_Hit_Mult[i] =
            new TH1I(Form("Layer%dSTT_Hit_Mult", i + 1),
                     Form("Layer%dSTT_Hit_Mult", i + 1), 10, -0, 10);
    }

    h_x = new TH1F("h_x", "h_x", 16, 0, 16);

    Int_t iev = (Int_t)tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl << endl;

    cout << STT_TRACK->stt_track_can.total_track_NTDCHits << endl;

    for (Int_t i = 0; i < iev; i++)
    {
        tree->GetEntry(i);
        if (i % 100 == 0) {
            cout << i << endl;
        }

        for (int n = 0; n < STT_TRACK->stt_track_can.total_track_NTDCHits; n++)
        {
            // char buff[200];
            // sprintf(buff, "can_%d", i);
            // TCanvas* c1 = new TCanvas(buff);
            // sprintf(buff, "h_%d", i);
            // int x_min = -10;
            // int x_max = 90;
            // int x_bins = (x_max - x_min) * 5;
            // int z_min = -10;
            // int z_max = 80;
            // int z_bins = (z_max - z_min) * 5;

            // TH2I* h = new TH2I(buff, "h;x,y [cm];z [cm]", x_bins, x_min,
            // x_max, z_bins, z_min, z_max);
            // c1->cd();
            // h->Draw();

            // h_track[i] = new TH2F (Form("Track No %d - Straw Vs
            // Layer",i),Form ("Track %d Straw VsLayer; Straw;Layer",i
            // ),18.0,-1.0,17.0,20,0.0,10.0);

            std::vector<SttHit*> vec_track_can;
            // SttHit* cal_hit  =
            // (SttHit*)STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt(n);
            SttTrackHit* track_hit =
                (SttTrackHit*)
                    STT_TRACK->stt_track_can.tdc_track_hits->ConstructedAt(n);
            vec_track_can = track_hit->vec_Track;

            std::map<int, int> layer_mult;
            for (int t = 0; t < vec_track_can.size(); t++)
            {
            	++layer_mult[vec_track_can[t]->layer - 1];
                h_STT_Hit_Diff->Fill(vec_track_can[t]->drifttime);
                L_STT_Hit_Diff[vec_track_can[t]->layer - 1]->Fill(
                    vec_track_can[t]->drifttime);
                // vec_o_test1.push_back(vec_track_can[t]->drifttime);
                // vec_o_test.push_back(vec_track_can[t]->drifttime);
                // vec_x.push_back(vec_track_can[t]->x);
                // vec_y.push_back(vec_track_can[t]->y);
                // vec_z.push_back(vec_track_can[t]->z);
                // vec_layer.push_back(vec_track_can[t]->layer);
                //          float uuu = 0.0;
                //          bool is_y = (vec_track_can[t]->layer % 2 == 0);
                //          if (is_y)
                //              uuu = vec_track_can[t]->y;
                //          else
                //              uuu = vec_track_can[t]->x;
                //          cout<<uuu<<"\t"<<vec_track_can[t]->z<<endl;
                //          TEllipse* straws = new TEllipse(uuu,
                //          vec_track_can[t]->z, 0.505, 0.505);
                //        //  box = new TBox(25,60,43,78);
                // //	box->SetFillColor(37);

                //          TMultiGraph* mg = new TMultiGraph();
                // mg->Draw();
                // TMultiGraph* mg1 = new TMultiGraph();

                // cout<<"see
                // :"<<a->drifttime<<"\t"<<a->layer<<"\t"<<a->module<<"\t"<<a->fee<<"\t"<<a->fee_channel<<"\t"<<a->x<<"\t"<<a->y<<"\t"<<a->z<<endl;
                /////////

                // cout << "\t" << vec_hits.size() << "\t" <<
                // vec_hits[n]->drifttime << "\t" << vec_hits[n]->layer << "\t"
                // << vec_hits[n]->module << "\t" << vec_hits[n]->fee << "\t" <<
                // vec_hits[n]->fee_channel << "\t" << vec_hits[n]->x << "\t" <<
                // vec_hits[n]->y << "\t" << vec_hits[n]->z << "\t" << endl;
                // straws->Draw("same");

                /*if (vec_track_can[t]->straw%2 ==0)
                {
                    h_track[i]->Fill((floor(vec_track_can[t]->straw)/2)-1,vec_track_can[t]->layer+0.5);
                }
                else
                    h_track[i]->Fill(vec_track_can[t]->straw/2,vec_track_can[t]->layer);

                if (vec_track_can[t]->layer == 1 || vec_track_can[t]->layer == 4
                || vec_track_can[t]->layer == 5 || vec_track_can[t]->layer == 8)
                    h_x->Fill(vec_track_can[t]->straw/2);
                else
                    printf("assasa\n");*/

                counterofdt++;
            }

            for (int l = 0; l < 8; ++l)
            L_STT_Hit_Mult[l]->Fill(layer_mult[l]);
            // DR_Tree->Fill();
            // vec_o_test.clear();
            // vec_x.clear();
            // vec_y.clear();
            // vec_z.clear();
        }
    }
    h_STT_Hit_Diff->Write();
    for (int h = 0; h < 8; h++)
    {
        L_STT_Hit_Diff[h]->Write();
        L_STT_Hit_Mult[h]->Write();
    }
    int layer[5] = {0, 0, 3, 4, 7};

    for (int a = 1; a < 5; a++)
    {
        cToT->cd(a);
        L_STT_Hit_Diff[layer[a]]->Draw();
        cDMult->cd(a);
        L_STT_Hit_Mult[layer[a]]->Draw();
    }

    // h_x->Write();

    cToT->Write();
    cDMult->Write();
    return kTRUE;

    ////////////////////
}
