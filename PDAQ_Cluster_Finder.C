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
#include "SttEvent.h"
#include "SttTrackEvent.h"

                          
#include "panda_subsystem.h"
#include "panda_subsystem_stt.h"
#include "panda_subsystem_sb.h"
#include "panda_stt_cal.h"
#include "panda_stt_track.h"

#pragma link C++ class SttRawHit+;
#pragma link C++ class SttHit+;
#pragma link C++ class SttEvent+;
#pragma link C++ class Stt_Track_Event+;

#pragma link C++ class PandaSubsystem+;
#pragma link C++ class PandaSubsystemSTT+;
#pragma link C++ class PandaSubsystemSB+;
#pragma link C++ class PandaSttCal+;
#pragma link C++ class PandaSttTrack+;



using namespace std;


bool f_sttHitCompareLeadTime(SttHit* a, SttHit* b)
{
    return (a->leadTime < b->leadTime);
}


bool f_sttHitCompareCell(SttHit* a, SttHit* b)
{
    return (a->cell < b->cell);
}


std::vector<std::vector<SttHit*>*>* clusterfinder(std::vector<SttHit*>* vec_flayer)
{
    std::vector<std::vector<SttHit*>*>* vec_Cl = new std::vector<std::vector<SttHit*>*>();
    std::vector<SttHit*>* clusterPointer;

    std::sort(vec_flayer->begin(), vec_flayer->end(), f_sttHitCompareCell);   
            
     if (vec_flayer->size() ==1)
    {
        clusterPointer = new std::vector<SttHit*>();
        clusterPointer->push_back(vec_flayer->at(0));
        vec_Cl->push_back(clusterPointer);
        return vec_Cl;
    }
    else 
    {
        for (Int_t aa=0; aa < vec_flayer->size(); aa++ )
        {
         if (aa < vec_flayer->size()-1)
            {
                 if(fabs(vec_flayer->at(aa)->cell - vec_flayer->at(aa+1)->cell) < 2)
                {
                    clusterPointer = new vector<SttHit*>();
                    clusterPointer->push_back(vec_flayer->at(aa));
                    clusterPointer->push_back(vec_flayer->at(aa+1));
                    vec_Cl->push_back(clusterPointer);                           
                } 
                                else if(aa==0)
                                {
                                    clusterPointer = new vector<SttHit*>();
                                    clusterPointer->push_back(vec_flayer->at(aa));
                                    vec_Cl->push_back(clusterPointer); 
                                }
                                else if(aa > 0)
                                {
                                    if(fabs(vec_flayer->at(aa-1)->cell - vec_flayer->at(aa)->cell) == 1) 
                                    {
                                        continue;
                                    } 
                                    
                                    else
                                   {
                                        clusterPointer = new vector<SttHit*>();
                                        clusterPointer->push_back(vec_flayer->at(aa));
                                        vec_Cl->push_back(clusterPointer); 

                                    }
                                }
            }

            else if (aa == vec_flayer->size()-1)
            {
                 if ((fabs(vec_flayer->at(aa-1)->cell - vec_flayer->at(aa)->cell) > 1))
                        {
                            clusterPointer = new vector<SttHit*>();
                            clusterPointer->push_back(vec_flayer->at(aa));
                            vec_Cl->push_back(clusterPointer);
                        } 
             }
        }
        return vec_Cl;
     }
}

std::vector<SttHit*>  GetPairs(std::vector<SttHit*> vec_get_pairs)
{
    std::vector<SttHit*> vec_fpair_clu;
    std::sort(vec_get_pairs.begin(), vec_get_pairs.end(), f_sttHitCompareCell);   


        for (Int_t sa=0; sa<vec_get_pairs.size(); sa++)
        {

            if ( sa < vec_get_pairs.size()-1)
            {
                if(fabs(vec_get_pairs[sa]->cell - vec_get_pairs[sa+1]->cell) == 1)
                {vec_fpair_clu.push_back(vec_get_pairs[sa]);
                cout<<"PAIR : "<<vec_get_pairs[sa]->layer<<"\t"<<vec_get_pairs[sa]->cell<<endl;}

                else if (sa >0)
                {
                    if(fabs(vec_get_pairs[sa-1]->cell - vec_get_pairs[sa]->cell) == 1)
                    vec_fpair_clu.push_back(vec_get_pairs[sa]);
                cout<<"PAIR : "<<vec_get_pairs[sa]->layer<<"\t"<<vec_get_pairs[sa]->cell<<endl;

                }

            }

            else if (sa == vec_get_pairs.size()-1)
            {
                if(fabs(vec_get_pairs[sa-1]->cell - vec_get_pairs[sa]->cell) == 1)
                vec_fpair_clu.push_back(vec_get_pairs[sa]);
                cout<<"PAIR : "<<vec_get_pairs[sa]->layer<<"\t"<<vec_get_pairs[sa]->cell<<endl;

            }

         }
            return vec_fpair_clu;

 }

Bool_t PDAQ_Cluster_Finder(void)
{

// PandaSttCal* CAL = new PandaSttCal();

// SttEvent* stt_event = &(CAL->stt_cal);

PandaSubsystemSTT* STT_RAW = 0;
PandaSubsystemSTT* CAL = 0;

//PandaSttCal* STT = 0;


PandaSttCal* STT_CAL = new PandaSttCal();



PandaSttTrack* STT_TRACKS = new PandaSttTrack();
Stt_Track_Event* stt_event = &(STT_TRACKS->stt_track_can);


    TFile file("PDAQ_Stt.root", "READ");
    TTree* tree = (TTree*)file.Get("PDAQ_EMC_STT_cluster_analysis");
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

    tree->SetBranchAddress("STT_CAL", &STT_CAL);
    tree->SetBranchAddress("STT_RAW", &STT_RAW);

    TFile* Ttree = new TFile("PDAQ_Stt.root", "UPDATE");
     TTree* PDAQ_EMC_STT_cluster_analysis =  new TTree("PDAQ_EMC_STT_cluster_analysis", "PDAQ_EMC_STT_cluster_analysis");

    //PDAQ_EMC_STT_cluster_analysis->Branch("STT_RAW", "PandaSubsystemSTT", &STT_RAW, 64000, 99);
    //PDAQ_EMC_STT_cluster_analysis->Branch("STT_CAL", "PandaSttCal", &STT_CAL, 64000, 99);
    PDAQ_EMC_STT_cluster_analysis->Branch("STT_TRACKS", "PandaSttTrack", &STT_TRACKS, 64000, 99);

    TH1F* h_STT_EMC_td1 = new TH1F("h_STT_EMC_td1", "h_STT_EMC_td1", 1000, 0, 1000);
    TH1F* h_STT_EMC_td2 = new TH1F("h_STT_EMC_td2", "h_STT_EMC_td2", 1000, 0, 1000);
    TH1F* h_STT_EMC_td3 = new TH1F("h_STT_EMC_td3", "h_STT_EMC_td3", 1000, 0, 1000);
    TH1F* h_STT_EMC_td4 = new TH1F("h_STT_EMC_td4", "h_STT_EMC_td4", 1000, 0, 1000);
    TH2F* h_XvsZ = new TH2F("h_XvsZ", "h_XvsZ", 100, -20, 80, 100, -20, 80);
    TH2F* h_YvsZ = new TH2F("h_YvsZ", "h_YvsZ", 100, 0, 100, 100, 0, 100);
    TH1F* h_STT_EMC_timeDiff = new TH1F("h_STT_EMC_timeDiff", "h_STT_EMC_timeDiff", 600, 100, 700);
    
    TH2F* h_Straw_DriftTime = new TH2F("h_Straw_DriftTime", "h_Straw_DriftTime", 296, 0, 296, 700, 0, 700);
    TH2F* h_Fee_DriftTime = new TH2F("h_Fee_DriftTime", "h_Fee_DriftTime", 16, 0, 16, 700, 0, 700);

    TH1F* h_straw_mean_straw = new TH1F("h_straw_mean_straw", "h_straw_mean_straw", 800, -100, 700);

    TH1F* h_Front = new TH1F("h_Front", "h_Front", 600, 0, 600);    
    TH1F* h_FrontNO = new TH1F("h_FrontNO", "h_FrontNO", 20, 0, 20);    

    TH1F* h_Fee[18];

    Double_t repeat =0;
    Double_t All_repeat =0;

    for (int f = 0; f < 18; f++) 
    {
        h_Fee[f] = new TH1F(TString::Format("h_Fee%f", f), "h_Fee", 600, 100, 700);
    }

    Int_t iev = (Int_t)tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl
         << endl;

//loop over all the vectors in the tree.
// sleep(5);

    for (Int_t i = 0; i < iev; i++) 
    {
        tree->GetEntry(i);
        cout << endl
             << endl;

        std::vector<SttHit*> vec_stthits;


        if (i%10000 ==0)
        {
            cout << "entry no. " << i << endl;
        }

        SttHit* hitOnLayer[4][500];
        memset(hitOnLayer, 0, 4*500*sizeof(SttHit*));
        int hitMultOnLayer[4];
        
        for (int h = 0; h < 4; h++)
        {
            hitMultOnLayer[h] = 0;
        }


    //Loop over the vector elements//////////////////////////////////////////////////////
        for (int n = 0; n < STT_CAL->stt_cal.total_cal_NTDCHits; n++)
        {

            SttHit* cal_hit  = (SttHit*)STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt(n); // retrieve particular hit


            int tdc_num = cal_hit->channel / 49;
            // hit on reference channel
            if (cal_hit->isRef == true)
            {
                 //cout<<"Reference Hit  -> "<<cal_hit->isRef<<"\t"<<cal_hit->layer<<"\t"<<cal_hit->module<<"\t"<<cal_hit->fee_channel<<"\t"<<cal_hit->cell<<endl;
            }
            else
            {

                if (cal_hit->layer == 0)
                {
                    printf("ERROR: cal_hit->layer %d %d %d\n", cal_hit->layer, cal_hit->channel, cal_hit->fee);
                    continue;
                }

                hitOnLayer[cal_hit->layer - 1][hitMultOnLayer[cal_hit->layer - 1]] = cal_hit;

                hitMultOnLayer[cal_hit->layer - 1]++;
            }

        }

        bool good_layers = true;
        for (int c = 0; c < 4; c++)
        {
            if (hitMultOnLayer[c] != 1)
            {
                 good_layers = false;
                 break;
            }
        }
        int layerCounter = 0;

        for (int m = 0; m < 4; m++)
        {
            if (hitMultOnLayer[m] > 0)
                layerCounter++;
        }

        if (layerCounter == 4)
        {
            std::vector<SttHit*> vec_leadTime;
            int filtercnt = 0;
            int fil_max = 250;

            for (int l = 0; l < 4; l++)
            {
                for (int h = 0; h < hitMultOnLayer[l]; h++)

                {
                    vec_leadTime.push_back(hitOnLayer[l][h]);

                }
            }

            for(Int_t je=0; je<vec_leadTime.size()-1;je++)
            {


                cout << "All entries  :  "<<vec_leadTime[je]->isRef<<"\t"<<" Lead Time : "<< vec_leadTime[je]->leadTime<<"\t" << " layer : "<< vec_leadTime[je]->layer<<"\t"<< "module : "<< vec_leadTime[je]->module<<"\t"<< " fee : "<< vec_leadTime[je]->fee<<"\t" << "FEE_Ch : "<< vec_leadTime[je]->fee_channel<<"\t"<<"Cell : "<<vec_leadTime[je]->cell<<endl;
                if ((vec_leadTime[je+1]->leadTime) == (vec_leadTime[je]->leadTime) && (vec_leadTime[je+1]->channel == vec_leadTime[je]->channel))
                {
                    repeat++;
                    //printf("{ LT,l,m,c,cf=%d,%d,%2d,%2d}{ LT1,l1,m1,c1,cf1=%d,%d,%2d,%2d}", vec_leadTime[je]->layer,vec_leadTime[je]->module,vec_leadTime[je]->fee,vec_leadTime[je]->fee_channel,vec_leadTime[je+1]->layer,vec_leadTime[je+1]->module,vec_leadTime[je+1]->fee,vec_leadTime[je+1]->fee_channel);
                    //cout<<vec_leadTime[je]->isRef<<"\t"<<vec_leadTime[je]->leadTime<<"\t"<<vec_leadTime[je]->tot<<"\t"<<vec_leadTime[je+1]->leadTime<<" "<<vec_leadTime[je]->tot<<endl;    
                 }
                 All_repeat++;
             }

            std::sort(vec_leadTime.begin(), vec_leadTime.end(), f_sttHitCompareLeadTime);
            vec_stthits.clear();


            if (vec_leadTime.size() > 7)
            {

                for (int v = 0; v < vec_leadTime.size(); v++)
                { // iterate over collected and sorted hits
                    vec_stthits.clear();
                    filtercnt = 1;

                    SttHit* h = vec_leadTime.at(v);
                    vec_stthits.push_back(h);


                    for (int vv = 0; vv < vec_leadTime.size(); vv++)
                    { // check each vs each if they fit into the window
                        if (vv == v)    continue;
                        if ((fabs(vec_leadTime[v]->leadTime - vec_leadTime[vv]->leadTime) < fil_max) )

                        {
                            SttHit* hh = vec_leadTime.at(vv);
                            vec_stthits.push_back(hh);
                            filtercnt++;
                        }
                        else
                        { // in case the hit is outside the window break and start next iteration
                            // with the next hit
                            break;
                            vec_stthits.clear();

                        }
                    }
                }
            }
        
            //double meanTime = 0;
            std::vector<SttHit*> vec_layer1;
            std::vector<SttHit*> vec_layer2;
            std::vector<SttHit*> vec_layer3;
            std::vector<SttHit*> vec_layer4;
            vec_layer1.clear();
            vec_layer2.clear();
            vec_layer3.clear();
            vec_layer4.clear();

                for (int w = 0; w < vec_stthits.size(); w++)
                {
                    cout<<"["<<w<<"]"<<" L, M, C, F, Ch, T {"<<vec_stthits[w]->layer<<", "<<vec_stthits[w]->module<<", "<<vec_stthits[w]->cell<<", "<<vec_stthits[w]->fee<<", "<<vec_stthits[w]->fee_channel<<", "<<vec_stthits[w]->channel<<"}  "<<" X , Y , Z :{"<<vec_stthits[w]->x<<", "<<vec_stthits[w]->y<<", "<<vec_stthits[w]->z<<"}   "<<" LT DT  {"<<vec_stthits[w]->leadTime<<",  "<<vec_stthits[w]->drifttime<<"}"<<endl;

                        if (vec_stthits.size()<40)
                        {
                            if (vec_stthits[w]->layer == 1)
                            {
                                vec_layer1.push_back(vec_stthits[w]);
                            }
                            else if (vec_stthits[w]->layer == 2)
                            {
                                vec_layer2.push_back(vec_stthits[w]);
                            }
                            else if (vec_stthits[w]->layer == 3)
                            {
                                vec_layer3.push_back(vec_stthits[w]);
                            }
                            else if (vec_stthits[w]->layer == 4)
                            {
                                vec_layer4.push_back(vec_stthits[w]);
                            }

                        }
            }

            std::vector< vector<SttHit*>* >* vec_Cl_L1;
            std::vector< vector<SttHit*>* >* vec_Cl_L2;
            std::vector< vector<SttHit*>* >* vec_Cl_L3;
            std::vector< vector<SttHit*>* >* vec_Cl_L4;

            std::vector<SttHit*> vec_Clusters;

            std::vector<double> vec_Chi2x;
            std::vector<double> vec_Chi2y;

            std::vector<double> vec_P0;
            std::vector<double> vec_P1;

            std::vector<double> vec_PP0;
            std::vector<double> vec_PP1;
// FILTER TO GET ONLY HITS WITH A PAIR////////////////////////////////////////////////////////////

            std::vector<SttHit*> vec_pair_clu;
            vec_pair_clu.clear();

            std::vector<SttHit*> vec_pair_layer11;
            std::vector<SttHit*> vec_pair_layer22;
            std::vector<SttHit*> vec_pair_layer33;
            std::vector<SttHit*> vec_pair_layer44;

            vec_pair_layer11.clear();
            vec_pair_layer22.clear();
            vec_pair_layer33.clear();
            vec_pair_layer44.clear();

            if (vec_layer1.size()>1 && vec_layer2.size()>1 && vec_layer3.size()>1 && vec_layer4.size()>1 )
            {
                vec_pair_layer11 = GetPairs(vec_layer1);
                vec_pair_layer22 = GetPairs(vec_layer2);
                vec_pair_layer33 = GetPairs(vec_layer3);
                vec_pair_layer44 = GetPairs(vec_layer4);
                
            }
            else continue;

for (Int_t a =0; a< vec_pair_layer11.size(); a++)
{
 
cout<<"COUT : "<<vec_pair_layer11[a]->cell<<endl;
 }

 // END OF FILTER TO GET ONLY HITS WITH A PAIR////////////////////////////////////////////////////////////
           
            if ( vec_pair_layer11.size() > 0 && vec_pair_layer22.size() > 0 && vec_pair_layer33.size() > 0 && vec_pair_layer44.size() > 0)
            {                
                vec_Cl_L1 = clusterfinder(&vec_pair_layer11);
                vec_Cl_L2 = clusterfinder(&vec_pair_layer22);
                vec_Cl_L3 = clusterfinder(&vec_pair_layer33);
                vec_Cl_L4 = clusterfinder(&vec_pair_layer44);                
            }
            else continue;

std::vector< vector<SttHit*> > vec_All_X;
std::vector< vector<SttHit*> > vec_All_Y;
vec_All_X.clear();
vec_All_Y.clear();

for (Int_t xa=0; xa<vec_Cl_L1->size(); xa++)
{

    for (Int_t xb=0; xb<vec_Cl_L2->size(); xb++)
    {
        for (Int_t xc=0; xc< vec_Cl_L3->size(); xc++)
        {
            for (Int_t xd=0; xd< vec_Cl_L4->size(); xd++)
            {
                vec_Clusters.clear();

                for (Int_t xaa =0; xaa< vec_Cl_L1->at(xa)->size(); xaa++)
                {
                    //cout<< "Vector  : "<<vec_Cl_L1->at(xa)->at(xaa)->layer <<"\t"<< vec_Cl_L1->at(xa)->at(xaa)->cell<<endl;
                    vec_Clusters.push_back(vec_Cl_L1->at(xa)->at(xaa));
                }               

                for (Int_t xbb=0; xbb< vec_Cl_L2->at(xb)->size(); xbb++ )
                {
                    //cout<< "Vector  : "<<vec_Cl_L2->at(xb)->at(xbb)->layer <<"\t"<< vec_Cl_L2->at(xb)->at(xbb)->cell<<endl;
                    vec_Clusters.push_back(vec_Cl_L2->at(xb)->at(xbb));
                }

                for (Int_t xcc=0; xcc< vec_Cl_L3->at(xc)->size(); xcc++ )
                {
                    //cout<< "Vector  : "<<vec_Cl_L3->at(xc)->at(xcc)->layer <<"\t"<< vec_Cl_L3->at(xc)->at(xcc)->cell<<endl;
                    vec_Clusters.push_back(vec_Cl_L3->at(xc)->at(xcc));
                }

                for (Int_t xdd=0; xdd< vec_Cl_L4->at(xd)->size(); xdd++ )
                {
                    //cout<< "Vector  : "<<vec_Cl_L4->at(xd)->at(xdd)->layer <<"\t"<< vec_Cl_L4->at(xd)->at(xdd)->cell<<endl;
                    vec_Clusters.push_back(vec_Cl_L4->at(xd)->at(xdd));

                }

                Double_t chi_valueX =0;
                Double_t chi_valueY =0;
                Int_t loneCounterX =0;
                Int_t loneCounterY =0;

                std::vector<SttHit*> vec_ClustersX;
                std::vector<SttHit*> vec_ClustersY;

                vec_ClustersX.clear();
                vec_ClustersY.clear();

                for (Int_t ya=0; ya<vec_Clusters.size(); ya++)
                {
                    cout<< " Cluster Vector :"<<vec_Clusters[ya]->layer<<"\t"<<vec_Clusters[ya]->cell<<"\t"<<vec_Clusters[ya]->x<<"\t"<<vec_Clusters[ya]->y<<endl;
                    if (vec_Clusters[ya]->layer == 1 || vec_Clusters[ya]->layer ==3)
                    {
                        vec_ClustersX.push_back(vec_Clusters[ya]);
                        loneCounterX++;
                    }
                    if (vec_Clusters[ya]->layer == 2 || vec_Clusters[ya]->layer ==4)
                    {
                        vec_ClustersY.push_back(vec_Clusters[ya]);
                        loneCounterY++;
                    }
                }

                vec_All_X.push_back(vec_ClustersX);
                vec_All_Y.push_back(vec_ClustersY);

                Double_t clusterArrayX[vec_ClustersX.size()];
                Double_t clusterArrayZx[vec_ClustersX.size()];              
                Double_t clusterArrayY[vec_ClustersY.size()];               
                Double_t clusterArrayZy[vec_ClustersY.size()];  


                for (Int_t yb=0; yb<vec_ClustersX.size(); yb++)
                {
                    clusterArrayX[yb]=vec_ClustersX[yb]->x;
                    clusterArrayZx[yb]=vec_ClustersX[yb]->z;
                }
    
                    
                for (Int_t yc=0; yc<vec_ClustersY.size(); yc++)
                {
                    clusterArrayY[yc]=vec_ClustersY[yc]->y;
                    clusterArrayZy[yc]=vec_ClustersY[yc]->z;
                }
            
                TF1* f1 = new TF1("f1", "pol1");
                TF1* f2 = new TF1("f2", "pol1");


                TGraph* chiX = new TGraph(vec_ClustersX.size(), clusterArrayX, clusterArrayZx);
                chiX->Fit(f1);
                chi_valueX = f1->GetChisquare();
                vec_Chi2x.push_back(chi_valueX);

                Double_t p0 = f1->GetParameter(0);
                Double_t p1 = f1->GetParameter(1);

                vec_P0.push_back(p0);
                vec_P1.push_back(p1);
                Double_t extrpX = ((65 - p0) / p1);

    
                TGraph* chiY = new TGraph(vec_ClustersY.size(), clusterArrayY, clusterArrayZy);
                chiY->Fit(f2,"QN");
                chi_valueY = f2->GetChisquare();
                vec_Chi2y.push_back(chi_valueY);
                Double_t pp0 = f2->GetParameter(0);
                Double_t pp1 = f2->GetParameter(1);

                vec_PP0.push_back(pp0);
                vec_PP1.push_back(pp1);
                Double_t extrpY = ((65 - pp0) / pp1);
                Double_t sumChi = (chi_valueX + chi_valueY);

                cout<<"X CHI :"<<chi_valueX<<"\t"<<"Y CHI"<<chi_valueY<<endl;

                cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Chi X : "<<chi_valueX<<"\t"<<"Chi Y : "<<chi_valueY<<"^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
            }
        }
    }   

} 

    std::vector<SttHit*> vec_tracks;
    vec_tracks.clear();

    Float_t smallestX = vec_Chi2x[0];
    Float_t smallestP0 = vec_P0[0];
    Float_t smallestP1 = vec_P1[0];

    Int_t chi_indexX = 0;

    for (Int_t ci = 0; ci < vec_Chi2x.size(); ci++) 
    {
        cout << "chiSquare X :" << vec_Chi2x[ci] << endl;

        if (smallestX > vec_Chi2x[ci])
         {
            smallestX = vec_Chi2x[ci];
            chi_indexX = ci;
        }
        if (smallestP0 > vec_P0[ci]) 
        {
            smallestP0 = vec_P0[ci];
        }
        if (smallestP1 > vec_P1[ci]) 
        {
            smallestP1 = vec_P1[ci];
        }
    }

    Float_t smallestY = vec_Chi2y[0];
    Float_t smallestPP0 = vec_PP0[0];
    Float_t smallestPP1 = vec_PP1[0];

    Int_t chi_indexY = 0;

    for (Int_t cj = 0; cj < vec_Chi2y.size(); cj++) 
    {
        cout << "chiSquare Y :" << vec_Chi2y[cj] << endl;

        if (smallestY > vec_Chi2y[cj]) 
        {
            smallestY = vec_Chi2y[cj];
            chi_indexY = cj;
        }
        if (smallestPP0 > vec_PP0[cj]) 
        {
            smallestPP0 = vec_PP0[cj];
        }
        if (smallestPP1 > vec_PP1[cj]) 
        {
            smallestPP1 = vec_PP1[cj];
        }
    }
    cout<<"INDEX X:"<<chi_indexX<<"\t"<<"INDEX Y:"<<chi_indexY<<endl;

    for (Int_t ck = 0; ck < vec_All_X.at(chi_indexX).size(); ck++)
    {
        cout<<"       BEST  X      "<< vec_All_X.at(chi_indexX).at(ck)->cell<<"\t"<<smallestX<<endl;
        vec_tracks.push_back(vec_All_X.at(chi_indexX).at(ck));
    }
    for (Int_t cl = 0; cl < vec_All_Y.at(chi_indexY).size(); cl++)
    {
        cout<<"       BEST  Y      "<< vec_All_Y.at(chi_indexY).at(cl)->cell<<"\t"<<smallestY<<endl;
        vec_tracks.push_back(vec_All_Y.at(chi_indexY).at(cl));
    }

    Double_t sum_leadtime = 0.0;
    Double_t mean_leadtime = 0.0;
    vec_driftTime.clear();


    double sumLeadTime = 0;
    double meanTime = 0;

    for (Int_t d = 0; d<vec_tracks.size(); d++)
    {
        sumLeadTime += vec_tracks.at(d)->leadTime;

     }

     meanTime = sumLeadTime/vec_tracks.size();
vec_event.push_back(i);

//Write Tracks
    stt_event->TrackClear();

    SttTrackHit* b = stt_event->AddTrackHit();
    b->vec_Track = vec_tracks;
    b->trackId = i;
    b->trackSize = vec_tracks.size();
    b->Px0 = smallestP0;
    b->Px1 = smallestP1;
    b->Py0 = smallestPP0;
    b->Py1 = smallestPP1;
    b->Chix = smallestX;
    b->Chiy = smallestY;


    for (Int_t tq=0; tq<vec_tracks.size(); tq++)
    {
        cout<<"TRACKS  : "<<vec_tracks[tq]->layer<<"\t"<<vec_tracks[tq]->cell<<"\t"<<vec_tracks[tq]->channel<<"\t"<<vec_tracks[tq]->leadTime<<"\t"<<meanTime<<"\t"<<(fabs((vec_tracks[tq]->leadTime)-meanTime))<<endl;

     }


    vec_Chi2x.clear();
    vec_Chi2y.clear();

    vec_P0.clear();
    vec_P1.clear();

    vec_PP0.clear();
    vec_PP1.clear();

    PDAQ_EMC_STT_cluster_analysis->Fill();

    }

  }//End of loop over events

cout<<"CGECH THE SIZE   :"<<vec_All_tracks.size()<<endl;

    PDAQ_EMC_STT_cluster_analysis->Write();
for (int ko=0; ko<vec_event.size(); ko++)
{
 
cout<<"EVENT : "<<vec_event[ko]<<endl;
 }


    return 0;
}

