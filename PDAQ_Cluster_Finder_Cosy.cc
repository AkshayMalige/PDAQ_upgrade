#include "PDAQ_Cluster_Finder_Cosy.h"
#include "math.h"
#include <TPaveText.h>
#include <cmath>


using namespace std;





std::vector<SciHit*> ex_Scint_pileup ( PandaSubsystemSCI* SCI_CAL )
{
    std::vector<SciHit*> B;

    int scint_diff_limit = 300;

    int hit =0;

    SciHit* l =0;
    SciHit* l2=0;
    SciHit* l3=0;
    bool scint =false;
    //cout<<"SIZE "<<SCI_CAL->sci_raw.totalNTDCHits<<endl;
    if ( SCI_CAL->sci_raw.totalNTDCHits>2 ) {
        for ( int a=0; a< SCI_CAL->sci_raw.totalNTDCHits-1; a++ ) {

            SciHit* l2 = ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( a );
            SciHit* l3 = ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( a +1 );
            if ( a >=1 ) {
                SciHit* l = ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( a -1 );
                if ( a==SCI_CAL->sci_raw.totalNTDCHits-2 ) {
//                     if ( ( l2->leadTime - l->leadTime ) >=scint_diff_limit && ( l3->leadTime - l2->leadTime ) >=scint_diff_limit ) {
                    if ( ( l2->leadTime - l->leadTime ) >=scint_diff_limit && ( l3->leadTime - l2->leadTime ) >=scint_diff_limit /*&& ( l2->leadTime - l->leadTime ) <=600 && ( l3->leadTime - l2->leadTime ) <=600*/ )  {
                        B.push_back ( l2 );
                        B.push_back ( l3 );
                    } else if ( ( l3->leadTime - l2->leadTime ) >=scint_diff_limit /*&& ( l3->leadTime - l2->leadTime ) <=600*/ ) {
                        B.push_back ( l3 );
                    }
                } else {
                    if ( ( l2->leadTime - l->leadTime ) >=scint_diff_limit && ( l3->leadTime - l2->leadTime ) >=scint_diff_limit /*&& ( l2->leadTime - l->leadTime ) <=600 && ( l3->leadTime - l2->leadTime ) <=600*/ ) {
                        B.push_back ( l2 );
                    }
                }

            } else if ( a==0 && l3->leadTime-l2->leadTime >=scint_diff_limit /*&& l3->leadTime-l2->leadTime <=600*/ ) {
                B.push_back ( l2 );
            }

        }
    } else {
        if ( SCI_CAL->sci_raw.totalNTDCHits==1 ) {
            B.push_back ( ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( 0 ) );
        }


        if ( SCI_CAL->sci_raw.totalNTDCHits==2 ) {
            SciHit* h1 = ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( 0 );
            SciHit* h2 = ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( 1 );
            if ( ( h2->leadTime - h1->leadTime ) >=scint_diff_limit ) {
                B.push_back ( h1 );
                B.push_back ( h2 );
            } else {
                //B.push_back ( h1 );
            }
        }
    }



    return B;

}

std::vector<SttHit> ex_Stt_double ( std::vector<std::vector<SttHit>> vec_layer_channel_hit )
{
    int diff_limit = 600;
    std::vector<SttHit> B;
    std::vector<SttHit> C;

    for ( int a=0; a<vec_layer_channel_hit.size(); a++ ) {
        B = vec_layer_channel_hit[a];
        if ( B.size() >2 ) {
            for ( int aa=0; aa< B.size()-1; aa++ ) {
                //printf("%lf %lf \n",B[aa+1]->leadTime - B[aa]->leadTime, B[aa]->leadTime - B[aa-1]->leadTime);
                if ( aa ==0 ) {
                    if ( ( B[aa+1].layer==B[aa].layer && B[aa+1].straw==B[aa].straw && fabs ( fabs ( B[aa+1].leadTime ) - fabs ( B[aa].leadTime ) ) <=diff_limit ) ||  B[aa+1].leadTime  ==  B[aa].leadTime ) {
                        B.erase ( B.begin() +aa+1 );
                        --aa;
                    }

                } else if ( aa >=1 ) {
                    if ( ( B[aa+1].layer==B[aa].layer && B[aa+1].straw==B[aa].straw && fabs ( fabs ( B[aa+1].leadTime ) - fabs ( B[aa].leadTime ) ) <=diff_limit ) || B[aa+1].leadTime  ==  B[aa].leadTime ) {
                        B.erase ( B.begin() +aa+1 );
                        --aa;
                    }
                }
            }
        } else {
            if ( B.size() ==2 ) {
                if ( ( B[1].layer==B[0].layer && B[1].straw==B[0].straw && fabs ( fabs ( B[1].leadTime ) - fabs ( B[0].leadTime ) ) <=diff_limit ) || B[1].leadTime  ==  B[0].leadTime ) {
                    B.erase ( B.begin() +1 );
                }
            }
        }
        //cout<<B.size()<<endl;
        for ( int d=0; d<B.size(); d++ ) {
            //cout<<B[d].layer<<"\t"<<B[d].channel<<"\t"<<B[d].leadTime<<endl;
            C.push_back ( B[d] );
        }
        B.clear();
    }

    std::sort ( C.begin(), C.end(), f_sttHitCompareLeadTime );
    return C;
}

int PDAQ_Cluster_Finder_Cosy ( char* intree, char* outtree, int maxEvents )
{

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;
    PandaSubsystemSCI* SCI_CAL = 0;

    histograms* h = new histograms();

    TCanvas * DT; //Canvas for ToT
    DT=new TCanvas ( "DT_Canvas","DT_Canvas" );
    DT->Divide ( 4,2 );

    TCanvas * TOT; //Canvas for ToT
    TOT=new TCanvas ( "TOT_Canvas","TOT_Canvas" );
    TOT->Divide ( 4,2 );

    TCanvas * pDTvsTOT; //Canvas for ToT
    pDTvsTOT=new TCanvas ( "pDTvsTOT_Canvas","pDTvsTOT_Canvas" );
    pDTvsTOT->Divide ( 4,2 );


    TCanvas * DTvsTOT; //Canvas for ToT
    DTvsTOT=new TCanvas ( "DTvsTOT_Canvas","DTvsTOT_Canvas" );
    DTvsTOT->Divide ( 4,2 );

    TCanvas * pChMult; //Canvas for ToT
    pChMult=new TCanvas ( "pChMult_Canvas","pChMult_Canvas" );
    pChMult->Divide ( 4,2 );

    TCanvas * ChMult; //Canvas for ToT
    ChMult=new TCanvas ( "ChMult_Canvas","ChMult_Canvas" );
    ChMult->Divide ( 4,2 );

    TCanvas * ChTOT; //Canvas for ToT
    ChTOT=new TCanvas ( "ChTOT_Canvas","ChTOT_Canvas" );
    ChTOT->Divide ( 4,2 );

    TCanvas * LtDiff; //Canvas for ToT
    LtDiff=new TCanvas ( "LtDiff_Canvas","LtDiff_Canvas" );
    LtDiff->Divide ( 4,2 );

    TCanvas * LtDiffLow; //Canvas for ToT
    LtDiffLow=new TCanvas ( "LtDiffLow_Canvas","LtDiffLow_Canvas" );
    LtDiffLow->Divide ( 4,2 );

    TCanvas * Eff; //Canvas for ToT
    Eff=new TCanvas ( "Eff_Canvas","Eff_Canvas" );
    Eff->Divide ( 4,2 );

    TCanvas * Eff2; //Canvas for ToT
    Eff2=new TCanvas ( "Eff_Canvas2","Eff_Canvas2" );
    Eff2->Divide ( 4,2 );

    TCanvas * ChDiff; //Canvas for ToT
    ChDiff=new TCanvas ( "ChDiff_Canvas","ChDiff_Canvas" );
    ChDiff->Divide ( 4,2 );

    TCanvas * PlaneVsStraw; //Canvas for ToT
    PlaneVsStraw=new TCanvas ( "PlaneVsStraw_Canvas","PlaneVsStraw_Canvas" );
    PlaneVsStraw->Divide ( 4,4 );

    TCanvas * PlaneDTTOT; //Canvas for ToT
    PlaneDTTOT=new TCanvas ( "PlaneDTTOT_Canvas","PlaneDTTOT_Canvas" );
    PlaneDTTOT->Divide ( 4,4 );

    TCanvas *Lay_ef = new TCanvas ( "Lay_ef","Lay_ef",200,10,500,300 );
    TCanvas *Lay_ef1 = new TCanvas ( "Lay_ef1","Lay_ef1",200,10,500,300 );

    printf ( "%s\n", outtree );
    PandaSttCal* STT_CAL = new PandaSttCal();
    PandaSubsystemSCI* SCI_TRACKS = new PandaSubsystemSCI();
    SciEvent* sci_event = & ( SCI_TRACKS->sci_raw );

    PandaSttTrack* STT_TRACKS = new PandaSttTrack();
    Stt_Track_Event* stt_event = & ( STT_TRACKS->stt_track_can );

    TFile inFile ( intree );
    TTree* tree = ( TTree* ) inFile.Get ( "PDAQ_tree" );
    if ( !tree ) {
        std::cerr << "Tree PDAQ_tree was not found\n";
        std::exit ( 1 );
    }


    std::vector<Double_t> vec_event;
    std::vector<Double_t> vec_test;
    std::vector<Double_t> vec_driftTime;
    std::vector<Double_t> vec_roundoff;
    std::vector<Double_t> vec_occurance;
    std::vector<Double_t> vec_pos_DT;
    std::vector<Double_t> vec_cumsum;
    std::vector<Double_t> vec_drift_radius;

    std::vector<SttHit*> vec_All_tracks;

    tree->SetBranchAddress ( "STT_CAL", &STT_CAL );
    tree->SetBranchAddress ( "SCI_CAL", &SCI_CAL );

    TFile* Ttree = new TFile ( outtree, "RECREATE" );
    TTree* PDAQ_tree = new TTree ( "PDAQ_tree", "PDAQ_tree" );
    PDAQ_tree->Branch ( "STT_TRACKS", "PandaSttTrack", &STT_TRACKS, 64000, 99 );


    TH1F* h_STT_Hit_Diff =  new TH1F ( "h_STT_Hit_Diff", "h_STT_Hit_Diff", 1000, -100, 900 );
    //TH1F* Lh_STT_Hit_Diff[7];

    TH1F* h_FrontNO = new TH1F ( "h_FrontNO", "h_FrontNO", 20, 0, 20 );

    TH1F* h_Fee[18];

    Int_t iev = ( Int_t ) tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl << endl;

    int event_counter = 0;
    int All_hit_counter = 0;
    int Good_hit_counter = 0;
    int Layer_eq_4_counter = 0;
    int vec_LT_g7_counter = 0;
    int vec_LT_l250_counter = 0;
    double final_counter = 0;

    float max_lead_time_diff = 0;
    int min_track_hits = 0;
    bool Scint = false;
    float refTime = 0;
    double total_particles=0;

    double plane_eff_cont =0;

    MParManager* a = MParManager::instance();
    MFTGeomPar* ftGeomPar = new MFTGeomPar();

    a->setParamSource ( "ftparams.txt" );
    a->parseSource();
    pm()->addParameterContainer ( "MFTPar", ftGeomPar );

    max_lead_time_diff = ftGeomPar->getLeadtimeWindow();
    min_track_hits = ftGeomPar->getTrackMinSize();
    int stations = ftGeomPar->getModules();
    int max_cluster_intake = ftGeomPar->getClusterLimit();
    int MAX_FT_TOTAL_LAYERS = 0;

    int High_ch[8]= {8,43,69,104,136,171,197,232};
    int Low_ch[8]= {0,32,64,96,128,160,192,224};
    int L0[2]= {7,11};
    int L1[2]= {10,14};
    int L2[2]= {4,8};
    int L3[2]= {7,11};
    int L4[2]= {7,11};
    int L5[2]= {10,14};
    int L6[2]= {4,8};
    int L7[2]= {7,11};

    double pl[8][3];
    double pl2[8][3];
    for ( int pp=0; pp<8; pp++ ) {
        for ( int p=0; p<3; p++ ) {
            pl[pp][p]=0;
            pl2[pp][p]=0;
        }
    }

    //int High_strw[8]= {9,12,6,9,9,12,6,9};

    TH1F* h_Ch_TOT[256];
    TH1F* h_PlaneMult_Straw[15];

    int cases[3] ;
    for ( int g=0; g<3; g++ ) {
        cases[g] =0;
    }


    for ( int mc = 0; mc < 5; mc++ ) {
        h->h_Cross_TOT[mc] = new TH1F ( Form ( "Case_%d_TOT", mc+1 ) , Form ( "Case_%d_TOT", mc+1 ), 4000, 0,4000 );
        h->h_Cross_DT[mc] = new TH1F ( Form ( "Case_%d_DT", mc+1 ) , Form ( "Case_%d_dT", mc+1 ), 500, 0,500 );
        h->h_Cross_DTvTOT[mc] = new TH2F ( Form ( "Cross_%d_DTvsTOT",mc+1 ),Form ( "Cross_%d_DTvsTOT;Drift Time [ns];Time Over Threshold [ns]",mc+1 ), 500, 0, 500,4000,0,4000 );
    }

    for ( int kc =0; kc<3 ; kc++ ) {
        h->h_CrossTalkCase[kc]=new TH1F ( Form ( "Cross_TlkCase%d", kc+1 ) , Form ( "Cross_TlkCase%d", kc+1 ), 3, 0,3 );
        h->h_CrossTalkDT_TOT[kc]=new TH2F ( Form ( "h_CrossTalkDT_TOT%d",kc+1 ),Form ( "h_CrossTalkDT_TOT%d;Drift Time [ns];Time Over Threshold [ns]",kc+1 ), 500, 0,500,4000,0,4000 );
        h->h_CrossTalkDT_DT[kc]= new TH2F ( Form ( "h_CrossTalkDT_DT%d",kc+1 ),Form ( "h_CrossTalkDT_DT%d;Drift Time [ns];Drift Time [ns]",kc+1 ), 500, 0, 500,500,0,500 );
    }


    for ( int mk=0; mk<3; mk++ ) {
        h->h_Cross_TOTvTOT[mk] = new TH2F ( Form ( "Cross_%d_TOTvsTOT",mk+1 ),Form ( "Cross_%d_TOTvsTOT1;Time Over Threshold [ns];Time Over Threshold [ns]",mk+1 ), 4000, 0, 4000,4000,0,4000 );
        h->h_Cross_DTvDT[mk] = new TH2F ( Form ( "Cross_%d_DTvsDT",mk+1 ),Form ( "Cross_%d_DTvsDT;Drift Time [ns];Drift Time [ns]",mk+1 ), 500, 0, 500,500,0,500 );
        h->h_Cross_DTsum[mk] = new TH1F ( Form ( "Cross_%d_DTSum", mk+1 ) , Form ( "Cross_%d_DTSum;Drift Time [ns]", mk+1 ), 600, -100,500 );
    }


    for ( int mh = 0; mh < 8; mh++ ) {
        h->h_pL_layerDT[mh] = new TH1F ( Form ( "Layer_%d_DTa", mh+1 ) , Form ( "Layer_%d_DTa;Drift Time [ns]", mh+1 ), 600, -100,500 );
        h->h_pL_dtvstot[mh] = new TH2F ( Form ( "Layer_%d_DTvsTOTa",mh+1 ),Form ( "Layer_%d_DTvsTOTa;Drift Time [ns];Time Over Threshold [ns]",mh+1 ), 500, 0, 500,700,0,700 );
        h->h_pL_TOT[mh] = new TH1F ( Form ( "Layer_%d_TOTa", mh+1 ) , Form ( "Layer_%d_TOTa;Time Over Threshold [ns]", mh+1 ), 600, -100, 500 );
        h->h_pL_channel_mult[mh] = new TH2F ( Form ( "Layer_%d_multa", mh+1 ) , Form ( "Layer_%d_multa;Channel;Multiplicity", mh+1 ), 33, 0, 33, 20, 0, 20 );

        h->h_L_layerDT[mh] = new TH1F ( Form ( "Layer_%d_DTb", mh+1 ) , Form ( "Layer_%d_DTb;Drift Time [ns]", mh+1 ), 600, -100, 500 );
        h->h_L_dtvstot[mh] = new TH2F ( Form ( "Layer_%d_DTvsTOTb",mh+1 ),Form ( "Layer_%d_DTvsTOTb;Drift Time [ns];Time Over Threshold [ns]",mh+1 ), 500, 0, 500,700,0,700 );
        h->h_L_TOT[mh] = new TH1F ( Form ( "Layer_%d_TOTb", mh+1 ) , Form ( "Layer_%d_TOTb;Time Over Threshold [ns]", mh+1 ), 600, -100, 500 );
        h->h_L_channel_mult[mh] = new TH2F ( Form ( "Layer_%d_multb", mh+1 ) , Form ( "Layer_%d_multb;Channel;Multiplicity", mh+1 ), 33, 0, 33, 20, 0, 20 );

        h->h_pChDiff[mh] = new TH1F ( Form ( "Layer_%d_pChDiff", mh+1 ) , Form ( "Layer_%d_pChDiff", mh+1 ), 5, 0, 5 );

        h->h_ChDiff[mh] = new TH1F ( Form ( "Layer_%d_ChDiff", mh+1 ) , Form ( "Layer_%d_ChDiff", mh+1 ), 6, 0, 6 );
        h->h_straw[mh] = new TH1F ( Form ( "Layer_%d_h_straw", mh+1 ) , Form ( "Layer_%d_h_straw", mh+1 ), 32, 0, 32 );


    }


    for ( int pl=0; pl<16; pl++ ) {
        h->h_PlaneMult_Straw[pl] = new TH2F ( Form ( "Plane_%d", pl+1 ) , Form ( "Plane_%d;Channel No.;Layer No.", pl+1 ), 33, 0, 33, 10, 0, 10 );
        h->h_PlaneMult_DT_TOT[pl] = new TH2F ( Form ( "PlaneDTvTOT_%d", pl+1 ) , Form ( "PlaneDTvTOT_%d;Drift Time [ns];Time Over Threshold [ns]", pl+1 ),500, 0, 500,700,0,700 );
    }

    for ( int ch=0; ch<256; ch++ ) {
        h->h_Ch_Dt[ch] = new TH1F ( Form ( "Ch_%d_Dt", ch+1 ) , Form ( "Ch_%d_Dt", ch+1 ), 700, -100,600 );
        h->h_pLT_Diff[ch] = new TH1F ( Form ( "Ch_%d_LT_diff", ch+1 ) , Form ( "Ch_%d_LT_diff", ch+1 ), 500, 1, 501 );
        h->h_Ch_TOT[ch] = new TH1F ( Form ( "Ch_%d_TOT", ch+1 ) , Form ( "Ch_%d_TOT", ch+1 ), 700, 0,700 );
    }



    if ( stations > 1 ) {
        for ( int a = 0; a < stations; a++ ) {
            MAX_FT_TOTAL_LAYERS += ftGeomPar->getLayers ( a );
        }
    }

    else {
        MAX_FT_TOTAL_LAYERS = ftGeomPar->getLayers ( 0 );
    }

    double repeat = 0;
    double All_repeat = 0;
    double scint_event =0;
    
    int lets_cont =0;
    int lets_notcont =0;
    

    SttHit* hitOnLayer[MAX_FT_TOTAL_LAYERS][5000];
    
    maxEvents = (maxEvents==100000000) ? iev :  maxEvents;


    for ( Int_t i = 0; i < iev; i++ ) {
//cout<<"check1"<<endl;
        event_counter++;

        tree->GetEntry ( i );

        std::vector<SttHit> vec_stthits;
        // cout<<i<<endl;

        bool scint = false;
        bool stt = false;
        
        bool Event_Marking = false;
        h->h_mark_ana->Fill(1);

        
        if ( i % 1000 == 0 ) {
            cout << "entry no. " << i <<":\t"<<(i*100)/maxEvents<<"%"<<endl;;
            //std::cout << '.';
            //std::flush(std::cout);

           // cout << i << " " << flush; 
          //  this_thread::sleep_for(chrono::seconds(1)); 
        }
       // cout<<maxEvents<<endl;
        if ( i >= maxEvents ) {
            break;
        }
        //cout<<"check0  "<<MAX_FT_TOTAL_LAYERS<<endl;
        memset ( hitOnLayer, 0, MAX_FT_TOTAL_LAYERS * 5000 * sizeof ( SttHit* ) );
        int hitMultOnLayer[MAX_FT_TOTAL_LAYERS];

        //cout<<"check1"<<endl;

        for ( int h = 0; h < MAX_FT_TOTAL_LAYERS; h++ ) {
            hitMultOnLayer[h] = 0;
        }
        // Loop over the vector elements//////////////////////////////////////////////////////
        h->h_hitmultiplicity0->Fill ( STT_CAL->stt_cal.total_cal_NTDCHits );
        h->h_hitmultiplicity1->Fill ( SCI_CAL->sci_raw.totalNTDCHits );

        std::vector<SciHit*> vec_scihits;

        
        (SCI_CAL->sci_raw.totalNTDCHits ==0) ? h->h_tracker_check->Fill(0) : h->h_tracker_check->Fill(1);
        if (SCI_CAL->sci_raw.totalNTDCHits >0 && STT_CAL->stt_cal.total_cal_NTDCHits >7 )  h->h_tracker_check->Fill(2);
        
        
        if ( SCI_CAL->sci_raw.totalNTDCHits >0 ) {
            h->h_scint_mult_b->Fill ( SCI_CAL->sci_raw.totalNTDCHits );
        }

        for ( int n = 0; n < STT_CAL->stt_cal.total_cal_NTDCHits; n++ ) {
            h->h_hitBlock->Fill ( 1 );
        }

        std::vector<SciHit*> B;
        B.clear();
        for ( int pn = 0; pn < SCI_CAL->sci_raw.totalNTDCHits; pn++ ) {
            SciHit* b = ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( pn );
            B.push_back ( b );
            h->h_scnt_plup_cnts->Fill ( 1 );
            h->h_ShitBlock->Fill ( 0 );
        }

        /////////////////////////////// scint diff ///////////////////////////////////////////////
        if ( B.size() >=2 ) {
            for ( int a=0; a< B.size()-1; a++ ) {
                h->h_scint_timediffa->Fill ( B[a+1]->leadTime - B[a]->leadTime );
            }
        }
        else if(B.size() ==1){h->h_scint_timediffa->Fill(0);}

        

        if ( STT_CAL->stt_cal.total_cal_NTDCHits >0 ) {
            SttHit* m = ( SttHit* ) STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt ( 0 );
            if(m->marking == true) {
                Event_Marking =true;
                h->h_mark_ana->Fill(3);
            }
            else {
                h->h_mark_ana->Fill(2);
                Event_Marking =false;
                lets_notcont++;
            }
        }
        else {
            Event_Marking =false;
            h->h_mark_ana->Fill(2);
            lets_notcont++;
        }

        h->track_mult=0;
        h->track_marked = 0;
        h->track_unmarked = 0;
        if ( SCI_CAL->sci_raw.totalNTDCHits >0 ) {


            vec_scihits =  ex_Scint_pileup ( SCI_CAL ) ;


            if ( vec_scihits.size() >=2 ) {
                for ( int a=0; a< vec_scihits.size()-1; a++ ) {
                    h->h_scint_timediffb->Fill ( vec_scihits[a+1]->leadTime - vec_scihits[a]->leadTime );
                }
            }

            h->h_scint_mult_a->Fill ( vec_scihits.size() );
          //  printf("\nSIZE : %i\n",vec_scihits.size());
            scint_event+=vec_scihits.size();
           // cout<<"**********"<<endl;

            for ( int sn = 0; sn < vec_scihits.size(); sn++ ) {

                SciHit* sh = vec_scihits[sn];
                h->h_scnt_plup_cnts->Fill ( 2 );
                h->h_ShitBlock->Fill ( 1 );
                if ( sh->tdcid ==0x6500 && sh->channel==1 && vec_scihits.size() >1 ) {
                    scint = true;
                } else {
                    scint = false;
                }
                
//                 cout<<"Bin : "<<int(fabs(sh->leadTime/2560))<<"\t"<<sh->leadTime<<endl;
                if(int(fabs(sh->leadTime/320))==0 || int(fabs(sh->leadTime/320))==1){
//                     cout<<"0\t"<<"0\t"<<endl;
                    h->h_scint_bin->Fill(0);
                    h->h_scint_bin->Fill(0);
                    
                }
                else if(int(fabs(sh->leadTime/320))==30 || int(fabs(sh->leadTime/320))==31){
//                     cout<<"14\t"<<"14\t"<<endl;
                    h->h_scint_bin->Fill(126);
                    h->h_scint_bin->Fill(126);
                }
                else{
//                     cout<<int(fabs(sh->leadTime/2560)/2)<<"\t"<<int(fabs(sh->leadTime/2560)/2)-1<<endl;
                    h->h_scint_bin->Fill(int(fabs(sh->leadTime/320)/2));
                    h->h_scint_bin->Fill(int(fabs(sh->leadTime/320)/2)-1);
                }
                
//                 
//DT = straw - scint

              //  h->h_hitmultiplicity1->Fill ( STT_CAL->stt_cal.total_cal_NTDCHits );
              //  printf ( "\n SCINT : %lf\n",sh->leadTime );
//cout<<  vec_scihits.size()<<"\t"<< STT_CAL->stt_cal.total_cal_NTDCHits<<endl;

                for ( int n = 0; n < STT_CAL->stt_cal.total_cal_NTDCHits; n++ ) {
                    SttHit* cal_hit = ( SttHit* ) STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt ( n ); // retrieve particular hit
                //   printf ( "XXXXXX:TDC: %x  , Layer : %i, Straw : %i, LT: %lf\n",cal_hit->tdcid,cal_hit->layer,cal_hit->straw,cal_hit->leadTime );
                    All_hit_counter++;
                    h->h_scint_timediff->Fill ( cal_hit->leadTime - sh->leadTime );
                    h->h_LTvsLayer0->Fill ( cal_hit->layer, cal_hit->leadTime );
                    h->h_raw_leadtimes->Fill ( cal_hit->leadTime );
                    h->h_tot0->Fill ( cal_hit->tot );

                    // hit on reference channel
                    if ( cal_hit->isRef == true ) {
                        // cout<<"Reference Hit  -> "<<cal_hit->isRef<<"\t"<<cal_hit->layer<<"\t"<< endl;
                    } else {
                        if ( cal_hit->layer == 0 ) {
                            continue;
                        }
                        //printf("%d %d\n", cal_hit->layer - 1, hitMultOnLayer[cal_hit->layer - 1]);

                        hitOnLayer[cal_hit->layer - 1][hitMultOnLayer[cal_hit->layer - 1]] = cal_hit;
                        hitMultOnLayer[cal_hit->layer - 1]++;
                        Good_hit_counter++;
                    }
                }


//                 bool good_layers = true;
//                 for ( int c = 0; c < MAX_FT_TOTAL_LAYERS; c++ ) {
//                     if ( hitMultOnLayer[c] != 1 ) {
//                         good_layers = false;
//                         break;
//                     }
//                 }
                int layerCounter = 0;

                for ( int m = 0; m < MAX_FT_TOTAL_LAYERS; m++ ) {
                    if ( hitMultOnLayer[m] > 0 ) {
                        layerCounter++;
                    }
                }
                h->h_0LMultiplicity->Fill ( layerCounter );

                if ( layerCounter == MAX_FT_TOTAL_LAYERS ) {
//                if ( layerCounter >0 ) {

                    Layer_eq_4_counter++;
                    stt = true;

                    std::vector<SttHit> vec_leadTime;
                    int filtercnt = 0;

                    for ( int l = 0; l < MAX_FT_TOTAL_LAYERS; l++ ) {
                        for ( int r = 0; r < hitMultOnLayer[l]; r++ )

                        {
                            if ( hitOnLayer[l][r]->tot >0 ) {
                                vec_leadTime.push_back ( *hitOnLayer[l][r] );
                                h->h_tot1->Fill ( hitOnLayer[l][r]->tot );
                                //cout<<"Initial : "<<hitOnLayer[l][r]->layer<<"\t"<<hitOnLayer[l][r]->channel<<"\t"<<hitOnLayer[l][r]->leadTime<<endl;
                                //printf("%lf  %lf  %lf \n",sh->leadTime,hitOnLayer[l][r]->leadTime,hitOnLayer[l][r]->leadTime - sh->leadTime);
                                //printf ( "YYYYY:TDC: %x  , Layer : %i, Straw : %i, LT: %lf\n",hitOnLayer[l][r]->tdcid,hitOnLayer[l][r]->layer,hitOnLayer[l][r]->straw,hitOnLayer[l][r]->leadTime );
                            }

                        }
                    }
//cout<<vec_leadTime.size()<<"\t";
                    std::sort ( vec_leadTime.begin(), vec_leadTime.end(), f_sttHitCompareLeadTime );
                    double doublehitdiff = 0;

////////////////////////////////////////////////////////// Multiple entry rejection ///////////////////
                    if ( vec_leadTime.size() >1 ) {
                        for ( Int_t je = 0; je < vec_leadTime.size()-1; je++ ) {
                           // cout<<vec_leadTime[je].layer<<"\t"<<vec_leadTime[je].channel<<"\t"<<vec_leadTime[je].leadTime<<endl;
                            if ( ( vec_leadTime[je + 1].layer ) == ( vec_leadTime[je].layer ) && ( vec_leadTime[je + 1].channel == vec_leadTime[je].channel ) && ( vec_leadTime[je + 1].leadTime == vec_leadTime[je].leadTime ) ) {

                                vec_leadTime.erase ( vec_leadTime.begin() +je+1 );
                                --je;
                              //  cout<<"REPEATED"<<endl;

                            }
                        }
                    }
//cout<<vec_leadTime.size()<<endl;
////////////////////////////////////////////////////////// DT 500 ns window /////////////////////////////////////////
                    std::vector<SttHit> vec_leadTime_e;
                    int MultLayer[8];
                    int ChHitsMult[8][32];
                    for ( int la = 0; la < 8; la++ ) {
                        MultLayer[la] = 0;
                        for ( int ma =0; ma< 32; ma++ ) {
                            ChHitsMult[la][ma] = 0;
                        }
                    }
                    bool high_tot=false;
                    int dt_cnt =0;
                    
                    for ( int jx=0; jx < vec_leadTime.size() ; jx++ ) {
//                   	 if(vec_leadTime[jx].trigger_no == 0x1174e2e4) printf("Lay :%i  Str:%i  DT : %lf \n",vec_leadTime[jx].layer,vec_leadTime[jx].straw,vec_leadTime[jx].leadTime - sh->leadTime);
//cout<<"mark : "<<vec_leadTime[jx].marking<<endl;
                        if ( vec_leadTime[jx].leadTime - sh->leadTime >-1 && vec_leadTime[jx].leadTime - sh->leadTime <500 ) {
                            h->h_pL_layerDT[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].leadTime - sh->leadTime );
                            h->h_pL_TOT[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].tot );
                            h->h_pL_dtvstot[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].leadTime - sh->leadTime,vec_leadTime[jx].tot );
                            h->h_p_drifttime->Fill ( vec_leadTime[jx].leadTime - sh->leadTime );
                            h->h_p_tot->Fill ( vec_leadTime[jx].tot );
                            MultLayer[vec_leadTime[jx].layer-1]++;
                            ChHitsMult[vec_leadTime[jx].layer-1][vec_leadTime[jx].straw-1]++;
                            vec_leadTime_e.push_back ( vec_leadTime[jx] );
                            h->h_hitBlock->Fill ( 2 );
                            h->h_p_drifttimevstot->Fill ( vec_leadTime[jx].leadTime - sh->leadTime,vec_leadTime[jx].tot );
                            if ( vec_leadTime[jx].tot>800 ) high_tot=true;
                            dt_cnt++;
                        //    printf ( "IN: %x  %i  %i  %lf  %lf  %lf\n",vec_leadTime[jx].tdcid,vec_leadTime[jx].layer,vec_leadTime[jx].straw,vec_leadTime[jx].leadTime , sh->leadTime,vec_leadTime[jx].leadTime - sh->leadTime);
                        } else {
                               // printf ( "OUT: %x  %i  %i  %lf  %lf  %lf \n",vec_leadTime[jx].tdcid,vec_leadTime[jx].layer,vec_leadTime[jx].straw,vec_leadTime[jx].leadTime , sh->leadTime,vec_leadTime[jx].leadTime - sh->leadTime );
                        }
                    }
                    if (dt_cnt >7 )  h->h_tracker_check->Fill(3);
                    //   if ( high_tot==true ) {
//                       cout<<"SIZE  : "<<vec_leadTime.size()<<endl;
//                       for (int t=0; t<vec_leadTime.size(); t++){
//                         cout<<vec_leadTime[t].tot<<"\t"<<vec_leadTime[t].layer<<"\t"<<vec_leadTime[t].straw<<endl;
//                       }
//                     }
//                     high_tot==false;
//                     cout<<"*************************************************\n"<<endl;

                    int mult = 0;
                    int layer_eff_count =0;
                    for ( int aa=0; aa<8; aa++ ) {
                        if ( MultLayer[aa]> 0 ) {
                            mult++;
                        }
                        h->h_pLayerMult->Fill ( aa,MultLayer[aa] );
                        //printf("dt multiplicity%i \n",dtHitsMult[a]);
                        for ( int ba=0; ba<32; ba++ ) {
                            if ( ChHitsMult[aa][ba] >0 ) {
                                h->h_pL_channel_mult[aa]->Fill ( ba,ChHitsMult[aa][ba] );
                                layer_eff_count++;
                                //cout<<aa+1<<"\t"<<ba+1<<"\t"<<ChHitsMult[aa][ba]<<endl;
                            }
                        }
                    }
                    if ( MultLayer[0]>0 && MultLayer[7]>0 ) {
                        h->h_pLayer_eff->Fill ( layer_eff_count );
                    }
                    h->h_pLMultiplicity->Fill ( mult );
                    if ( mult==8 ) {
                        total_particles++;
                    }

                    //cout<<"First : "<<vec_leadTime_e.size() <<endl;
                    std::vector<SttHit> vec_channel_hit;
                    std::vector<std::vector<SttHit>> vec_layer_channel_hit;
                    std::vector<SttHit> vec_p_L5S9;
                    for ( Int_t jc = 0; jc < vec_leadTime_e.size(); jc++ ) {

                        if ( vec_leadTime_e[jc].layer == 5 && vec_leadTime_e[jc].straw==9 ) {
                            vec_p_L5S9.push_back ( vec_leadTime_e[jc] );
                        }

                    }
                    if ( vec_p_L5S9.size() >1 ) {
                        for ( int L=0; L<vec_p_L5S9.size()-1; L++ ) {
                            h->h_pL5_S9LTdiff->Fill ( vec_p_L5S9[L+1].leadTime - vec_p_L5S9[L].leadTime );
                            //cout<<vec_L5S9.size() <<"\t"<<vec_L5S9[L+1].leadTime<<"\t"<<vec_L5S9[L].leadTime<<endl;
                        }
                    }

                    for ( int l=0; l<8; l++ ) {

                        for ( int st=0; st<32; st++ ) {
                            vec_channel_hit.clear();
                            for ( Int_t jc = 0; jc < vec_leadTime_e.size(); jc++ ) {
                                if ( vec_leadTime_e[jc].layer-1 == l && vec_leadTime_e[jc].straw-1 == st ) {
                                    vec_channel_hit.push_back ( vec_leadTime_e[jc] );
                                }
                            }
                            vec_layer_channel_hit.push_back ( vec_channel_hit );
                        }
                    }

                    std::vector<SttHit> B_lt;
                    int ch_sq=0;

                    for ( int ae=0; ae<vec_layer_channel_hit.size(); ae++ ) {
                        B_lt = vec_layer_channel_hit[ae];
                        if ( B_lt.size() >1 ) {
                            for ( int aee=0; aee< B_lt.size()-1; aee++ ) {
                                ch_sq = ( ( B_lt[aee].layer-1 ) * 32 ) +B_lt[aee].straw-1;
                                h->h_pLT_Diff[ch_sq-1]->Fill ( B_lt[aee+1].leadTime - B_lt[aee].leadTime );

                                if ( ( B_lt[aee+1].leadTime - B_lt[aee].leadTime ) >0 && ( B_lt[aee+1].leadTime - B_lt[aee].leadTime ) <20 ) {
                                    //printf ( "Event : %i, Lay:  %i , Ch  :%i, Lt1 : %lf, Lt2 : %lf ,diff %f \n ",i,B_lt[aee].layer,B_lt[aee].straw,B_lt[aee].leadTime,B_lt[aee+1].leadTime, ( B_lt[aee+1].leadTime - B_lt[aee].leadTime ) );
                                    //exit;

                                }
                                h->h_Ch_TOT[ch_sq-1]->Fill ( B_lt[aee+1].tot );
                            }
                        }
                    }

                    std::vector<SttHit> vec_leadTime_f;
                    vec_leadTime_f = ex_Stt_double ( vec_layer_channel_hit );
                    //cout<<"Next : "<<vec_leadTime_e.size() <<endl;

                    h->h_hitmultiplicity3->Fill ( vec_leadTime_f.size() );

                    int MultLayer2[8];
                    int ChHitsMult2[8][32];
                    for ( int lb = 0; lb < 8; lb++ ) {
                        MultLayer2[lb] = 0;
                        for ( int mb =0; mb< 32; mb++ ) {
                            ChHitsMult2[lb][mb] = 0;
                        }
                    }

                    std::vector<SttHit> vec_L5S9;
                    vec_stthits.clear();
                    int sq_ch=0;
                    
                //    cout<<"SIZE  : "<<vec_stthits.size()<<endl;
                    for ( Int_t jc = 0; jc < vec_leadTime_f.size(); jc++ ) {
                     //   if(vec_leadTime_f[jc].layer == 1 ||vec_leadTime_f[jc].layer == 4||vec_leadTime_f[jc].layer == 5||vec_leadTime_f[jc].layer == 8){
                            vec_leadTime_f[jc].drifttime = vec_leadTime_f[jc].leadTime - sh->leadTime;
                            vec_stthits.push_back ( vec_leadTime_f[jc] );
                            h->h_L_layerDT[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime );
                            h->h_drifttime->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime );
                            h->h_tot->Fill ( vec_leadTime_f[jc].tot );
                            h->h_hitBlock->Fill ( 3 );

                            sq_ch = ( ( vec_leadTime_f[jc].layer-1 ) * 32 ) +vec_leadTime_f[jc].straw-1;
                            h->h_Ch_Dt[sq_ch]->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime );
                            h->h_sq_ch->Fill ( sq_ch );
                            h->h_straw[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].straw );
                            
                             //if(vec_leadTime_f[jc].trigger_no == 0x1174e2e4)
                            //cout<<vec_leadTime_f[jc].layer<<"\t"<<vec_leadTime_f[jc].straw<<"\t"<<vec_leadTime_f[jc].drifttime<<"\t"<<vec_leadTime_f[jc].marking<<endl;
                    //    }

                    }//cout<<"*******\n\n"<<endl;

                    if ( vec_stthits.size() >= 8 && vec_stthits.size() <= max_cluster_intake) {
                            PDAQ_Event_Finder ( vec_stthits, i, PDAQ_tree, stt_event, ftGeomPar, SCI_CAL, h );
                            final_counter++;                       
                        
                    }                    
  
                    
                    // }
                }
            }
                     

        }
        
            if(h->track_mult>0){
                h->h_track_mult->Fill(h->track_mult);
            }
            if(h->track_marked >0){
                h->h_mark_ana->Fill(4);
            }
            if(h->track_unmarked >0){
                h->h_mark_ana->Fill(5);
            }
            if(Event_Marking==true && h->track_mult==0){
                h->h_mark_ana->Fill(6);
            }
            if(Event_Marking==false && h->track_mult>0){
                    h->h_mark_ana->Fill(7);
                    
            }
            if(Event_Marking==false && h->track_unmarked==0){
                    h->h_mark_ana->Fill(8);
                    lets_cont++;
            }
// if(h->track_mult>0 && h->track_unmarked >0)cout<<h->track_mult <<"="<<h->track_marked<<"+"<<h->track_unmarked<<endl;
//             if(Event_Marking==false) cout<<"NOT MARKED "<<lets_notcont<<"\t"<<h->track_mult <<"="<<h->track_marked<<"+"<<h->track_unmarked;
//             if(Event_Marking==false && h->track_mult==0) cout<<"**********NOT MARKED  NO TRACKS "<<h->track_mult <<"="<<h->track_marked<<"+"<<h->track_unmarked<<"//"<<lets_cont<<endl;

    } // End of loop over events

//    cout<<"CASE COUNT  "<<cases[0]<<endl;
    double maxbin =0;
    double DTmin =0;
    double dt_at_10=0;
    Int_t maximum = h->h_drifttime->GetBinContent ( h->h_drifttime->GetMaximumBin() );
    DTmin = h->h_drifttime->FindFirstBinAbove ( maximum/10,1 );
    dt_at_10 = h->h_drifttime->GetXaxis()->GetBinCenter ( h->h_drifttime->GetXaxis()->FindBin ( DTmin ) );

//     cout<< dt_at_10 - 100<<endl;
//
//    cout<<"Plane Mult "<<pl2[0][0]<<"\t"<<pl2[0][1]<<"\t"<<pl2[7][0]<<"\t"<<pl2[7][1]<<endl;

    h_STT_Hit_Diff->Write();
    h->h_cluster_size->Write();
    h->h_LTvsLayer0->Write();
    h->h_LTvsLayer1->Write();
    h->h_hitmultiplicity0->Write();
    
    h->HitMultiplicity->cd();
    h->h_hitmultiplicity0->Draw();
    h->h_hitmultiplicity1->Draw ( "same" );
    h->h_hitmultiplicity2->Draw ( "same" );
    h->h_hitmultiplicity3->Draw ( "same" );
    h->h_hitmultiplicity4->Draw ( "same" );
    h->h_hitmultiplicity0->SetLineColor ( kBlack );
    h->h_hitmultiplicity1->SetLineColor ( kMagenta );
    h->h_hitmultiplicity2->SetLineColor ( kTeal-8 );
    h->h_hitmultiplicity3->SetLineColor ( kGreen );
    h->h_hitmultiplicity4->SetLineColor ( kOrange+1 );
    h->h_hitmultiplicity0->SetLineWidth ( 3 );
    h->h_hitmultiplicity1->SetLineWidth ( 3 );
    h->h_hitmultiplicity2->SetLineWidth ( 3 );
    h->h_hitmultiplicity3->SetLineWidth ( 3 );
    h->h_hitmultiplicity4->SetLineWidth ( 3 );
    TLegend* leg = new TLegend ( 0.5,0.7,0.9,0.9 );
    leg->SetHeader ( "Multiplicity for hits in" );
    leg->SetFillColor ( 1 );
    leg->AddEntry ( h->h_hitmultiplicity0,"In an STT event","lep" );
    leg->AddEntry ( h->h_hitmultiplicity1,"In an Sint event","lep" );
    leg->AddEntry ( h->h_hitmultiplicity2,"Events with hit in all layers ","lep" );
    leg->AddEntry ( h->h_hitmultiplicity3,"After rejecting double hits","lep" );
    leg->AddEntry ( h->h_hitmultiplicity4,"After event-cluster finding","lep" );
    leg->SetFillStyle ( 0 );
    leg->Draw();
    h->HitMultiplicity->SetLogy();
    h->HitMultiplicity->Write();

    h->h_drifttimevstot->Write();
    h->h_drifttimevstot->GetXaxis()->SetLabelSize ( 0.055 );
    h->h_drifttimevstot->GetXaxis()->SetNdivisions ( 310 );
    h->h_drifttimevstot->GetYaxis()->SetLabelSize ( 0.055 );
    h->h_drifttimevstot->GetYaxis()->SetNdivisions ( 310 );

    h->TimeOverThreshold->cd();
    h->h_tot0->Draw();
    h->h_tot1->Draw ( "same" );
    h->h_tot2->Draw ( "same" );
    h->h_tot3->Draw ( "same" );
    h->h_tot4->Draw ( "same" );
    h->h_tot0->SetLineColor ( kOrange+1 );
    h->h_tot1->SetLineColor ( kMagenta );
    h->h_tot2->SetLineColor ( kTeal-8 );
    h->h_tot3->SetLineColor ( kGreen );
    h->h_tot4->SetLineColor ( kBlack );
    TLegend* leg1 = new TLegend ( 0.5,0.7,0.9,0.9 );
    leg1->SetHeader ( "Time Over Threshold for hits" );
    leg1->SetFillColor ( 1 );
    leg1->AddEntry ( h->h_tot0,"In an event","lep1" );
    leg1->AddEntry ( h->h_tot1,"Events with hit in all layers","lep1" );
    leg1->AddEntry ( h->h_tot2,"After rejecting double hits ","lep1" );
    leg1->AddEntry ( h->h_tot3,"After event cluster finding","lep1" );
    leg1->AddEntry ( h->h_tot4,"Track candiadtes","lep1" );
    leg1->SetFillStyle ( 0 );

    h->h_drifttimevsLayer->Write();


    h->h_drifttime->GetXaxis()->SetNdivisions ( 310 );
    h->h_drifttime->GetYaxis()->SetLabelSize ( 0.055 );
    h->h_drifttime->GetXaxis()->SetTitleSize ( 0.050 );
    h->h_drifttime->GetYaxis()->SetNdivisions ( 310 );
    h->h_drifttime->Write();

    h->h_tot->GetXaxis()->SetNdivisions ( 310 );
    h->h_tot->GetYaxis()->SetNdivisions ( 310 );
    h->h_tot->GetXaxis()->SetLabelSize ( 0.055 );
    h->h_tot->GetYaxis()->SetLabelSize ( 0.055 );
    h->h_tot->GetXaxis()->SetTitleSize ( 0.050 );
    h->h_tot->Write();


    h->h_scint_mult_b->Write();
    h->h_scint_mult_a->Write();
    h->h_scnt_plup_cnts->Write();

    h->h_scint_timediff->Write();
    leg1->Draw();
    h->TimeOverThreshold->Write();

    h->h_raw_leadtimes->Write();

    h->h_drifttimeTRB1->Write();
    h->h_drifttimeTRB2->Write();
    h->h_p_drifttimevstot->Write();
    h->h_p_drifttime->Write();
    h->h_drifttime->Write();
    h->h_p_tot->Write();


    h->h_scint_timediffa->SetLineColor ( kOrange+1 );
    h->h_scint_timediffa->GetXaxis()->SetLabelSize ( 0.055 );
    h->h_scint_timediffa->GetYaxis()->SetLabelSize ( 0.055 );
    h->h_scint_timediffa->GetXaxis()->SetTitleSize ( 0.050 );
    h->h_scint_timediffa->Write();

    h->h_scint_timediffb->SetLineColor ( kOrange+1 );
    h->h_scint_timediffb->GetXaxis()->SetLabelSize ( 0.055 );
    h->h_scint_timediffb->GetYaxis()->SetLabelSize ( 0.055 );
    h->h_scint_timediffb->GetXaxis()->SetTitleSize ( 0.050 );
    h->h_scint_timediffb->Write();

//     h->h_L5_S9LTdiff->Write();
//     h->h_pL5_S9LTdiff->Write();
//     h->h_pLayerMult->Write();
//     h->h_LayerMult->Write();
//     h->h_hitBlock->Write();
//     h->h_ShitBlock->Write();
//     h->h_ShitBlock->Write();
//     h->h_0LMultiplicity->Write();
//     h->h_pLMultiplicity->Write();
//     h->h_LMultiplicity->Write();
// 
//     h->h_L_layerEff->Write();
// 
//     double norm = h->h_Layer_eff->GetEntries();
//     h->h_Layer_eff->Scale ( 1/norm );
//     h->h_Layer_eff->Write();
// 
//     double norm2 = h->h_Layer_effPlane->GetEntries();
//     h->h_Layer_effPlane->Scale ( 1/norm2 );
//     h->h_Layer_effPlane->Write();
// 
//     double norm4 = h->h_Layer_effStraw->GetEntries();
//     h->h_Layer_effStraw->Scale ( 1/norm4 );
//     h->h_Layer_effStraw->Write();
// 
//     h->h_pLayer_eff->Write();
//     h->h_sq_ch->Write();


    /*for ( int hh = 0; hh < 8; hh++ ) {


        h->h_pL_layerDT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pL_layerDT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_pL_layerDT[hh]->GetXaxis()->SetTitleSize ( 0.045 );
        h->h_pL_layerDT[hh]->SetLineColor ( kOrange+1 );
        h->h_pL_layerDT[hh]->SetLineWidth ( 2 );
        h->h_pL_layerDT[hh]->Write();

        h->h_L_layerDT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_L_layerDT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_L_layerDT[hh]->GetXaxis()->SetTitleSize ( 0.045 );
        h->h_L_layerDT[hh]->SetLineColor ( kTeal-8 );
        h->h_L_layerDT[hh]->SetLineWidth ( 2 );
        h->h_L_layerDT[hh]->Write();
//------------------------------------------------------------------------

        h->h_pL_TOT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pL_TOT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_pL_TOT[hh]->GetXaxis()->SetTitleSize ( 0.05 );
        h->h_pL_TOT[hh]->SetLineColor ( kOrange+1 );
        h->h_pL_TOT[hh]->SetLineWidth ( 2 );
        h->h_pL_TOT[hh]->Write();



        h->h_L_TOT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_L_TOT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_L_TOT[hh]->GetXaxis()->SetTitleSize ( 0.045 );
        h->h_L_TOT[hh]->SetLineColor ( kTeal-8 );
        h->h_L_TOT[hh]->SetLineWidth ( 2 );
        h->h_L_TOT[hh]->Write();
//------------------------------------------------------------------------
        h->h_pL_dtvstot[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pL_dtvstot[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_pL_dtvstot[hh]->GetXaxis()->SetTitleSize ( 0.05 );
        h->h_pL_dtvstot[hh]->GetYaxis()->SetTitleSize ( 0.05 );
        h->h_pL_dtvstot[hh]->GetXaxis()->SetNdivisions ( 5 );
        h->h_pL_dtvstot[hh]->Write();

        h->h_L_dtvstot[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_L_dtvstot[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_L_dtvstot[hh]->GetXaxis()->SetTitleSize ( 0.05 );
        h->h_L_dtvstot[hh]->GetYaxis()->SetTitleSize ( 0.05 );
        h->h_L_dtvstot[hh]->GetXaxis()->SetNdivisions ( 5 );
        h->h_L_dtvstot[hh]->Write();
//------------------------------------------------------------------------

        h->h_pL_channel_mult[hh]->Write();



        h->h_L_channel_mult[hh]->Write();
        h->h_L_channel_mult[hh]->GetXaxis()->SetLabelSize ( 0.045 );

        ChDiff->cd ( hh+1 );
        h->h_ChDiff[hh]->Draw();
        Double_t norm1 = h->h_pChDiff[hh]->GetEntries();
        h->h_pChDiff[hh]->Scale ( 1/norm1 );
        h->h_pChDiff[hh]->Write();
        h->h_pChDiff[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pChDiff[hh]->GetYaxis()->SetLabelSize ( 0.045 );


//         for ( int a=1; a<8; a++ ) {
//             int b = h->h_ChDiff[hh]->GetBinContent ( a );
//             stringstream ss;
//             ss << b;
//             TString str = ss.str();
//             h->h_ChDiff[hh]->GetXaxis()->SetBinLabel ( a,str );
//         }
        Double_t norm2 = h->h_ChDiff[hh]->GetEntries();
        h->h_ChDiff[hh]->Scale ( 1/norm2 );
        h->h_ChDiff[hh]->SetLineWidth ( 3 );
        h->h_ChDiff[hh]->SetLineColor ( kTeal-8 );
        h->h_ChDiff[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_ChDiff[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_ChDiff[hh]->GetXaxis()->SetTitleSize ( 0.045 );
        h->h_ChDiff[hh]->GetYaxis()->SetTitleSize ( 0.045 );
        h->h_ChDiff[hh]->SetMarkerSize ( 3 );
        h->h_ChDiff[hh]->SetMarkerColor ( kOrange+1 );
        h->h_ChDiff[hh]->Draw ( "h,TEXT60" );
        h->h_ChDiff[hh]->Write();


        h->h_straw[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_straw[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_straw[hh]->Write();

        DT->cd ( hh+1 );
        h->h_pL_layerDT[hh]->Draw();
        h->h_pL_layerDT[hh]->SetLineWidth ( 2 );
        h->h_pL_layerDT[hh]->SetLineColor ( kOrange+1 );

        h->h_L_layerDT[hh]->Draw ( "same" );
        h->h_L_layerDT[hh]->SetLineWidth ( 2 );
        h->h_L_layerDT[hh]->SetLineColor ( kTeal-8 );

        TOT->cd ( hh+1 );
        h->h_pL_TOT[hh]->Draw();
        h->h_pL_TOT[hh]->SetLineWidth ( 2 );
        h->h_pL_TOT[hh]->SetLineColor ( kOrange+1 );

        h->h_L_TOT[hh]->Draw ( "same" );
        h->h_L_TOT[hh]->SetLineWidth ( 2 );
        h->h_L_TOT[hh]->SetLineColor ( kTeal-8 );

        pDTvsTOT->cd ( hh+1 );
        h->h_pL_dtvstot[hh]->Draw ( "colz" );

        DTvsTOT->cd ( hh+1 );
        h->h_L_dtvstot[hh]->Draw ( "colz" );

        pChMult->cd ( hh+1 );
        h->h_pL_channel_mult[hh]->Draw ( "colz" );

        ChMult->cd ( hh+1 );
        h->h_L_channel_mult[hh]->Draw ( "colz" );

        ChTOT->cd ( hh+1 );
        h->h_Ch_TOT[High_ch[hh]]->Draw ();

        LtDiff->cd ( hh+1 );
        h->h_pLT_Diff[High_ch[hh]]->Draw ();
        gPad->SetLogy();

        LtDiffLow->cd ( hh+1 );
        h->h_pLT_Diff[Low_ch[hh]]->Draw ();
        gPad->SetLogy();

    }*/


//     std::vector<double> vec_DT_start;
//     cout<<norm<<"\t"<<total_particles<<endl;
//     float sc=1/total_particles;

//     for ( int cha=0; cha<256; cha++ ) {
//         h->h_Ch_Dt[cha]->Scale ( 1/total_particles );
//         h->h_Ch_Dt[cha]->Write();
//         h->h_Ch_Dt[cha]->GetXaxis()->SetLabelSize ( 0.045 );
//         h->h_Ch_Dt[cha]->GetYaxis()->SetLabelSize ( 0.045 );
// 
//         h->h_pLT_Diff[cha]->Scale ( 1/total_particles );
//         h->h_pLT_Diff[cha]->GetXaxis()->SetLabelSize ( 0.045 );
//         h->h_pLT_Diff[cha]->GetYaxis()->SetLabelSize ( 0.045 );
//         h->h_pLT_Diff[cha]->GetXaxis()->SetTitle ( "Time Diff [ns]" );
//         h->h_pLT_Diff[cha]->GetXaxis()->SetTitleSize ( 0.045 );
//         h->h_pLT_Diff[cha]->GetXaxis()->SetNdivisions ( 5 );
//         h->h_pLT_Diff[cha]->Write();
// 
//         h->h_Ch_TOT[cha]->Scale ( 1/total_particles );
//         h->h_Ch_TOT[cha]->GetXaxis()->SetLabelSize ( 0.045 );
//         h->h_Ch_TOT[cha]->GetYaxis()->SetLabelSize ( 0.045 );
//         h->h_Ch_TOT[cha]->GetXaxis()->SetTitle ( "Time Over Threshold [ns]" );
//         h->h_Ch_TOT[cha]->GetXaxis()->SetTitleSize ( 0.045 );
//         h->h_Ch_TOT[cha]->GetXaxis()->SetNdivisions ( 5 );
//         h->h_Ch_TOT[cha]->Write();
//     }

//     for(int dot=0; dot<vec_DT_start.size(); dot++){
//       printf("Channel : %i   DT_crr: %d \n",dot,vec_DT_start[dot]-dt_at_10);
//     }
//     for ( int x=0; x<16; x++ ) {
// 
//         h->h_PlaneMult_Straw[x]->GetXaxis()->SetLabelSize ( 0.045 );
//         h->h_PlaneMult_Straw[x]->GetYaxis()->SetLabelSize ( 0.045 );
//         h->h_PlaneMult_Straw[x]->GetXaxis()->SetTitleSize ( 0.05 );
//         h->h_PlaneMult_Straw[x]->GetYaxis()->SetTitleSize ( 0.05 );
//         h->h_PlaneMult_Straw[x]->GetXaxis()->SetLabelSize ( 0.045 );
//         h->h_PlaneMult_DT_TOT[x]->GetYaxis()->SetLabelSize ( 0.045 );
//         h->h_PlaneMult_DT_TOT[x]->GetXaxis()->SetTitleSize ( 0.05 );
//         h->h_PlaneMult_DT_TOT[x]->GetYaxis()->SetTitleSize ( 0.05 );
//         h->h_PlaneMult_DT_TOT[x]->Write();
//         h->h_PlaneMult_DT_TOT[x]->Write();
//         PlaneVsStraw->cd ( x+1 );
//         h->h_PlaneMult_Straw[x]->Draw ( "colz" );
//         PlaneDTTOT->cd ( x+1 );
//         h->h_PlaneMult_DT_TOT[x]->Draw ( "colz" );
//     }


    DT->Write();
    TOT->Write();
    pDTvsTOT->Write();
    DTvsTOT->Write();
    pChMult->Write();
    ChMult->Write();
    LtDiff->Write();
    LtDiffLow->Write();
    ChTOT->Write();
    ChDiff->Write();

    PlaneVsStraw->Write();
    PlaneDTTOT->Write();


    Double_t norm3 = h->h_CrossCase->GetEntries();
    //h->h_CrossCase->Scale ( 1/norm3 );
    h->h_CrossCase->Scale ( 1/final_counter );
    h->h_CrossCase->SetLineWidth ( 3 );
    h->h_CrossCase->SetLineColor ( kOrange+1 );
    h->h_CrossCase->SetFillColor ( kOrange+1 );
    h->h_CrossCase->SetMarkerSize ( 3 );
    h->h_CrossCase->SetMarkerColor ( kOrange+1 );
    h->h_CrossCase->Draw ( "h,TEXT" );
    h->h_CrossCase->Write();
    
    h->h_CrossMaxTOT->Scale ( 1/h->h_CrossMaxTOT->GetEntries() );
    h->h_CrossMaxTOT->SetLineWidth ( 3 );
    h->h_CrossMaxTOT->SetLineColor ( kTeal-8 );
    h->h_CrossMaxTOT->SetFillColor ( kTeal-8 );
    h->h_CrossMaxTOT->SetMarkerSize ( 3 );
    h->h_CrossMaxTOT->SetMarkerColor ( kTeal-8 );
    h->h_CrossMaxTOT->GetYaxis()->SetRangeUser ( 0,1 );
    h->h_CrossMaxTOT->Draw ( "h,TEXT" );
    h->h_CrossMaxTOT->Write();

    double lay_ef[8];
    double lay_ef2[8];
    double lay_index[8];

    for ( int u=0; u<8; u++ ) {
//         lay_ef[u]=pl2[u][2]/pl2[u][1];
//         lay_ef2[u]=pl2[u][2]/pl2[u][0];
        lay_index[u]=u+1;

        lay_ef[u]=pl2[u][2]/plane_eff_cont;
        lay_ef2[u]=pl2[u][2]/plane_eff_cont;
    }

    for ( int xc=0; xc<5; xc++ ) {

        h->h_Cross_DT[xc]->Write();
        h->h_Cross_TOT[xc]->Write();
        h->h_Cross_DTvTOT[xc]->Write();
    }

    for ( int xd=0; xd<3; xd++ ) {

        h->h_Cross_TOTvTOT[xd]->Write();
        h->h_Cross_DTvDT[xd]->Write();
        h->h_Cross_DTsum[xd]->Write();
    }

    Lay_ef->cd();
    TGraph* gLayerEff = new TGraph ( 8, lay_index, lay_ef );
    gLayerEff->Draw ( "AB" );
    gLayerEff->SetTitle ( "Plane_0" );
    gLayerEff->GetXaxis()->SetTitle ( "Layer" );
    gLayerEff->Write();
    gLayerEff->SetFillColor ( kOrange+1 );

    Lay_ef1->cd();
    TGraph* gLayerEff2 = new TGraph ( 8, lay_index, lay_ef2 );
    gLayerEff2->Write();
    gLayerEff2->Draw ( "AB" );
    gLayerEff2->SetTitle ( "Plane_1" );
    gLayerEff2->GetXaxis()->SetTitle ( "Layer" );
    gLayerEff2->SetFillColor ( kOrange+1 );

    Lay_ef->Write();
    Lay_ef1->Write();

    h->h_Corr_Layer_Mult->Scale ( 1/h->h_Corr_Layer_Mult->GetEntries() );
    h->h_Corr_Layer_Mult->Write();
    h->h_Corr_Plane_Mult->Scale ( 1/h->h_Corr_Plane_Mult->GetEntries() );
    h->h_Corr_Plane_Mult->Write();
    h->h_Corr_Channel_Mult->Scale ( 1/h->h_Corr_Channel_Mult->GetEntries() );
    h->h_Corr_Channel_Mult->Write();

    for ( int y=0; y<2; y++ ) {
        h->h_CrossTalkCase[y]->Scale ( 1/final_counter );
        h->h_CrossTalkCase[y]->SetLineWidth(3);
        h->h_CrossTalkCase[y]->SetLineColor(kOrange+1);
        h->h_CrossTalkCase[y]->SetFillColor(kOrange+1);
        h->h_CrossTalkCase[y]->Write();
        h->h_CrossTalkDT_TOT[y]->Write();
        h->h_CrossTalkDT_DT[y]->Write();
    }
    
    h->h_High_TOT_Layer->Scale(1/h->h_High_TOT_Layer->GetEntries());
	h->h_High_TOT_Layer->SetLineWidth(3);
	h->h_High_TOT_Layer->SetLineColor(kOrange+1);
	h->h_High_TOT_Layer->SetFillColor(kOrange+1);
    h->h_High_TOT_Layer->Write();
    
    h->h_scint_bin->Write();
    h->h_pairs_in_lay->Write();
    
    //h->h_tracker_check->Scale(1/h->h_tracker_check->GetBinContent(1));
    h->h_tracker_check->SetLineWidth(3);
    h->h_tracker_check->Write();
    
  //  cout<<"tracks : "<<h->h_marker_check->GetBinContent(1)<<endl;
 //       cout<<"Total no. events : "<<(float)maxEvents<<" total tracks : "<<h->track_count<<endl;

    Float_t tot_events = (float)iev;
   // h->h_marker_check->Scale(1/tot_events);
    h->h_track_mult->Write();
    h->h_marker_check->Write();
    h->h_mark_ana->Write();
    
    PDAQ_tree->Write();
    Ttree->Close();
    //inFile.Close();
    printf("\n\nEvents :\n Total events :\t\t\t\t\t%f\t100%\n Total marked events :\t\t\t\t%f\t%f%\n Total marked events with tracks :\t\t%f\t%f%\n Total marked events with no tracks :\t\t%f\t%f%\n Total unmarked events :\t\t\t%f\t%f%\n Total unmarked events with tracks :\t\t%f\t%f%\n Total unmarked events with no tracks :\t\t%f\t%f%\n\nTracks :\n Total tracks :\t\t\t\t\t%f\t100%\n Total tracks in marked events :\t\t%f\t%f%\n Total tracks in unmarked events :\t\t%f\t%f%\n\n",                        h->h_mark_ana->GetBinContent(2),   /*Total events*/                                                                                                                                      h->h_mark_ana->GetBinContent(4),h->h_mark_ana->GetBinContent(4)*100/h->h_mark_ana->GetBinContent(2),   /*Total marked events*/                                                                  h->h_mark_ana->GetBinContent(5),h->h_mark_ana->GetBinContent(5)*100/h->h_mark_ana->GetBinContent(2),   /*Total marked events with tracks*/                                                                  h->h_mark_ana->GetBinContent(7),h->h_mark_ana->GetBinContent(7)*100/h->h_mark_ana->GetBinContent(2),   /*Total marked events with no tracks*/                                                                  h->h_mark_ana->GetBinContent(3),h->h_mark_ana->GetBinContent(3)*100/h->h_mark_ana->GetBinContent(2),   /*Total unmarked events*/                                                                  h->h_mark_ana->GetBinContent(8),h->h_mark_ana->GetBinContent(8)*100/h->h_mark_ana->GetBinContent(2),   /*Total unmarked events with tracks*/                                                                                    h->h_mark_ana->GetBinContent(9),h->h_mark_ana->GetBinContent(9)*100/h->h_mark_ana->GetBinContent(2),   /*Total unmarked events with no tracks*/                                                                  h->track_count,                                                                                        /*Total tracks*/                                                                  h->h_marker_check->GetBinContent(3),h->h_marker_check->GetBinContent(3)*100/h->track_count,            /*Total tracks in marked events*/                                                                  h->h_marker_check->GetBinContent(4), h->h_marker_check->GetBinContent(4)*100/h->track_count);          /*Total tracks in unmarked events*/
    
    cout<<"lets_cont "<<lets_cont<<endl;
        
//     printf ( "Total Hits processed : %f Repeated Hits in an Event : %f\n\n", All_repeat, repeat );
    printf ( "In_File: %s   Out_File:  %s\n", intree, outtree );
    return 0;
}

int main ( int argc, char** argv )
{

    if ( argc >= 3 )

        if ( !argv[3] ) {
            printf ( "\n\nNote : One Million Events will be processed. To change "
                     "add the number of events to be processed after the ouput "
                     "file name.\n" );
            // atoi(argv[3]) == 1000;
            sleep ( 2 );
            -         PDAQ_Cluster_Finder_Cosy ( argv[1], argv[2], 100000000 );
        } else {
            PDAQ_Cluster_Finder_Cosy ( argv[1], argv[2], atoi ( argv[3] ) );
        }

    else {
        return 1;
    }

    return 0;
}

//go4/trb3
//source trb3login
//hldprint -filename -tdc 0x64ff -skip 1000 -num 100 | more


















