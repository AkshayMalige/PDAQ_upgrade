#include "PDAQ_Cluster_Finder_Cosy.h"
#include "math.h"

using namespace std;

long fact ( int n )
{
    unsigned long long factorial = 1;
    for ( int i = 1; i <=n; ++i ) {
        factorial *= i;
    }
    return factorial;
}

std::vector<float> eff_job ( double ax )
{

    std::vector<float> A;
    double probability =0;
    double sum =0;
    for ( int a=0; a<17; a++ ) {
        probability = ( ( fact ( 16 )  / ( fact ( a ) * fact ( 16-a ) ) )    * pow ( ax,a ) * pow ( ( 1-ax ),16-a ) );
        sum += probability;
        cout<<a<<"\t"<<probability*100<<endl;
        A.push_back ( probability );


    }

    cout<<"Sum "<<sum<<endl;
    return A;
}

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
                    if ( ( l2->leadTime - l->leadTime ) >=scint_diff_limit && ( l3->leadTime - l2->leadTime ) >=scint_diff_limit ) {
                        B.push_back ( l2 );
                        B.push_back ( l3 );
                    } else if ( ( l3->leadTime - l2->leadTime ) >=scint_diff_limit ) {
                        B.push_back ( l3 );
                    }
                } else {
                    if ( ( l2->leadTime - l->leadTime ) >=scint_diff_limit && ( l3->leadTime - l2->leadTime ) >=scint_diff_limit ) {
                        B.push_back ( l2 );
                    }
                }

            } else if ( a==0 && l3->leadTime-l2->leadTime >=scint_diff_limit ) {
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
    int diff_limit = 500;
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

    TCanvas * LtDiff; //Canvas for ToT
    LtDiff=new TCanvas ( "LtDiff_Canvas","LtDiff_Canvas" );
    LtDiff->Divide ( 4,2 );

    TCanvas * LtDiffLow; //Canvas for ToT
    LtDiffLow=new TCanvas ( "LtDiffLow_Canvas","LtDiffLow_Canvas" );
    LtDiffLow->Divide ( 4,2 );

    TCanvas * Eff; //Canvas for ToT
    Eff=new TCanvas ( "Eff_Canvas","Eff_Canvas" );
    Eff->Divide ( 4,2 );

    TCanvas * ChDiff; //Canvas for ToT
    ChDiff=new TCanvas ( "ChDiff_Canvas","ChDiff_Canvas" );
    ChDiff->Divide ( 4,2 );

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
    int final_counter = 0;

    float max_lead_time_diff = 0;
    int min_track_hits = 0;
    bool Scint = false;
    float refTime = 0;
    double total_particles=0;

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

    //int High_strw[8]= {9,12,6,9,9,12,6,9};

    for ( int mh = 0; mh < 8; mh++ ) {
        h->h_pL_layerDT[mh] = new TH1F ( Form ( "Layer_%d_DTa", mh+1 ) , Form ( "Layer_%d_DTa", mh+1 ), 600, -100, 500 );
        h->h_pL_dtvstot[mh] = new TH2F ( Form ( "Layer_%d_DTvsTOTa",mh+1 ),Form ( "Layer_%d_DTvsTOTa;Drift Time;Time Over Threshold",mh+1 ), 500, 0, 500,700,0,700 );
        h->h_pL_TOT[mh] = new TH1F ( Form ( "Layer_%d_TOTa", mh+1 ) , Form ( "Layer_%d_TOTa", mh+1 ), 600, -100, 500 );
        h->h_pL_channel_mult[mh] = new TH2F ( Form ( "Layer_%d_multa", mh+1 ) , Form ( "Layer_%d_multa;Channel;Multiplicity", mh+1 ), 33, 0, 33, 20, 0, 20 );

        h->h_L_layerDT[mh] = new TH1F ( Form ( "Layer_%d_DTb", mh+1 ) , Form ( "Layer_%d_DTb", mh+1 ), 600, -100, 500 );
        h->h_L_dtvstot[mh] = new TH2F ( Form ( "Layer_%d_DTvsTOTb",mh+1 ),Form ( "Layer_%d_DTvsTOTb;Drift Time;Time Over Threshold",mh+1 ), 500, 0, 500,700,0,700 );
        h->h_L_TOT[mh] = new TH1F ( Form ( "Layer_%d_TOTb", mh+1 ) , Form ( "Layer_%d_TOTb", mh+1 ), 600, -100, 500 );
        h->h_L_channel_mult[mh] = new TH2F ( Form ( "Layer_%d_multb", mh+1 ) , Form ( "Layer_%d_multb;Channel;Multiplicity", mh+1 ), 33, 0, 33, 20, 0, 20 );

        h->h_pChDiff[mh] = new TH1F ( Form ( "Layer_%d_pChDiff", mh+1 ) , Form ( "Layer_%d_pChDiff", mh+1 ), 5, 0, 5 );
        h->h_straw[mh] = new TH1F ( Form ( "Layer_%d_h_straw", mh+1 ) , Form ( "Layer_%d_h_straw", mh+1 ), 32, 0, 32 );
    }

    for ( int ch=0; ch<256; ch++ ) {
        h->h_Ch_Dt[ch] = new TH1F ( Form ( "Ch_%d_Dt", ch+1 ) , Form ( "Ch_%d_Dt", ch+1 ),  500, 1,501 );
        h->h_pLT_Diff[ch] = new TH1F ( Form ( "Ch_%d_LT_diff", ch+1 ) , Form ( "Ch_%d_LT_diff", ch+1 ), 500, 1, 501 );
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
    int scint_event =0;

    SttHit* hitOnLayer[MAX_FT_TOTAL_LAYERS][5000];

    for ( Int_t i = 0; i < iev; i++ ) {
//cout<<"check1"<<endl;
        event_counter++;

        tree->GetEntry ( i );

        std::vector<SttHit> vec_stthits;


        bool scint = false;
        bool stt = false;

        if ( scint_event >= maxEvents ) {
            break;
        }

        if ( i % 1000 == 0 ) {
            cout << "entry no. " << i << endl;
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

        //if (STT_CAL->SCI_CAL->sci_raw.totalNTDCHits>0){


        if ( SCI_CAL->sci_raw.totalNTDCHits >0 ) {


            vec_scihits =  ex_Scint_pileup ( SCI_CAL ) ;


            if ( vec_scihits.size() >=2 ) {
                for ( int a=0; a< vec_scihits.size()-1; a++ ) {
                    //if
                    h->h_scint_timediffb->Fill ( vec_scihits[a+1]->leadTime - vec_scihits[a]->leadTime );
                }
            }


            h->h_scint_mult_a->Fill ( vec_scihits.size() );
            //printf("\nSIZE : %i\n",vec_scihits.size());
            for ( int sn = 0; sn < vec_scihits.size(); sn++ ) {
                scint_event++;
                SciHit* sh = vec_scihits[sn];
                h->h_scnt_plup_cnts->Fill ( 2 );
                h->h_ShitBlock->Fill ( 1 );
                if ( sh->tdcid ==0x6500 && sh->channel==1 && vec_scihits.size() >1 ) {
                    scint = true;
                } else {
                    scint = false;
                }
//DT = straw - scint

                h->h_hitmultiplicity1->Fill ( STT_CAL->stt_cal.total_cal_NTDCHits );
                //printf("\nSCINT : %lf\n",sh->leadTime);

                for ( int n = 0; n < STT_CAL->stt_cal.total_cal_NTDCHits; n++ ) {
                    SttHit* cal_hit = ( SttHit* ) STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt ( n ); // retrieve particular hit
                    //printf("\nStraw : %lf",cal_hit->leadTime);
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

//cout<<"\n\n"<<endl;

                bool good_layers = true;
                for ( int c = 0; c < MAX_FT_TOTAL_LAYERS; c++ ) {
                    if ( hitMultOnLayer[c] != 1 ) {
                        good_layers = false;
                        break;
                    }
                }
                int layerCounter = 0;
// cout<<"check5"<<endl;
                for ( int m = 0; m < MAX_FT_TOTAL_LAYERS; m++ ) {
                    if ( hitMultOnLayer[m] > 0 ) {
                        layerCounter++;
                    }
                }
                h->h_0LMultiplicity->Fill ( layerCounter );
//                cout<<"check6"<<endl;

                if ( layerCounter >0/*== MAX_FT_TOTAL_LAYERS*/ ) {

                    Layer_eq_4_counter++;
                    stt = true;

                    std::vector<SttHit> vec_leadTime;
                    int filtercnt = 0;

                    for ( int l = 0; l < MAX_FT_TOTAL_LAYERS; l++ ) {
                        for ( int r = 0; r < hitMultOnLayer[l]; r++ )

                        {
                            vec_leadTime.push_back ( *hitOnLayer[l][r] );
                            h->h_tot1->Fill ( hitOnLayer[l][r]->tot );
                            //cout<<"Initial : "<<hitOnLayer[l][r]->layer<<"\t"<<hitOnLayer[l][r]->channel<<"\t"<<hitOnLayer[l][r]->leadTime<<endl;
                            //printf("%lf  %lf  %lf \n",sh->leadTime,hitOnLayer[l][r]->leadTime,hitOnLayer[l][r]->leadTime - sh->leadTime);

                        }
                    }

                    std::sort ( vec_leadTime.begin(), vec_leadTime.end(), f_sttHitCompareLeadTime );
                    double doublehitdiff = 0;

////////////////////////////////////////////////////////// Multiple entry rejection ///////////////////
                    if ( vec_leadTime.size() >1 ) {
                        for ( Int_t je = 0; je < vec_leadTime.size()-1; je++ ) {
                            if ( ( vec_leadTime[je + 1].layer ) == ( vec_leadTime[je].layer ) && ( vec_leadTime[je + 1].channel == vec_leadTime[je].channel ) && ( vec_leadTime[je + 1].leadTime == vec_leadTime[je].leadTime ) ) {

                                vec_leadTime.erase ( vec_leadTime.begin() +je+1 );
                                --je;
                            }
                        }
                    }

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
                    for ( int jx=0; jx < vec_leadTime.size() ; jx++ ) {

                        if ( vec_leadTime[jx].leadTime - sh->leadTime >-100 && vec_leadTime[jx].leadTime - sh->leadTime <500 ) {
                            h->h_pL_layerDT[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].leadTime - sh->leadTime );
                            h->h_pL_TOT[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].tot );
                            h->h_pL_dtvstot[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].leadTime - sh->leadTime,vec_leadTime[jx].tot );
                            h->h_p_drifttime->Fill ( vec_leadTime[jx].leadTime - sh->leadTime );
                            h->h_p_tot->Fill ( vec_leadTime[jx].tot );
                            MultLayer[vec_leadTime[jx].layer-1]++;
                            ChHitsMult[vec_leadTime[jx].layer-1][vec_leadTime[jx].straw-1]++;
                            vec_leadTime_e.push_back ( vec_leadTime[jx] );
                            h->h_hitBlock->Fill ( 2 );
                            //printf ( "STRAW: %i  %i  %lf \n",vec_leadTime_d[jx].layer,vec_leadTime_d[jx].straw,vec_leadTime_d[jx].leadTime );
                        } else {
                            //vec_leadTime_d.erase ( vec_leadTime_d.begin() + jx );
                        }
                    }

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
                    //cout<<vec_leadTime_e.size() <<"\t"<<vec_leadTime_f.size() <<endl;
                    for ( Int_t jc = 0; jc < vec_leadTime_f.size(); jc++ ) {
                        // printf ( "TWO : %i  %i  %lf\n",vec_leadTime_f[jc].layer,vec_leadTime_f[jc].channel,vec_leadTime_f[jc].leadTime );
                        vec_leadTime_f[jc].drifttime = vec_leadTime_f[jc].leadTime - sh->leadTime;
                        vec_stthits.push_back ( vec_leadTime_f[jc] );
                        h->h_L_layerDT[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime );
                        h->h_L_TOT[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].tot );
                        h->h_L_dtvstot[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime,vec_leadTime_f[jc].tot );
                        h->h_drifttime->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime );
                        h->h_tot->Fill ( vec_leadTime_f[jc].tot );
                        h->h_hitBlock->Fill ( 3 );

                        sq_ch = ( ( vec_leadTime_f[jc].layer-1 ) * 32 ) +vec_leadTime_f[jc].straw-1;
                        h->h_Ch_Dt[sq_ch]->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime );

                        MultLayer2[vec_leadTime_f[jc].layer-1]++;
                        ChHitsMult2[vec_leadTime_f[jc].layer-1][vec_leadTime_f[jc].straw-1]++;
                        if ( vec_leadTime_f[jc].layer == 5 && vec_leadTime_f[jc].straw==9 ) {
                            vec_L5S9.push_back ( vec_leadTime_f[jc] );
                        }
                        h->h_sq_ch->Fill ( sq_ch );
                        h->h_straw[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].straw );

                    }
                    int mult2 = 0;
                    int layer_eff_count2 =0;
                    for ( int ab=0; ab<8; ab++ ) {
                        if ( MultLayer2[ab]> 0 ) {
                            mult2++;
                        }
                        h->h_LayerMult->Fill ( ab,MultLayer2[ab] );
                        //printf("dt multiplicity%i \n",MultLayer2[ab]);
                        for ( int bb=0; bb<32; bb++ ) {
                            if ( ChHitsMult2[ab][bb] >0 ) {
                                h->h_L_channel_mult[ab]->Fill ( bb,ChHitsMult2[ab][bb] );    //Fill(b,dtHitsMult[a][b]);
                                layer_eff_count2++;
                            }
                        }
                    }
                    if ( MultLayer2[0]>0 && MultLayer2[7]>0 ) {
                        h->h_Layer_eff->Fill ( layer_eff_count2 );
                        //cout<<"ONE "<<layer_eff_count<<endl;
                    }
                    h->h_LMultiplicity->Fill ( mult2 );

                    if ( vec_L5S9.size() >1 ) {
                        for ( int L=0; L<vec_L5S9.size()-1; L++ ) {
                            h->h_L5_S9LTdiff->Fill ( vec_L5S9[L+1].leadTime - vec_L5S9[L].leadTime );
                            //cout<<vec_L5S9.size() <<"\t"<<vec_L5S9[L+1].leadTime<<"\t"<<vec_L5S9[L].leadTime<<endl;
                        }
                    }

                    int sqq_ch=0;
                    bool ch_hi =false;
                    bool diff_1=false;
                    bool diff_2=false;
                    for ( int hc=0; hc<8; hc++ ) {
                        for ( int fc=0; fc<vec_stthits.size(); fc++ ) {
                            sqq_ch = ( ( vec_stthits[fc].layer-1 ) * 32 ) +vec_stthits[fc].straw-1;
                            if ( sqq_ch == High_ch[hc]+1 ) {
                                ch_hi = true;
                            }
                            if ( sqq_ch == High_ch[hc]+2 || sqq_ch == High_ch[hc] ) {
                                diff_1 = true;
                            }
                            if ( sqq_ch == High_ch[hc]+3 || sqq_ch == High_ch[hc]-1 ) {
                                diff_2 = true;
                            }
                        }
                        //cout<<ch_hi<<"\t"<<diff_1<<"\t"<<diff_2<<endl;
                        if ( ch_hi==true && diff_2==true && diff_1==true ) {
                            h->h_pChDiff[hc]->Fill ( 4 );
                        } else if ( ch_hi==true && diff_2==false && diff_1==true ) {
                            h->h_pChDiff[hc]->Fill ( 2 );
                        } else if ( ch_hi==true && diff_1==false && diff_2==true ) {
                            h->h_pChDiff[hc]->Fill ( 3 );
                        } else if ( ch_hi==true  && diff_1==false && diff_2==false ) {
                            h->h_pChDiff[hc]->Fill ( 1 );
                        }

                        ch_hi=false;
                        diff_1 = false;
                        diff_2 = false;


                    }


///////////////////////////////////////////////////////////////////////////////////
                    std::vector<SttHit> vec_corridor;
                    vec_corridor.clear();
                    int MultLayer3[8];
                    int ChHitsMult3[8][32];

                    for ( int lb = 0; lb < 8; lb++ ) {
                        MultLayer3[lb] = 0;
                        for ( int mb =0; mb< 32; mb++ ) {
                            ChHitsMult3[lb][mb] = 0;
                        }
                    }

                    for ( int cr=0; cr<vec_stthits.size(); cr++ ) {
                        if ( ( vec_stthits[cr].layer==1 || vec_stthits[cr].layer==4 || vec_stthits[cr].layer==5 || vec_stthits[cr].layer==8 ) && ( vec_stthits[cr].straw>=7 && vec_stthits[cr].straw<=11 ) ) {
                            vec_corridor.push_back ( vec_stthits[cr] );
                        } else if ( ( vec_stthits[cr].layer==2 || vec_stthits[cr].layer==6 ) && ( vec_stthits[cr].straw>=10 && vec_stthits[cr].straw<=14 ) ) {
                            vec_corridor.push_back ( vec_stthits[cr] );
                        } else if ( ( vec_stthits[cr].layer==3 || vec_stthits[cr].layer== 7 ) && ( vec_stthits[cr].straw>=4 && vec_stthits[cr].straw<=8 ) ) {
                            vec_corridor.push_back ( vec_stthits[cr] );
                        }
                    }

                    for ( int cr=0; cr<vec_corridor.size(); cr++ ) {
                        MultLayer3[vec_corridor[cr].layer-1]++;
                        ChHitsMult3[vec_corridor[cr].layer-1][vec_corridor[cr].straw-1]++;

                    }

                    int mult3 = 0;
                    int layer_eff_count3 =0;

                    for ( int kb=0; kb<8; kb++ ) {
                        if ( MultLayer3[kb]> 0 ) {
                            mult3++;
                            //cout<<"mult 3  "<<kb<<"\t"<<MultLayer3[kb]<<endl;
                        }

                        for ( int bk=0; bk<32; bk++ ) {
                            if ( ChHitsMult3[kb][bk] >0 ) {
                                layer_eff_count3++;
                            }
                        }
                    }
                    if ( MultLayer3[0]>0 && MultLayer3[7]>0 ) {
                        h->h_Layer_eff2->Fill ( layer_eff_count3 );
                        //cout<<"Two "<<layer_eff_count3<<endl;
                    }


                    //cout<<"***************"<<endl;

                    // for (int sw=0; sw<vec_scihits.size(); sw++){
                    //cout<<vec_stthits[sw].leadTime<<"\t"<<vec_stthits[sw].layer<<"\t"<<vec_stthits[sw].straw<<endl;
                    //}


                    if ( vec_stthits.size() >= min_track_hits && vec_stthits.size() <= max_cluster_intake && mult2==8 ) {
                        // PDAQ_Event_Finder ( vec_stthits, i, PDAQ_tree, stt_event, ftGeomPar, SCI_CAL, h );

                    }

                }
            }
        }

    } // End of loop over events

    double maxbin =0;
    double DTmin =0;
    double dt_at_10=0;
    Int_t maximum = h->h_drifttime->GetBinContent ( h->h_drifttime->GetMaximumBin() );
    DTmin = h->h_drifttime->FindFirstBinAbove ( maximum/10,1 );
    dt_at_10 = h->h_drifttime->GetXaxis()->GetBinCenter ( h->h_drifttime->GetXaxis()->FindBin ( DTmin ) );

//     cout<< dt_at_10 - 100<<endl;
//
//     cout<<"New  "<<maximum<<"\t"<<DTmin <<"\t"<<dt_at_10<<endl;

    h_STT_Hit_Diff->Write();
    h->h_cluster_size->Write();
    h->h_LTvsLayer0->Write();
    h->h_LTvsLayer1->Write();

    h->HitMultiplicity->cd();
    h->h_hitmultiplicity0->Draw();
    h->h_hitmultiplicity1->Draw ( "same" );
    h->h_hitmultiplicity2->Draw ( "same" );
    h->h_hitmultiplicity3->Draw ( "same" );
    h->h_hitmultiplicity4->Draw ( "same" );
    h->h_hitmultiplicity0->SetLineColor ( kBlack );
    h->h_hitmultiplicity1->SetLineColor ( kMagenta );
    h->h_hitmultiplicity2->SetLineColor ( kBlue );
    h->h_hitmultiplicity3->SetLineColor ( kGreen );
    h->h_hitmultiplicity4->SetLineColor ( kRed );
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

    h->TimeOverThreshold->cd();
    h->h_tot0->Draw();
    h->h_tot1->Draw ( "same" );
    h->h_tot2->Draw ( "same" );
    h->h_tot3->Draw ( "same" );
    h->h_tot4->Draw ( "same" );
    h->h_tot0->SetLineColor ( kRed );
    h->h_tot1->SetLineColor ( kMagenta );
    h->h_tot2->SetLineColor ( kBlue );
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

    h->h_drifttime->Write();
    h->h_drifttime->GetXaxis()->SetLabelSize ( 0.050 );

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
    h->h_tot->Write();
    h->h_scint_timediffa->Write();
    h->h_scint_timediffb->Write();
    h->h_L5_S9LTdiff->Write();
    h->h_pL5_S9LTdiff->Write();
    h->h_pLayerMult->Write();
    h->h_LayerMult->Write();
    h->h_hitBlock->Write();
    h->h_ShitBlock->Write();
    h->h_ShitBlock->Write();
    h->h_0LMultiplicity->Write();
    h->h_pLMultiplicity->Write();
    h->h_LMultiplicity->Write();

    Double_t norm = h->h_Layer_eff->GetEntries();
    h->h_Layer_eff->Scale ( 1/norm );
    h->h_Layer_eff->Write();

    Double_t norm2 = h->h_Layer_eff2->GetEntries();
    h->h_Layer_eff2->Scale ( 1/norm2 );
    h->h_Layer_eff2->Write();

    h->h_pLayer_eff->Write();
    h->h_sq_ch->Write();

    double from_data[30];
    double from_corridor[30];
    double index[30];
    std::vector<float> B;
    std::vector<float> C;
    B = eff_job ( 0.90 );
    C = eff_job ( 0.95 );

    float from_90per[B.size()];
    float from_95per[C.size()];
    float from_90index[B.size()];
    for ( int jh=0; jh<B.size(); jh++ ) {
        from_90per[jh] = B[jh];
        from_95per[jh] = C[jh];
        from_90index[jh] = jh;
    }


    for ( int r=0; r<30; r++ ) {
        from_data[r]= h->h_Layer_eff->GetBinContent ( r );
        from_corridor[r]= h->h_Layer_eff2->GetBinContent ( r );
        index[r] =r-1;
    }

    Eff->cd();
    TGraph* gEffeciency = new TGraph ( 30, index , from_data );
    TGraph* gEffeciency2 = new TGraph ( 30, index , from_corridor );
    TGraph* gFrom90Data = new TGraph ( B.size(), from_90index , from_90per );
    TGraph* gFrom95Data = new TGraph ( C.size(), from_90index , from_95per );
    
    gEffeciency->Write();
    gEffeciency->Draw ( "AB" );
    gEffeciency->SetFillColor ( kRed );
    gEffeciency->SetFillStyle ( 3007 );
    gEffeciency->SetMarkerStyle ( 3 );
    gEffeciency->SetMarkerColor ( kRed );
    

    gEffeciency2->Write();
    gEffeciency2->Draw ( "same,P" );
    //gEffeciency2->SetFillColor ( kBlack );
    //gEffeciency2->SetFillStyle ( 3006 );
    gEffeciency2->SetMarkerStyle ( 3 );
    gEffeciency2->SetMarkerColor ( kBlack );

    //gFrom90Data->Draw("AP");

    gFrom90Data->SetMarkerStyle ( 3 );
    gFrom90Data->SetMarkerColor ( kBlue );
    gFrom90Data->Draw ( "same,P" );
    gFrom90Data->Write();

    gFrom95Data->SetMarkerStyle ( 3 );
    gFrom95Data->SetMarkerColor ( kGreen );
    gFrom95Data->Draw ( "same,P" );
    gFrom95Data->Write();



    for ( int hh = 0; hh < 8; hh++ ) {
        h->h_pL_layerDT[hh]->Write();
        h->h_pL_layerDT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pL_layerDT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_pL_dtvstot[hh]->Write();
        h->h_pL_layerDT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pL_layerDT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_pL_TOT[hh]->Write();
        h->h_pL_TOT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pL_TOT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_pL_channel_mult[hh]->Write();
        h->h_pL_TOT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pL_TOT[hh]->GetYaxis()->SetLabelSize ( 0.045 );

        h->h_L_layerDT[hh]->Write();
        h->h_L_layerDT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_L_layerDT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_L_dtvstot[hh]->Write();
        h->h_L_dtvstot[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_L_dtvstot[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_L_TOT[hh]->Write();
        h->h_L_TOT[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_L_TOT[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_L_channel_mult[hh]->Write();
        h->h_L_channel_mult[hh]->GetXaxis()->SetLabelSize ( 0.045 );

        Double_t norm1 = h->h_pChDiff[hh]->GetEntries();
        h->h_pChDiff[hh]->Scale ( 1/norm1 );
        h->h_pChDiff[hh]->Write();
        h->h_pChDiff[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pChDiff[hh]->GetYaxis()->SetLabelSize ( 0.045 );
        h->h_straw[hh]->Write();
        h->h_straw[hh]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_straw[hh]->GetYaxis()->SetLabelSize ( 0.045 );

        DT->cd ( hh+1 );
        h->h_pL_layerDT[hh]->Draw();
        h->h_pL_layerDT[hh]->SetLineWidth ( 2 );
        h->h_L_layerDT[hh]->Draw ( "same" );
        h->h_L_layerDT[hh]->SetLineWidth ( 2 );
        h->h_L_layerDT[hh]->SetLineColor ( kRed );

        TOT->cd ( hh+1 );
        h->h_pL_TOT[hh]->Draw();
        h->h_pL_TOT[hh]->SetLineWidth ( 2 );
        h->h_L_TOT[hh]->Draw ( "same" );
        h->h_L_TOT[hh]->SetLineWidth ( 2 );
        h->h_L_TOT[hh]->SetLineColor ( kRed );

        pDTvsTOT->cd ( hh+1 );
        h->h_pL_dtvstot[hh]->Draw ( "colz" );

        DTvsTOT->cd ( hh+1 );
        h->h_L_dtvstot[hh]->Draw ( "colz" );

        pChMult->cd ( hh+1 );
        h->h_pL_channel_mult[hh]->Draw ( "colz" );

        ChMult->cd ( hh+1 );
        h->h_L_channel_mult[hh]->Draw ( "colz" );

        LtDiff->cd ( hh+1 );
        h->h_pLT_Diff[High_ch[hh]]->Draw ();
        gPad->SetLogy();

        LtDiffLow->cd ( hh+1 );
        h->h_pLT_Diff[Low_ch[hh]]->Draw ();
        gPad->SetLogy();

        ChDiff->cd ( hh+1 );
        h->h_pChDiff[hh]->Draw();
        //gPad->SetLogy();

    }



//     cout<< dt_at_10 - 100<<endl;
//
//     cout<<"New  "<<maximum<<"\t"<<DTmin <<"\t"<<dt_at_10<<endl;




    std::vector<double> vec_DT_start;
    cout<<norm<<"\t"<<total_particles<<endl;
    float sc=1/total_particles;

    for ( int cha=0; cha<256; cha++ ) {
        //h->h_Ch_Dt[cha]->Scale ( 1/total_particles );
        h->h_Ch_Dt[cha]->Write();
        h->h_Ch_Dt[cha]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_Ch_Dt[cha]->GetYaxis()->SetLabelSize ( 0.045 );

        h->h_pLT_Diff[cha]->Scale ( 1/total_particles );
        h->h_pLT_Diff[cha]->Write();
        h->h_pLT_Diff[cha]->GetXaxis()->SetLabelSize ( 0.045 );
        h->h_pLT_Diff[cha]->GetYaxis()->SetLabelSize ( 0.045 );
    }

//     for(int dot=0; dot<vec_DT_start.size(); dot++){
//       printf("Channel : %i   DT_crr: %d \n",dot,vec_DT_start[dot]-dt_at_10);
//     }



    DT->Write();
    TOT->Write();
    pDTvsTOT->Write();
    DTvsTOT->Write();
    pChMult->Write();
    ChMult->Write();
    LtDiff->Write();
    LtDiffLow->Write();
    ChDiff->Write();
    Eff->Write();


    PDAQ_tree->Write();
    //Ttree->Close();
    //inFile.Close();
    printf ( "Total Hits processed : %f Repeated Hits in an Event : %f\n\n", All_repeat, repeat );
    printf ( "In_File: %s 	Out_File:  %s\n", intree, outtree );
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
























