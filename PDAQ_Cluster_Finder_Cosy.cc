#include "PDAQ_Cluster_Finder_Cosy.h"

using namespace std;

std::vector<SciHit*> ex_Scint_pileup ( std::vector<SciHit*> B )
{
    //std::vector<SciHit*> B;
    std::vector<SciHit*> C;
    int scint_diff_limit = 200;

    bool scint =false;

    if ( B.size() >2 ) {
        for ( int aa=0; aa< B.size()-1; aa++ ) {
            //printf("%lf %lf \n",B[aa+1]->leadTime - B[aa]->leadTime, B[aa]->leadTime - B[aa-1]->leadTime);
            if ( aa ==0 ) {
                if ( B[aa+1]->leadTime - B[aa]->leadTime >=scint_diff_limit ) {}
                else {
                    B.erase ( B.begin() +aa+1 );
                    //cout<<B.size() <<endl;
                    --aa;
                }
            } else if ( aa >=1 ) {
                if ( B[aa+1]->leadTime - B[aa]->leadTime >=scint_diff_limit && B[aa]->leadTime - B[aa-1]->leadTime >=scint_diff_limit ) {}
                else {
                    B.erase ( B.begin() +aa+1 );
                    //cout<<B.size() <<endl;
                    --aa;
                }
            }
        }
    } else {
        if ( B.size() ==2 ) {
            if ( B[1]->leadTime - B[0]->leadTime >=scint_diff_limit ) {}
            else {
                B.erase ( B.begin() +1 );
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
                    if ( B[aa+1].layer==B[aa].layer && B[aa+1].straw==B[aa].straw && fabs ( fabs ( B[aa+1].leadTime ) - fabs ( B[aa].leadTime ) ) <=diff_limit ||  B[aa+1].leadTime  ==  B[aa].leadTime ) {
                        B.erase ( B.begin() +aa+1 );
                        --aa;
                    }

                } else if ( aa >=1 ) {
                    if ( B[aa+1].layer==B[aa].layer && B[aa+1].straw==B[aa].straw && fabs ( fabs ( B[aa+1].leadTime ) - fabs ( B[aa].leadTime ) ) <=diff_limit || B[aa+1].leadTime  ==  B[aa].leadTime ) {
                        B.erase ( B.begin() +aa+1 );
                        --aa;
                    }
                }
            }
        } else {
            if ( B.size() ==2 ) {
                if ( B[1].layer==B[0].layer && B[1].straw==B[0].straw && fabs ( fabs ( B[1].leadTime ) - fabs ( B[0].leadTime ) ) <=diff_limit  || B[1].leadTime  ==  B[0].leadTime ) {
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
    // int MAX_FT_TOTAL_LAYERS1 =0;

//     TH1F* h_p_layerDT[8];
//     TH2F* h_pL_dtvstot[8];
//     TH1F* h_p_TOT[8];
//     TH1F* h_pL_channel_mult[8];

    for ( int mh = 0; mh < 8; mh++ ) {
        h->h_pL_layerDT[mh] = new TH1F ( Form ( "Layer_%d_DTa", mh+1 ) , Form ( "Layer_%d_DTa", mh+1 ), 600, -100, 500 );
        h->h_pL_dtvstot[mh] = new TH2F ( Form ( "Layer_%d_DTvsTOTa",mh+1 ),Form ( "Layer_%d_DTvsTOTa;Drift Time;Time Over Threshold",mh+1 ), 560, -10, 550,550,0,550 );
        h->h_pL_TOT[mh] = new TH1F ( Form ( "Layer_%d_TOTa", mh+1 ) , Form ( "Layer_%d_TOTa", mh+1 ), 600, -100, 500 );
        h->h_pL_channel_mult[mh] = new TH2F ( Form ( "Layer_%d_multa", mh+1 ) , Form ( "Layer_%d_multa;Channel;Multiplicity", mh+1 ), 33, 0, 33, 20, 0, 20 );

        h->h_L_layerDT[mh] = new TH1F ( Form ( "Layer_%d_DTb", mh+1 ) , Form ( "Layer_%d_DTb", mh+1 ), 600, -100, 500 );
        h->h_L_dtvstot[mh] = new TH2F ( Form ( "Layer_%d_DTvsTOTb",mh+1 ),Form ( "Layer_%d_DTvsTOTb;Drift Time;Time Over Threshold",mh+1 ), 560, -10, 550,550,0,550 );
        h->h_L_TOT[mh] = new TH1F ( Form ( "Layer_%d_TOTb", mh+1 ) , Form ( "Layer_%d_TOTb", mh+1 ), 600, -100, 500 );
        h->h_L_channel_mult[mh] = new TH2F ( Form ( "Layer_%d_multb", mh+1 ) , Form ( "Layer_%d_multb;Channel;Multiplicity", mh+1 ), 33, 0, 33, 20, 0, 20 );

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

    SttHit* hitOnLayer[MAX_FT_TOTAL_LAYERS][5000];

    for ( Int_t i = 0; i < iev; i++ ) {
//cout<<"check1"<<endl;
        event_counter++;

        tree->GetEntry ( i );

        std::vector<SttHit> vec_stthits;
        std::vector<SttHit> vec_stthits2;

        bool scint = false;
        bool stt = false;

        if ( i == maxEvents ) {
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

        std::vector<SciHit*> B;
        B.clear();
        for ( int pn = 0; pn < SCI_CAL->sci_raw.totalNTDCHits; pn++ ) {
            SciHit* b = ( SciHit* ) SCI_CAL->sci_raw.adc_hits->ConstructedAt ( pn );
            B.push_back ( b );
            h->h_scnt_plup_cnts->Fill ( 1 );

        }

        /////////////////////////////// scint diff ///////////////////////////////////////////////
        if ( B.size() >=2 ) {
            for ( int a=0; a< B.size()-1; a++ ) {
                h->h_scint_timediffa->Fill ( B[a+1]->leadTime - B[a]->leadTime );
            }
        }

        //if (STT_CAL->SCI_CAL->sci_raw.totalNTDCHits>0){


        if ( B.size() >0 ) {


            vec_scihits =  ex_Scint_pileup ( B ) ;


            if ( vec_scihits.size() >=2 ) {
                for ( int a=0; a< vec_scihits.size()-1; a++ ) {
                    h->h_scint_timediffb->Fill ( vec_scihits[a+1]->leadTime - vec_scihits[a]->leadTime );
                }
            }


            h->h_scint_mult_a->Fill ( vec_scihits.size() );
            //printf("\nSIZE : %i\n",vec_scihits.size());
            for ( int sn = 0; sn < vec_scihits.size(); sn++ ) {
                SciHit* sh = vec_scihits[sn];
                h->h_scnt_plup_cnts->Fill ( 2 );
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
                    vec_stthits.clear();
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
                    for ( int jx=0; jx<vec_leadTime.size() ; jx++ ) {
                        if ( vec_leadTime[jx].leadTime - sh->leadTime >0 && vec_leadTime[jx].leadTime - sh->leadTime <500 ) {
                            h->h_pL_layerDT[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].leadTime - sh->leadTime );
                            h->h_pL_TOT[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].tot );
                            h->h_pL_dtvstot[vec_leadTime[jx].layer-1]->Fill ( vec_leadTime[jx].leadTime - sh->leadTime,vec_leadTime[jx].tot );
                            MultLayer[vec_leadTime[jx].layer-1]++;
                            ChHitsMult[vec_leadTime[jx].layer-1][vec_leadTime[jx].straw-1]++;
                            vec_leadTime_e.push_back ( vec_leadTime[jx] );
                            //printf ( "STRAW: %i  %i  %lf \n",vec_leadTime_d[jx].layer,vec_leadTime_d[jx].straw,vec_leadTime_d[jx].leadTime );
                        } else {
                            //vec_leadTime_d.erase ( vec_leadTime_d.begin() + jx );
                        }
                    }

                    int mult = 0;
                    for ( int aa=0; aa<8; aa++ ) {
                        if ( MultLayer[aa]> 0 ) {
                            mult++;
                        }
                        //printf("dt multiplicity%i \n",dtHitsMult[a]);
                        for ( int ba=0; ba<32; ba++ ) {
                            if ( ChHitsMult[aa][ba] >0 ) {
                                h->h_pL_channel_mult[aa]->Fill ( ba,ChHitsMult[aa][ba] );    //Fill(b,dtHitsMult[a][b]);
                                //cout<<aa+1<<"\t"<<ba+1<<"\t"<<ChHitsMult[aa][ba]<<endl;
                            }
                        }
                    }

                    //cout<<"First : "<<vec_leadTime_e.size() <<endl;
                    std::vector<SttHit> vec_channel_hit;
                    std::vector<std::vector<SttHit>> vec_layer_channel_hit;
                    std::vector<SttHit> vec_p_L5S9;
                    for ( Int_t jc = 0; jc < vec_leadTime_e.size(); jc++ ) {
                        //printf ( "ONE : %i  %i  %lf\n",vec_leadTime_e[jc].layer,vec_leadTime_e[jc].channel,vec_leadTime_e[jc].leadTime );
                        h->h_L_layerDT[vec_leadTime_e[jc].layer-1]->Fill ( vec_leadTime_e[jc].leadTime - sh->leadTime );
                        h->h_L_TOT[vec_leadTime_e[jc].layer-1]->Fill ( vec_leadTime_e[jc].tot );
                        h->h_L_dtvstot[vec_leadTime_e[jc].layer-1]->Fill ( vec_leadTime_e[jc].leadTime - sh->leadTime,vec_leadTime_e[jc].tot );
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


                    std::vector<SttHit> vec_leadTime_f;
                    vec_leadTime_f=ex_Stt_double ( vec_layer_channel_hit );
                    //cout<<"Next : "<<vec_leadTime_e.size() <<endl<<endl<<endl;;

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
                    //cout<<vec_leadTime_e.size() <<"\t"<<vec_leadTime_f.size() <<endl;
                    for ( Int_t jc = 0; jc < vec_leadTime_f.size(); jc++ ) {
                        // printf ( "TWO : %i  %i  %lf\n",vec_leadTime_f[jc].layer,vec_leadTime_f[jc].channel,vec_leadTime_f[jc].leadTime );
                        h->h_L_layerDT[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime );
                        h->h_L_TOT[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].tot );
                        h->h_L_dtvstot[vec_leadTime_f[jc].layer-1]->Fill ( vec_leadTime_f[jc].leadTime - sh->leadTime,vec_leadTime_f[jc].tot );

                        MultLayer2[vec_leadTime_f[jc].layer-1]++;
                        ChHitsMult2[vec_leadTime_f[jc].layer-1][vec_leadTime_f[jc].straw-1]++;
                        if ( vec_leadTime_f[jc].layer == 5 && vec_leadTime_f[jc].straw==9 ) {
                            vec_L5S9.push_back ( vec_leadTime_f[jc] );
                        }

                    }
                    int mult2 = 0;
                    for ( int ab=0; ab<8; ab++ ) {
                        if ( MultLayer2[ab]> 0 ) {
                            mult2++;
                        }
                        //printf("dt multiplicity%i \n",dtHitsMult[a]);
                        for ( int bb=0; bb<32; bb++ ) {
                            if ( ChHitsMult2[ab][bb] >0 ) {
                                h->h_L_channel_mult[ab]->Fill ( bb,ChHitsMult2[ab][bb] );    //Fill(b,dtHitsMult[a][b]);
                                //cout<<aa+1<<"\t"<<ba+1<<"\t"<<ChHitsMult[aa][ba]<<endl;
                            }
                        }
                    }

                    if ( vec_L5S9.size() >1 ) {
                        for ( int L=0; L<vec_L5S9.size()-1; L++ ) {
                            h->h_L5_S9LTdiff->Fill ( vec_L5S9[L+1].leadTime - vec_L5S9[L].leadTime );
                            //cout<<vec_L5S9.size() <<"\t"<<vec_L5S9[L+1].leadTime<<"\t"<<vec_L5S9[L].leadTime<<endl;
                        }
                    }
                    //cout<<endl<<endl;
                    //cout<<stt<<scint<<endl;
                    /*
                                        if ( stt == true && scint==true ) {

                                            if ( vec_leadTime.size() >= min_track_hits && vec_leadTime.size() < max_cluster_intake ) {
                                                vec_stthits.clear();

                                                for ( int e = 0; e < vec_leadTime.size(); e++ ) {
                    //                             printf ( "diff %lf\n",vec_leadTime[e].leadTime - sh->leadTime );
                                                    //printf ( "Track hit time: %lf Scint: %lf diff: %lf (%x %i %i %i| %i)\n",vec_leadTime[e].leadTime,sh->leadTime,vec_leadTime[e].leadTime - sh->leadTime,vec_leadTime[e].tdcid,vec_leadTime[e].layer,vec_leadTime[e].channel,vec_leadTime[e].straw,sh->channel );

                                                    if ( ( vec_leadTime[e].leadTime - sh->leadTime >0 ) && ( vec_leadTime[e].leadTime - sh->leadTime <220 ) ) {

                                                        //printf ( "Track hit time: %lf Scint: %lf diff: %lf (%x %i %i %i| %i)\n",vec_leadTime[e].leadTime,sh->leadTime,vec_leadTime[e].leadTime - sh->leadTime,vec_leadTime[e].tdcid,vec_leadTime[e].layer,vec_leadTime[e].channel,vec_leadTime[e].straw,sh->channel );
                                                        SttHit& s = vec_leadTime.at ( e );
                                                        s.drifttime = ( vec_leadTime[e].leadTime - sh->leadTime );
                                                        vec_stthits.push_back ( s );
                                                    }
                                                }

                                            }
                                        }
                                        //cout<<"\n\n";
                                        h->h_cluster_size->Fill ( vec_stthits.size() );


                                        if ( vec_stthits.size() >= min_track_hits && vec_stthits.size() < max_cluster_intake ) {
                                            //PDAQ_Event_Finder ( vec_stthits, i, PDAQ_tree, stt_event, ftGeomPar, SCI_CAL, h );
                                            h->h_hitmultiplicity4->Fill ( vec_stthits.size() );
                                        }*/

                    /* if ( vec_leadTime.size() >= min_track_hits && vec_leadTime.size() < max_cluster_intake ) {
                         const int minNumber = min_track_hits;
                         const int maxDifference = max_lead_time_diff - 1;
                         int currentNumber = 0;
                         int currentSum = 0;
                         int numGroups = 0;
                         double last = vec_leadTime[0]->leadTime - maxDifference - 1;
                         vec_stthits.clear();

                         for ( int e = 0; e < vec_leadTime.size(); e++ ) {
                             if ( vec_leadTime[e]->leadTime - last <= maxDifference ) { // Continue with current group
                                 currentNumber++;
                                 SttHit* h = vec_leadTime.at ( e );
                                 vec_stthits.push_back ( h );

                                 }
                             else {   // Finish previous group and start anew
                                 if ( currentNumber >= minNumber ) { // Previous was a valid group
                                     numGroups++;
                                     PDAQ_Event_Finder ( vec_stthits, i, PDAQ_tree, stt_event, ftGeomPar, SCI_CAL, h );
                                     h->h_cluster_size->Fill ( currentNumber );
                                     h->h_hitmultiplicity3->Fill ( currentNumber );
                                     }
                                 vec_stthits.clear();
                                 SttHit* h = vec_leadTime.at ( e );
                                 vec_stthits.push_back ( h ); // Start afresh
                                 currentNumber = 1;
                                 }
                             last = vec_leadTime[e]->leadTime;
                             }
                         // Deal with leftovers at the end
                         if ( currentNumber >= minNumber ) { // Previous was a valid group
                             numGroups++;
                             PDAQ_Event_Finder ( vec_stthits, i, PDAQ_tree, stt_event, ftGeomPar, SCI_CAL, h );
                             h->h_cluster_size->Fill ( currentNumber );
                             h->h_hitmultiplicity3->Fill ( currentNumber );
                             }
                         }*/
                }
            }
        }

    } // End of loop over events
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
    h->h_scint_timediffa->Write();
    h->h_scint_timediffb->Write();
    h->h_L5_S9LTdiff->Write();
    h->h_pL5_S9LTdiff->Write();

    for ( int hh = 0; hh < 8; hh++ ) {
        h->h_pL_layerDT[hh]->Write();
        h->h_pL_dtvstot[hh]->Write();
        h->h_pL_TOT[hh]->Write();
        h->h_pL_channel_mult[hh]->Write();

        h->h_L_layerDT[hh]->Write();
        h->h_L_dtvstot[hh]->Write();
        h->h_L_TOT[hh]->Write();
        h->h_L_channel_mult[hh]->Write();

    }


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










