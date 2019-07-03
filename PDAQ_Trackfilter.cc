#include "PDAQ_Trackfilter.h"

bool f_sttHitCompareLeadTime ( SttHit a, SttHit b )
{
    return ( a.leadTime < b.leadTime );
}

bool f_sttHitCompareCell ( SttHit a, SttHit b )
{
    return ( a.straw < b.straw );
}

VecSttHit GetPairs ( VecSttHit vec_get_pairs ,histograms* h )
{
    VecSttHit vec_fpair_clu;
    std::sort ( vec_get_pairs.begin(), vec_get_pairs.end(), f_sttHitCompareCell );

    for ( int sa = 0; sa < vec_get_pairs.size(); sa++ ) {
        if ( sa < vec_get_pairs.size() - 1 ) {
            if ( fabs ( vec_get_pairs[sa].straw - vec_get_pairs[sa + 1].straw ) ==
                    1 ) {
                vec_fpair_clu.push_back ( vec_get_pairs[sa] );
            }

            else if ( sa > 0 ) {
                if ( fabs ( vec_get_pairs[sa - 1].straw -
                            vec_get_pairs[sa].straw ) == 1 ) {
                    vec_fpair_clu.push_back ( vec_get_pairs[sa] );
                } else {
                  //  h->h_Cross_TOT[2]->Fill ( vec_get_pairs[sa].tot );
                   // h->h_Cross_DT[2]->Fill ( vec_get_pairs[sa].tot );
                }
            }
        } else if ( sa == vec_get_pairs.size() - 1 ) {
            if ( fabs ( vec_get_pairs[sa - 1].straw - vec_get_pairs[sa].straw ) ==
                    1 ) {
                vec_fpair_clu.push_back ( vec_get_pairs[sa] );
            } else {
              //  h->h_Cross_TOT[2]->Fill ( vec_get_pairs[sa].tot );
              //  h->h_Cross_DT[2]->Fill ( vec_get_pairs[sa].tot );
            }
        }
    }
    return vec_fpair_clu;
}

std::vector<VecSttHit> clusterfinder ( VecSttHit vec_flayer );


void add_tuples ( const std::vector<std::vector<VecSttHit>>& vectors,
                  std::size_t pos, std::vector<VecSttHit> prefix,
                  std::vector<std::vector<VecSttHit>>& result )
{
    if ( pos == vectors.size() ) {
        result.push_back ( prefix );
    } else if ( vectors[pos].empty() ) {
        add_tuples ( vectors, pos + 1, prefix,
                     result ); // note: skip empty vectors
    } else {
        for ( VecSttHit v : vectors[pos] ) {
            prefix.push_back ( v );
            add_tuples ( vectors, pos + 1, prefix, result );
            prefix.pop_back();
        }
    }
}

std::vector<VecSttHit> bestPair ( std::vector<VecSttHit> vec_pair ,histograms* h )
{
//     for ( int b=0; b<vec_pair.size(); b++ ) {
//         printf ( "Pair %i\n",b+1 );
//         for ( int c=0; c<vec_pair[b].size(); c++ ) {
//             printf ( "%i %f %f %f\n",vec_pair.at ( b ).at ( c ).straw,vec_pair.at ( b ).at ( c ).tot,vec_pair.at ( b ).at ( c ).drifttime,vec_pair.at ( b ).at ( c ).leadTime );
//         }
//         printf ( "\n\n" );
//     }

    std::vector<int> index;
    int No_pairs = vec_pair.size();
    for ( int a =0; a<No_pairs-1; a++ ) {
        if ( vec_pair.at ( a ).at ( 1 ).straw == vec_pair.at ( a+1 ).at ( 0 ).straw ) {
            h->h_Cross_TOT[1]->Fill ( vec_pair.at ( a ).at ( 0 ).tot );
            h->h_Cross_DT[1]->Fill ( vec_pair.at ( a ).at ( 0 ).drifttime );

            h->h_Cross_TOT[2]->Fill ( vec_pair.at ( a ).at ( 1 ).tot );
            h->h_Cross_DT[2]->Fill ( vec_pair.at ( a ).at ( 1 ).drifttime );

            h->h_Cross_TOT[3]->Fill ( vec_pair.at ( a +1).at ( 1 ).tot );
            h->h_Cross_DT[3]->Fill ( vec_pair.at ( a +1).at ( 1 ).drifttime );
            
            h->h_Cross_TOTvTOT[0]->Fill(vec_pair.at ( a +1).at ( 0 ).tot,vec_pair.at ( a +1).at ( 1 ).tot);
            h->h_Cross_TOTvTOT[1]->Fill(vec_pair.at ( a +1).at ( 0 ).tot,vec_pair.at ( a +1).at ( 1 ).tot);
            
            h->h_Cross_DTvDT[0]->Fill(vec_pair.at ( a +1).at ( 0 ).drifttime,vec_pair.at ( a +1).at ( 1 ).drifttime);
            h->h_Cross_DTvDT[1]->Fill(vec_pair.at ( a +1).at ( 0 ).drifttime,vec_pair.at ( a +1).at ( 1 ).drifttime);
            
            h->h_Cross_DTsum[0]->Fill(vec_pair.at ( a ).at ( 0 ).drifttime + vec_pair.at ( a ).at ( 1 ).drifttime);
            h->h_Cross_DTsum[1]->Fill(vec_pair.at ( a +1).at ( 0 ).drifttime + vec_pair.at ( a +1).at ( 1 ).drifttime);
            
            h->h_Cross_DTvTOT[1]->Fill(vec_pair.at ( a ).at ( 0 ).drifttime,vec_pair.at ( a).at ( 0 ).tot);
            h->h_Cross_DTvTOT[2]->Fill(vec_pair.at ( a ).at ( 1 ).drifttime,vec_pair.at ( a).at ( 1 ).tot);
            h->h_Cross_DTvTOT[3]->Fill(vec_pair.at ( a +1).at ( 0 ).drifttime,vec_pair.at ( a +1).at ( 0 ).tot);
            
            if ( vec_pair.at ( a ).at ( 0 ).tot > vec_pair.at ( a+1 ).at ( 1 ).tot ) {
                index.push_back ( a );
            } else {
                index.push_back ( a+1 );
            }
        }
    }
    int offset =0;
    for ( int b=0; b<index.size(); b++ ) {
        //printf ( "Indexs : %i \t",index[b] );
        vec_pair.erase ( vec_pair.begin() +index[b]-offset );
        offset++;
    }
//     printf ( "\n" );
// 
//     for ( int b=0; b<vec_pair.size(); b++ ) {
//         printf ( "Pair %i\n",b+1 );
//         for ( int c=0; c<vec_pair[b].size(); c++ ) {
//             printf ( "%i \n",vec_pair.at ( b ).at ( c ).straw );
//         }
//         printf ( "\n\n" );
//     }
//     printf ( "**************************************************\n" );
    return vec_pair;
}


std::vector<std::vector<VecSttHit>>
                                 make_tuples ( const std::vector<std::vector<VecSttHit>>& vectors )
{
    std::vector<std::vector<VecSttHit>> result;
    add_tuples ( vectors, 0, {}, result );
    return result;
}

std::vector<VecSttHit> clusterfinder ( VecSttHit vec_flayer )
{
    // clusterfinder ( VecSttHit vec_flayer ) {
    std::vector<VecSttHit> vec_Cl;
    VecSttHit clusterPointer;

    std::sort ( vec_flayer.begin(), vec_flayer.end(), f_sttHitCompareCell );

    if ( vec_flayer.size() == 1 ) {
        clusterPointer.clear();
        clusterPointer.push_back ( vec_flayer[0] );
        vec_Cl.push_back ( clusterPointer );
        return vec_Cl;
    } else {
        for ( Int_t aa = 0; aa < vec_flayer.size(); aa++ ) {
            if ( aa < vec_flayer.size() - 1 ) {
                if ( fabs ( vec_flayer[aa].straw - vec_flayer[aa + 1].straw ) < 2 ) {
                    clusterPointer.clear();
                    clusterPointer.push_back ( vec_flayer[aa] );
                    clusterPointer.push_back ( vec_flayer[aa + 1] );
                    vec_Cl.push_back ( clusterPointer );
                } else if ( aa == 0 ) {
                    clusterPointer.clear();
                    clusterPointer.push_back ( vec_flayer[aa] );
                    vec_Cl.push_back ( clusterPointer );
                } else if ( aa > 0 ) {
                    if ( fabs ( vec_flayer[aa - 1].straw -
                    vec_flayer[aa].straw ) == 1 ) {
                        continue;
                    }

                    else {
                        clusterPointer.clear();
                        clusterPointer.push_back ( vec_flayer[aa] );
                        vec_Cl.push_back ( clusterPointer );
                    }
                }
            }

            else if ( aa == vec_flayer.size() - 1 ) {
                if ( ( fabs ( vec_flayer[aa - 1].straw - vec_flayer[aa].straw ) >
                1 ) ) {
                    clusterPointer.clear();
                    clusterPointer.push_back ( vec_flayer[aa] );
                    vec_Cl.push_back ( clusterPointer );
                }
            }
        }
        return vec_Cl;
    }
}

bool PDAQ_Event_Finder ( VecSttHit vec_stthits, int i,
                         TTree* PDAQ_tree, Stt_Track_Event* stt_event,
                         MFTGeomPar* ftGeomPar, PandaSubsystemSCI* SCI_CAL,  histograms* h )
{
    //    printf("*** Event %d\n", i);
    // for ( Int_t tq=0; tq<vec_stthits.size(); tq++ ) {
    // vec_tracks[tq].drifttime =  max_dt_offset+(meanTime - (
    // vec_tracks[tq].leadTime ) ) ;
    for ( int ac=0; ac< vec_stthits.size(); ac++ ) {
        //h->h_tot3->Fill ( vec_stthits[ac].tot );
    }
    // 	printf("\n*********************\n");
    //}

    int max_cluster_intake = 0;

    // max_cluster_intake = ftGeomPar->getClusterLimit();
    //     float max_dt_offset = ftGeomPar->getDTOffset();
    int stations = ftGeomPar->getModules();
    int MAX_FT_TOTAL_LAYERS = 0;

    if ( stations > 1 ) {
        for ( int a = 0; a < stations; a++ ) {
            MAX_FT_TOTAL_LAYERS += ftGeomPar->getLayers ( a );
        }
    }

    else {
        MAX_FT_TOTAL_LAYERS = ftGeomPar->getLayers ( 0 );
    }

    std::vector<VecSttHit> vec_layer;
    VecSttHit vec_hitlayer;
    vec_layer.clear();
    vec_hitlayer.clear();

    //Check to filter out hits with no layer info in any//////
    for ( int s = 1; s < MAX_FT_TOTAL_LAYERS + 1; s++ ) {
        for ( int ss = 0; ss < vec_stthits.size(); ss++ ) {
            if ( vec_stthits[ss].layer == s ) {
                vec_hitlayer.push_back ( vec_stthits[ss] );
            }
        }
        vec_layer.push_back ( vec_hitlayer );
        vec_hitlayer.clear();
    }

    VecSttHit vec_Clusters;

    std::vector<double> vec_Chi2x;
    std::vector<double> vec_Chi2y;

    std::vector<double> vec_P0;
    std::vector<double> vec_P1;

    std::vector<double> vec_PP0;
    std::vector<double> vec_PP1;

    // // FILTER TO GET ONLY HITS WITH A PAIR ////////////////////////////////////////////////////////////

    VecSttHit vec_pair_clu;
    vec_pair_clu.clear();

    std::vector<VecSttHit> vec_pair_layer;
    vec_pair_layer.clear();
    VecSttHit vec_player;
    vec_player.clear();

    //////////////////////////////////////////////////////////////////////////

    for ( int j = 0; j < vec_layer.size(); j++ ) {
        if ( vec_layer[j].size() ==1 ) {
            h->h_Cross_TOT[0]->Fill ( vec_layer.at ( j ).at ( 0 ).tot );
            h->h_Cross_DT[0]->Fill ( vec_layer.at ( j ).at ( 0 ).drifttime );
            h->h_Cross_DTvTOT[0]->Fill(vec_layer.at ( j ).at ( 0 ).drifttime,vec_layer.at ( j ).at ( 0 ).tot);
        } else if ( vec_layer[j].size() > 1 ) {
            vec_player = GetPairs ( vec_layer[j],h );
        }
        vec_pair_layer.push_back ( vec_player );
        vec_player.clear();
    }

    std::vector<vector<VecSttHit>> vec_cluster_layer;
    vec_cluster_layer.clear();
    std::vector<VecSttHit> vec_clayer;
    vec_clayer.clear();

    for ( int j = 0; j < vec_pair_layer.size(); j++ ) {
        if ( vec_pair_layer[j].size() > 1 ) {
            VecSttHit vec_imi = vec_pair_layer[j];
           // cout<<"Layer "<< j+1<<endl;
            vec_clayer = bestPair ( clusterfinder ( vec_pair_layer[j] ), h );
        }
        vec_cluster_layer.push_back ( vec_clayer );
        vec_clayer.clear();
    }
  //  cout<<"############################################################################"<<endl;


    int track_sanity = 0;
    for ( int da = 0; da < vec_cluster_layer.size(); da++ ) {
        if ( vec_cluster_layer[da].size() > 0 ) {
            track_sanity++;
        }
    }
    if ( track_sanity >= MAX_FT_TOTAL_LAYERS ) {
        const std::vector<std::vector<VecSttHit>> vectors =
            vec_cluster_layer;

        const auto tuples = make_tuples ( vectors );
        std::vector<VecSttHit> vec_All_X;
        std::vector<VecSttHit> vec_All_Y;
        vec_All_X.clear();
        vec_All_Y.clear();

        std::size_t ntups = 0;
        for ( const auto& tup : tuples ) {
            vec_Clusters.clear();
            for ( int a = 0; a < tup.size(); a++ ) {
                for ( int b = 0; b < tup[a].size(); b++ ) {
                    vec_Clusters.push_back ( tup.at ( a ).at ( b ) );
                }
            }
//cout<<"SIZE : "<<vec_Clusters.size()<<endl;

            Double_t chi_valueX = 0;
            Double_t chi_valueY = 0;
            Int_t loneCounterX = 0;
            Int_t loneCounterY = 0;

            VecSttHit vec_ClustersX;
            VecSttHit vec_ClustersY;

            vec_ClustersX.clear();
            vec_ClustersY.clear();

            for ( Int_t ya = 0; ya < vec_Clusters.size(); ya++ ) {
                if ( vec_Clusters[ya].layer == 2 ||
                        vec_Clusters[ya].layer == 3 ||
                        vec_Clusters[ya].layer == 6 ||
                vec_Clusters[ya].layer == 7 ) {
                    vec_ClustersY.push_back ( vec_Clusters[ya] );
                    loneCounterY++;
//                     printf ( "X : TDC :%x , Layer -%d , Straw -%d \n",vec_Clusters[ya]->tdcid,vec_Clusters[ya]->layer,vec_Clusters[ya]->straw );
                } else {
                    vec_ClustersX.push_back ( vec_Clusters[ya] );
                    loneCounterX++;
                    // printf("X : TDC :%x , Layer : %d , Straw : %d , X :
                    // %lf\n",vec_Clusters[ya]->tdcid,vec_Clusters[ya]->layer,vec_Clusters[ya]->straw,vec_Clusters[ya]->x);
                }
                // printf("*****************\n\n");
            }

            vec_All_X.push_back ( vec_ClustersX );
            vec_All_Y.push_back ( vec_ClustersY );

            Double_t clusterArrayX[vec_ClustersX.size()];
            Double_t clusterArrayZx[vec_ClustersX.size()];
            Double_t clusterArrayY[vec_ClustersY.size()];
            Double_t clusterArrayZy[vec_ClustersY.size()];

            for ( Int_t yb = 0; yb < vec_ClustersX.size(); yb++ ) {
                clusterArrayX[yb] = vec_ClustersX[yb].x;
                clusterArrayZx[yb] = vec_ClustersX[yb].z;
                // cout<<yb<<"\t"<<vec_ClustersX[yb]->x<<endl;
            }

            for ( Int_t yc = 0; yc < vec_ClustersY.size(); yc++ ) {
                clusterArrayY[yc] = vec_ClustersY[yc].y;
                clusterArrayZy[yc] = vec_ClustersY[yc].z;

            }

            TF1* f1 = new TF1 ( "f1", "pol1" );
            TF1* f2 = new TF1 ( "f2", "pol1" );
            TGraph* chiX = new TGraph ( vec_ClustersX.size(), clusterArrayX, clusterArrayZx );
            chiX->Fit ( f1, "q" );
            chi_valueX = f1->GetChisquare();
            vec_Chi2x.push_back ( chi_valueX );

            Double_t p0 = f1->GetParameter ( 0 );
            Double_t p1 = f1->GetParameter ( 1 );

            vec_P0.push_back ( p0 );
            vec_P1.push_back ( p1 );

            TGraph* chiY =
                new TGraph ( vec_ClustersY.size(), clusterArrayY, clusterArrayZy );
            chiY->Fit ( f2, "q" );
            chi_valueY = f2->GetChisquare();
            vec_Chi2y.push_back ( chi_valueY );
            Double_t pp0 = f2->GetParameter ( 0 );
            Double_t pp1 = f2->GetParameter ( 1 );

            vec_PP0.push_back ( pp0 );
            vec_PP1.push_back ( pp1 );

            delete f1;
            delete f2;
            delete chiX;
            delete chiY;
        }

        ////////////////////////////////////////////////////////////////////////////////////////

        std::vector<SttHit> vec_tracks;

        Float_t smallestX = vec_Chi2x[0];
        Float_t smallestP0 = vec_P0[0];
        Float_t smallestP1 = vec_P1[0];

        Int_t chi_indexX = 0;

        for ( Int_t ci = 0; ci < vec_Chi2x.size(); ci++ ) {

            if ( smallestX > vec_Chi2x[ci] ) {
                smallestX = vec_Chi2x[ci];
                chi_indexX = ci;
            }
            if ( smallestP0 > vec_P0[ci] ) {
                smallestP0 = vec_P0[ci];
            }
            if ( smallestP1 > vec_P1[ci] ) {
                smallestP1 = vec_P1[ci];
            }
        }
        // printf("smallest chix :%2.3f\n", smallestX);

        Float_t smallestY = vec_Chi2y[0];
        Float_t smallestPP0 = vec_PP0[0];
        Float_t smallestPP1 = vec_PP1[0];

        Int_t chi_indexY = 0;

        for ( Int_t cj = 0; cj < vec_Chi2y.size(); cj++ ) {

            if ( smallestY > vec_Chi2y[cj] ) {
                smallestY = vec_Chi2y[cj];
                chi_indexY = cj;
            }
            if ( smallestPP0 > vec_PP0[cj] ) {
                smallestPP0 = vec_PP0[cj];
            }
            if ( smallestPP1 > vec_PP1[cj] ) {
                smallestPP1 = vec_PP1[cj];
            }
        }

        for ( Int_t ck = 0; ck < vec_All_X.at ( chi_indexX ).size(); ck++ ) {
            vec_tracks.push_back ( vec_All_X.at ( chi_indexX ).at ( ck ) );
            // printf("Chi best : TDC :%x , Layer -%d , Straw -%d
            // \n",vec_All_X.at ( chi_indexX ).at ( ck )->tdcid, vec_All_X.at (
            // chi_indexX ).at ( ck )->layer,vec_All_X.at ( chi_indexX ).at ( ck
            // )->straw);
        }
        // printf("\n$$$$$$$$$$$$$$$$$$$$$$$$$\n");

        for ( Int_t cl = 0; cl < vec_All_Y.at ( chi_indexY ).size(); cl++ ) {
            vec_tracks.push_back ( vec_All_Y.at ( chi_indexY ).at ( cl ) );
        }
        double sumLeadTime = 0;
        double meanTime = 0;

        for ( Int_t d = 0; d < vec_tracks.size(); d++ ) {
            sumLeadTime += vec_tracks.at ( d ).leadTime;
            //printf("tdc: %x ch: %i LT: %lf DT: %lf \n",vec_tracks[d]->tdcid,vec_tracks[d]->straw,vec_tracks[d]->leadTime,vec_tracks[d]->drifttime);
        }
        // printf("\n*********************\n");

        meanTime = sumLeadTime / vec_tracks.size();
        // printf("mean = %f\n", meanTime);
        // Write Tracks
        stt_event->TrackClear();

        //     b->Py0 = smallestPP0;
        //     b->Py1 = smallestPP1;
        //     b->Chix = smallestX;
        //     b->Chiy = smallestY;
        // printf(" Mean time %f track size %d \n", meanTime,vec_tracks.size());

        SttTrackHit& b = stt_event->AddTrackHit();
        b.vec_Track = vec_tracks;
        b.trackId = i;
        b.trackSize = vec_tracks.size();
        b.Px0 = smallestP0;
        b.Px1 = smallestP1;


        for ( Int_t tq = 0; tq < vec_tracks.size(); tq++ ) {
            h->h_tot4->Fill ( vec_tracks[tq].tot );
            h->h_L_dtvstot[vec_tracks[tq].layer-1]->Fill ( vec_tracks[tq].drifttime,vec_tracks[tq].tot );
            h->h_drifttimevsLayer->Fill ( vec_tracks[tq].drifttime,vec_tracks[tq].layer );
            h->h_drifttime->Fill ( vec_tracks[tq].drifttime );
            if ( vec_tracks[tq].tdcid==0x6400||vec_tracks[tq].tdcid==0x6410||vec_tracks[tq].tdcid==0x6411 ) {
                h->h_drifttimeTRB1->Fill ( vec_tracks[tq].drifttime );
            } else {
                h->h_drifttimeTRB2->Fill ( vec_tracks[tq].drifttime );
            }
            h->h_hitBlock->Fill ( 4 );
        }

        vec_Chi2x.clear();
        vec_Chi2y.clear();

        vec_P0.clear();
        vec_P1.clear();

        vec_PP0.clear();
        vec_PP1.clear();

        PDAQ_tree->Fill();

        return true;
    }
    return false;
}
