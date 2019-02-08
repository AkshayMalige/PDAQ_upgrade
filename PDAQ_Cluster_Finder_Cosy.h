#ifndef PDAQ_CLUSTER_FINDER_COSY_H
#define PDAQ_CLUSTER_FINDER_COSY_H

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

#include <cctype>
#include <fstream>
#include <sstream>
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>
#include "FTGeo.h"
#include <algorithm>

using namespace std;


template< typename T >
struct delete_pointer_element
{
    void operator()( T element ) const
    {
        delete element;
    }
};

bool f_sttHitCompareLeadTime ( SttHit* a, SttHit* b ) {
    return ( a->leadTime < b->leadTime );
}


bool f_sttHitCompareCell ( SttHit* a, SttHit* b ) {
    return ( a->straw < b->straw );
}


std::vector<SttHit*>  GetPairs ( std::vector<SttHit*> vec_get_pairs ) {
    std::vector<SttHit*> vec_fpair_clu;
    std::sort ( vec_get_pairs.begin(), vec_get_pairs.end(), f_sttHitCompareCell );

    for ( Int_t sa=0; sa<vec_get_pairs.size(); sa++ ) {
        if ( sa < vec_get_pairs.size()-1 ) {
            if ( fabs ( vec_get_pairs[sa]->straw - vec_get_pairs[sa+1]->straw ) == 1 ) {
                vec_fpair_clu.push_back ( vec_get_pairs[sa] );
            }

            else if ( sa >0 ) {
                if ( fabs ( vec_get_pairs[sa-1]->straw - vec_get_pairs[sa]->straw ) == 1 )
                    vec_fpair_clu.push_back ( vec_get_pairs[sa] );
            }

        }
        else if ( sa == vec_get_pairs.size()-1 ) {
            if ( fabs ( vec_get_pairs[sa-1]->straw - vec_get_pairs[sa]->straw ) == 1 )
                vec_fpair_clu.push_back ( vec_get_pairs[sa] );
        }

    }
    return vec_fpair_clu;
}

std::vector<std::vector<SttHit*>> clusterfinder ( std::vector<SttHit*> vec_flayer ) {
    std::vector<std::vector<SttHit*>> vec_Cl;
    std::vector<SttHit*> clusterPointer;

    std::sort ( vec_flayer.begin(), vec_flayer.end(), f_sttHitCompareCell );

    if ( vec_flayer.size() ==1 ) {
        clusterPointer.clear();
        clusterPointer.push_back ( vec_flayer[ 0 ] );
        vec_Cl.push_back ( clusterPointer );
        return vec_Cl;
    } else {
        for ( Int_t aa=0; aa < vec_flayer.size(); aa++ ) {
            if ( aa < vec_flayer.size()-1 ) {
                if ( fabs ( vec_flayer[ aa ]->straw - vec_flayer[ aa+1]->straw ) < 2 ) {
                    clusterPointer.clear();
                    clusterPointer.push_back ( vec_flayer[ aa ] );
                    clusterPointer.push_back ( vec_flayer[ aa+1 ] );
                    vec_Cl.push_back ( clusterPointer );
                } else if ( aa==0 ) {
                    clusterPointer.clear();
                    clusterPointer.push_back ( vec_flayer[ aa ] );
                    vec_Cl.push_back ( clusterPointer );
                } else if ( aa > 0 ) {
                    if ( fabs ( vec_flayer[ aa-1]->straw - vec_flayer[ aa]->straw ) == 1 ) {
                        continue;
                    }

                    else {
                        clusterPointer.clear();
                        clusterPointer.push_back ( vec_flayer[ aa ] );
                        vec_Cl.push_back ( clusterPointer );

                    }
                }
            }

            else if ( aa == vec_flayer.size()-1 ) {
                if ( ( fabs ( vec_flayer[ aa-1 ]->straw - vec_flayer[ aa ]->straw ) > 1 ) ) {
                    clusterPointer.clear();
                    clusterPointer.push_back ( vec_flayer[ aa ] );
                    vec_Cl.push_back ( clusterPointer );
                }
            }
        }
        return vec_Cl;
    }
}

void add_tuples( const std::vector< std::vector<std::vector<SttHit*> >>& vectors, std::size_t pos,
                 std::vector<std::vector<SttHit*>> prefix, std::vector< std::vector<std::vector<SttHit*> >>& result )
{
    if( pos == vectors.size() ) result.push_back(prefix) ;
    else if( vectors[pos].empty() ) add_tuples( vectors, pos+1, prefix, result ) ; // note: skip empty vectors
    else
    {
        for( std::vector<SttHit*> v : vectors[pos] )
        {
            prefix.push_back(v) ;
            add_tuples( vectors, pos+1, prefix, result ) ;
            prefix.pop_back() ;
        }
    }
}

std::vector< std::vector<std::vector<SttHit*> >> make_tuples( const std::vector<std::vector< std::vector<SttHit*> >>& vectors )
{
    std::vector< std::vector<std::vector<SttHit*> >> result ;
    add_tuples( vectors, 0, {}, result ) ;
    return result ;
}

int PDAQ_Cluster_Finder_Cosy(char* intree, int maxEvents = 1000000000, char* outtree="PDAQ_cluster_output.root");

bool PDAQ_Event_Finder(std::vector<SttHit*> vec_stthits, int i,TTree* PDAQ_tree,Stt_Track_Event* stt_event)
{
      //for ( Int_t tq=0; tq<vec_stthits.size(); tq++ ) {
	//vec_tracks[tq]->drifttime =  max_dt_offset+(meanTime - ( vec_tracks[tq]->leadTime ) ) ;
// 	for (int ac=0; ac< vec_stthits.size(); ac++)
// 	{
// 	  printf("TDC :%x , Layer -%d , Straw -%d  \n",vec_stthits[ac]->tdcid, vec_stthits[ac]->layer,vec_stthits[ac]->straw);
// 	}
// 	printf("\n*********************\n");
    //}
    int max_cluster_intake=0;
    float max_dt_offset = 0;     

    MParManager* a = MParManager::instance();
    MFTGeomPar* ftGeomPar = new MFTGeomPar();
    MPar * d;
	
    a->setParamSource("ftparams.txt");
    a->parseSource();
    pm()->addParameterContainer("MFTPar", ftGeomPar);	
 
    //max_cluster_intake = ftGeomPar->getClusterLimit();
    max_dt_offset = ftGeomPar->getDTOffset();
    int stations = ftGeomPar->getModules();    
    int MAX_FT_TOTAL_LAYERS =0;
    
    if (stations > 1)
    {
      for (int a=0; a<stations; a++)
      {
	MAX_FT_TOTAL_LAYERS +=  ftGeomPar->getLayers(a);
      }          
    }
    
    else
    {
      MAX_FT_TOTAL_LAYERS = ftGeomPar->getLayers(0);
    }

    std::vector< vector<SttHit*> > vec_layer;
    std::vector<SttHit*> vec_hitlayer;
    vec_layer.clear();
    vec_hitlayer.clear();

    for (int s =1; s<MAX_FT_TOTAL_LAYERS+1; s++)
      {
      for (int ss=0; ss<vec_stthits.size(); ss++)
	{
	  if (vec_stthits[ss]->layer == s)
	  {
	    vec_hitlayer.push_back(vec_stthits[ss]);
	}
      }
      vec_layer.push_back(vec_hitlayer);
      vec_hitlayer.clear();      
      }

    std::vector<SttHit*> vec_Clusters;

    std::vector<double> vec_Chi2x;
    std::vector<double> vec_Chi2y;

    std::vector<double> vec_P0;
    std::vector<double> vec_P1;

    std::vector<double> vec_PP0;
    std::vector<double> vec_PP1;
    
// // FILTER TO GET ONLY HITS WITH A PAIR////////////////////////////////////////////////////////////

    std::vector<SttHit*> vec_pair_clu;
    vec_pair_clu.clear();
   
    std::vector< vector<SttHit*> > vec_pair_layer;
    vec_pair_layer.clear();
    std::vector<SttHit*> vec_player;
    vec_player.clear();

//////////////////////////////////////////////////////////////////////////
    
    for (int j=0; j<vec_layer.size();j++)
    {
      if (vec_layer[j].size()>1)
      {
	vec_player = GetPairs (vec_layer[j]);
      }
      vec_pair_layer.push_back(vec_player);
      vec_player.clear();
    }

    std::vector<vector<vector<SttHit*>>> vec_cluster_layer;
    vec_cluster_layer.clear();
    std::vector< vector<SttHit*> > vec_clayer;
    vec_clayer.clear();

    for (int j=0; j<vec_pair_layer.size();j++)
    {
      if (vec_pair_layer[j].size()>1)
      {
	std::vector<SttHit*> vec_imi = vec_pair_layer[j];
	vec_clayer = clusterfinder (vec_pair_layer[j]);
      }
      vec_cluster_layer.push_back(vec_clayer);
      vec_clayer.clear();
    }   
    int track_sanity = 0;
     for (int da=0; da< vec_cluster_layer.size();da++)
     {
       if (vec_cluster_layer[da].size()>0)
       {
	 track_sanity++ ; 
      }
     }
     if (track_sanity >= MAX_FT_TOTAL_LAYERS)
     {
      const std::vector< std::vector<std::vector<SttHit*> > >vectors =  vec_cluster_layer ;
 
      const auto tuples = make_tuples(vectors) ;
      std::vector< vector<SttHit*> > vec_All_X;
      std::vector< vector<SttHit*> > vec_All_Y;
      vec_All_X.clear();
      vec_All_Y.clear();
      
      std::size_t ntups = 0 ;
      for( const auto& tup : tuples )
      {
	  vec_Clusters.clear();	    
	  for (int a =0; a<tup.size(); a++)
	  {
	    for (int b=0; b<tup[a].size(); b++)
	    {
	      vec_Clusters.push_back(tup.at(a).at(b));
	    }	    
	  }
	    
	  Double_t chi_valueX =0;
	  Double_t chi_valueY =0;
	  Int_t loneCounterX =0;
	  Int_t loneCounterY =0;

	  std::vector<SttHit*> vec_ClustersX;
	  std::vector<SttHit*> vec_ClustersY;

	  vec_ClustersX.clear();
	  vec_ClustersY.clear();

	  for ( Int_t ya=0; ya<vec_Clusters.size(); ya++ ) {
	      
	    if ( vec_Clusters[ya]->layer ==2 || vec_Clusters[ya]->layer ==3 || vec_Clusters[ya]->layer ==6 || vec_Clusters[ya]->layer ==7) {
		  vec_ClustersY.push_back ( vec_Clusters[ya] );
		  loneCounterY++;
		  //printf("X : TDC :%x , Layer -%d , Straw -%d  \n",vec_Clusters[ya]->tdcid, vec_Clusters[ya]->layer,vec_Clusters[ya]->straw);

	      }
	      else {
		  vec_ClustersX.push_back ( vec_Clusters[ya] );
		  loneCounterX++;
		  //printf("Y : TDC :%x , Layer -%d , Straw -%d  \n",vec_Clusters[ya]->tdcid, vec_Clusters[ya]->layer,vec_Clusters[ya]->straw);

	      }
	      //printf("*****************\n\n");
	  }

	  vec_All_X.push_back ( vec_ClustersX );
	  vec_All_Y.push_back ( vec_ClustersY );

	  Double_t clusterArrayX[vec_ClustersX.size()];
	  Double_t clusterArrayZx[vec_ClustersX.size()];
	  Double_t clusterArrayY[vec_ClustersY.size()];
	  Double_t clusterArrayZy[vec_ClustersY.size()];


	  for ( Int_t yb=0; yb<vec_ClustersX.size(); yb++ ) {
	      clusterArrayX[yb]=vec_ClustersX[yb]->x;
	      clusterArrayZx[yb]=vec_ClustersX[yb]->z;
	  }


	  for ( Int_t yc=0; yc<vec_ClustersY.size(); yc++ ) {
	      clusterArrayY[yc]=vec_ClustersY[yc]->y;
	      clusterArrayZy[yc]=vec_ClustersY[yc]->z;
	  }


	  TF1* f1 = new TF1 ( "f1", "pol1" );
	  //TF1* f2 = new TF1 ( "f2", "pol1" );	 
	  TGraph* chiX = new TGraph ( vec_ClustersX.size(), clusterArrayX, clusterArrayZx );
	  chiX->Fit ( f1,"q" );
	  chi_valueX = f1->GetChisquare();
	  vec_Chi2x.push_back ( chi_valueX );

	  Double_t p0 = f1->GetParameter ( 0 );
	  Double_t p1 = f1->GetParameter ( 1 );

	  vec_P0.push_back ( p0 );
	  vec_P1.push_back ( p1 );

// 	  TGraph* chiY = new TGraph ( vec_ClustersY.size(), clusterArrayY, clusterArrayZy );
// 	  chiY->Fit ( f2,"q" );
// 	  chi_valueY = f2->GetChisquare();
// 	  vec_Chi2y.push_back ( chi_valueY );
// 	  Double_t pp0 = f2->GetParameter ( 0 );
// 	  Double_t pp1 = f2->GetParameter ( 1 );

// 	  vec_PP0.push_back ( pp0 );
// 	  vec_PP1.push_back ( pp1 );

	  delete f1;
//	  delete f2;
	  delete chiX;
// 	  delete chiY;
 	  
      }
   
////////////////////////////////////////////////////////////////////////////////////////
	    
    std::vector<SttHit*> vec_tracks;
    vec_tracks.clear();

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
    	    //printf("smallest chix :%2.3f\n", smallestX);

//     Float_t smallestY = vec_Chi2y[0];
//     Float_t smallestPP0 = vec_PP0[0];
//     Float_t smallestPP1 = vec_PP1[0];
// 
//     Int_t chi_indexY = 0;
// 
//     for ( Int_t cj = 0; cj < vec_Chi2y.size(); cj++ ) {
// 
// 	if ( smallestY > vec_Chi2y[cj] ) {
// 	    smallestY = vec_Chi2y[cj];
// 	    chi_indexY = cj;
// 	}
// 	if ( smallestPP0 > vec_PP0[cj] ) {
// 	    smallestPP0 = vec_PP0[cj];
// 	}
// 	if ( smallestPP1 > vec_PP1[cj] ) {
// 	    smallestPP1 = vec_PP1[cj];
// 	}
//     }

    for ( Int_t ck = 0; ck < vec_All_X.at ( chi_indexX ).size(); ck++ ) {
	vec_tracks.push_back ( vec_All_X.at ( chi_indexX ).at ( ck ) );
	//printf("Chi best : TDC :%x , Layer -%d , Straw -%d  \n",vec_All_X.at ( chi_indexX ).at ( ck )->tdcid, vec_All_X.at ( chi_indexX ).at ( ck )->layer,vec_All_X.at ( chi_indexX ).at ( ck )->straw);

    }
    	//printf("\n$$$$$$$$$$$$$$$$$$$$$$$$$\n");

//     for ( Int_t cl = 0; cl < vec_All_Y.at ( chi_indexY ).size(); cl++ ) {
// 	vec_tracks.push_back ( vec_All_Y.at ( chi_indexY ).at ( cl ) );
//     }

    double sumLeadTime = 0;
    double meanTime = 0;

    for ( Int_t d = 0; d<vec_tracks.size(); d++ ) {
	sumLeadTime += vec_tracks.at ( d )->leadTime;
	//printf("track size : %d", vec_tracks.size());
// 	for (int ac=0; ac< 8; ac++)
// 	  {
	   // printf("TDC :%x , Layer -%d , Straw -%d  \n",vec_tracks[d]->tdcid, vec_tracks[d]->layer,vec_tracks[d]->straw);
	  //}
    }
	//printf("\n*********************\n");

    meanTime = sumLeadTime/vec_tracks.size();

//Write Tracks
    stt_event->TrackClear();

    SttTrackHit* b = stt_event->AddTrackHit();
    b->vec_Track = vec_tracks;
    b->trackId = i;
    b->trackSize = vec_tracks.size();
    b->Px0 = smallestP0;
    b->Px1 = smallestP1;
//     b->Py0 = smallestPP0;
//     b->Py1 = smallestPP1;
//     b->Chix = smallestX;
//     b->Chiy = smallestY;


    for ( Int_t tq=0; tq<vec_tracks.size(); tq++ ) {
	vec_tracks[tq]->drifttime =  max_dt_offset+(meanTime - ( vec_tracks[tq]->leadTime ) ) ;
// 	for (int ac=0; ac< vec_tracks.size(); ac++)
// 	{
// 	  printf("TDC :%x , Layer -%d , Straw -%d  \n",vec_tracks[tq][ac].tdcid, vec_tracks[tq][ac].layer,vec_tracks[tq][ac].straw);
// 	}
// 	printf("\n*********************\n");
    }

    
    vec_Chi2x.clear();
    vec_Chi2y.clear();

    vec_P0.clear();
    vec_P1.clear();

//     vec_PP0.clear();
//     vec_PP1.clear();

    PDAQ_tree->Fill();

return true;
     }  
}

#endif