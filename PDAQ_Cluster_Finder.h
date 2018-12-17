#ifndef PDAQ_CLUSTER_FINDER_H
#define PDAQ_CLUSTER_FINDER_H

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
//#include "FT_vector_def.h"


#include <cctype>
#include <fstream>
#include <sstream>
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>
#include "FTGeo.h"

using namespace std;

bool f_sttHitCompareLeadTime ( SttHit* a, SttHit* b ) {
    return ( a->leadTime < b->leadTime );
}


bool f_sttHitCompareCell ( SttHit* a, SttHit* b ) {
    return ( a->straw < b->straw );
}

int a =0;

std::vector<std::vector<SttHit*>*>* clusterfinder ( std::vector<SttHit*>* vec_flayer ) {
    std::vector<std::vector<SttHit*>*>* vec_Cl = new std::vector<std::vector<SttHit*>*>();
    std::vector<SttHit*>* clusterPointer;

    std::sort ( vec_flayer->begin(), vec_flayer->end(), f_sttHitCompareCell );

    if ( vec_flayer->size() ==1 ) {
        clusterPointer = new std::vector<SttHit*>();
        clusterPointer->push_back ( vec_flayer->at ( 0 ) );
        vec_Cl->push_back ( clusterPointer );
        return vec_Cl;
    } else {
        for ( Int_t aa=0; aa < vec_flayer->size(); aa++ ) {
            if ( aa < vec_flayer->size()-1 ) {
                if ( fabs ( vec_flayer->at ( aa )->straw - vec_flayer->at ( aa+1 )->straw ) < 2 ) {
                    clusterPointer = new vector<SttHit*>();
                    clusterPointer->push_back ( vec_flayer->at ( aa ) );
                    clusterPointer->push_back ( vec_flayer->at ( aa+1 ) );
                    vec_Cl->push_back ( clusterPointer );
                } else if ( aa==0 ) {
                    clusterPointer = new vector<SttHit*>();
                    clusterPointer->push_back ( vec_flayer->at ( aa ) );
                    vec_Cl->push_back ( clusterPointer );
                } else if ( aa > 0 ) {
                    if ( fabs ( vec_flayer->at ( aa-1 )->straw - vec_flayer->at ( aa )->straw ) == 1 ) {
                        continue;
                    }

                    else {
                        clusterPointer = new vector<SttHit*>();
                        clusterPointer->push_back ( vec_flayer->at ( aa ) );
                        vec_Cl->push_back ( clusterPointer );

                    }
                }
            }

            else if ( aa == vec_flayer->size()-1 ) {
                if ( ( fabs ( vec_flayer->at ( aa-1 )->straw - vec_flayer->at ( aa )->straw ) > 1 ) ) {
                    clusterPointer = new vector<SttHit*>();
                    clusterPointer->push_back ( vec_flayer->at ( aa ) );
                    vec_Cl->push_back ( clusterPointer );
                }
            }
        }
        return vec_Cl;
    }
}

std::vector<SttHit*>  GetPairs ( std::vector<SttHit*> vec_get_pairs ) {
    std::vector<SttHit*> vec_fpair_clu;
    std::sort ( vec_get_pairs.begin(), vec_get_pairs.end(), f_sttHitCompareCell );


    for ( Int_t sa=0; sa<vec_get_pairs.size(); sa++ ) {

        if ( sa < vec_get_pairs.size()-1 ) {
            if ( fabs ( vec_get_pairs[sa]->straw - vec_get_pairs[sa+1]->straw ) == 1 ) {
                vec_fpair_clu.push_back ( vec_get_pairs[sa] );
                //cout<<"PAIR : "<<vec_get_pairs[sa]->layer<<"\t"<<vec_get_pairs[sa]->straw<<endl;

            }

            else if ( sa >0 ) {
                if ( fabs ( vec_get_pairs[sa-1]->straw - vec_get_pairs[sa]->straw ) == 1 )
                    vec_fpair_clu.push_back ( vec_get_pairs[sa] );
                //cout<<"PAIR : "<<vec_get_pairs[sa]->layer<<"\t"<<vec_get_pairs[sa]->straw<<endl;

            }

        }

        else if ( sa == vec_get_pairs.size()-1 ) {
            if ( fabs ( vec_get_pairs[sa-1]->straw - vec_get_pairs[sa]->straw ) == 1 )
                vec_fpair_clu.push_back ( vec_get_pairs[sa] );
            //cout<<"PAIR : "<<vec_get_pairs[sa]->layer<<"\t"<<vec_get_pairs[sa]->straw<<endl;

        }

    }
    return vec_fpair_clu;

}

//class vector_parameter::;

//Bool_t PDAQ_Cluster_Finder(void);
int PDAQ_Cluster_Finder(char* intree, int maxEvents = 1000000000, char* outtree="PDAQ_cluster_output.root");

bool PDAQ_Event_Finder(std::vector<SttHit*> vec_stthits, int i,TTree* PDAQ_tree,Stt_Track_Event* stt_event)
{
//for (int ak =0; ak<1; ak++){
//     PandaSttTrack* STT_TRACKS = new PandaSttTrack();
//     Stt_Track_Event* stt_event = & ( STT_TRACKS->stt_track_can ); 

//double meanTime = 0;

    
    //printf("VECTOR PARAMETERS : %i, %f, %i, %f", max_cluster_intake,max_lead_time_diff,min_track_hits,max_cluster);
  
    int max_cluster_intake=0;
    float max_dt_offset = 0;
     

    MParManager* a = MParManager::instance();
    MFTGeomPar* ftGeomPar = new MFTGeomPar();
    MPar * d;
	
    a->setParamSource("ftparams.txt");
    a->parseSource();
    pm()->addParameterContainer("MFTPar", ftGeomPar);
	
 
    max_cluster_intake = ftGeomPar->getClusterLimit();
    max_dt_offset = ftGeomPar->getDTOffset();
  
    std::vector<SttHit*> vec_layer1;
    std::vector<SttHit*> vec_layer2;
    std::vector<SttHit*> vec_layer3;
    std::vector<SttHit*> vec_layer4;
	    
    vec_layer1.clear();
    vec_layer2.clear();
    vec_layer3.clear();
    vec_layer4.clear();
	    
  
      //cout<<"vec_stthits ssize : "<<vec_stthits.size()<<endl;

    for ( int w = 0; w < vec_stthits.size(); w++ ) {
	  
    //if (vec_stthits.size()>4){final_counter++;}

	if ( vec_stthits.size() <max_cluster_intake ) {
	    if ( vec_stthits[w]->layer == 1 ) {
		vec_layer1.push_back ( vec_stthits[w] );
	    } else if ( vec_stthits[w]->layer == 2 ) {
		vec_layer2.push_back ( vec_stthits[w] );
	    } else if ( vec_stthits[w]->layer == 3 ) {
		vec_layer3.push_back ( vec_stthits[w] );
	    } else if ( vec_stthits[w]->layer == 4 ) {
		vec_layer4.push_back ( vec_stthits[w] );
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
// // FILTER TO GET ONLY HITS WITH A PAIR////////////////////////////////////////////////////////////

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
    
    if ( vec_layer1.size() >1 && vec_layer2.size() >1 && vec_layer3.size() >1 && vec_layer4.size() >1 ) {
	vec_pair_layer11 = GetPairs ( vec_layer1 );
	vec_pair_layer22 = GetPairs ( vec_layer2 );
	vec_pair_layer33 = GetPairs ( vec_layer3 );
	vec_pair_layer44 = GetPairs ( vec_layer4 );

    } else return false;
    
    //cout<<"Pair SIZES "<<vec_pair_layer11.size()<<"\t"<<vec_pair_layer22.size()<<"\t"<<vec_pair_layer33.size()<<"\t"<<vec_pair_layer44.size()<<"\t"<<endl;
    

    if ( vec_pair_layer11.size() > 1 && vec_pair_layer22.size() > 1 && vec_pair_layer33.size() > 1 && vec_pair_layer44.size() > 1 ) {
	vec_Cl_L1 = clusterfinder ( &vec_pair_layer11 );
	vec_Cl_L2 = clusterfinder ( &vec_pair_layer22 );
	vec_Cl_L3 = clusterfinder ( &vec_pair_layer33 );
	vec_Cl_L4 = clusterfinder ( &vec_pair_layer44 );
    } else return false;
//cout<<"Cluster SIZES "<<vec_Cl_L1->size()<<"\t"<<vec_Cl_L2->size()<<"\t"<<vec_Cl_L3->size()<<"\t"<<vec_Cl_L4->size()<<endl;

    
    std::vector< vector<SttHit*> > vec_All_X;
    std::vector< vector<SttHit*> > vec_All_Y;
    vec_All_X.clear();
    vec_All_Y.clear();

    

    for ( Int_t xa=0; xa<vec_Cl_L1->size(); xa++ ) {

	for ( Int_t xb=0; xb<vec_Cl_L2->size(); xb++ ) {
	    for ( Int_t xc=0; xc< vec_Cl_L3->size(); xc++ ) {
		for ( Int_t xd=0; xd< vec_Cl_L4->size(); xd++ ) {
 		    vec_Clusters.clear();

		    for ( Int_t xaa =0; xaa< vec_Cl_L1->at( xa )->size(); xaa++ ) {
			vec_Clusters.push_back ( vec_Cl_L1->at ( xa )->at ( xaa ) );
		    }

		    for ( Int_t xbb=0; xbb< vec_Cl_L2->at ( xb )->size(); xbb++ ) {
			//cout<< "Vector  : "<<vec_Cl_L2->at(xb)->at(xbb)->layer <<"\t"<< vec_Cl_L2->at(xb)->at(xbb)->straw<<endl;
			vec_Clusters.push_back ( vec_Cl_L2->at ( xb )->at ( xbb ) );
		    }

		    for ( Int_t xcc=0; xcc< vec_Cl_L3->at ( xc )->size(); xcc++ ) {
			//cout<< "Vector  : "<<vec_Cl_L3->at(xc)->at(xcc)->layer <<"\t"<< vec_Cl_L3->at(xc)->at(xcc)->straw<<endl;
			vec_Clusters.push_back ( vec_Cl_L3->at ( xc )->at ( xcc ) );
		    }

		    for ( Int_t xdd=0; xdd< vec_Cl_L4->at ( xd )->size(); xdd++ ) {
			//cout<< "Vector  : "<<vec_Cl_L4->at(xd)->at(xdd)->layer <<"\t"<< vec_Cl_L4->at(xd)->at(xdd)->straw<<endl;
			vec_Clusters.push_back ( vec_Cl_L4->at ( xd )->at ( xdd ) );

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
			//cout<< " Cluster Vector :"<<vec_Clusters[ya]->layer<<"\t"<<vec_Clusters[ya]->straw<<"\t"<<vec_Clusters[ya]->x<<"\t"<<vec_Clusters[ya]->y<<endl;
			if ( vec_Clusters[ya]->layer == 1 || vec_Clusters[ya]->layer ==3 ) {
			    vec_ClustersX.push_back ( vec_Clusters[ya] );
			    loneCounterX++;
			}
			if ( vec_Clusters[ya]->layer == 2 || vec_Clusters[ya]->layer ==4 ) {
			    vec_ClustersY.push_back ( vec_Clusters[ya] );
			    loneCounterY++;
			}
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
		    TF1* f2 = new TF1 ( "f2", "pol1" );	 
		    TGraph* chiX = new TGraph ( vec_ClustersX.size(), clusterArrayX, clusterArrayZx );
		    chiX->Fit ( f1,"q" );
		    chi_valueX = f1->GetChisquare();
		    vec_Chi2x.push_back ( chi_valueX );

		    Double_t p0 = f1->GetParameter ( 0 );
		    Double_t p1 = f1->GetParameter ( 1 );

		    vec_P0.push_back ( p0 );
		    vec_P1.push_back ( p1 );
		    //Double_t extrpX = ((65 - p0) / p1);


		    TGraph* chiY = new TGraph ( vec_ClustersY.size(), clusterArrayY, clusterArrayZy );
		    chiY->Fit ( f2,"q" );
		    chi_valueY = f2->GetChisquare();
		    vec_Chi2y.push_back ( chi_valueY );
		    Double_t pp0 = f2->GetParameter ( 0 );
		    Double_t pp1 = f2->GetParameter ( 1 );

		    vec_PP0.push_back ( pp0 );
		    vec_PP1.push_back ( pp1 );
		    //Double_t extrpY = ((65 - pp0) / pp1);
		    //Double_t sumChi = (chi_valueX + chi_valueY);

		    //cout<<"X CHI :"<<chi_valueX<<"\t"<<"Y CHI"<<chi_valueY<<endl;

		    //cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Chi X : "<<chi_valueX<<"\t"<<"Chi Y : "<<chi_valueY<<"^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
		    delete f1;
		    delete f2;
		    delete chiX;
		    delete chiY;
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

    for ( Int_t ci = 0; ci < vec_Chi2x.size(); ci++ ) {
	//cout << "chiSquare X :" << vec_Chi2x[ci] << endl;

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

    Float_t smallestY = vec_Chi2y[0];
    Float_t smallestPP0 = vec_PP0[0];
    Float_t smallestPP1 = vec_PP1[0];

    Int_t chi_indexY = 0;

    for ( Int_t cj = 0; cj < vec_Chi2y.size(); cj++ ) {
	//cout << "chiSquare Y :" << vec_Chi2y[cj] << endl;

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
    //cout<<"INDEX X:"<<chi_indexX<<"\t"<<"INDEX Y:"<<chi_indexY<<endl;

    for ( Int_t ck = 0; ck < vec_All_X.at ( chi_indexX ).size(); ck++ ) {
	//cout<<"       BEST  X      "<< vec_All_X.at(chi_indexX).at(ck)->straw<<"\t"<<smallestX<<endl;
	vec_tracks.push_back ( vec_All_X.at ( chi_indexX ).at ( ck ) );
    }
    for ( Int_t cl = 0; cl < vec_All_Y.at ( chi_indexY ).size(); cl++ ) {
	//cout<<"       BEST  Y      "<< vec_All_Y.at(chi_indexY).at(cl)->straw<<"\t"<<smallestY<<endl;
	vec_tracks.push_back ( vec_All_Y.at ( chi_indexY ).at ( cl ) );
    }

    //vec_driftTime.clear();


    double sumLeadTime = 0;
    double meanTime = 0;

    for ( Int_t d = 0; d<vec_tracks.size(); d++ ) {
	sumLeadTime += vec_tracks.at ( d )->leadTime;

    }

    meanTime = sumLeadTime/vec_tracks.size();
    //vec_event.push_back ( i );

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


    for ( Int_t tq=0; tq<vec_tracks.size(); tq++ ) {
	//cout<<"TRACKS  : "<<vec_tracks[tq]->layer<<"\t"<<vec_tracks[tq]->straw<<"\t"<<vec_tracks[tq]->channel<<"\t"<<vec_tracks[tq]->leadTime<<"\t"<<meanTime<<"\t"<<(fabs((vec_tracks[tq]->leadTime)-meanTime))<<endl;

	vec_tracks[tq]->drifttime =  max_dt_offset+(meanTime - ( vec_tracks[tq]->leadTime ) ) ;
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
//}

#endif