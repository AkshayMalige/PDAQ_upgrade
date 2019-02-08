#include "PDAQ_Cluster_Finder_Cosy.h"

using namespace std;

int PDAQ_Cluster_Finder_Cosy ( char* intree, int maxEvents, char* outtree) {

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;

    printf("%s\n",outtree);
    PandaSttCal* STT_CAL = new PandaSttCal();



    PandaSttTrack* STT_TRACKS = new PandaSttTrack();
    Stt_Track_Event* stt_event = & ( STT_TRACKS->stt_track_can );

    TFile inFile(intree);
    TTree* tree = (TTree*)inFile.Get("PDAQ_tree");
    if (!tree)
    {
        std::cerr << "Tree PDAQ_tree was not found\n";
        std::exit(1);
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
    TFile* Ttree = new TFile (outtree, "RECREATE" );
    TTree* PDAQ_tree =  new TTree ( "PDAQ_tree", "PDAQ_tree" );
    PDAQ_tree->Branch ( "STT_TRACKS", "PandaSttTrack", &STT_TRACKS, 64000, 99 );

    TH1F* h_STT_EMC_td1 = new TH1F ( "h_STT_EMC_td1", "h_STT_EMC_td1", 1000, 0, 1000 );
    TH1F* h_STT_EMC_td2 = new TH1F ( "h_STT_EMC_td2", "h_STT_EMC_td2", 1000, 0, 1000 );
    TH1F* h_STT_EMC_td3 = new TH1F ( "h_STT_EMC_td3", "h_STT_EMC_td3", 1000, 0, 1000 );
    TH1F* h_STT_EMC_td4 = new TH1F ( "h_STT_EMC_td4", "h_STT_EMC_td4", 1000, 0, 1000 );
    TH2F* h_XvsZ = new TH2F ( "h_XvsZ", "h_XvsZ", 100, -20, 80, 100, -20, 80 );
    TH2F* h_YvsZ = new TH2F ( "h_YvsZ", "h_YvsZ", 100, 0, 100, 100, 0, 100 );
    TH1F* h_STT_EMC_timeDiff = new TH1F ( "h_STT_EMC_timeDiff", "h_STT_EMC_timeDiff", 600, 100, 700 );

    TH2F* h_Straw_DriftTime = new TH2F ( "h_Straw_DriftTime", "h_Straw_DriftTime", 296, 0, 296, 700, 0, 700 );
    TH2F* h_Fee_DriftTime = new TH2F ( "h_Fee_DriftTime", "h_Fee_DriftTime", 16, 0, 16, 700, 0, 700 );

    TH1F* h_straw_mean_straw = new TH1F ( "h_straw_mean_straw", "h_straw_mean_straw", 800, -100, 700 );

    TH1F* h_Front = new TH1F ( "h_Front", "h_Front", 600, 0, 600 );
    TH1F* h_FrontNO = new TH1F ( "h_FrontNO", "h_FrontNO", 20, 0, 20 );

    TH1F* h_Fee[18];

    Double_t repeat =0;
    Double_t All_repeat =0;

    Int_t iev = ( Int_t ) tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl
    << endl;

    int event_counter =0;
    int All_hit_counter =0;
    int Good_hit_counter =0;
    int Layer_eq_4_counter =0;
    int vec_LT_g7_counter =0;
    int vec_LT_l250_counter =0;
    int final_counter =0; 
    
    float max_lead_time_diff = 0;
    int min_track_hits = 0;
    
     

    MParManager* a = MParManager::instance();
    MFTGeomPar* ftGeomPar = new MFTGeomPar();
    MPar * d;
	
    a->setParamSource("ftparams.txt");
    a->parseSource();
    pm()->addParameterContainer("MFTPar", ftGeomPar);
	
    max_lead_time_diff = ftGeomPar->getLeadtimeWindow();
    min_track_hits = ftGeomPar->getTrackMinSize();
    int stations = ftGeomPar->getModules();
    int max_cluster_intake = ftGeomPar->getClusterLimit();
    int MAX_FT_TOTAL_LAYERS =0;
    //int MAX_FT_TOTAL_LAYERS1 =0;
    
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
    
        
    SttHit* hitOnLayer[MAX_FT_TOTAL_LAYERS][500];

    for ( Int_t i = 0; i < iev; i++ ) {
           
	event_counter++;

        tree->GetEntry ( i );

        std::vector<SttHit*> vec_stthits;
        std::vector<SttHit*> vec_stthits2;

	if (i == maxEvents)
            break;

        if ( i%2 ==0 ) {
            cout << "entry no. " << i << endl;
        }

        memset ( hitOnLayer, 0, MAX_FT_TOTAL_LAYERS*500*sizeof ( SttHit* ) );
        int hitMultOnLayer[MAX_FT_TOTAL_LAYERS];

        for ( int h = 0; h < MAX_FT_TOTAL_LAYERS; h++ ) {
            hitMultOnLayer[h] = 0;
        }


        //Loop over the vector elements//////////////////////////////////////////////////////
        for ( int n = 0; n < STT_CAL->stt_cal.total_cal_NTDCHits; n++ ) 
	{
            SttHit* cal_hit  = ( SttHit* ) STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt ( n ); // retrieve particular hit
            All_hit_counter++;
            // hit on reference channel
            if ( cal_hit->isRef == true ) 
	    {
                //cout<<"Reference Hit  -> "<<cal_hit->isRef<<"\t"<<cal_hit->layer<<"\t"<<cal_hit->module<<"\t"<<cal_hit->fee_channel<<"\t"<<cal_hit->straw<<endl;
            } 
            else 
	    {
                if ( cal_hit->layer == 0 ) 
		{
                    continue;
		}
		
	      hitOnLayer[cal_hit->layer - 1][hitMultOnLayer[cal_hit->layer - 1]] = cal_hit;

	      hitMultOnLayer[cal_hit->layer - 1]++;
	      Good_hit_counter++;	      
	      }
        }

        bool good_layers = true;
        for ( int c = 0; c < MAX_FT_TOTAL_LAYERS; c++ ) {
            if ( hitMultOnLayer[c] != 1 ) {
                good_layers = false;
                break;
            }
        }
        int layerCounter = 0;

        for ( int m = 0; m < MAX_FT_TOTAL_LAYERS; m++ ) {
            if ( hitMultOnLayer[m] > 0 )
                layerCounter++;
        }

        if ( layerCounter == MAX_FT_TOTAL_LAYERS ) {
	  
              Layer_eq_4_counter++;

            std::vector<SttHit*> vec_leadTime;
            int filtercnt = 0;
            int fil_max = max_lead_time_diff;


            for ( int l = 0; l < MAX_FT_TOTAL_LAYERS; l++ ) {
                for ( int h = 0; h < hitMultOnLayer[l]; h++ )

                {
                    vec_leadTime.push_back ( hitOnLayer[l][h] );

                }
            }

            for ( Int_t je=0; je<vec_leadTime.size()-1; je++ ) {

                if ( ( vec_leadTime[je+1]->leadTime ) == ( vec_leadTime[je]->leadTime ) && ( vec_leadTime[je+1]->channel == vec_leadTime[je]->channel ) ) {
                    repeat++;
                }
                All_repeat++;
            }

            std::sort ( vec_leadTime.begin(), vec_leadTime.end(), f_sttHitCompareLeadTime );
            vec_stthits.clear();
   
	    if ( vec_leadTime.size() >= min_track_hits && vec_leadTime.size() < max_cluster_intake)
	    {
	      const int minNumber = min_track_hits;
	      const int maxDifference = max_lead_time_diff - 1;
	      int currentNumber = 0;
	      int currentSum    = 0;
	      int numGroups     = 0;
	      double last          = vec_leadTime[0]->leadTime - maxDifference - 1;
	      vec_stthits.clear();

	      for ( int e =0 ; e<vec_leadTime.size();e++ )
	      {
		  if ( vec_leadTime[e]->leadTime - last <= maxDifference )         // Continue with current group
		  {
		    currentNumber++;
		    SttHit* h = vec_leadTime.at(e);
		    vec_stthits.push_back(h);
		  }
		  else                                     // Finish previous group and start anew
		  {
		    if ( currentNumber >= minNumber )     // Previous was a valid group
		    {
			numGroups++;
			PDAQ_Event_Finder(vec_stthits, i, PDAQ_tree,stt_event);
		    }
		    vec_stthits.clear();
		    SttHit* h = vec_leadTime.at(e);
		    vec_stthits.push_back(h);        // Start afresh
		    currentNumber = 1;
		  }
		  last = vec_leadTime[e]->leadTime;
	      }

	      // Deal with leftovers at the end
	      if ( currentNumber >= minNumber )           // Previous was a valid group
	      {
		  numGroups++;
		  PDAQ_Event_Finder(vec_stthits, i, PDAQ_tree,stt_event);
	      }
	    }

	  }

      }//End of loop over events

      PDAQ_tree->Write();
      printf("In_File: %s 	Out_File:  %s\n",intree,outtree);
      return 0;
}

int main ( int argc, char ** argv ) {

    if ( argc >= 3 )
        PDAQ_Cluster_Finder_Cosy ( argv[1], atoi(argv[2]) , argv[3] );
    else return 1;

    return 0;
}

