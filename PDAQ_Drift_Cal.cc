#include "PDAQ_Drift_Cal.h"

using namespace std;

Bool_t PDAQ_Drift_Cal ( char* intree, char* outtree, int maxEvents )
{
    cout << "Open" << endl;

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;
    PandaSttCal* STT_CAL = 0;
    PandaSttTrack* STT_TRACK = new PandaSttTrack();

    PandaSttTrack* DT_TRACKS = new PandaSttTrack();
    Stt_Track_Event* stt_event = & ( DT_TRACKS->stt_track_can );

    std::vector<double> vec_o_test;
    std::vector<double> vec_x;
    std::vector<double> vec_y;
    std::vector<double> vec_z;
    std::vector<double> vec_layer;
    std::vector<double> vec_module;
    std::vector<double> vec_straw;
    std::vector<double> vec_fee_ch;
    std::vector<double> vec_tdc_ch;

    TFile inFile ( intree );
    TTree* tree = ( TTree* ) inFile.Get ( "PDAQ_tree" );
    if ( !tree ) {
        std::cerr << "Tree PDAQ_tree was not found\n";
        std::exit ( 1 );
    }


    //TFile* file = TFile::Open ( "c.root","READ" );
    TH1F* DR = new TH1F ( "DR", "DRX", 1000,-0.1, 0.6 );
    //TTree* tree = 0;
    //file->GetObject ( "PDAQ_tree", tree );
    if ( !tree ) {
        std::cerr << "Tree doesn't exists" << std::endl;
        return 1;
    }
    tree->Print();
    printf ( "%s\n", outtree );
    tree->SetBranchAddress ( "STT_TRACKS", &STT_TRACK );
    TFile* Ttree = new TFile ( outtree, "RECREATE" );
    TTree* PDAQ_tree = new TTree ( "PDAQ_tree", "PDAQ_tree" );
    PDAQ_tree->Branch ( "DT_TRACKS", "PandaSttTrack", &DT_TRACKS, 64000, 99 );

    Int_t iev = ( Int_t ) tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl << endl;

    cout << STT_TRACK->stt_track_can.total_track_NTDCHits << endl;


    PDAQ_tree->Branch ( "vec_Drifttime", &vec_o_test );
    PDAQ_tree->Branch ( "vec_x", &vec_x );
    PDAQ_tree->Branch ( "vec_y", &vec_y );
    PDAQ_tree->Branch ( "vec_z", &vec_z );
    PDAQ_tree->Branch ( "vec_layer", &vec_layer );
    PDAQ_tree->Branch ( "vec_straw", &vec_straw );

    TH1F* h_Ch_Dt[256];
    TH1F* h_drifttime = new TH1F ( "h_drifttime", "h_drifttime;Drift Time [ns]", 800, -100,700 );
    TH1F* h_Cal_drifttime = new TH1F ( "h_Cal_drifttime", "h_Cal_drifttime;Drift Time [ns]", 800, -100,700 );
    TH1F* h_Cal_Ch_Dt[256];
    TH1F* h_Dt[256];

    for ( int ch=0; ch<256; ch++ ) {
        h_Cal_Ch_Dt[ch] = new TH1F ( Form ( "Cal_Ch_%d_Dt", ch+1 ) , Form ( "Cal_Ch_%d_Dt", ch+1 ),  800, -100,700 );
        h_Dt[ch] = new TH1F ( Form ( "h_%d_Dt", ch+1 ) , Form ( "h_%d_Dt", ch+1 ),  800, -100,700 );

    }

    h_drifttime = ( TH1F* ) inFile.Get ( "h_drifttime" );

    double maxbin =0;
    double DTmin =0;

    Int_t maximum = h_drifttime->GetBinContent ( h_drifttime->GetMaximumBin() );
    DTmin = h_drifttime->FindFirstBinAbove ( maximum/10,1 );


    cout<<"New  "<<maximum<<"\t"<<DTmin <<endl;
    std::vector<double> vec_DT_start;

    for ( int chh=0; chh<256; chh++ ) {
        h_Ch_Dt[chh] = ( TH1F* ) inFile.Get ( Form ( "Ch_%d_Dt",chh+1 ) );

        if ( h_Ch_Dt[chh]->GetBinContent ( h_Ch_Dt[chh]->GetMaximumBin() ) >=20 ) {
            vec_DT_start.push_back ( DTmin - ( h_Ch_Dt[chh]->FindFirstBinAbove ( ( h_Ch_Dt[chh]->GetBinContent ( h_Ch_Dt[chh]->GetMaximumBin() ) ) /10,1 ) ) );
        } else {
            vec_DT_start.push_back ( 0 );
        }

    }

//     for ( int as=0; as<vec_DT_start.size(); as++ ) {
//
//         cout<<vec_DT_start[as]<<endl;
//     }




    for ( Int_t i = 0; i < iev; i++ ) {
        tree->GetEntry ( i );


        if ( i == maxEvents ) {
            break;
        }
        if ( i % 1000 == 0 ) {
            cout << i << endl;
        }

        for ( int n = 0; n < STT_TRACK->stt_track_can.total_track_NTDCHits; n++ ) {

            SttTrackHit& track_hit = STT_TRACK->stt_track_can.tdc_track_hits[n];
            const std::vector<SttHit> & vec_track_can = track_hit.vec_Track;

            SttHit a;
            double dt_crr =0;
            int sq_ch=0;
            std::vector<SttHit> vec_tracks;
            for ( int t = 0; t < vec_track_can.size(); t++ ) {
                sq_ch = ( ( vec_track_can[t].layer-1 ) * 32 ) +vec_track_can[t].straw-1;
                dt_crr = vec_track_can[t].drifttime + vec_DT_start[sq_ch];
                if ( dt_crr>0.0 && dt_crr<=200.0 ) {
                    a = vec_track_can[t];
                    a.drifttime = dt_crr;
                    vec_tracks.push_back ( a );
                    h_Cal_Ch_Dt[sq_ch]->Fill ( a.drifttime );
                    h_Dt[sq_ch]->Fill ( vec_track_can[t].drifttime );

                    h_Cal_drifttime->Fill ( dt_crr );


                    vec_o_test.push_back ( a.drifttime );
                    vec_x.push_back ( vec_track_can[t].x );
                    vec_y.push_back ( vec_track_can[t].y );
                    vec_z.push_back ( vec_track_can[t].z );
                    vec_layer.push_back ( vec_track_can[t].layer );
                    vec_straw.push_back ( vec_track_can[t].straw );
                }
                //cout<<"Channel "<<sq_ch<<"\t"<<vec_track_can[t].drifttime<<"\t"<<a.drifttime<<endl;
            }
            
            
            //cout<<vec_tracks.size()<<endl;
            SttTrackHit& b = stt_event->AddTrackHit();
            b.vec_Track = vec_tracks;
            PDAQ_tree->Fill();
            stt_event->TrackClear();
            vec_o_test.clear();
            vec_x.clear();
            vec_y.clear();
            vec_z.clear();
            vec_straw.clear();

        }
    }

    //h_Ch_Dt[1]->Draw();
    // h_Ch_Dta->Draw();
    h_drifttime->Write();
    h_Cal_drifttime->Write();

    double C =0;
    int xmin=0;
    int xmax=0;
    std::vector<double> vec_drift_radius;
    double a1[201];
    double b1[201];

    C = ( h_Cal_drifttime->GetEntries() ) /0.505;

    for ( int x=0; x<=200; x++ ) {
        TAxis *axis = h_Cal_drifttime->GetXaxis();
        int bmin = axis->FindBin ( 0.0 );
        int bmax = axis->FindBin ( x );
        double drift_radius = ( h_Cal_drifttime->Integral ( bmin,bmax ) ) /C;
        //vec_drift_radius.push_back ( drift_radius );
        cout<<x<<"\t"<<drift_radius<<endl;
        a1[x]=x;
        b1[x]=drift_radius;
    }

    TGraph* gDriftRadius = new TGraph ( 201, a1, b1 );
    gDriftRadius->SetName ( "PDAQ_DR" );

    gDriftRadius->Write();



    for ( int chh=0; chh<256; chh++ ) {

        h_Cal_Ch_Dt[chh]->Write();
        h_Dt[chh]->Write();

    }
    PDAQ_tree->Write();

    return kTRUE;

    ////////////////////
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
            -         PDAQ_Drift_Cal ( argv[1], argv[2], 100000000 );
        } else {
            PDAQ_Drift_Cal ( argv[1], argv[2], atoi ( argv[3] ) );
        }

    else {
        return 1;
    }

    return 0;
}





