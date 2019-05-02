#include "PDAQ_Drift_Cal.h"

using namespace std;

Bool_t PDAQ_Drift_Cal ( void )
{
    cout << "Open" << endl;

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;
    PandaSttCal* STT_CAL = 0;
    PandaSttTrack* STT_TRACK = new PandaSttTrack();

    PandaSttTrack* DT_TRACKS = new PandaSttTrack();
    Stt_Track_Event* stt_event = & ( DT_TRACKS->stt_track_can );

    TFile* file = TFile::Open ( "c.root","READ" );
    TH1F* DR = new TH1F ( "DR", "DRX", 1000,-0.1, 0.6 );
    TTree* tree = 0;
    file->GetObject ( "PDAQ_tree", tree );
    if ( !tree ) {
        std::cerr << "Tree doesn't exists" << std::endl;
        return 1;
    }
    tree->Print();

    tree->SetBranchAddress ( "STT_TRACKS", &STT_TRACK );
    TFile* Ttree = new TFile ( "d.root", "RECREATE" );
    TTree* PDAQ_tree = new TTree ( "PDAQ_tree", "PDAQ_tree" );
    PDAQ_tree->Branch ( "DT_TRACKS", "PandaSttTrack", &DT_TRACKS, 64000, 99 );

    Int_t iev = ( Int_t ) tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl << endl;

    cout << STT_TRACK->stt_track_can.total_track_NTDCHits << endl;

    TH1F* h_Ch_Dt[256];
    TH1F* h_drifttime = new TH1F ( "h_drifttime", "h_drifttime;Drift Time [ns]", 800, -100,700 );
    TH1F* h_Cal_Ch_Dt[256];
    TH1F* h_Dt[256];

    for ( int ch=0; ch<256; ch++ ) {
        h_Cal_Ch_Dt[ch] = new TH1F ( Form ( "Cal_Ch_%d_Dt", ch+1 ) , Form ( "Cal_Ch_%d_Dt", ch+1 ),  800, -100,700 );
        h_Dt[ch] = new TH1F ( Form ( "h_%d_Dt", ch+1 ) , Form ( "h_%d_Dt", ch+1 ),  800, -100,700 );

    }

    h_drifttime = ( TH1F* ) file->Get ( "h_drifttime" );

    double maxbin =0;
    double DTmin =0;

    Int_t maximum = h_drifttime->GetBinContent ( h_drifttime->GetMaximumBin() );
    DTmin = h_drifttime->FindFirstBinAbove ( maximum/10,1 );


    cout<<"New  "<<maximum<<"\t"<<DTmin <<endl;
    std::vector<double> vec_DT_start;

    for ( int chh=0; chh<256; chh++ ) {
        h_Ch_Dt[chh] = ( TH1F* ) file->Get ( Form ( "Ch_%d_Dt",chh+1 ) );

        if ( h_Ch_Dt[chh]->GetBinContent ( h_Ch_Dt[chh]->GetMaximumBin() ) >=20 ) {
            vec_DT_start.push_back ( DTmin - ( h_Ch_Dt[chh]->FindFirstBinAbove ( ( h_Ch_Dt[chh]->GetBinContent ( h_Ch_Dt[chh]->GetMaximumBin() ) ) /10,1 ) ) );
        } else {
            vec_DT_start.push_back ( 0 );
        }

    }

    for(int as=0; as<vec_DT_start.size();as++){

      cout<<vec_DT_start[as]<<endl;
    }




    for ( Int_t i = 0; i < iev; i++ ) {
        tree->GetEntry ( i );
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
                dt_crr = vec_track_can[t].drifttime + vec_DT_start[sq_ch-1];
                a = vec_track_can[t];
                a.drifttime = dt_crr;
                vec_tracks.push_back ( a );
                h_Cal_Ch_Dt[sq_ch]->Fill(a.drifttime);
                h_Dt[sq_ch]->Fill(vec_track_can[t].drifttime);
                //cout<<"Channel "<<sq_ch<<"\t"<<vec_track_can[t].drifttime<<"\t"<<a.drifttime<<endl;
            }
            //cout<<vec_tracks.size()<<endl;
            SttTrackHit& b = stt_event->AddTrackHit();
            b.vec_Track = vec_tracks;
            PDAQ_tree->Fill();
            stt_event->TrackClear();

        }
    }

    //h_Ch_Dt[1]->Draw();
    // h_Ch_Dta->Draw();
    h_drifttime->Write();
    for ( int chh=0; chh<256; chh++ ) {

        h_Cal_Ch_Dt[chh]->Write();
        h_Dt[chh]->Write();

    }
    PDAQ_tree->Write();

    return kTRUE;

    ////////////////////
}

int main()
{
    return PDAQ_Drift_Cal();
}
