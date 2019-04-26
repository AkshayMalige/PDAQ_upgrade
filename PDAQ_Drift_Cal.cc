#include "PDAQ_Drift_Cal.h"

using namespace std;

Bool_t PDAQ_Drift_Cal(void)
{
    cout << "Open" << endl;

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;
    PandaSttCal* STT_CAL = 0;
    PandaSttTrack* STT_TRACK = new PandaSttTrack();

    TFile* file = TFile::Open("1750_Tp20_Th20c.root","READ");
    TH1F* DR = new TH1F ( "DR", "DRX", 1000,-0.1, 0.6 );
    TTree* tree = 0;
    file->GetObject("PDAQ_tree", tree);
    if (!tree) {
        std::cerr << "Tree doesn't exists" << std::endl;
        return 1;
    }
    tree->Print();

    tree->SetBranchAddress("STT_TRACKS", &STT_TRACK);
    TFile* ftree = new TFile("DT1750.root", "RECREATE");
    TTree* DR_Tree = new TTree("DR_Tree", "DR_Tree");

    Int_t iev = (Int_t)tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl << endl;

    cout << STT_TRACK->stt_track_can.total_track_NTDCHits << endl;
    
    SttHit Array[256][100000];
    
    for (Int_t i = 0; i < iev; i++)
    {
        tree->GetEntry(i);
        if (i % 100 == 0) {
            cout << i << endl;
        }

        for (int n = 0; n < STT_TRACK->stt_track_can.total_track_NTDCHits; n++)
        {

            SttTrackHit& track_hit = STT_TRACK->stt_track_can.tdc_track_hits[n];
            const std::vector<SttHit> & vec_track_can = track_hit.vec_Track;
          

            for (int t = 0; t < vec_track_can.size(); t++)
            {
                
            }
          
        }
    }
   
    return kTRUE;

    ////////////////////
}

int main() { return PDAQ_Drift_Cal(); }
