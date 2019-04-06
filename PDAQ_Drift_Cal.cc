#include "PDAQ_Drift_Cal.h"

using namespace std;

Bool_t PDAQ_Drift_Cal(void)
{
    cout << "Open" << endl;

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;
    PandaSttCal* STT_CAL = 0;
    PandaSttTrack* STT_TRACK = new PandaSttTrack();

    TFile* file = TFile::Open("c.root","READ");
  TH1F* DR = new TH1F ( "DR", "DRX", 1000,-0.1, 0.6 );
    TTree* tree = 0;
    file->GetObject("PDAQ_tree", tree);
    if (!tree) {
        std::cerr << "Tree doesn't exists" << std::endl;
        return 1;
    }
    tree->Print();

    std::vector<Double_t> vec_event;
    std::vector<Double_t> vec_test;
    std::vector<Double_t> vec_driftTime;
    std::vector<Double_t> vec_roundoff;
    std::vector<Double_t> vec_occurance;
    std::vector<Double_t> vec_pos_DT;
    std::vector<Double_t> vec_cumsum;
    std::vector<Double_t> vec_drift_radius;
    std::vector<SttHit*> vec_All_tracks;

    std::vector<double> vec_o_test;
    std::vector<double> vec_o_test1;

    std::vector<double> vec_x;
    std::vector<double> vec_y;
    std::vector<double> vec_z;
    std::vector<double> vec_layer;
    std::vector<double> vec_module;
    std::vector<double> vec_straw;
    std::vector<double> vec_fee_ch;
    std::vector<double> vec_tdc_ch;

    int counterofdt = 0;
    int driftTimeCounter2 = 0;

    tree->SetBranchAddress("STT_TRACKS", &STT_TRACK);
    TFile* ftree = new TFile("DT1750.root", "RECREATE");
    TTree* DR_Tree = new TTree("DR_Tree", "DR_Tree");

    DR_Tree->Branch("vec_Drifttime", &vec_o_test);
    DR_Tree->Branch("vec_x", &vec_x);
    DR_Tree->Branch("vec_y", &vec_y);
    DR_Tree->Branch("vec_z", &vec_z);
    DR_Tree->Branch("vec_layer", &vec_layer);
    DR_Tree->Branch("vec_straw", &vec_straw);

    Int_t iev = (Int_t)tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl << endl;

    cout << STT_TRACK->stt_track_can.total_track_NTDCHits << endl;

    for (Int_t i = 0; i < 25000; i++)
    {
        tree->GetEntry(i);
        if (i % 100 == 0) {
            cout << i << endl;
        }

        for (int n = 0; n < STT_TRACK->stt_track_can.total_track_NTDCHits; n++)
        {

//             std::vector<SttHit> vec_track_can;
            // SttHit* cal_hit  =
            // (SttHit*)STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt(n);
            SttTrackHit& track_hit =
                    STT_TRACK->stt_track_can.tdc_track_hits[n];
            const std::vector<SttHit> & vec_track_can = track_hit.vec_Track;

            for (int t = 0; t < vec_track_can.size(); t++)
            {
                vec_o_test1.push_back(vec_track_can[t].drifttime);
                vec_o_test.push_back(vec_track_can[t].drifttime);
                vec_x.push_back(vec_track_can[t].x);
                vec_y.push_back(vec_track_can[t].y);
                vec_z.push_back(vec_track_can[t].z);
                vec_layer.push_back(vec_track_can[t].layer);
                vec_straw.push_back(vec_track_can[t].straw);

                counterofdt++;
            }

            DR_Tree->Fill();
            vec_o_test.clear();
            vec_x.clear();
            vec_y.clear();
            vec_z.clear();
            vec_straw.clear();

          
}
    }
    cout << "Number of drifttimes :  " << counterofdt << endl;

    int d = 0;
    int sum = 0;

    // Only under 450ns.

    for (int r = 0; r < vec_o_test1.size(); r++)
    {
        if (vec_o_test1[r] <= 220) {
            vec_test.push_back(vec_o_test1[r]);
            driftTimeCounter2++;
        }
    }

    std::sort(vec_test.begin(), vec_test.end());
    // cout << "counttts:  " << driftTimeCounter2 << endl;

    // Ignore the decimals of the drifttimes.

    for (int t = 0; t < vec_test.size(); t++)
    {
        // cout << seeme[t] << endl;
        int x = (vec_test[t] * 100) / 100;
        vec_roundoff.push_back(x);
    }

    // Calculate the occurances of each drift time.

    for (int j = 0; j < vec_test.size(); j++)
    {
        // cout << vec1[j] << endl;
        int occ_count = 1;
        int limit = (vec_test.size() - 1);

        while (j < limit && vec_roundoff[j] == vec_roundoff[j + 1])
        {
            occ_count++;
            j++;
        }
        vec_occurance.push_back(occ_count);
        vec_pos_DT.push_back(vec_roundoff[j]);
        // cout<< vec_roundoff[j] << "\t" << occ_count <<  endl;
    }

    // Calculate the cummulative sum (integral) of the occurances.

    for (int m = 0; m < vec_occurance.size(); m++)
    {
        int sum = 0;

        for (int k = 0; k < m + 1; k++)
        {
            sum += vec_occurance[k];
        }

        vec_cumsum.push_back(sum);
    }

    // Calculate the drift radius.

    int max_dr = driftTimeCounter2;
    int dt_range = vec_cumsum.size();

    double a1[dt_range];
    double b1[dt_range];
    double C = 0;
    double R = 0.505;
    double drift_radius = 0;

    C = max_dr / R;

    for (int e = 0; e < dt_range; e++)
    {

        drift_radius = vec_cumsum[e] / C;
        vec_drift_radius.push_back(drift_radius);
    }

    for (int l = 0; l < vec_pos_DT.size(); l++)
    {
         //cout << "Drift Time: "<< vec_pos_DT[l] << "\t" << vec_drift_radius[l];
        // << endl;
  DR->Fill(vec_drift_radius[l]);
        a1[l] = vec_pos_DT[l];
        b1[l] = vec_drift_radius[l];
    }

    TGraph* gDriftRadius = new TGraph(dt_range, a1, b1);
    gDriftRadius->SetName("PDAQ_DR");

    gDriftRadius->Write();
    DR_Tree->Write();
  DR->Write();
    return kTRUE;

    ////////////////////
}

int main() { return PDAQ_Drift_Cal(); }