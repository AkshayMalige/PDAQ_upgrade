#include "PDAQ_Cluster_Finder_Cosy.h"

using namespace std;

int PDAQ_Cluster_Finder_Cosy(char* intree, char* outtree, int maxEvents)
{

    PandaSubsystemSTT* STT_RAW = 0;
    PandaSubsystemSTT* CAL = 0;
    PandaSubsystemSCI* SCI_CAL = 0;

    printf("%s\n", outtree);
    PandaSttCal* STT_CAL = new PandaSttCal();
    PandaSubsystemSCI* SCI_TRACKS = new PandaSubsystemSCI();
    SciEvent* sci_event = &(SCI_TRACKS->sci_raw);

    PandaSttTrack* STT_TRACKS = new PandaSttTrack();
    Stt_Track_Event* stt_event = &(STT_TRACKS->stt_track_can);

    TFile inFile(intree);
    TTree* tree = (TTree*)inFile.Get("PDAQ_tree");
    if (!tree) {
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

    tree->SetBranchAddress("STT_CAL", &STT_CAL);
    tree->SetBranchAddress("SCI_CAL", &SCI_CAL);

    TFile* Ttree = new TFile(outtree, "RECREATE");
    TTree* PDAQ_tree = new TTree("PDAQ_tree", "PDAQ_tree");
    PDAQ_tree->Branch("STT_TRACKS", "PandaSttTrack", &STT_TRACKS, 64000, 99);

    TH1F* h_STT_EMC_td1 =
        new TH1F("h_STT_EMC_td1", "h_STT_EMC_td1", 1000, 0, 1000);
    TH1F* h_STT_EMC_td2 =
        new TH1F("h_STT_EMC_td2", "h_STT_EMC_td2", 1000, 0, 1000);
    TH1F* h_STT_EMC_td3 =
        new TH1F("h_STT_EMC_td3", "h_STT_EMC_td3", 1000, 0, 1000);
    TH1F* h_STT_EMC_td4 =
        new TH1F("h_STT_EMC_td4", "h_STT_EMC_td4", 1000, 0, 1000);
    TH2F* h_XvsZ = new TH2F("h_XvsZ", "h_XvsZ", 100, -20, 80, 100, -20, 80);
    TH2F* h_YvsZ = new TH2F("h_YvsZ", "h_YvsZ", 100, 0, 100, 100, 0, 100);
    TH1F* h_STT_EMC_timeDiff =
        new TH1F("h_STT_EMC_timeDiff", "h_STT_EMC_timeDiff", 600, 100, 700);

    TH2F* h_Straw_DriftTime = new TH2F("h_Straw_DriftTime", "h_Straw_DriftTime",
                                       296, 0, 296, 700, 0, 700);
    TH2F* h_Fee_DriftTime =
        new TH2F("h_Fee_DriftTime", "h_Fee_DriftTime", 16, 0, 16, 700, 0, 700);

    TH1F* h_straw_mean_straw =
        new TH1F("h_straw_mean_straw", "h_straw_mean_straw", 800, -100, 700);

    TH1F* h_STT_Hit_Diff =
        new TH1F("h_STT_Hit_Diff", "h_STT_Hit_Diff", 10000, 0, 10000);
    TH1F* Lh_STT_Hit_Diff[7];

    //     for (int i = 0; i < 8; i++) {
    // 	    Lh_STT_Hit_Diff[i] = new TH1F(Form("Layer%dh_STT_Hit_Diff", i + 1),
    // Form("Layer%dh_STT_Hit_Diff", i + 1), 1000, 0, 1000);
    //
    //     }

    TH1F* h_FrontNO = new TH1F("h_FrontNO", "h_FrontNO", 20, 0, 20);

    TH1F* h_Fee[18];

    Int_t iev = (Int_t)tree->GetEntries();
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

    a->setParamSource("ftparams.txt");
    a->parseSource();
    pm()->addParameterContainer("MFTPar", ftGeomPar);

    max_lead_time_diff = ftGeomPar->getLeadtimeWindow();
    min_track_hits = ftGeomPar->getTrackMinSize();
    int stations = ftGeomPar->getModules();
    int max_cluster_intake = ftGeomPar->getClusterLimit();
    int MAX_FT_TOTAL_LAYERS = 0;
    // int MAX_FT_TOTAL_LAYERS1 =0;

    if (stations > 1) {
        for (int a = 0; a < stations; a++)
        {
            MAX_FT_TOTAL_LAYERS += ftGeomPar->getLayers(a);
        }
    }

    else
    {
        MAX_FT_TOTAL_LAYERS = ftGeomPar->getLayers(0);
    }

    double repeat = 0;
    double All_repeat = 0;

    SttHit* hitOnLayer[MAX_FT_TOTAL_LAYERS][500];

    for (Int_t i = 0; i < iev; i++)
    {

        event_counter++;

        tree->GetEntry(i);

        std::vector<SttHit*> vec_stthits;
        std::vector<SttHit*> vec_stthits2;

        if (i == maxEvents) break;

        if (i % 1000 == 0) {
            cout << "entry no. " << i << endl;
        }

        memset(hitOnLayer, 0, MAX_FT_TOTAL_LAYERS * 500 * sizeof(SttHit*));
        int hitMultOnLayer[MAX_FT_TOTAL_LAYERS];

        for (int h = 0; h < MAX_FT_TOTAL_LAYERS; h++)
        {
            hitMultOnLayer[h] = 0;
        }
        // cout<<"Sizes
        // "<<STT_CAL->stt_cal.total_cal_NTDCHits<<"\t"<<SCI_CAL->sci_raw.totalNTDCHits<<"\t"<<t->channel<<endl;

        // Loop over the vector
        // elements//////////////////////////////////////////////////////
        for (int n = 0; n < STT_CAL->stt_cal.total_cal_NTDCHits; n++)
        {
            SttHit* cal_hit =
                (SttHit*)STT_CAL->stt_cal.tdc_cal_hits->ConstructedAt(
                    n); // retrieve particular hit
            All_hit_counter++;

            // hit on reference channel
            if (cal_hit->isRef == true) {
                // cout<<"Reference Hit  ->
                // "<<cal_hit->isRef<<"\t"<<cal_hit->layer<<"\t"<<cal_hit->module<<"\t"<<cal_hit->fee_channel<<"\t"<<cal_hit->straw<<endl;
            }
            else
            {
                if (cal_hit->layer == 0) {
                    continue;
                }

                hitOnLayer[cal_hit->layer - 1]
                          [hitMultOnLayer[cal_hit->layer - 1]] = cal_hit;

                hitMultOnLayer[cal_hit->layer - 1]++;
                Good_hit_counter++;
            }
        }

        bool good_layers = true;
        for (int c = 0; c < MAX_FT_TOTAL_LAYERS; c++)
        {
            if (hitMultOnLayer[c] != 1) {
                good_layers = false;
                break;
            }
        }
        int layerCounter = 0;

        for (int m = 0; m < MAX_FT_TOTAL_LAYERS; m++)
        {
            if (hitMultOnLayer[m] > 0) layerCounter++;
        }

        if (layerCounter == MAX_FT_TOTAL_LAYERS) {

            Layer_eq_4_counter++;

            std::vector<SttHit*> vec_leadTime;
            int filtercnt = 0;
            int fil_max = max_lead_time_diff;

            for (int l = 0; l < MAX_FT_TOTAL_LAYERS; l++)
            {
                for (int h = 0; h < hitMultOnLayer[l]; h++)

                {
                    vec_leadTime.push_back(hitOnLayer[l][h]);

                    // cout<<"Initial :
                    // "<<hitOnLayer[l][h]->layer<<"\t"<<hitOnLayer[l][h]->channel<<"\t"<<hitOnLayer[l][h]->leadTime<<endl;
                }
            }
            // cout<<"Initial Size: "<<vec_leadTime.size()<<"\t";

            // cout<<"Erased Size: "<<vec_leadTime.size()<<endl;

            // cout<<repeat <<"\t"<<All_repeat<<endl;
            std::sort(vec_leadTime.begin(), vec_leadTime.end(),
                      f_sttHitCompareLeadTime);
            vec_stthits.clear();

            double doublehitdiff = 0;
//             for (Int_t je = 0; je < vec_leadTime.size() - 1; je++)
//             {
//                 if ((vec_leadTime[je + 1]->layer) ==
//                         (vec_leadTime[je]->layer) &&
//                     (vec_leadTime[je + 1]->channel ==
//                      vec_leadTime[je]->channel))
//                 {
//                     doublehitdiff = fabs(vec_leadTime[je]->leadTime -
//                                          vec_leadTime[je + 1]->leadTime);
//                     h_STT_Hit_Diff->Fill(doublehitdiff);
//                     if (doublehitdiff <= 500) {
//                         vec_leadTime.erase(vec_leadTime.begin() + je + 1);
// 			
//                         if (doublehitdiff > 0) {
//                             // cout<<"double "<<endl;
//                             repeat++;
//                         }
//                         // sleep(1);
//                     }
//                 }
//                 All_repeat++;
//             }
            for (Int_t je = 0; je < vec_leadTime.size() - 1; je++)
            {
                if ((vec_leadTime[je + 1]->layer) ==
                        (vec_leadTime[je]->layer) &&
                    (vec_leadTime[je + 1]->channel ==
                     vec_leadTime[je]->channel))
                {
                    doublehitdiff = fabs(vec_leadTime[je]->leadTime -
                                         vec_leadTime[je + 1]->leadTime);
                    h_STT_Hit_Diff->Fill(doublehitdiff);
                    if (doublehitdiff <= 500) {
                        vec_leadTime.erase(vec_leadTime.begin() + je + 1);
			--je;
                        if (doublehitdiff > 0) {
                            // cout<<"double "<<endl;
                            repeat++;
                        }
                        // sleep(1);
                    }
                }
                All_repeat++;
            }
// 	    for (Int_t je = vec_leadTime.size() - 2; je >= 0; --je)
//             {
//                 if ((vec_leadTime[je + 1]->layer) ==
//                         (vec_leadTime[je]->layer) &&
//                     (vec_leadTime[je + 1]->channel ==
//                      vec_leadTime[je]->channel))
//                 {
//                     doublehitdiff = fabs(vec_leadTime[je]->leadTime -
//                                          vec_leadTime[je + 1]->leadTime);
//                     h_STT_Hit_Diff->Fill(doublehitdiff);
//                     if (doublehitdiff <= 500) {
//                         vec_leadTime.erase(vec_leadTime.begin() + je + 1);
// 			
//                         if (doublehitdiff > 0) {
//                             // cout<<"double "<<endl;
//                             repeat++;
//                         }
//                         // sleep(1);
//                     }
//                 }
//                 All_repeat++;
//             }
            for (int ex = 0; ex < vec_leadTime.size(); ex++)
            {
                // cout<<"Final :
                // "<<vec_leadTime[ex]->layer<<"\t"<<vec_leadTime[ex]->channel<<"\t"<<vec_leadTime[ex]->leadTime<<endl;
            }
            if (vec_leadTime.size() >= min_track_hits &&
                vec_leadTime.size() < max_cluster_intake)
            {
                const int minNumber = min_track_hits;
                const int maxDifference = max_lead_time_diff - 1;
                int currentNumber = 0;
                int currentSum = 0;
                int numGroups = 0;
                double last = vec_leadTime[0]->leadTime - maxDifference - 1;
                vec_stthits.clear();

                for (int e = 0; e < vec_leadTime.size(); e++)
                {
                    if (vec_leadTime[e]->leadTime - last <=
                        maxDifference) // Continue with current group
                    {
                        currentNumber++;
                        SttHit* h = vec_leadTime.at(e);
                        vec_stthits.push_back(h);
                    }
                    else // Finish previous group and start anew
                    {
                        if (currentNumber >=
                            minNumber) // Previous was a valid group
                        {
                            numGroups++;
                            PDAQ_Event_Finder(vec_stthits, i, PDAQ_tree,
                                              stt_event, ftGeomPar, SCI_CAL);
                        }
                        vec_stthits.clear();
                        SttHit* h = vec_leadTime.at(e);
                        vec_stthits.push_back(h); // Start afresh
                        currentNumber = 1;
                    }
                    last = vec_leadTime[e]->leadTime;
                }

                // Deal with leftovers at the end
                if (currentNumber >= minNumber) // Previous was a valid group
                {
                    numGroups++;
                    PDAQ_Event_Finder(vec_stthits, i, PDAQ_tree, stt_event,
                                      ftGeomPar, SCI_CAL);
                }
            }
        }

    } // End of loop over events
    h_STT_Hit_Diff->Write();
    PDAQ_tree->Write();
    printf("Total Hits processed : %f       Repeated Hits in an Event : %f\n\n",
           All_repeat, repeat);
    printf("In_File: %s 	Out_File:  %s\n", intree, outtree);
    return 0;
}

int main(int argc, char** argv)
{

    if (argc >= 3)

        if (!argv[3]) {
            printf("\n\nNote : One Million Events will be processed. To change "
                   "add the number of events to be processed after the ouput "
                   "file name.\n");
            // atoi(argv[3]) == 1000;
            sleep(2);
            PDAQ_Cluster_Finder_Cosy(argv[1], argv[2], 1000000);
        }
        else
            PDAQ_Cluster_Finder_Cosy(argv[1], argv[2], atoi(argv[3]));

    else
        return 1;

    return 0;
}
