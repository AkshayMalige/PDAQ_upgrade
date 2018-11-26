#include "PDAQ_Stt_Calibirator.h"
#include <MLookup.h>
#include <MLookupContainer.h>
#include <MLookupManager.h>
#include <MLookupManager.cc>
                           
using namespace std;

class TestChannel : public MLookupChannel
{
public:
    uint mod, lay, cell, str;

    void setAddress(const char * address) {
        sscanf(address, "%*s %*s %d %d %d\n", &mod, &lay, &cell, &str);
    }

    void print(const char * prefix) {
        printf("%s %d  %d  %d\n", prefix, mod, lay, cell, str);
    }

};

class TestLookupTable : public MLookupTable {
public:
    TestLookupTable(const std::string & container, UInt_t addr_min, UInt_t addr_max, UInt_t channels) :
        MLookupTable(container, addr_min, addr_max, channels) {}

    MLookupChannel * initial() { return new TestChannel(); }
};



bool f_sttHitCompareLeadTime(SttHit* a, SttHit* b)
{
    return (a->leadTime < b->leadTime);
}

int PDAQ_Stt_Calibirator(void)
{
  
    MLookupManager* lm = MLookupManager::instance();
    lm->setSource("lookup.txt");
    lm->parseSource();        
    MLookupTable * t = (MLookupTable*) new TestLookupTable("TestLookup", 0xe100, 0xe202, 49);
 
    MParManager* a = MParManager::instance();
    MFTGeomPar* ftGeomPar = new MFTGeomPar();
    MPar * d;
	
    a->setParamSource("ftparams.txt");
    a->parseSource();
    pm()->addParameterContainer("MFTPar", ftGeomPar);
    d = pm()->getParameterContainer("MFTPar");

    if (!ftGeomPar)
    {
	std::cerr << "Parameter container 'PFTGeomPar' was not obtained!" << std::endl;
    }
    else
    { 
	printf("*ftGeomPar:%i\n", ftGeomPar);
    }

    ftGeomPar->print();
    
    cout<<ftGeomPar->getStrawRadius(1)<<endl;
  
    PandaSubsystemSTT* STT = 0;
    PandaSttCal* CAL = new PandaSttCal();
    Stt_Cal_Event* stt_event = &(CAL->stt_cal);
    SttDetector stt_det;


    TFile file("PDAQ_SttRAW.root", "READ");
    TTree* tree = (TTree*)file.Get("PDAQ_EMC_STT_cluster_analysis");
    if (!tree) {
        std::cerr << "Tree doesn't exists" << std::endl;
        return 1;
    }

    //tree->Print();
    cout<<"tree loading: "<<tree->SetBranchAddress("STT", &STT)<<endl;


    TH1F* h_STT_layers = new TH1F("h_STT_layers", "h_STT_layers",  10, 0, 10);
    TH1F* h_STT_channelsWithoutRef = new TH1F("h_STT_channelsWithoutRef", "h_STT_channelsWithoutRef", 350, 0, 350);
    TH1F* h_STT_tot = new TH1F("h_STT_tot", "h_STT_tot", 500, 0, 500);
    TH1F* h_STT_channels_4layers = new TH1F("h_STT_channels_4layers", "h_STT_channels_4layers", 500, 0, 500);
    TH1F* h_STT_layerhits = new TH1F("h_STT_layerhits", "h_STT_layerhits", 10, 0, 10);
    TH1F* h_STT_layerMult = new TH1F("h_STT_layerMult", "h_STT_layerMult", 17, 0, 17);
    TH1F* h_stt_cosmic = new TH1F("h_stt_cosmic", "h_stt_cosmic", 5, 0, 5);
    TH1F* h_STT_layerhitsmult = new TH1F("h_STT_layerhitsmult", "h_STT_layerhitsmult", 10, 0, 10);
    TH1F* h_STT_counts = new TH1F("h_STT_counts", "h_STT_counts", 10, 0, 10);
    TH2F* h_STT_TOTvsCh = new TH2F("h_STT_TOTvsCh", "h_STT_TOTvsCh", 300, 0, 300, 500, 0, 500);

    Double_t repeat =0;
    Double_t All_repeat =0;
    int good_counter =0;

    Int_t nentries = (Int_t)tree->GetEntries();
    std::cout << nentries << "\n";

    TFile* Ttree = new TFile("PDAQ_Stt_CAL.root", "RECREATE");
    TTree* PDAQ_EMC_STT_cluster_analysis = new TTree("PDAQ_EMC_STT_cluster_analysis", "PDAQ_EMC_STT_cluster_analysis");    
    PDAQ_EMC_STT_cluster_analysis->Branch("STT_CAL", "PandaSttCal", &CAL, 64000, 99);

    Long_t global_cnt = 0;

     for (Long_t e = 0; e < nentries; e++)
    {
    	tree->GetEntry(e);
        int hitsInEvent = 0;
        std::vector<SttHit*> vec_filterLeadTime;
        std::vector<SttHit*> vec_L1;
        std::vector<SttHit*> vec_L2;
        std::vector<SttHit*> vec_L3;
        std::vector<SttHit*> vec_L4;

        if (e % 10000 == 0)
        {
            printf("%ld\n", e);
        }

        SttHit* hitOnLayer[4][500];
        memset(hitOnLayer, 0, 4*500*sizeof(SttHit*));
        int hitMultOnLayer[4];
        int layerCounterr = 0;

        for (int i = 0; i < 4; i++)
        {
            hitMultOnLayer[i] = 0;
        }
	PDAQ_EMC_STT_cluster_analysis->Fill();
	stt_event->CalClear();
        
        //cout<<"^^^^^^^^^totalNTDCHits : "<<STT->stt_raw.totalNTDCHits<<endl<<endl;

        for (int i = 0; i < STT->stt_raw.totalNTDCHits; i++)
        {
	  
            SttRawHit* hit  = (SttRawHit*)STT->stt_raw.tdc_hits->ConstructedAt(i); // retrieve particular hit 
	    TestChannel *tc = (TestChannel*)t->getAddress(hit->tdcid, hit->new_channel);       
   
            SttHit* cal_hit = stt_event->AddCalHit(hit->new_channel);  
	    cal_hit->tdcid = hit->tdcid;
            cal_hit->leadTime = hit->leadTime;
            cal_hit->trailTime = hit->trailTime;
            cal_hit->tot = hit->tot;
            cal_hit->isRef = hit->isRef;       
            cal_hit->layer = tc->lay;
            cal_hit->cell = tc->cell;
	    cal_hit->station= tc->mod;
	    if (tc->cell%2==0){cal_hit->plane=1;}
	    else   {cal_hit->plane = 0;} 
	    
	    
	    if (hit->tot > 0)
	    {
	    good_counter++;
	    }

            // hit on reference channel
            if (cal_hit->isRef == true){}
            else
            {
                h_STT_channelsWithoutRef->Fill(cal_hit->channel);
                h_STT_tot->Fill(cal_hit->tot);
                h_STT_layers->Fill(cal_hit->layer);
                h_STT_counts->Fill(0);
                hitOnLayer[cal_hit->layer - 1][hitMultOnLayer[cal_hit->layer - 1]] = cal_hit;

                hitMultOnLayer[cal_hit->layer - 1]++;
            }

          hitsInEvent++;
        } // end of loop over hits

        bool good_layers = true;
        for (int c = 0; c < 4; c++)
        {
            if (hitMultOnLayer[c] != 1)
            {
                 good_layers = false;
                 break;
            }
        }
        if(good_layers)
            h_STT_layerMult->Fill(1);
        else
            h_STT_layerMult->Fill(0);

        // find the number of layers hit
        int layerCounter = 0;

        for (int i = 0; i < 4; i++)
        {
            if (hitMultOnLayer[i] > 0)
                layerCounter++;
        }
        h_STT_layerhits->Fill(layerCounter);
        h_STT_layerhits->Fill(layerCounter);

        // condition where all 4 layers were hit
        if (layerCounter == 4)
        {
            h_STT_channels_4layers->Fill(hitMultOnLayer[0] + hitMultOnLayer[1] + hitMultOnLayer[2] + hitMultOnLayer[3]);
            h_STT_counts->Fill(1);

            double leadTimesOfHits[2000];
            std::vector<SttHit*> vec_leadTime;
            std::vector<double> sumT;
            int hitsCtr = 0;
            int filtercnt = 0;
            int fil_max = 300;
            int s = 8;
            double ldiff = 0;

            for (int l = 0; l < 4; l++)
            {
                for (int h = 0; h < hitMultOnLayer[l]; h++)

                {
                    vec_leadTime.push_back(hitOnLayer[l][h]);

                    hitsCtr++;

                    if (hitMultOnLayer[l] > 1)
                    {
                        h_STT_layerhitsmult->Fill(h);
                    }
                }
            }

//Check for repeated entries of a hit

            for(Int_t je=0; je<vec_leadTime.size()-1;je++)
            {
                if ((vec_leadTime[je+1]->leadTime) == (vec_leadTime[je]->leadTime) && (vec_leadTime[je+1]->channel == vec_leadTime[je]->channel))
                {
                    repeat++;
                 }
                 All_repeat++;
             }

            std::sort(vec_leadTime.begin(), vec_leadTime.end(), f_sttHitCompareLeadTime);
            vec_filterLeadTime.clear();

            if (vec_leadTime.size() > 7)
            {
                h_STT_counts->Fill(2);

                for (int v = 0; v < vec_leadTime.size(); v++)
                { // iterate over collected and sorted hits
                    vec_filterLeadTime.clear();
                    filtercnt = 1;

                    SttHit* h = vec_leadTime.at(v);
                    vec_filterLeadTime.push_back(h);


                    for (int vv = 0; vv < vec_leadTime.size(); vv++)
                    { // check each vs each if they fit into the window
                        if (vv == v)    continue;
                        if ((fabs(vec_leadTime[v]->leadTime - vec_leadTime[vv]->leadTime) < fil_max) )

                        {
                            SttHit* hh = vec_leadTime.at(vv);
                            vec_filterLeadTime.push_back(hh);
                            filtercnt++;
                        }
                        else
                        { // in case the hit is outside the window break and start next iteration
                            // with the next hit
                            break;
                            vec_filterLeadTime.clear();

                        }
                    }
                }
            }
            //  in case we have at least 8 hits within timewindow go with tracking
            if (vec_filterLeadTime.size() > 4)
             {

                for (Int_t h =0; h< vec_filterLeadTime.size(); h++)
                {
                    SttHit* a = stt_event->AddCalHit(vec_filterLeadTime[h]->channel);
                    a->leadTime = vec_filterLeadTime[h]->leadTime;
                    a->trailTime = vec_filterLeadTime[h]->trailTime;
                    a->tot = vec_filterLeadTime[h]->tot;
		    a->tdcid = vec_filterLeadTime[h]->tdcid;
                    a->layer = vec_filterLeadTime[h]->layer;
                    a->cell = vec_filterLeadTime[h]->cell;
		    a->station = vec_filterLeadTime[h]->station;
		    a->plane = vec_filterLeadTime[h]->plane; 
		    
		    double pit =ftGeomPar->getStrawPitch(a->station-1);
		    double rad =ftGeomPar->getStrawRadius(a->station-1);
		    double offX = ftGeomPar->getOffsetX(a->station-1,a->layer-1,a->plane);
		    double offY =ftGeomPar->getOffsetY(a->station-1,a->layer-1,a->plane);
		    double offZ =ftGeomPar->getOffsetZ(a->station-1,a->layer-1,a->plane);
  
		    if (a->layer%2==0){a->x=0;}
		    else { a->x=offX + (a->cell*rad);}    
		    if(a->layer%2==0){a->y = offY+ (a->cell*rad);}
		    else {a->y=0;}
		    a->z = offZ; 

                }

            }

        } 


    }// over events

 
    PDAQ_EMC_STT_cluster_analysis->Write();
    cout << "Repeated entries  :"<< repeat<<"/"<<All_repeat<<endl;  
    cout << "Good Hits : "<<good_counter<<endl;
    // if (fp)
    //     fclose(fp);


return 0;
}
