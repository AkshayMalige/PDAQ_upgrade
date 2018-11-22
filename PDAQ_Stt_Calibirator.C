#include "PDAQ_Stt_Calibirator.h"

                           
using namespace std;

bool f_sttHitCompareLeadTime(SttHit* a, SttHit* b)
{
    return (a->leadTime < b->leadTime);
}

int PDAQ_Stt_Calibirator(void)
{

PandaSubsystemSTT* STT = 0;
//PandaSubsystemSTT* raw_stt = new PandaSubsystemSTT();
//PandaSubsystemSTT* stt = new PandaSubsystemSTT();


PandaSttCal* CAL = new PandaSttCal();

Stt_Cal_Event* stt_event = &(CAL->stt_cal);
//SttEvent* raw_event = &(stt->stt_raw);



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
	
    //tree=(TTree*)Ttree->Get("PDAQ_Events");

    std::vector<Double_t> vec_Stage_DT;
    std::vector<Double_t> vec_Stage_x;
    std::vector<Double_t> vec_Stage_y;
    std::vector<Double_t> vec_Stage_z;
    std::vector<Double_t> vec_Stage_layer;
    std::vector<Double_t> vec_Stage_module;
    std::vector<Double_t> vec_Stage_fee;
    std::vector<Double_t> vec_Stage_cell;
    std::vector<Double_t> vec_Stage_fee_ch;
    std::vector<Double_t> vec_Stage_tdc_ch;
    std::vector<Double_t> vec_Stage_leadtime;
    std::vector<Double_t> vec_Stage_trailtime;

        double x_offset = 13.105;
        double y2_offset = 71.705;
        double y4_offset = 55.305;
        double mod_width = 8.08;
        double ss = 0.505;


   //  TFile* Ttree = new TFile("Drift_Radius_test.root", "RECREATE");
   //  TTree* DR_Tree = new TTree("DR_Tree", "DR_Tree");

     TFile* Ttree = new TFile("PDAQ_Stt_CAL.root", "RECREATE");
     TTree* PDAQ_EMC_STT_cluster_analysis = new TTree("PDAQ_EMC_STT_cluster_analysis", "PDAQ_EMC_STT_cluster_analysis");
    
    PDAQ_EMC_STT_cluster_analysis->Branch("STT_CAL", "PandaSttCal", &CAL, 64000, 99);
    //PDAQ_EMC_STT_cluster_analysis->Branch("STT_RAW", "PandaSubsystemSTT", &STT, 64000, 99);

    Long_t global_cnt = 0;

     for (Long_t e = 0; e < nentries; e++)
    {
    	tree->GetEntry(e);

    	//cout<<STT->stt_raw.totalNTDCHits<<endl;

    	// for (int i = 0; i < STT->stt_raw.totalNTDCHits; i++) {
    	// 	SttHit* hit = (SttHit*)STT->stt_raw.tdc_hits->ConstructedAt(i);
		//cout<<hit->channel<<endl;

        int hitsInEvent = 0;

        //bool emcCosmic = false;
        bool sttCosmic = false;
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
            //raw_event->clear();
        
        //cout<<"^^^^^^^^^Check 1 : "<<STT->stt_raw.totalNTDCHits<<endl<<endl;

        for (int i = 0; i < STT->stt_raw.totalNTDCHits; i++)
        {
            SttRawHit* hit  = (SttRawHit*)STT->stt_raw.tdc_hits->ConstructedAt(i); // retrieve particular hit
         
         cout<<"^^^^^^^^^Check 2 : "<<hit->channel<<"\t"<<hit->leadTime<<"\t"<<hit->trailTime<<"\t"<<hit->tot<<"\t"<<hit->isRef<<endl;
   
            SttHit* cal_hit = stt_event->AddCalHit(hit->channel);
            
            
            cal_hit->channel = hit->channel;
            cal_hit->leadTime = hit->leadTime;
            cal_hit->trailTime = hit->trailTime;
            cal_hit->tot = hit->tot;
            cal_hit->isRef = hit->isRef;

if (hit->tot > 0)
{
 good_counter++;
 }

            detLoc l = stt_det.GetDetectorLocFromTDCChannel(hit->channel);
       
            cal_hit->layer = l.layer;
            cal_hit->module = l.module;
            cal_hit->fee = l.fee;
            cal_hit->fee_channel = l.channel_no;
            cal_hit->cell = (32 * (l.module -1)) + (16 * (l.fee -1)) + (16-(l.channel_no-1));
            
            //cout<<"^^^^^^^^^Check 3"<<endl;

            int tdc_num = hit->channel / 49;
            // hit on reference channel
            if (cal_hit->isRef == true)
            {
            }
            else
            {
                h_STT_channelsWithoutRef->Fill(cal_hit->channel);
                h_STT_tot->Fill(cal_hit->tot);
                h_STT_layers->Fill(cal_hit->layer);
                h_STT_counts->Fill(0);

                if (cal_hit->layer == 0)
                {
                    printf("ERROR: cal_hit->layer %d %d %d\n", cal_hit->layer, cal_hit->channel, cal_hit->fee);
                    continue;
                }

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


               // cout << "All entries  :  " <<"\t"<<" Lead Time : "<< vec_leadTime[je]->leadTime<<"\t" << " layer : "<< vec_leadTime[je]->layer<<"\t"<< "module : "<< vec_leadTime[je]->module<<"\t"<< " fee : "<< vec_leadTime[je]->fee<<"\t" << "FEE_Ch : "<< vec_leadTime[je]->fee_channel<<"\t"<<endl;
                if ((vec_leadTime[je+1]->leadTime) == (vec_leadTime[je]->leadTime) && (vec_leadTime[je+1]->channel == vec_leadTime[je]->channel))
                {
                    repeat++;
                    printf("{ LT,l,m,c,cf=%d,%d,%2d,%2d}{ LT1,l1,m1,c1,cf1=%d,%d,%2d,%2d}", vec_leadTime[je]->layer,vec_leadTime[je]->module,vec_leadTime[je]->fee,vec_leadTime[je]->fee_channel,vec_leadTime[je+1]->layer,vec_leadTime[je+1]->module,vec_leadTime[je+1]->fee,vec_leadTime[je+1]->fee_channel);
                    cout<<vec_leadTime[je]->leadTime<<" "<<vec_leadTime[je]->tot<<"\t"<<vec_leadTime[je+1]->leadTime<<" "<<vec_leadTime[je]->tot<<endl;    
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

                sttCosmic = true;

                for (Int_t h =0; h< vec_filterLeadTime.size(); h++)
                {
                    SttHit* a = stt_event->AddCalHit(vec_filterLeadTime[h]->channel);
                    a->leadTime = vec_filterLeadTime[h]->leadTime;
                    a->trailTime = vec_filterLeadTime[h]->trailTime;
                    a->tot = vec_filterLeadTime[h]->tot;
                   
                    a->layer = vec_filterLeadTime[h]->layer;
                    a->module = vec_filterLeadTime[h]->module;
                    a->fee = vec_filterLeadTime[h]->fee;
                    a->fee_channel = vec_filterLeadTime[h]->fee_channel;
                    a->cell = vec_filterLeadTime[h]->cell;

                    cout<<a->leadTime<<"\t"<<a->trailTime<<"\t"<<a->tot<<"\t"<<a->layer<<"\t"<<a->module<<"\t"<<a->fee<<"\t"<<a->fee_channel<<"\t"<<a->cell<<endl;

                    if (vec_filterLeadTime[h]->layer ==1){
                     a->x = x_offset + (((vec_filterLeadTime[h]->module)-1) * mod_width) + ((vec_filterLeadTime[h]->fee) * mod_width) - ((vec_filterLeadTime[h]->fee_channel)*ss);
                     a->y = 0;

                         if ((vec_filterLeadTime[h]->fee_channel)%2==0){
                             a->z = 1.01;
                         }

                         else {
                               a->z = 0;
                         }

                    }


                    else if (vec_filterLeadTime[h]->layer ==2){
                     a->y = y2_offset + (((vec_filterLeadTime[h]->module)-1) * mod_width) + ((vec_filterLeadTime[h]->fee) * mod_width) - ((vec_filterLeadTime[h]->fee_channel)*ss);
                     a->x = 0;

                     if ((vec_filterLeadTime[h]->fee_channel)%2==0){
                         a->z = 12.0;
                     }

                     else {
                           a->z = 13.01;
                     }

                    }


                    else if (vec_filterLeadTime[h]->layer ==3){
                     a->x = x_offset + (((vec_filterLeadTime[h]->module)-1) * mod_width) + ((vec_filterLeadTime[h]->fee) * mod_width) - ((vec_filterLeadTime[h]->fee_channel)*ss);
                     a->y = 0;

                         if ((vec_filterLeadTime[h]->fee_channel)%2==0){
                             a->z = 28.51;
                         }

                         else {
                               a->z = 27.5;
                         }
                    }

                    else if (vec_filterLeadTime[h]->layer ==4){
                     a->y = y4_offset + (((vec_filterLeadTime[h]->module)-1) * mod_width) + ((vec_filterLeadTime[h]->fee) * mod_width) - ((vec_filterLeadTime[h]->fee_channel)*ss);
                     a->x = 0;

                     if ((vec_filterLeadTime[h]->fee_channel)%2==0){
                         a->z = 39.59;
                     }

                     else {
                           a->z = 40.61;
                     }

                    }

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

// int main(int argc, char **argv) {
//     TApplication theApp("App", &argc, argv);


//     MParManager* a = MParManager::instance();
//     MFTGeomPar* ftGeomPar = new MFTGeomPar();
//     MPar * d;
    
// //Unpacker * u = new Unpacker();

//     a->setParamSource("ftparams.txt");
//     a->parseSource();
//     pm()->addParameterContainer("MFTPar", ftGeomPar);
//     d = pm()->getParameterContainer("MFTPar");

//     if (!ftGeomPar)
//     {
//         std::cerr << "Parameter container 'PFTGeomPar' was not obtained!" << std::endl;
//         //exit(EXIT_FAILURE);
//     }
//     else{ 
//     printf("*ftGeomPar:%i\n", ftGeomPar);
//     }
//     //a->getParameterContainer("MFibersStackGeomPar");
//     //a->print(); 

//     //printf("\n\n Layers from ->getLayers() : %d\n",ftGeomPar->getStraws(1)); 
//     //ftGeomPar->print();
//     //printf("AAAA:\n");

//     ftGeomPar->print();

//     PDAQ_Stt_Calibirator();

//     cout<<"Run Finished"<<endl;
//     theApp.Run();
//     return 0;
// }