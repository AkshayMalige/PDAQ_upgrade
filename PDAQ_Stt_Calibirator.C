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

int get_mod(Int_t a)
{

  int b =0;
  if ((a)%32 ==0){b = (floor((a)/32));}
  else {b = (floor((a)/32)+1);}
  return b;

}

int get_fee(Int_t a)
{
  int b =0;
  int c=0;
  if (a % 49 ==0){ c =0;}
  else
  {
    if ((a)%16 ==0){b = (floor((a)/16));}
    else {b = (floor((a)/16)+1);}
      
    if (b%2 == 0) {c =2;}
    else {c = 1;}
  }
    return c;

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



     TFile* Ttree = new TFile("PDAQ_Stt_CAL.root", "RECREATE");
     TTree* PDAQ_EMC_STT_cluster_analysis = new TTree("PDAQ_EMC_STT_cluster_analysis", "PDAQ_EMC_STT_cluster_analysis");
    
    PDAQ_EMC_STT_cluster_analysis->Branch("STT_CAL", "PandaSttCal", &CAL, 64000, 99);
    //PDAQ_EMC_STT_cluster_analysis->Branch("STT_RAW", "PandaSubsystemSTT", &STT, 64000, 99);

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
        
        //cout<<"^^^^^^^^^Check 1 : "<<STT->stt_raw.totalNTDCHits<<endl<<endl;

        for (int i = 0; i < STT->stt_raw.totalNTDCHits; i++)
        {
	  
            SttRawHit* hit  = (SttRawHit*)STT->stt_raw.tdc_hits->ConstructedAt(i); // retrieve particular hit 
	    TestChannel *tc = (TestChannel*)t->getAddress(hit->tdcid, hit->new_channel);       
	    detLoc l = stt_det.GetDetectorLocFromTDCChannel(hit->channel);
	    
	    int a = get_mod(tc->cell);
	    int b = get_fee(tc->cell);
	    
// 	    int fe = 0;
// 	    if (a <=16){fe =1;}
// 	    else if (a >16 && a<33){fe=2;}
// 	    else if (a>32 && a<49){fe=1;}
// 	    else if (a>48 && a<65){fe=2;}
// 	    else {fe=0;}
	    
	    
	    int str =0;
	    
	    if ( tc->cell%49 ==0) {str =0;}
	    else {
	      str = (16- (tc->cell - ((32*(a -1))+(16*(b-1)))))+1;
	    }
    
/*	    if (hit->channel%49){
	    printf("New : Layer %i, Module %i, Cell %i |||| Old : Layer %i, Module %i, Cell %i\n",
	    tc->lay,tc->mod,tc->cell,la.layer,la.module,(32 * (la.module -1)) + (16 * (la.fee -1)) + (16-(la.channel_no-1)));}	*/    
            SttHit* cal_hit = stt_event->AddCalHit(hit->new_channel);            
            cal_hit->leadTime = hit->leadTime;
            cal_hit->trailTime = hit->trailTime;
            cal_hit->tot = hit->tot;
            cal_hit->isRef = hit->isRef;       
            cal_hit->layer = tc->lay;
            cal_hit->module = get_mod(tc->cell);
            cal_hit->fee = get_fee(tc->cell);
            cal_hit->fee_channel = str;
            cal_hit->cell = tc->cell;

// 	    cal_hit->layer = l.layer;
//             cal_hit->module = l.module;
//             cal_hit->fee = l.fee;
//             cal_hit->fee_channel = l.channel_no;
//             cal_hit->cell = (32 * (l.module -1)) + (16 * (l.fee -1)) + (16-(l.channel_no-1));

// 	      if (hit->channel%49){
// 	     printf("old fee : {%i, %i, %i %i}  new fee : {%i, %i, %i %i %i}\n",l.layer, l.module, l.fee, l.channel_no, tc->lay, get_mod(tc->cell), get_fee(tc->cell),str,str2 );
// 		//cout<<l.module<<"\t"<<get_mod(tc->cell)<<"\t"<<l.fee<<"\t"<<get_fee(tc->cell)<<endl;
// 	    }
	    
	    
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
                    //printf("{ LT,l,m,c,cf=%d,%d,%2d,%2d}{ LT1,l1,m1,c1,cf1=%d,%d,%2d,%2d}", vec_leadTime[je]->layer,vec_leadTime[je]->module,vec_leadTime[je]->fee,vec_leadTime[je]->fee_channel,vec_leadTime[je+1]->layer,vec_leadTime[je+1]->module,vec_leadTime[je+1]->fee,vec_leadTime[je+1]->fee_channel);
                    //cout<<vec_leadTime[je]->leadTime<<" "<<vec_leadTime[je]->tot<<"\t"<<vec_leadTime[je+1]->leadTime<<" "<<vec_leadTime[je]->tot<<endl;    
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
                    a->module = vec_filterLeadTime[h]->module;
                    a->fee = vec_filterLeadTime[h]->fee;
                    a->fee_channel = vec_filterLeadTime[h]->fee_channel;
                    a->cell = vec_filterLeadTime[h]->cell;
		    
		   // TestChannel *tx = (TestChannel*)t->getAddress(a->tdcid, a->fee_channel);     
		    
		    int plane =0;
		    int stn =0;
		    int str_act =0;
		    
		    if ((vec_filterLeadTime[h]->fee_channel)%2==0)
		    {plane =1;}		
		      else {plane =2;}
		    if (a->layer <=2){stn =1;}
		    else {stn =2;}
	    
		    if ( a->cell%49 ==0) {str_act =0;}
		    else {
		      str_act = (16- (a->cell - ((32*(a->module -1))+(16*(a->fee-1)))))+1;
		    }	     

  double pit =ftGeomPar->getStrawPitch(stn);
  double rad =ftGeomPar->getStrawRadius(1);
  double offX = ftGeomPar->getOffsetX(stn,a->layer,plane);
  double offY =ftGeomPar->getOffsetY(stn,a->layer,plane);
  double offZ =ftGeomPar->getOffsetZ(stn,a->layer,plane);
//cout<<str_act<<"\t"<<rad<<endl;

                    //cout<<a->leadTime<<"\t"<<a->trailTime<<"\t"<<a->tot<<"\t"<<a->layer<<"\t"<<a->module<<"\t"<<a->fee<<"\t"<<a->fee_channel<<"\t"<<a->cell<<endl;

                    if (a->layer ==1){
//                      a->x = x_offset + (((vec_filterLeadTime[h]->module)-1) * mod_width) + ((vec_filterLeadTime[h]->fee) * mod_width) - ((vec_filterLeadTime[h]->fee_channel)*ss);
		      a->x =offX + (str_act*rad*2);
		      a->y = 0;

                         if ((a->fee_channel)%2==0){
                             a->z = 1.01;
			     cout<<a->module<<"\t"<<a->fee<<"\t"<<a->fee_channel<<"\t"<<str_act<<"\t"<<a->cell<<"\t"<<plane<<endl;
			     printf("Old : {X %0.3f, Y %0.3f, Z %0.3f }  ||| New : {offX %0.3f, offY %0.3f, offZ %0.3f}\n",a->x,a->y,a->z,offX,offY,offZ);

                         }

                         else {
                               a->z = 0;
			      //printf("Old : {X %0.3f, Y %0.3f, Z %0.3f }  ||| New : {offX %0.3f, offY %0.3f, offZ %0.3f}\n",a->x,a->y,a->z,offX,offY,offZ);

                         }

                    }


                    else if (vec_filterLeadTime[h]->layer ==2){
                     a->y = y2_offset + (((vec_filterLeadTime[h]->module)-1) * mod_width) + ((vec_filterLeadTime[h]->fee) * mod_width) - ((vec_filterLeadTime[h]->fee_channel)*ss);
                     a->x = 0;

                     if ((vec_filterLeadTime[h]->fee_channel)%2==0){
                         a->z = 12.0;
			cout<<a->module<<"\t"<<a->fee<<"\t"<<a->fee_channel<<"\t"<<str_act<<"\t"<<a->cell<<"\t"<<plane<<endl;

			printf("Old : {X %0.3f, Y %0.3f, Z %0.3f }  ||| New : {offX %0.3f, offY %0.3f, offZ %0.3f}\n",a->x,a->y,a->z,offX,offY,offZ);

                     }

                     else {
                           a->z = 13.01;
			   			      printf("Old : {X %0.3f, Y %0.3f, Z %0.3f }  ||| New : {offX %0.3f, offY %0.3f, offZ %0.3f}\n",a->x,a->y,a->z,offX,offY,offZ);

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