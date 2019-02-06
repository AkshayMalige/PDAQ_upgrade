#include "PDAQ_Stt_Calibirator.h"
#include <MLookup.h>
#include <MLookupContainer.h>
#include <MLookupManager.h>
#include <MLookupManager.cc>
                           
using namespace std;

class TestChannel : public MLookupChannel
{
public:
    uint mod, lay, straw, cell;

    void setAddress(const char * address) {
        sscanf(address, "%*s %*s %d %d %d\n", &mod, &lay, &straw, &cell);
    }

    void print(const char * prefix) {
        printf("%s %d  %d  %d\n", prefix, mod, lay, straw, cell);
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

int PDAQ_Stt_Calibirator(char* intree , char* outtree, int maxEvents)
{
  
    MLookupManager* lm = MLookupManager::instance();
    lm->setSource("lookup.txt");
    lm->parseSource();        
    MLookupTable * t = (MLookupTable*) new TestLookupTable("TestLookup", 0x6000, 0xe203, 49);
 
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
    
    //cout<<ftGeomPar->getStrawRadius(1)<<endl;
  
    PandaSubsystemSTT* STT = 0;
    PandaSttCal* CAL = new PandaSttCal();
    Stt_Cal_Event* stt_event = &(CAL->stt_cal);
    SttDetector stt_det;


//     TFile file("Raw3b.root", "READ");
//     TTree* tree = (TTree*)file.Get("PDAQ_tree");
//     if (!tree) {
//         std::cerr << "Tree doesn't exists" << std::endl;
//         return 1;
//     }
    
    
    TFile inFile(intree);
    TTree* tree = (TTree*)inFile.Get("PDAQ_tree");
    if (!tree)
    {
        std::cerr << "Tree PDAQ_tree was not found\n";
        std::exit(1);
    }

    //tree->Print();
    cout<<"tree loading: "<<tree->SetBranchAddress("STT", &STT)<<endl;


    Double_t repeat =0;
    Double_t All_repeat =0;
    int good_counter =0;
    double percentage =0.0;

    Int_t nentries = (Int_t)tree->GetEntries();
    std::cout << nentries << "\n";

    TFile* Ttree = new TFile(outtree, "RECREATE");
    TTree* PDAQ_tree = new TTree("PDAQ_tree", "PDAQ_tree");    
    PDAQ_tree->Branch("STT_CAL", "PandaSttCal", &CAL, 64000, 99);

    Long_t global_cnt = 0;

     for (Long_t e = 0; e < nentries; e++)
    {
    	tree->GetEntry(e);
        int hitsInEvent = 0;
        std::vector<SttHit*> vec_filterLeadTime;

        //percentage = (e*100)/nentries;

        if (e % 10000 == 0)
        {
            printf("%d\n",e);
        }
                if (e == maxEvents)
            break;
	    

	PDAQ_tree->Fill();
	stt_event->CalClear();
        
        //cout<<"^^^^^^^^^totalNTDCHits : "<<STT->stt_raw.totalNTDCHits<<endl<<endl;

        for (int i = 0; i < STT->stt_raw.totalNTDCHits; i++)
        { //cout<<endl<<endl;
	//cout<<"check0 "<<endl; 
	//printf("i %i  NTDCHits %i\n",i,STT->stt_raw.totalNTDCHits);
	SttRawHit* hit  = (SttRawHit*)STT->stt_raw.tdc_hits->ConstructedAt(i); // retrieve particular hit 
	printf("tdc : %x ch: %i\n",hit->tdcid,hit->channel);
	
// 	if ((TestChannel*)t->getAddress(hit->tdcid, hit->new_channel)){continue;}
// 	else {
// 	  cout<<"Bad TDC Address"<<endl;}
	 
	//if (hit->tdcid!=0x8200 && hit->tdcid!=0x8100){
//cout<<"Before"<<endl;	  
	    TestChannel *tc = 0;
	    tc = (TestChannel*)t->getAddress(hit->tdcid, hit->channel);   
	    //cout<<"After"<<endl;
	    
	    if (tc == 0) {
	    }
	    else {
	    
	    
	    //tc->print("   address");
            SttHit* cal_hit = stt_event->AddCalHit(hit->channel);  
	    cal_hit->tdcid = hit->tdcid;
            cal_hit->leadTime = hit->leadTime;
            cal_hit->trailTime = hit->trailTime;
            cal_hit->tot = hit->tot;
            cal_hit->isRef = hit->isRef;       
            cal_hit->layer = tc->lay;
            cal_hit->straw = tc->straw;
	    cal_hit->station= tc->mod;
	    if (tc->straw%2==0){cal_hit->plane=1;}
	    else   {cal_hit->plane = 0;} 
	    
	    
	    if (hit->tot > 0)
	    {
	    good_counter++;
	    }
	    	     	   		//cout<<"check 1"<<endl;    

	    if (cal_hit->isRef==false && cal_hit->tdcid != 0xe103)
	    {
	    	     	   		//cout<<"check 2"<<endl;    

	     // printf("TDC: %x Ch: %i Stn: %i Lay: %i Cell: %i Pln: %i X: %.3f Y: %.3f Z: %.3f\n",
	    //  cal_hit->tdcid,cal_hit->channel,cal_hit->station,cal_hit->layer,cal_hit->straw,cal_hit->plane,cal_hit->x,cal_hit->y,cal_hit->z);

	      double pit =ftGeomPar->getStrawPitch(cal_hit->station-1);
	      double rad =ftGeomPar->getStrawRadius(cal_hit->station-1);
	      double offX = ftGeomPar->getOffsetX(cal_hit->station-1,cal_hit->layer-1,cal_hit->plane);
	      double offY =ftGeomPar->getOffsetY(cal_hit->station-1,cal_hit->layer-1,cal_hit->plane);
	      double offZ =ftGeomPar->getOffsetZ(cal_hit->station-1,cal_hit->layer-1,cal_hit->plane);

 
	      if (offX == 0)
	      cal_hit->x= 0;
	      else
		cal_hit->x = (offX + (cal_hit->straw*rad)); 
	      if (offY ==0)
		cal_hit->y=0;
	      else
	      {cal_hit->y = (offY+ (cal_hit->straw*rad));}
	      
	      cal_hit->z = offZ; 
	      
 	      //printf("TDC: %x Ch: %i Stn: %i Lay: %i Cell: %i Pln: %i X: %.3f Y: %.3f Z: %.3f\n",
 	      //cal_hit->tdcid,cal_hit->channel,cal_hit->station,cal_hit->layer,cal_hit->straw,cal_hit->plane,cal_hit->x,cal_hit->y,cal_hit->z);

	    }
	    
	    else 
	    {
	      cal_hit->x=0;   
	      cal_hit->y =0;
	      cal_hit->z =0;
	    }
	    }//}
 	  
	}

    }// over events

 
    PDAQ_tree->Write();
    cout << "Repeated entries  :"<< repeat<<"/"<<All_repeat<<endl;  
    cout << "Good Hits : "<<good_counter<<endl;
    // if (fp)
    //     fclose(fp);
printf("In_File: %s 	Out_File:  %s\n",intree,outtree);

return 0;
}
