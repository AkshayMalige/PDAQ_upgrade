#include "PDAQ_Stt_Calibirator.h"
#include "TH2F.h"
#include "panda_subsystem_sci.h"
#include <MLookup.h>
#include <MLookupContainer.h>
#include <MLookupManager.cc>
#include <MLookupManager.h>

using namespace std;

class TestChannel : public MLookupChannel
{
public:
    uint mod, lay, straw, cell;

    void setAddress ( const char* address ) {
        sscanf ( address, "%*s %*s %d %d %d\n", &mod, &lay, &straw, &cell );
    }

    void print ( const char* prefix ) {
        printf ( "%s %d  %d  %d\n", prefix, mod, lay, straw, cell );
    }
};

class TestLookupTable : public MLookupTable
{
public:
    TestLookupTable ( const std::string& container, UInt_t addr_min,
                      UInt_t addr_max, UInt_t channels )
        : MLookupTable ( container, addr_min, addr_max, channels ) {
    }

    MLookupChannel* initial() {
        return new TestChannel();
    }
};

bool f_sttHitCompareLeadTime ( SttHit* a, SttHit* b )
{
    return ( a->leadTime < b->leadTime );
}

//////////////////////////////////         TDC CALIBRATION FORMULA     /////////////////////////////////////////////////////
/////  ( ( ( hit->leadTime- ( tdc_ref[index]-tdc_ref[0] ) ) +tdc_ref[index] ) +trb_diff )- ( tdc_ref[0]+trb_diff );       //
/////                                                                                                                     //
/////   ( ( ( hit->leadTime- ( tdc_ref[index]-tdc_ref[0] ) ) +tdc_ref[index] ) +trb_diff ) =  A                           //
/////   ( tdc_ref[0]+trb_diff ) = B                                                                                       //
/////   A - B = hit time  - ref time                                                                                      //
/////   ( tdc_ref[index]-tdc_ref[0] ) -> To shift the TDC ref time to the TRB ref time                                    //
/////   tdc_ref[index] -> To get hit time from leadtime                                                                   //
/////   +trb_diff -> To shift the hit time to the reference TRB time                                                      //
/////   ( tdc_ref[0]+trb_diff ) -> Diff btwn the TRB's                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int PDAQ_Stt_Calibirator ( char* intree, char* outtree, int maxEvents )
{

    MLookupManager* lm = MLookupManager::instance();
    lm->setSource ( "lookup.txt" );
    lm->parseSource();
    MLookupTable* t =
        ( MLookupTable* ) new TestLookupTable ( "TestLookup", 0x6000, 0x6500, 49 );

    MParManager* a = MParManager::instance();
    MFTGeomPar* ftGeomPar = new MFTGeomPar();
    MPar* d;

    a->setParamSource ( "ftparams.txt" );
    a->parseSource();
    pm()->addParameterContainer ( "MFTPar", ftGeomPar );
    d = pm()->getParameterContainer ( "MFTPar" );

    if ( !ftGeomPar ) {
        std::cerr << "Parameter container 'PFTGeomPar' was not obtained!"
                  << std::endl;
    } else {
        printf ( "*ftGeomPar:%i\n", ftGeomPar );
    }

    ftGeomPar->print();
    //     double ax=ftGeomPar->getOffsetX ( 0, 0,  0 );
    //     double tx = ftGeomPar->getStrawRadius ( 0 );
    //     double bx=ftGeomPar->getOffsetX ( 0, 0,  1 );
    //     double cx=ftGeomPar->getOffsetX ( 0, 3,  1 );
    //     for ( int tr=1; tr<33; tr++ ) {
    //         double ay = ( ax + ( tr * tx ) );
    //         cout<<ax<<"\t"<<ay <<endl;
    //
    //     }
    // double offX = ftGeomPar->getOffsetX ( 0, 0,  0 );

    PandaSubsystemSTT* STT = 0;
    PandaSttCal* CAL = new PandaSttCal();
    Stt_Cal_Event* stt_event = & ( CAL->stt_cal );
    SttDetector stt_det;
    PandaSubsystemSCI* SCI = 0;
    PandaSubsystemSCI* SCI_CAL = new PandaSubsystemSCI();
    SciEvent* sci_event = & ( SCI_CAL->sci_raw );

    //     TFile file("Raw3b.root", "READ");
    //     TTree* tree = (TTree*)file.Get("PDAQ_tree");
    //     if (!tree) {
    //         std::cerr << "Tree doesn't exists" << std::endl;
    //         return 1;
    //     }

    TFile inFile ( intree );
    TTree* tree = ( TTree* ) inFile.Get ( "PDAQ_tree" );
    if ( !tree ) {
        std::cerr << "Tree PDAQ_tree was not found\n";
        std::exit ( 1 );
    }

    // tree->Print();
    cout << "tree loading: " << tree->SetBranchAddress ( "STT", &STT ) << endl;
    cout << "tree loading: " << tree->SetBranchAddress ( "SCI", &SCI ) << endl;

    Double_t repeat = 0;
    Double_t All_repeat = 0;
    int good_counter = 0;
    double percentage = 0.0;

    Int_t nentries = ( Int_t ) tree->GetEntries();
    std::cout << nentries << "\n";

    TFile* Ttree = new TFile ( outtree, "RECREATE" );
    TTree* PDAQ_tree = new TTree ( "PDAQ_tree", "PDAQ_tree" );
    PDAQ_tree->Branch ( "STT_CAL", "PandaSttCal", &CAL, 64000, 99 );
    PDAQ_tree->Branch ( "SCI_CAL", "PandaSubsystemSCI", &SCI_CAL, 64000, 99 );

    TH1F* h_TRB_ref_diff =
        new TH1F ( "h_TRB_ref_diff", "h_TRB_ref_diff;Time diff [ns]", 600000,
                   -100000, 1000000 );
    TH1F* h_LMult = new TH1F("h_LMult", "h_LMult;Multiplicity", 10, 0, 10);


    Long_t global_cnt = 0;
    // TH2F* h_TotVsChannel = new TH2F("TotVsChannel", "TotVsChannel", 1000,
    // -100, 700, 300, 0, 300);

    Float_t scint_offset = ftGeomPar->getDTOffset();
    double rad = ftGeomPar->getStrawRadius ( 0 );

    for ( Long_t e = 0; e < nentries; e++ ) {
        tree->GetEntry ( e );
        int hitsInEvent = 0;

        double trb_diff = 0;
        double tdc0 = 0;
        double tdc1 = 0;
        double tdc2 = 0;
        double tdc3 = 0;
        double tdc4 = 0;
        double tdc5 = 0;
        double tdc6 = 0;

        UInt_t tdc[7] = {0x6400, 0x6410, 0x6411, 0x6420, 0x6430, 0x6431, 0x6500 };
        double tdc_ref[7];
        for ( int a = 0; a < 7; a++ ) {
            tdc_ref[a] = 0;
        }

        std::vector<SttHit*> vec_filterLeadTime;

        // percentage = (e*100)/nentries;

        if ( e % 2 == 0 ) {
            printf ( "%d\n", e );
        }
        if ( e == maxEvents ) {
            break;
        }

       // printf("Event ::   %i\n",e);
        PDAQ_tree->Fill();
        stt_event->CalClear();
        sci_event->Clear();

///Loop to get TDC Ref's
        for ( int s = 0; s < STT->stt_raw.totalNTDCHits; s++ ) {
            SttRawHit* s_hit = ( SttRawHit* ) STT->stt_raw.tdc_hits->ConstructedAt ( s );
            int tdc_index = 0;

            for ( int t = 0; t < 7; t++ ) {
                if ( s_hit->tdcid == tdc[t] && s_hit->channel == 0 ) {
                    tdc_index = t;
                    tdc_ref[tdc_index] = s_hit->leadTime;
                    //printf ( "%i  %x  %x\n",t,tdc[t],s_hit->tdcid );
                }
            }
        }

        //printf ( "%lf  %lf  %lf\n",tdc_ref[0],tdc_ref[3],tdc_ref[0]- tdc_ref[3] );



//Loop over Scintillator data
        for ( int s = 0; s < SCI->sci_raw.totalNTDCHits; s++ ) {
            SciHit* scihit = ( SciHit* ) SCI->sci_raw.adc_hits->ConstructedAt ( s );
            if ( scihit->channel==0 ) {
                tdc_ref[6]=scihit->leadTime;
            }
        }

        trb_diff = tdc_ref[6]- tdc_ref[3];
        h_TRB_ref_diff->Fill ( trb_diff );

        for ( int s = 0; s < SCI->sci_raw.totalNTDCHits; s++ ) {
            SciHit* scihit = ( SciHit* ) SCI->sci_raw.adc_hits->ConstructedAt ( s );
            if ( !scihit->channel==0 ) {
                SciHit* shit = sci_event->AddSciHit();
                shit->tdcid = scihit->tdcid;
                shit->channel = scihit->channel;
                shit->isRef = scihit->isRef;
                shit->leadTime = ( ( scihit->leadTime - tdc_ref[6] ) - scint_offset );
                shit->trailTime = ( ( scihit->trailTime- tdc_ref[6] ) - scint_offset );
//                 printf("%x , %lf \n",shit->tdcid,shit->leadTime );
            }

        }

        //printf ( "ref 1:%lf  ref2: %lf  diff: %lf \n",tdc0,tdc3,trb_diff );

        int layers_hit[8] = { 0, 0, 0, 0, 0, 0, 0, 0};
        
//Loop over straw hits
        for ( int i = 0; i < STT->stt_raw.totalNTDCHits; i++ ) {
            SttRawHit* hit = ( SttRawHit* ) STT->stt_raw.tdc_hits->ConstructedAt ( i ); // retrieve particular hit
            // printf("%x \n",hit->tdcid);
//printf ( "tdcid %x, channel %i\n",hit->tdcid,hit->channel );

            if ( hit->channel == 0 || hit->tdcid==0xe103 ) {

            } else {
                
                
                //printf ( "tdcid %x, channel %i\n",hit->tdcid,hit->channel );
                TestChannel* tc = 0;
                int index = 0;
                double ref_diff=0;
                tc = ( TestChannel* ) t->getAddress ( hit->tdcid, hit->channel );
                if ( tc == 0 ) {
                } else if ( tc->mod == 1 ) {
                    for ( int j = 0; j < 7; j++ ) {
                        if ( hit->tdcid == tdc[j] ) {
                            index = j;
                        }
                    }
                    SttHit* cal_hit = stt_event->AddCalHit ( hit->channel );
                    cal_hit->tdcid = hit->tdcid;
                    cal_hit->tot = hit->tot;
                    cal_hit->isRef = hit->isRef;
                    cal_hit->layer = tc->lay;
                    cal_hit->straw = tc->straw;
                    cal_hit->station = tc->mod;
                    cal_hit->trigger_no = hit->trigger_no;
                    
                    if(hit->trigger_no == 0x3ae8c9cc) printf("lay :%i str :%i\n",cal_hit->layer,cal_hit->straw);
                    
             //      if(hit->trigger_no == 0x3ae8a294) printf("Trig : %x\n",hit->trigger_no);
// if (hit->trigger_no == 0x7c6c78ad)  81077175
//  if (hit->trigger_no == 0x81077175)
// {
//     printf("Evnt : %li  , TDC: %x , lay : %i, ch : %i,str:%i \n",e, cal_hit->tdcid,cal_hit->layer,hit->channel,cal_hit->straw);
// }
                    

                    if ( index < 3 ) {
                        ref_diff = ( tdc_ref[6] - tdc_ref[index] ) - ( tdc_ref[3] - tdc_ref[0] );
                    } else {
                        ref_diff = tdc_ref[6] - tdc_ref[index];
                    }
                    cal_hit->leadTime = ( hit->leadTime - tdc_ref[index] ) + ref_diff;
                    cal_hit->trailTime = ( hit->trailTime - tdc_ref[index] ) + ref_diff;

                    //printf("index %i, ref_diff %lf  %lf  %lf  %lf\n",index,ref_diff,hit->leadTime, tdc_ref[index], hit->leadTime - tdc_ref[index] );

                    //printf ( "TRB1ref: %lf, TRB2ref:%lf refdiff:%lf  CLT:%lf  id:%x  L: %i  S: %i\n",tdc_ref[0],tdc_ref[3],tdc_ref[0]-tdc_ref[3],cal_hit->leadTime,cal_hit->tdcid,cal_hit->layer,cal_hit->straw );
                    //cout<<e<<endl;
                    //if (cal_hit->tdcid==0x6411 && cal_hit->layer==4 && cal_hit->straw==11 && cal_hit->leadTime > -6230 && cal_hit->leadTime < -6218){
                    
                 //     printf("Evnt : %li  , TDC: %x , lay : %i, ch : %i,str:%i, lt:%lf \n",e, cal_hit->tdcid,cal_hit->layer,hit->channel,cal_hit->straw,cal_hit->leadTime);
                      //exit;
                    //}
                    
                    layers_hit[tc->lay - 1]++;

                    if ( tc->lay ==1 ||tc->lay ==3||tc->lay ==5||tc->lay ==7 ) {
                        if ( tc->straw%2 ==0 ) {
                            cal_hit->plane =0;
                        } else {
                            cal_hit->plane =1;
                        }
                    } else if ( tc->lay ==2 ||tc->lay ==4||tc->lay ==6||tc->lay ==8 ) {
                        if ( tc->straw%2 ==0 ) {
                            cal_hit->plane =1;
                        } else {
                            cal_hit->plane =0;
                        }

                    }

                    if ( hit->tot > 0 ) {
                        good_counter++;
                    }
                    // cout<<"check 1"<<endl;

                    if ( cal_hit->isRef == false && cal_hit->tdcid != 0xe103 ) {


                        double offX = ftGeomPar->getOffsetX ( cal_hit->station - 1, cal_hit->layer - 1, 0 );
                        double offY = ftGeomPar->getOffsetY ( cal_hit->station - 1, cal_hit->layer - 1,
                                                              cal_hit->plane );
                        double offZ = ftGeomPar->getOffsetZ ( cal_hit->station - 1, cal_hit->layer - 1,
                                                              cal_hit->plane );

                        if ( offX == 0 ) {
                            cal_hit->x = 0;
                        } else {
                            cal_hit->x = ( offX + ( cal_hit->straw * rad ) );
                        }
                        if ( offY == 0 ) {
                            cal_hit->y = 0;
                        } else {
                            cal_hit->y = ( offY + ( cal_hit->straw * rad ) );
                        }

                        cal_hit->z = offZ;

                        // printf("TDC: %x Ch: %i Stn: %i Lay: %i Straw: %i Pln:
                        // %i X: %.3f Y: %.3f Z: %.3f \n",
                        // cal_hit->tdcid,cal_hit->channel,cal_hit->station,cal_hit->layer,cal_hit->straw,cal_hit->plane,cal_hit->x,cal_hit->y,cal_hit->z);
                    }

                    else {
                        cal_hit->x = 0;
                        cal_hit->y = 0;
                        cal_hit->z = 0;
                    }
                }
            }
        }
        
       // printf("%d: ", e);
        int mult =0;
        for (int ll =0 ;ll < 8; ll++)
        {
            if (layers_hit[ll] > 0)
            {
                mult++;
            }
           // printf("%d\t", layers_hit[ll]);
        }
        if (mult > 0){
        h_LMult->Fill(mult);}
//if (mult > 0 && mult < 8)cout<<e<<endl;       
//         printf("\n");
        // cout<<"***********\n\n"<<endl;

    } // over events
    h_TRB_ref_diff->Write();
    h_LMult->Write();
    PDAQ_tree->Write();
    Ttree->Close();
    cout << "Repeated entries  :" << repeat << "/" << All_repeat << endl;
    cout << "Good Hits : " << good_counter << endl;
    printf ( "In_File: %s   Out_File:  %s\n", intree, outtree );

    return 0;
}



