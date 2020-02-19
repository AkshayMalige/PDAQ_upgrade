#include "PDAQ_RawDecoder_EMC_STT.h"
#include "panda_subsystem_sci.h"
#include "TH2F.h"


// std::stringstream sstream;
//===================================================================
// Histograms

TH1F* h_stt_tdc_event_sizes;
TH1F* h_stt_tdc_channels;
TH1F* h_stt_tdc_leadTimes;

//===================================================================
// Constants
const int stt_Header_system_id = 0xabcd; // 0xFBBB;
std::map<UInt_t, int> stt_channel_offsets;

//
//===================================================================
// Event structure
UShort_t Header_status;
UInt_t SB_number;
//
//-------------------------------------------------------------------
// STT

std::vector<UInt_t> stt_tdc_event_sizes;
std::vector<SttEvent> stt_tdc_hits;

SttDetector stt_det;

//-------------------------------------------------------------------
PandaSubsystemSB* sb = new PandaSubsystemSB();
PandaSubsystemSTT* stt = new PandaSubsystemSTT();
PandaSubsystemEMC* emc = new PandaSubsystemEMC();
PandaSubsystemSCI* sci = new PandaSubsystemSCI();

//===================================================================
// Zero event
void pd_Zero_Event()
{
    // STT
    stt_tdc_event_sizes.clear();
    stt_tdc_hits.clear();
}

void pd_Reserve_Event ( int Nsize )
{
    // STT
    stt_tdc_event_sizes.reserve ( Nsize );
    stt_tdc_hits.reserve ( Nsize );
}

//===================================================================
#define MAKE_HST_IND(hstname, bins, minx, maxx, ind, tstr)                     \
    {                                                                          \
        std::string tmpstr;                                                    \
        tmpstr = #hstname;                                                     \
        tmpstr += "_%2i";                                                      \
        if (hstname[ind])                                                      \
            hstname[ind]->Reset();                                             \
        else                                                                   \
        {                                                                      \
            sprintf(tstr, tmpstr.c_str(), i + 1);                              \
            hstname[ind] = new TH1F(tstr, tstr, bins, minx, maxx);             \
        }                                                                      \
    }
#define MAKE_HST2_IND(hstname, binsX, minx, maxx, binsY, miny, maxy, ind,      \
                      tstr)                                                    \
    {                                                                          \
        std::string tmpstr;                                                    \
        tmpstr = #hstname;                                                     \
        tmpstr += "_%2i";                                                      \
        if (hstname[ind])                                                      \
            hstname[ind]->Reset();                                             \
        else                                                                   \
        {                                                                      \
            sprintf(tsrt, tmpstr.c_str(), i + 1);                              \
            hstname[ind] =                                                     \
                new TH2F(tsrt, tsrt, binsX, minx, maxx, binsY, miny, maxy);    \
        }                                                                      \
    }
#define MAKE_HST(hstname, bins, minx, maxx)                                    \
    {                                                                          \
        if (hstname)                                                           \
            hstname->Reset();                                                  \
        else                                                                   \
            hstname = new TH1F(#hstname, #hstname, bins, minx, maxx);          \
    }
#define MAKE_HST2(hstname, bins, minx, maxx, binsy, miny, maxy)                \
    {                                                                          \
        if (hstname)                                                           \
            hstname->Reset();                                                  \
        else                                                                   \
            hstname = new TH2F(#hstname, #hstname, bins, minx, maxx, binsy,    \
                               miny, maxy);                                    \
    }

//===================================================================
void pd_init_hst()
{
    MAKE_HST ( h_stt_tdc_event_sizes, 100, 0, 100 );

    MAKE_HST ( h_stt_tdc_channels, 350, 0, 350 );

    MAKE_HST ( h_stt_tdc_leadTimes, 40000, 0, 40000 );
}
#undef MAKE_HST
#undef MAKE_HST2
#undef MAKE_HST_IND
#undef MAKE_HST2_IND

UInt_t readWord ( std::ifstream* in_file )
{
    UInt_t data4 = 0;
    UInt_t word = 0;
    in_file->read ( ( char* ) &data4, 4 );
    word = ( ( data4 & 0xff ) << 24 ) | ( ( ( data4 >> 8 ) & 0xff ) << 16 ) |
           ( ( ( data4 >> 16 ) & 0xff ) << 8 ) | ( ( ( data4 >> 24 ) & 0xff ) );
    return word;
}

//===================================================================
// Reading and decoding input file

void PDAQ_RawDecoder_HADES ( char* in_file_name, char* out_file_name = 0,
                             int maxEvents = 100 )
{
    TH1F* h_tdc_ref = new TH1F ( "h_tdc_ref", "h_tdc_ref;Number of ref hits", 10, 0, 10 );
    TH1F* h_refTimeTRB1 = new TH1F ( "h_refTimeTRB1", "h_refTimeTRB1;Time diff [ns]", 500, -25, 25 );
    TH1F* h_refTimeTRB2 = new TH1F ( "h_refTimeTRB2", "h_refTimeTRB2;Time diff [ns]", 500, -25, 25 );
    TH1F* h_currpt = new TH1F ( "h_currpt", "h_currpt", 3, 0, 3 );

    TH1F* h_refTimeTDC[6];
    double RefTime[8];

    for ( int r =0; r< 6; r++ ) {
        h_refTimeTDC[r] = new TH1F ( Form ( "ref_Time_TDC%d",r ), Form ( "ref_Time_TDC%d",r ),500,-2500,2500 );
    }

    pd_init_hst();
    //---------------------------------------------------------------
    // Open input file
    std::ifstream in_file ( in_file_name );
    if ( !in_file ) {
        std::cerr << "Can not open input file!\n";
        return;
        // std::exit(1);
    }
    //---------------------------------------------------------------

    cout << "In File: "
         << "\t" << in_file_name << "\t"
         << "Out File: "
         << "\t" << out_file_name << endl;

    SttEvent* stt_event = & ( stt->stt_raw );

    EmcEvent* emc_event = & ( emc->emc_raw );

    SciEvent* sci_event = & ( sci->sci_raw );

    // Open output file
    bool use_tree_output = false;
    TFile* ofile;
    TTree* tree;
    if ( out_file_name ) {
        use_tree_output = true;
        //
        ofile = new TFile ( out_file_name, "recreate" );
        tree = new TTree ( "PDAQ_tree", "PDAQ_tree" );
        //.......................

        //.......................
        tree->Branch ( "SB", "PandaSubsystemSB", &sb, 64000, 2 );
        tree->Branch ( "STT", "PandaSubsystemSTT", &stt, 64000, 99 );
        tree->Branch ( "EMC", "PandaSubsystemEMC", &emc, 64000, 99 );
        tree->Branch ( "SCI", "PandaSubsystemSCI", &sci, 64000, 99 );
    } else {
        abort();
    }

    //---------------------------------------------------------------

    // skip first file headers from the EB
    in_file.ignore ( 32 );

    UInt_t N_events = 0;
    UInt_t word = 0;
    UInt_t queue_size = 0;
    UInt_t sub_size = 0;
    UInt_t tdc_size = 0;
    UInt_t sub_id = 0;
    UInt_t tdc_id = 0;
    UInt_t queue_words = 0;
    UInt_t tdc_words = 0;
    bool end_of_queue = false;

    UInt_t trigger_nr = 0;
    UInt_t decoding = 0;
    UInt_t channel_nr = 0;
    UInt_t edge = 0;
    UInt_t fine = 0;
    UInt_t coarse = 0;
    UInt_t epoch = 0;


    double doubleCntr = 0;
    double doubleCntr1 = 0;
    
    int prev_time=0;

    while ( !in_file.eof() ) {

        if ( N_events == maxEvents ) {
            break;
        }

        // skip queue headers
        in_file.read ( ( char* ) &word, 4 );
        queue_size = word / 4;
        //printf("\t Queue: size: %d\n", queue_size);

        if ( queue_size < 20 ) { // skipping empty queues
            in_file.ignore ( ( queue_size - 1 ) * 4 );
        }

        else {
            in_file.ignore ( 28 );

            queue_words = 8;

            end_of_queue = false;


            int tdc_ptr=0;

            while ( !in_file.eof() ) {
                double refTime = 0;
                double lastRise = 0;

                // decode sub headers
                word = readWord ( &in_file ); // sub_size
                sub_size = word / 4;
                // in_file.ignore(4);  // decoding
                word = readWord ( &in_file );
                decoding = word;
                word = readWord ( &in_file ); // sub_id
                sub_id = word & 0xffff;
                // in_file.ignore(4);  // trigger nr
                word = readWord ( &in_file );
                trigger_nr = word;
                //printf("Subevent: id: %x size: %d trg:%x\n", sub_id,sub_size, trigger_nr);

                queue_words += 4;

                // loop over TDCs
                while ( !in_file.eof() ) {
                    
                    word = readWord ( &in_file ); // tdc headers
                    // printf("%x\n", word);
                    tdc_size = ( word >> 16 ) & 0xffff; // tdc size
                    tdc_id = word & 0xffff;
                   // printf("\tTDC: id: %x size: %d\n", tdc_id, tdc_size);

                    queue_words++;
                    tdc_words = 0;


                    if ( tdc_id == 0x5555 ) { // skip control subsub
                        in_file.ignore ( 4 );
                        queue_words++;
                        if ( sub_size % 2 != 0 ) { // skip padding word
                            in_file.ignore ( 4 );
                            queue_words++;
                        }
                       // printf("trailer: queue:%d queue_words:%d\n", queue_size, queue_words);
                        if ( queue_size == queue_words ) {
                            end_of_queue = true;
                        }
                        break;
                    
                    } else {

                        if ( (tdc_id >> 8) != 0x8b) {
                            tdc_ptr++;
                        }

                        // loop over TDC data
                        while ( !in_file.eof() ) {
                            word = readWord ( &in_file );
                            tdc_words++;
                            queue_words++;
                            //printf("%x\n", word);

                            if ( ( ( word >> 28 ) & 0xf ) == 0x8 ) { // hit marker
                                channel_nr = ( ( word >> 22 ) & 0x7f );
                                edge = ( ( word >> 11 ) & 0x1 );
                                fine = ( ( word >> 12 ) & 0x3ff );
                                coarse = ( word & 0x7ff );

                                h_stt_tdc_channels->Fill ( channel_nr );

                                double time =
                                    ( double ) ( ( ( ( unsigned ) epoch ) << 11 ) * 5.0 );
                                time += ( ( coarse * 5. ) - ( fine / 100.0 ) );
                                 //printf ( "%lf %x\n", time, tdc_id );

                                if ( channel_nr == 0 ) { // ref time
                                    refTime = time;
                                    if ( tdc_id == 0x6400 ||tdc_id == 0x6410 ||tdc_id == 0x6411||tdc_id == 0x6420||tdc_id == 0x6430||tdc_id == 0x6431 ) {
                                        SttRawHit* a = stt_event->AddHit ( channel_nr );
                                        a->tdcid = tdc_id;
                                        a->leadTime = time;
                                        a->trailTime = 0;
                                        a->isRef = true;
                                        RefTime[tdc_ptr-1] = refTime;
                                    }


                                    //printf("%d %x\n", tdc_ptr-1, tdc_id);


                                    else if ( tdc_id == 0x6500 ) {
                                        SciHit* s = sci_event->AddSciHit();
                                        s->tdcid = tdc_id;
                                        s->channel = channel_nr;
                                        s->leadTime = time;
                                        s->trailTime = ( time );
                                        s->isRef = true;
                                        RefTime[tdc_ptr-1] = refTime;

                                        //h_tdc_ref->Fill ( 7 );
                                        //printf("\tRef R: %f on channel %d on %x on %x tdcptr: %d\n", s->leadTime, channel_nr,  tdc_id, sub_id, tdc_ptr - 1);
                                    
                                    }
                                    
                                } else {
                                    if ( edge == 1 && (fabs(time-refTime)<100000)) { // rising edge

                                        lastRise = time;
                                        
                                         // if (prev_time- lastRise ==0){
                                          
                                          prev_time = time;

                                        if ( lastRise-refTime <=100000 ) {
                                            h_currpt->Fill ( 1 );
                                          // printf("Good %lf %lf %lf\n",lastRise,refTime,lastRise-refTime);
                                        } else {
                                            h_currpt->Fill ( 2 );
                                           // printf("corrupt %lf  %lf %lf\n",lastRise,refTime,lastRise-refTime);
                                        }
//                                                                                 printf("%lf %lf %lf %x\n", time, refTime ,lastRise, tdc_id);

                                        h_stt_tdc_leadTimes->Fill ( time -                                                     refTime );
                                        if ( tdc_id == 0x6500 && channel_nr==1 ) {
                                            SciHit* s = sci_event->AddSciHit();
                                            s->tdcid = tdc_id;
                                            s->channel = channel_nr;
                                            s->leadTime = lastRise;
                                            s->trailTime = ( time );
                                            s->isRef = false;
                                            //printf("0x6500:  %i\n",channel_nr);
                                        }
                                    } else {
                                        // falling edge
  ///////////////////////(fabs(time-refTime)<100000) -> To reject the corrupted entries from the epoch counter/////////////
                                        if ( lastRise != 0 && (fabs(time-refTime)<100000)) { // only in case
                                            // there was a
                                            // rising to pair
                                            bool doubleHit = false;

                                            /* for ( Int_t ui=0;
                                             ui<stt_event->totalNTDCHits; ui++ )
                                             {
                                                 //cout<<"totalNTDCHits
                                             "<<stt_event->totalNTDCHits<<"\t"<<endl;
                                                 if ( ( ( SttRawHit* )
                                             stt_event->tdc_hits->ConstructedAt
                                             ( ui ) )->channel == channel_nr +
                                             stt_channel_offsets[tdc_id] ) {
                                                     //if ( ( ( SttRawHit* )
                                             stt_event->tdc_hits->ConstructedAt
                                             ( ui ) )->leadTime == lastRise ) {
                                                         //cout<<"Double hit
                                             :"<<endl;
                                                         doubleHit = true;
                                                         doubleCntr++;
                                                     }
                                                 //}

                                                 doubleCntr1++;
                                             }*/
                                            /* if (tdc_id == 0x6500 &&
                                             channel_nr ==1)
                                             {
                                               SciHit* s =
                                             sci_event->AddSciHit();
                                               s->tdcid = tdc_id;
                                               s->channel = channel_nr;
                                                                   s->leadTime =
                                             lastRise;
                                                                   s->trailTime
                                             = ( refTime - time );
                                               s->isRef = false;
                                               //printf("0x6500");
                                             }*/
                                            
    for (Int_t ui = 0; ui < stt_event->totalNTDCHits; ui++) 
    {
  //cout<<"totalNTDCHits"<<stt_event->totalNTDCHits<<"\t"<<endl;
        if (((SttRawHit * ) stt_event->tdc_hits->ConstructedAt(ui))->channel == channel_nr +
        stt_channel_offsets[tdc_id]) 
        {
            if (((SttRawHit * ) stt_event->tdc_hits->ConstructedAt(ui))->leadTime == lastRise) 
            {
            //cout<<"Double hit :"<<endl;
                doubleHit = true;
                doubleCntr++;
            }
        }
        doubleCntr1++;
    }

                                            SttRawHit* a =
                                            stt_event->AddHit ( channel_nr );
                                            a->tdcid = tdc_id;
                                            a->leadTime = lastRise;
                                            a->trailTime = ( time );
                                            a->tot = - ( a->leadTime - a->trailTime );
                                            a->isRef = false;

                                            lastRise = 0;

                                             //printf("\tHit: %f on channel %d // on %x on %x\n", a->leadTime,channel_nr,tdc_id, sub_id);
                                        }
                                    }
                                }
                            } else if ( ( ( word >> 28 ) & 0xf ) == 0x6 ) {
                                // epoch counter
                                epoch = word & 0xffffff;
                                // printf("\tEpoch: %x on %x on %x\n", epoch,
                                // tdc_id, sub_id);
                            }

                            if ( tdc_words == tdc_size ) {
                                 //printf("\t Sizes: queue: %d tdc:%d current queue:%d tdc:%d\n", queue_size, tdc_size, queue_words, tdc_words);
                                break;
                            }
                        } // end of tdc data loop
                    }     // tdc select if

                } // end of loop over tdcs
                
                
                //dummy comment

                h_refTimeTRB1->Fill ( RefTime[0] - RefTime[1] );
                h_refTimeTRB1->Fill ( RefTime[0] - RefTime[2] );
                //h_refTimeTRB1->Fill ( RefTime[1] - RefTime[3] );
                h_refTimeTRB2->Fill ( RefTime[4] - RefTime[5] );
                h_refTimeTRB2->Fill ( RefTime[4] - RefTime[6] );
                h_refTimeTRB2->Fill ( RefTime[4] - RefTime[7] );

                //printf ( " 0:%lf  1:%lf  2: %lf  4: %lf  5: %lf 6: %lf\n", RefTime[0], RefTime[1], RefTime[2], RefTime[4], RefTime[5],RefTime[6] );

                h_refTimeTDC[0]->Fill ( RefTime[0] - RefTime[1] );
                h_refTimeTDC[1]->Fill ( RefTime[0] - RefTime[2] );
                //h_refTimeTDC[2]->Fill ( RefTime[1] - RefTime[4] );
                h_refTimeTDC[3]->Fill ( RefTime[4] - RefTime[5] );
                h_refTimeTDC[4]->Fill ( RefTime[4] - RefTime[6] );
                h_refTimeTDC[5]->Fill ( RefTime[4] - RefTime[7] );

//                                 printf("0-1:%lf  0-2:%lf  4-5:%lf  4-6:%lf  4-7:%lf  0-4:%lf\n",RefTime[0] - RefTime[1],RefTime[0] - RefTime[2],RefTime[4] - RefTime[5],RefTime[4] - RefTime[6],RefTime[4] - RefTime[7],RefTime[0] - RefTime[4]);
//                                h_refTimeTDC[6]->Fill ( RefTime[0] - RefTime[4] );


                if ( end_of_queue == true ) {

                    // printf("CLOSING EVENT\n");

                    tree->Fill();
                    stt_event->Clear();
                    emc_event->Clear();
                    sci_event->Clear();

                    N_events++;

                    if ( N_events % 1000 == 0 ) {
                        printf ( "%d\n", N_events );
                    }
                    break;
                    //}
                }

            } // end of sub loop

        } // end of queue size if
    }

    cout << "Repeated hits :" << doubleCntr << "/" << doubleCntr1 << endl;

    in_file.close();
    if ( use_tree_output ) {
        printf ( "writing file\n" );
        h_tdc_ref->Write();
        h_refTimeTRB1->Write();
        h_refTimeTRB2->Write();
        h_currpt->Write();
        for ( int h=0; h< 6; h++ ) {
            h_refTimeTDC[h]->Write();
        }
        tree->Write();
        ofile->Close();
    }
    std::cout << "Total number of processed events: " << N_events << std::endl;
    cout << "In File: "
         << "\t" << in_file_name << "\t"
         << "Out File: "
         << "\t" << out_file_name << endl;
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
            PDAQ_RawDecoder_HADES ( argv[1], argv[2], 100000000 );
        } else {
            PDAQ_RawDecoder_HADES ( argv[1], argv[2], atoi ( argv[3] ) );
        }

    else {
        return 1;
    }

    return 0;
}





