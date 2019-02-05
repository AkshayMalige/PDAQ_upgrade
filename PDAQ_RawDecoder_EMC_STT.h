#ifndef PDAQ_RAWDECODER_EMC_STT_H
#define PDAQ_RAWDECODER_EMC_STT_H

#include <fstream>
#include <sstream>
#include <iomanip>
#include <TH1F.h>
#include <TF1.h>
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <TObject.h>
#include <TClonesArray.h>
#include "TClassEdit.h"
#include <TBranch.h>
#include <TBranchElement.h>

#include "SttRawHit.h"
#include "SttEvent.h"
#include "EmcEvent.h"
#include "EmcHit.h"



#include "SttDetector.h"

#include "panda_subsystem.h"
#include "panda_subsystem_stt.h"
#include "panda_subsystem_sb.h"
#include "panda_subsystem_emc.h"

using namespace std;
// class SttEvent;
// class SttRawHit;
// // class SttHit;

// class EmcHit;
// class EmcEvent;


// class PandaSubsystem;
// class PandaSubsystemSTT;
// class PandaSubsystemSB;

//class Unpacker : public TObject

//{

//public:
    //Unpacker() {}

void pd_Zero_Event();

void pd_Reserve_Event(int Nsize);

void pd_init_hst();
void pd_Event_Processor();//{
void PDAQ_RawDecoder_EMC_STT(char *in_file_name,char *out_file_name);//{

    //~Unpacker() {}

//ClassDef(Unpacker, 1)

//};

#endif

