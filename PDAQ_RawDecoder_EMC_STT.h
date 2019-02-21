#ifndef PDAQ_RAWDECODER_EMC_STT_H
#define PDAQ_RAWDECODER_EMC_STT_H

#include "TClassEdit.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include <TBranch.h>
#include <TBranchElement.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TH1F.h>
#include <TObject.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "EmcEvent.h"
#include "EmcHit.h"
#include "SttEvent.h"
#include "SttRawHit.h"

#include "SttDetector.h"

#include "panda_subsystem.h"
#include "panda_subsystem_emc.h"
#include "panda_subsystem_sb.h"
#include "panda_subsystem_stt.h"

using namespace std;
// class SttEvent;
// class SttRawHit;
// // class SttHit;

// class EmcHit;
// class EmcEvent;

// class PandaSubsystem;
// class PandaSubsystemSTT;
// class PandaSubsystemSB;

// class Unpacker : public TObject

//{

// public:
// Unpacker() {}

void pd_Zero_Event();

void pd_Reserve_Event(int Nsize);

void pd_init_hst();
void pd_Event_Processor();                                             //{
void PDAQ_RawDecoder_EMC_STT(char* in_file_name, char* out_file_name); //{

//~Unpacker() {}

// ClassDef(Unpacker, 1)

//};

#endif
