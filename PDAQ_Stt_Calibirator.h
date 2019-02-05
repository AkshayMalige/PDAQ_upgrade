
#ifndef PDAQ_STT_CALIBIRATOR_H
#define PDAQ_STT_CALIBIRATOR_H

#include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <TLinearFitter.h>
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <TGraph.h>
#include <math.h>
#include <cstdlib>
#include <TLinearFitter.h>

#include "panda_subsystem.h"
#include "panda_subsystem_stt.h"
#include "panda_subsystem_sb.h"
#include "panda_stt_cal.h"
#include "SttDetector.h"


#include "TRandom.h"
#include "SttRawHit.h"
#include "SttHit.h"
#include "SttEvent.h"
#include "Stt_Cal_Event.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>
//#include <MFTGeomPar.h>
//#include <MParManager.cc>
#include "FTGeo.h"
                           
using namespace std;


int PDAQ_Stt_Calibirator(char* intree, int maxEvents=1000000, char* outtree="PDAQ_calibrator_tree.root");

#endif