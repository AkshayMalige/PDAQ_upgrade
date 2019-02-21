
#ifndef PDAQ_STT_CALIBIRATOR_H
#define PDAQ_STT_CALIBIRATOR_H

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include <TF1.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLinearFitter.h>
#include <TLinearFitter.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "SttDetector.h"
#include "panda_stt_cal.h"
#include "panda_subsystem.h"
#include "panda_subsystem_sb.h"
#include "panda_subsystem_stt.h"

#include "SttEvent.h"
#include "SttHit.h"
#include "SttRawHit.h"
#include "Stt_Cal_Event.h"
#include "TRandom.h"

#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>
#include <cctype>
#include <fstream>
#include <sstream>
//#include <MFTGeomPar.h>
//#include <MParManager.cc>
#include "FTGeo.h"

using namespace std;

int PDAQ_Stt_Calibirator(char* intree, char* outtree, int maxEvents);

#endif