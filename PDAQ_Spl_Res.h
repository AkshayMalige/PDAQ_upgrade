#ifndef PDAQ_SPL_RES_H
#define PDAQ_SPL_RES_H

#include "TFile.h"
#include <TGraph.h>
#include <TH1F.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "SttEvent.h"
#include "SttHit.h"
#include "TF1.h"
#include "TRandom.h"
#include "TTree.h"
#include <TH2D.h>
#include <TLinearFitter.h>
#include <TNamed.h>
#include <TObject.h>
#include <cstdlib>
#include <math.h>

#include "TCanvas.h"
#include <TChain.h>
#include <TROOT.h>
#include "TObject.h"
#include "TEllipse.h"
#include "TLegend.h"

using namespace std;

//Bool_t PDAQ_Spl_Res(void);

Bool_t PDAQ_Spl_Res(char* intree, char* outtree, int maxEvents);


#endif
