#ifndef PDAQ_TRACK_FILTER_H
#define PDAQ_TRACK_FILTER_H

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include <TF1.h>
#include <TGraph.h>
#include <TLinearFitter.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <cctype>
#include <sstream>

#include "SciHit.h"
#include "SttEvent.h"
#include "SttHit.h"
#include "SttRawHit.h"
#include "SttTrackEvent.h"

#include "panda_stt_cal.h"
#include "panda_stt_track.h"
#include "panda_subsystem.h"
#include "panda_subsystem_sb.h"
#include "panda_subsystem_sci.h"
#include "panda_subsystem_stt.h"

#include "FTGeo.h"
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>


using namespace std;

struct histograms{
    TH1F* h_filtered_cluster_size = new TH1F(
        "h_filtered_cluster_size", "h_filtered_cluster_size", 110, -10, 100);
};

bool f_sttHitCompareLeadTime(SttHit* a, SttHit* b);

bool f_sttHitCompareCell(SttHit* a, SttHit* b);

std::vector<SttHit*> GetPairs(std::vector<SttHit*> vec_get_pairs);

std::vector<std::vector<SttHit*>> clusterfinder(std::vector<SttHit*> vec_flayer);


void add_tuples(const std::vector<std::vector<std::vector<SttHit*>>>& vectors,
                std::size_t pos, std::vector<std::vector<SttHit*>> prefix,
                std::vector<std::vector<std::vector<SttHit*>>>& result);

std::vector<std::vector<std::vector<SttHit*>>>
make_tuples(const std::vector<std::vector<std::vector<SttHit*>>>& vectors);


bool PDAQ_Event_Finder(std::vector<SttHit*> vec_stthits, int i,
                       TTree* PDAQ_tree, Stt_Track_Event* stt_event,
                       MFTGeomPar* ftGeomPar, PandaSubsystemSCI* SCI_CAL,  histograms* h);

#endif