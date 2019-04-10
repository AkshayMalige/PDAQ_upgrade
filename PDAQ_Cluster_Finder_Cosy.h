#ifndef PDAQ_CLUSTER_FINDER_COSY_H
#define PDAQ_CLUSTER_FINDER_COSY_H
#include "PDAQ_Trackfilter.h"

int PDAQ_Cluster_Finder_Cosy(char* intree, char* outtree, int maxEvents);

std::vector<SciHit*> ex_Scint_pileup(PandaSubsystemSCI* SCI_CAL);

std::vector<SttHit> ex_STT_double(PandaSttCal* STT_CAL);

#endif
