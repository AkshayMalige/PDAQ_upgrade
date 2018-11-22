#ifndef __PANDA_SUBSYSTEM_EMC_H__
#define __PANDA_SUBSYSTEM_EMC_H__

#include "panda_subsystem.h"
#include "EmcEvent.h"
#include "string.h"

class PandaSubsystemEMC : public PandaSubsystem {

    public:

    EmcEvent emc_raw;
    
    PandaSubsystemEMC();
    ~PandaSubsystemEMC();

	ClassDef(PandaSubsystemEMC, 1)
};

#endif
