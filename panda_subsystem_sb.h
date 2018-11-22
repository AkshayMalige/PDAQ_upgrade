#ifndef __PANDA_SUBSYSTEM_SB_H__
#define __PANDA_SUBSYSTEM_SB_H__

#include "panda_subsystem.h"
#include "string.h"

class PandaSubsystemSB : public PandaSubsystem {

    public:

    int SB_number;
    int Header_Status;
    
    PandaSubsystemSB();
    ~PandaSubsystemSB();

	ClassDef(PandaSubsystemSB, 1)
};

#endif