#ifndef __PANDA_SUBSYSTEM_STT_H__
#define __PANDA_SUBSYSTEM_STT_H__

#include "panda_subsystem.h"
#include "SttEvent.h"
#include "string.h"

class PandaSubsystemSTT : public PandaSubsystem {

    public:

    SttEvent stt_raw;

    
    PandaSubsystemSTT();
    ~PandaSubsystemSTT();

	ClassDef(PandaSubsystemSTT, 1)
};

#endif
