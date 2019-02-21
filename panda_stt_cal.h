#ifndef __PANDA_STT_CAL_H__
#define __PANDA_STT_CAL_H__

#include "panda_subsystem.h"
//#include "SttEvent.h"
#include "Stt_Cal_Event.h"
#include "string.h"

class PandaSttCal : public PandaSubsystem
{

  public:
    // SttEvent stt_track_can;

    Stt_Cal_Event stt_cal;

    PandaSttCal();
    ~PandaSttCal();

    ClassDef(PandaSttCal, 1)
};

#endif
