#ifndef __PANDA_STT_TRACK_H__
#define __PANDA_STT_TRACK_H__

#include "SttTrackEvent.h"
#include "panda_subsystem.h"
#include "string.h"

class PandaSttTrack : public PandaSubsystem
{

  public:
    Stt_Track_Event stt_track_can;
    // SttEvent stt_drift_radius;

    PandaSttTrack();
    ~PandaSttTrack();

    ClassDef(PandaSttTrack, 1)
};

#endif
