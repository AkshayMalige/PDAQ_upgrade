#ifndef H_STT_Track_EVENT
#define H_STT_Track_EVENT

#include "SttHit.h"
#include "SttRawHit.h"
#include "SttTrackHit.h"
#include <TClonesArray.h>
#include <TObject.h>

class Stt_Track_Event : public TObject
{
  public:
//     TClonesArray* tdc_track_hits;
    std::vector<SttTrackHit> tdc_track_hits;
    int total_track_NTDCHits;

    // TClonesArray* tdc_events;
    // int totalNTDCEvents;

    // TClonesArray* tdc_raw;
    // int rawNTDCEvents;

    Stt_Track_Event();
    virtual ~Stt_Track_Event();

    void TrackClear(void);

    SttTrackHit& AddTrackHit();
    // SttHit* event_size(int stt_tdc_event_sizes);

    // SttRawHit* AddRawHit(int channel);

    // SttRawHit* AddHit(int channel);
    // SttRawHit* event_size(int stt_tdc_event_sizes);

    ClassDef(Stt_Track_Event, 1)
};

#endif
