#ifndef H_STT_EVENT
#define H_STT_EVENT

#include "SttHit.h"
#include "SttRawHit.h"
#include <TClonesArray.h>
#include <TNamed.h>

class SttEvent : public TObject
{
  public:
    TClonesArray* tdc_hits;
    int totalNTDCHits;

    // TClonesArray* tdc_events;
    // int totalNTDCEvents;

    // TClonesArray* tdc_raw;
    // int rawNTDCEvents;

    SttEvent();
    virtual ~SttEvent() { Clear(); delete tdc_hits; }

    void Clear(Option_t * t = "");

    SttRawHit* AddHit(int channel);
    // SttHit* event_size(int stt_tdc_event_sizes);

    // SttRawHit* AddRawHit(int channel);

    // SttRawHit* AddHit(int channel);
    // SttRawHit* event_size(int stt_tdc_event_sizes);

    ClassDef(SttEvent, 1)
};

#endif