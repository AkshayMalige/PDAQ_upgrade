

#include "SttRawHit.h"

ClassImp(SttRawHit)

    SttRawHit::SttRawHit()
{

    stt_tdc_event_sizes = 0;
    tdcid = 0;
    trigger_no = 0;
    channel = 0;
    leadTime = 0;
    trailTime = 0;
    tot = 0;
    marking = false;

    isRef = false;
}
