#include "Stt_Cal_Event.h"
#include "SttHit.h"

ClassImp(Stt_Cal_Event)

    Stt_Cal_Event::Stt_Cal_Event()
{
    tdc_cal_hits = new TClonesArray("SttHit", 1000);
    // tdc_events = new TClonesArray("SttHit", 1000);
    // tdc_raw = new TClonesArray("SttRawHit", 1000);

    total_cal_NTDCHits = 0;
    // totalNTDCEvents = 0;
    // rawNTDCEvents =0;
}

SttHit* Stt_Cal_Event::AddCalHit(int channel)
{
    TClonesArray& thits = *tdc_cal_hits;
    SttHit* cal_hit = new (thits[total_cal_NTDCHits++]) SttHit();
    cal_hit->SetChannel(channel);

    return cal_hit;
}

// SttHit* SttEvent::event_size(int stt_tdc_event_sizes) {
// 	TClonesArray& tevents = *tdc_events;
// 	SttHit* event = new (tevents[totalNTDCEvents++]) SttHit();
// 	//hit->SetChannel(channel);

// 	return event;
// }

// SttRawHit* SttEvent::AddRawHit(int channel) {
// 	TClonesArray& trawevents = *tdc_raw;
// 	SttRawHit* raw_hit = new (trawevents[rawNTDCEvents++]) SttRawHit();
// 	raw_hit->SetChannel(channel);

// 	return trawevents;
// }

void Stt_Cal_Event::CalClear(Option_t*)
{
    tdc_cal_hits->Clear("C");

    // tdc_events->Clear("C");
    // tdc_raw->Clear("C");

    //	delete tdc_hits;
    //	tdc_hits = new TClonesArray("SttHit", 1000);

    total_cal_NTDCHits = 0;
    // totalNTDCEvents = 0;
    // rawNTDCEvents = 0;
}
