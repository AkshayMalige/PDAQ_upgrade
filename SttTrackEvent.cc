#include "SttTrackEvent.h"
#include <TBuffer.h>

ClassImp(Stt_Track_Event)

    Stt_Track_Event::Stt_Track_Event() : TObject()
{
//     std::vector<SttTrackHit *> tdc_track_hits;// = new TClonesArray("SttTrackHit", 1000);
    // tdc_events = new TClonesArray("SttHit", 1000);
    // tdc_raw = new TClonesArray("SttRawHit", 1000);

    total_track_NTDCHits = 0;
    // totalNTDCEvents = 0;
    // rawNTDCEvents =0;
}

SttTrackHit& Stt_Track_Event::AddTrackHit()
{
//     TClonesArray& thits = *tdc_track_hits;
//     SttTrackHit* track_hit = new (thits[total_track_NTDCHits++]) SttTrackHit();
    // track_hit->SetChannel(trackId);
      SttTrackHit track_hit;
    tdc_track_hits.push_back(track_hit);
    total_track_NTDCHits = tdc_track_hits.size();

    return tdc_track_hits[total_track_NTDCHits-1];
}

Stt_Track_Event::~Stt_Track_Event()
{
//   TrackClear();
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

void Stt_Track_Event::TrackClear(void)
{
//   for (uint i = 0; i < tdc_track_hits.size(); ++i)
//     delete tdc_track_hits[i];
  tdc_track_hits.clear();

//     tdc_track_hits->Clear("C");

    // tdc_events->Clear("C");
    // tdc_raw->Clear("C");

    //	delete tdc_hits;
    //	tdc_hits = new TClonesArray("SttHit", 1000);

    total_track_NTDCHits = 0;
    // totalNTDCEvents = 0;
    // rawNTDCEvents = 0;
}



