#ifndef H_STT_EVENT
#define H_STT_EVENT

#include <TNamed.h>
#include <TClonesArray.h>
#include "SttHit.h"
#include "SttRawHit.h"


class SttEvent : public TNamed {
public:

	TClonesArray* tdc_hits; 
	int totalNTDCHits;
	
	// TClonesArray* tdc_events; 
	// int totalNTDCEvents;

	// TClonesArray* tdc_raw; 
	// int rawNTDCEvents;
	
	SttEvent();
	virtual ~SttEvent() { Clear(); }

	void Clear(void);

	SttRawHit* AddHit(int channel);
	//SttHit* event_size(int stt_tdc_event_sizes);

	//SttRawHit* AddRawHit(int channel);

	// SttRawHit* AddHit(int channel);
	// SttRawHit* event_size(int stt_tdc_event_sizes);

	ClassDef(SttEvent, 1)
};

#endif