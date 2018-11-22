#ifndef H_EMC_EVENT
#define H_EMC_EVENT

#include <TObject.h>
#include <TClonesArray.h>
#include "EmcHit.h"

class EmcEvent : public TObject {
public:

	TClonesArray* adc_hits; 
	TClonesArray* adc_parms; 

	int totalNTDCHits;
	int totalNADCHits;
	
	EmcEvent();
	virtual ~EmcEvent() { Clear(); } 

	void Clear(Option_t * opt = 0);

	EmcHit* AddHit();

	EmcHit* AddParm();

	ClassDef(EmcEvent, 1)
};

#endif
