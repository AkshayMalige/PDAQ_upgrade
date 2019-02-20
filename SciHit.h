#ifndef H_SCI_HIT
#define H_SCI_HIT

#include <TObject.h>

class SciHit : public TObject {
public:
	
	double leadTime;
	double trailTime;
	int channel;
	UInt_t tdcid;

	UInt_t sci_tdc_event_sizes;


	bool isRef;

	SciHit();

	virtual ~SciHit() { }

	ClassDef(SciHit, 1)

};
#endif


