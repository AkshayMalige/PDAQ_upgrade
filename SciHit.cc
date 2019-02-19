#include "SciHit.h"

ClassImp(SciHit)

SciHit::SciHit() : TObject()

    {
	leadTime =0;
	trailTime =0;
	channel =0;
	tdcid =0;

	sci_tdc_event_sizes =0;

	isRef = false;
	}

