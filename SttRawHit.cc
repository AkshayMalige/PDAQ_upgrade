

#include "SttRawHit.h"

ClassImp(SttRawHit)

SttRawHit::SttRawHit() {

			stt_tdc_event_sizes =0;
			tdcid=0;
			channel = 0;
			//new_channel=0;
			leadTime = 0;
			trailTime = 0;
			tot = 0;
			//emc_time =0;

			isRef = false;
	}
