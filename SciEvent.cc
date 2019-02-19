#include "SciEvent.h"

ClassImp(SciEvent)

SciEvent::SciEvent() {
	adc_hits = new TClonesArray("SciHit", 1000);
	//adc_parms = new TClonesArray("EmcHit", 1000);

	totalNTDCHits = 0;

	//totalNADCHits = 0;

}

SciHit* SciEvent::AddSciHit() {
	TClonesArray& thits = *adc_hits;
	SciHit* hit = new (thits[totalNTDCHits++]) SciHit();
	//hit->SetChannel(emc_Hits_ADC_channel);

	return hit;

}

// EmcHit* EmcEvent::AddParm() {

// 	TClonesArray& tparms = *adc_parms;
// 	EmcHit* parm = new (tparms[totalNADCHits++]) EmcHit();
// 	//parm->SetDiameter(emc_Cluster_diameter);

// 	return parm;


// }

void SciEvent::Clear(Option_t *) {
	adc_hits->Clear("C");
	//adc_parms->Clear("C");

	////delete adc_hits;
	//adc_hits = new TClonesArray("EmcHit", 1000);

	totalNTDCHits = 0;
	//totalNADCHits = 0;

}


