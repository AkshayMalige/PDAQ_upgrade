#include "EmcEvent.h"

ClassImp(EmcEvent)

    EmcEvent::EmcEvent()
{
    adc_hits = new TClonesArray("EmcHit", 1000);
    adc_parms = new TClonesArray("EmcHit", 1000);

    totalNTDCHits = 0;

    totalNADCHits = 0;
}

EmcHit* EmcEvent::AddHit()
{
    TClonesArray& thits = *adc_hits;
    EmcHit* hit = new (thits[totalNTDCHits++]) EmcHit();
    // hit->SetChannel(emc_Hits_ADC_channel);

    return hit;
}

EmcHit* EmcEvent::AddParm()
{

    TClonesArray& tparms = *adc_parms;
    EmcHit* parm = new (tparms[totalNADCHits++]) EmcHit();
    // parm->SetDiameter(emc_Cluster_diameter);

    return parm;
}

void EmcEvent::Clear(Option_t*)
{
    adc_hits->Clear("C");
    adc_parms->Clear("C");

    ////delete adc_hits;
    // adc_hits = new TClonesArray("EmcHit", 1000);

    totalNTDCHits = 0;
    totalNADCHits = 0;
}
