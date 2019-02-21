#ifndef H_SCI_EVENT
#define H_SCI_EVENT

#include "SciHit.h"
#include <TClonesArray.h>
#include <TObject.h>

class SciEvent : public TObject
{
  public:
    TClonesArray* adc_hits;
    // TClonesArray* adc_parms;

    int totalNTDCHits;
    // int totalNADCHits;

    SciEvent();
    virtual ~SciEvent() { Clear(); }

    void Clear(Option_t* opt = 0);

    SciHit* AddSciHit();

    // EmcHit* AddParm();

    ClassDef(SciEvent, 1)
};

#endif
