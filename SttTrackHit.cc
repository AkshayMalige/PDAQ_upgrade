

#include "SttTrackHit.h"

ClassImp(SttTrackHit)

    SttTrackHit::SttTrackHit() : TObject()
{

    // stt_tdc_event_sizes =0;
    // channel = 0;
    // leadTime = 0;
    // trailTime = 0;
    // tot = 0;

    vec_Track.clear();

    trackId = 0;
    trackSize = 0;
    Px0 = 0;
    Px1 = 0;
    Chix = 0;

    Py0 = 0;
    Py1 = 0;
    Chiy = 0;
    DriftT = 0;
    scint_time_diff = 0;

    // drifttime =0;
    // DriftRadius =0;
    // layer = 0;
    // module = 0;
    // fee = 0;
    // fee_channel = 0;
    // cell = 0;

    // //cell2 = 0;
    // x = 0;
    // y = 0;
    // z = 0;

    // isRef = false;
}

SttTrackHit::~SttTrackHit()
{
//   for (uint i = 0; i < vec_Track.size(); ++i)
//     delete vec_Track[i];
  vec_Track.clear();

}
