#ifndef H_STT_TRACK_HIT
#define H_STT_TRACK_HIT

#include <TObject.h>
#include <SttHit.h>

class SttTrackHit : public TObject {
public:
	// int channel;
	// int stt_tdc_event_sizes;
	// double leadTime;
	// double trailTime;
	// double tot;
	// double drifttime;
	// double DriftRadius;

	// short layer;
	// short module;
	// short fee;
	// short fee_channel;
	// short cell;

	// double x;
	// double y;
	// double z;

     std::vector<SttHit*> vec_Track;

	 double trackId;
     double trackSize;
     double Px0;
     double Px1;
     double Chix;

     double Py0;
     double Py1;
     double Chiy;
	//short cell2;

	//bool isRef;

	SttTrackHit();

	~SttTrackHit() { }

	//void SetChannel(int c) { trackId = c; } 


	ClassDef(SttTrackHit, 1)

};
#endif


