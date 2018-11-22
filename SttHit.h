#ifndef H_STT_HIT
#define H_STT_HIT

#include <TObject.h>
#include <SttRawHit.h>


class SttHit : public SttRawHit {
public:
	// int channel;
	// int stt_tdc_event_sizes;
	// double leadTime;
	// double trailTime;
	// double tot;
	double drifttime;
	double DriftRadius;

	short layer;
	short module;
	short fee;
	short fee_channel;
	short cell;

	double x;
	double y;
	double z;

	//short cell2;

	//bool isRef;

	SttHit();

	~SttHit() { }

	void SetChannel(int c) { channel = c; } 


	ClassDef(SttHit, 1)

};
#endif


