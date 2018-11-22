#ifndef H_EMC_HIT
#define H_EMC_HIT

#include <TObject.h>

class EmcHit : public TObject {
public:
	
	UInt_t emc_Cluster_local_time;
	UShort_t emc_Cluster_diameter;
	UShort_t emc_Cluster_2x;
	UShort_t emc_Cluster_2y;
	UShort_t emc_Cluster_N_hits;
	UInt_t emc_Hits_local_time;
	UChar_t emc_Hits_status;
	UShort_t emc_Hits_ADC_channel;
	UShort_t emc_Hits_Energy;

	UInt_t stt_tdc_event_sizes;

	// double emc_Cluster_local_time;
	// int emc_Cluster_diameter;
	// int emc_Cluster_2x;
	// int emc_Cluster_2y;
	// double emc_Cluster_N_hits;
	// double emc_Hits_local_time;
	// int emc_Hits_status;
	// int emc_Hits_ADC_channel;
	// double emc_Hits_Energy;

	// double stt_tdc_event_sizes;

	bool isRef;

	EmcHit();

	virtual ~EmcHit() { }

	//void SetChannel(int c) { emc_Hits_ADC_channel = c; } 

	//void SetDiameter(short d) { emc_Cluster_diameter = d; } 

	ClassDef(EmcHit, 1)

};
#endif


