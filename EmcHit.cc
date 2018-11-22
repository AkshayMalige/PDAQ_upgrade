

#include "EmcHit.h"

ClassImp(EmcHit)

EmcHit::EmcHit() : TObject()

{
	 emc_Cluster_local_time =0;
	 emc_Cluster_diameter=0;
	 emc_Cluster_2x=0;
	 emc_Cluster_2y=0;
	 emc_Cluster_N_hits=0;
	 emc_Hits_local_time=0;
	 emc_Hits_status=0;
	 emc_Hits_ADC_channel=0;
	 emc_Hits_Energy=0;

	 stt_tdc_event_sizes=0;

		isRef = false;
	}

