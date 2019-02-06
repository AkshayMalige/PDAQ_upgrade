#include "PDAQ_Stt_Calibirator.h"

// int main() {
// 	return PDAQ_Stt_Calibirator();
// 	
// 
// }

int main ( int argc, char ** argv ) {

    if ( argc >= 3 )
      
      	if (!argv[3])
	{
	  printf("\n\nNote : One Million Events will be processed. To change add the number of events to be processed after the ouput file name.\n");
	  //atoi(argv[3]) == 1000;
	  sleep(2);
	  PDAQ_Stt_Calibirator ( argv[1],argv[2],1000000 );
	}
	else 
        PDAQ_Stt_Calibirator ( argv[1] , argv[2],atoi(argv[3]) );
    else return 1;

    return 0;
}