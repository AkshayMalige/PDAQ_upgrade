#include "PDAQ_Stt_Calibirator.h"

// int main() {
// 	return PDAQ_Stt_Calibirator();
// 	
// 
// }

int main ( int argc, char ** argv ) {

    if ( argc >= 3 )
        PDAQ_Stt_Calibirator ( argv[1], atoi(argv[2]) , argv[3] );
    else return 1;

    return 0;
}