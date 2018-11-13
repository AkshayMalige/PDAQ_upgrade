#include <iostream>
#include "TApplication.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>
#include <MParManager.cc>
//#include "MFTGeomPar.h"
//#include "pFTGeomPar.h"
#include "FTGeo.h"

// Header file for the classes stored in the TTree if any.
using  namespace std;
int main(int argc, char **argv) {
    TApplication theApp("App", &argc, argv);


	MParManager* a = MParManager::instance();
	MFTGeomPar* ftGeomPar = new MFTGeomPar();
	MPar * d;

	a->setParamSource("ftparams.txt");
	a->parseSource();
	pm()->addParameterContainer("MFTPar", ftGeomPar);
	d = pm()->getParameterContainer("MFTPar");

	if (!ftGeomPar)
    {
        std::cerr << "Parameter container 'PFTGeomPar' was not obtained!" << std::endl;
        //exit(EXIT_FAILURE);
    }
    else{ 
	printf("D: %d\n", ftGeomPar);
	}
	//a->getParameterContainer("MFibersStackGeomPar");
	a->print(); 

	printf("\n\nGet Modules : %d\n",ftGeomPar->getModules()); 
	//ftGeomPar->print();
	//printf("AAAA:\n");

	//b->print();
    cout<<"Run Finished"<<endl;
    theApp.Run();
    return 0;
}
