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
	//a->print(); 

	//printf("\n\n Layers from ->getLayers() : %d\n",ftGeomPar->getStraws(1)); 
	//ftGeomPar->print();
	//printf("AAAA:\n");

	ftGeomPar->print();


	// for (int i =0; i<2 ; i++)
	// {
	// 	for (int j =0; j<4; j++)
	// 	{

	// 		for (int k =0; k<2; k++)
	// 		{
	// 			cout<<"offset x : "<<ftGeomPar->getOffsetX(i,j,k)<<endl;

	// 			cout<<"offset z : "<<ftGeomPar->getOffsetZ(i,j,k)<<endl;
	// 		}
	// 	}
	// }
	
    cout<<"Run Finished"<<endl;
    theApp.Run();
    return 0;
}
