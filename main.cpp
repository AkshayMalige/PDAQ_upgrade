#include <iostream>
#include "TApplication.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>
#include <MParManager.cc>
#include "MFTGeomPar.h"

// Header file for the classes stored in the TTree if any.
using  namespace std;
int main(int argc, char **argv) {
    TApplication theApp("App", &argc, argv);


	MParManager* a = MParManager::instance();
	MFTGeomPar* b;
	MParContainer* c;
	MPar * d;
	a->setParamSource("params.txt");
	//pm()->addParameterContainer("MFibersStackGeomPar",  new MFTGeomPar());
	pm()->getParameterContainer("MFibersStackGeomPar");

	//b->getParams();
	//c->MParContainer();
	a->parseSource();

	//a->getParameterContainer("MFibersStackGeomPar");
	a->print();  
    cout<<"Run Finished"<<endl;
    theApp.Run();
    return 0;
}
