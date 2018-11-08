#include <iostream>
#include "TApplication.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <MPar.h>
#include <MParContainer.h>
#include <MParManager.h>
#include <MParManager.cc>


// Header file for the classes stored in the TTree if any.
using  namespace std;
int main(int argc, char **argv) {
    TApplication theApp("App", &argc, argv);


MParManager* a = MParManager::instance();

a->setParamSource("parms.txt");

// a->addParameterContainer("dummy");
 a->parseSource();
a->print();  
    cout<<"Run Finished"<<endl;
    theApp.Run();
    return 0;
}
