#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TH1F.h>
#include "TFile.h"
#include <TGraph.h>
#include <TLinearFitter.h>
#include <TH2D.h>
#include "TTree.h"
#include <TNamed.h>
#include <TObject.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include "TMultiGraph.h"

using namespace std;


bool hard_data()
{


	TFile* outfile = new TFile ( "48_hard.root", "RECREATE" );

	TFile hard_file ( "48tracking_out.root" );
	TTree* hard_tree = ( TTree* ) hard_file.Get ( "t1" );
	hard_tree->Print();

	UInt_t tr[1];
	Float_t sp[1];
	Float_t cs[1];

	hard_tree->SetBranchAddress ( "tr", &tr );
	hard_tree->SetBranchAddress ( "sp", &sp );
	hard_tree->SetBranchAddress ( "cs", &cs );


	int hard_iev = (int)hard_tree->GetEntries() ;
	
	cout << "number of hardware tracks:" << hard_iev << endl << endl;

    	TH1F* h_hard_slope = new TH1F ( "h_hard_slope", "h_hard_slope", 100, -5, 5 );
    	TH1F* h_hard_const = new TH1F ( "h_hard_const", "h_hard_const", 500, 0, 50 );
    	TH1F* h_hard_size = new TH1F ( "h_hard_size", "h_hard_size", 50, 0, 50 );

   	ofstream myfile1;
    	myfile1.open ("48hard_dump.txt");

	for(int i=0; i< hard_iev ; i++){
		hard_tree->GetEntry(i);
		//myfile1 << hex << tr[0]<<endl;
		h_hard_const->Fill(cs[0]);
		h_hard_slope->Fill(sp[0]);
		myfile1<< hex << tr[0]<<"\t"<<sp[0]<<"\t"<<cs[0]<<endl;

	}

	outfile->cd();
	h_hard_const->Write();
	h_hard_slope->Write();
	outfile->Close();
	myfile1.close();

return true;
}
