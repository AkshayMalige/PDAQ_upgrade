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


bool check_tracks()
{


	TFile* outfile = new TFile ( "compare_track.root", "RECREATE" );

	TFile hard_file ( "tracking_out.root" );
	TFile soft_file ( "48Trig_f.root" );

	TTree* hard_tree = ( TTree* ) hard_file.Get ( "t1" );
	TTree* soft_tree = ( TTree* ) soft_file.Get ( "TRIG_tree" );
	
	hard_tree->Print();
	soft_tree->Print();

	std::vector<uint> vec_h_trigger;
	std::vector<UInt_t> *vec_trigger =0;
	UInt_t tr[1];
	Float_t sp[1];
	Float_t cs[1];

	hard_tree->SetBranchAddress ( "tr", &tr );
	hard_tree->SetBranchAddress ( "sp", &sp );
	hard_tree->SetBranchAddress ( "cs", &cs );

	soft_tree->SetBranchAddress("vec_trigger",&vec_trigger);

	int hard_iev = (int)hard_tree->GetEntries() ;
	int soft_iev = (int)soft_tree->GetEntries() ;
	
	cout << "number of hardware tracks:" << hard_iev << endl << endl;
	cout << "number of software tracks:" << soft_iev << endl << endl;


	TH1D * match	= new TH1D ("match","match",5,0,5);

    TH1F* h_hard_slope = new TH1F ( "h_hard_slope", "h_hard_slope", 100, -5, 5 );
    TH1F* h_hard_const = new TH1F ( "h_hard_const", "h_hard_const", 500, 0, 50 );

	bool fo=false;
	bool fe=false;
	
    //	ofstream myfile0;
    //	myfile0.open ("soft_dump.txt");
   	ofstream myfile1;
    	myfile1.open ("48hard_dump.txt");

	for(int i=0; i< soft_iev; i++){
		soft_tree->GetEntry(i);		
		for(int j=0; j< hard_iev; j++){
			hard_tree->GetEntry(j);
			if(vec_trigger->at(0) == tr[0]){
				//match->Fill(1);
				fo = true;
			}
					
		}
		fo == true ? match->Fill(1) : match->Fill(10);
		fo = false;
	}

	for(int i=0; i< hard_iev ; i++){
		hard_tree->GetEntry(i);
		myfile1 << hex << tr[0]<<"\t"<<sp[0]<<"\t"<<cs[0]<<endl;
		h_hard_const->Fill(cs[0]);
		h_hard_slope->Fill(sp[0]);
		//myfile1<<hex<<tr[0]<<"\t"<<"slope :"<<sp[0]<<"\t"<<cs[0]<<endl;
	}

	outfile->cd();
	match->Write();
	h_hard_const->Write();
	h_hard_slope->Write();
	outfile->Close();
	myfile1.close();

return true;
}
