#include <iostream>
#include <fstream>
#include <vector>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include "TFile.h"
#include <string>
#include <TMultiGraph.h>
#include <TLegend.h>

using namespace std;


class track {
   public:
	    	string trigger;
		double s_slope;
		double s_const;
		double h_slope;
		double h_const;
		double s_diff;
		double c_diff;
};
//Int_t Fee_scan_analysis(const char *  	inputFile, const char* outputFile)
Int_t tk_diff()
{



  


	TFile* scan_results = new TFile("tk_diff.root", "RECREATE");


    	ifstream iFile("pandas.txt");  

    	track tk_obj;

    	vector<track> vec_data; 
 

	    while (!iFile.eof())
	    {
	    	string trigger;
		double s_slope;
		double s_const;
		double h_slope;
		double h_const;
		double s_diff;
		double c_diff;

		iFile >> trigger >>s_slope >> s_const >> h_slope >> h_const >> s_diff >> c_diff;

		tk_obj.trigger = trigger;
		tk_obj.s_slope = s_slope;
		tk_obj.s_const = s_const;
		tk_obj.h_slope = h_slope;
		tk_obj.h_const = h_const;
		tk_obj.s_diff = s_diff;
		tk_obj.c_diff = c_diff;

		vec_data.push_back(tk_obj);

	    }
	    
	    
	    TH1F * s_diff = new TH1F("s_diff","s_diff",300,-3,3);
	    TH1F * c_diff = new TH1F("c_diff","c_diff",300,-3,3);
	    TH2F * sc_diff = new TH2F("sc_diff","sc_diff",12,-3,3,12,-3,3);
	    
	    
	    for(int i=0; i< vec_data.size(); i++){
	    	s_diff->Fill(vec_data[i].s_diff);
	    	c_diff->Fill(vec_data[i].c_diff - 1);
	    	sc_diff->Fill(vec_data[i].s_diff , vec_data[i].c_diff - 1);
	    //	printf("%f\n",vec_data.at(i).s_diff);
	    }
	    
	    s_diff->Scale(1/s_diff->GetEntries());
	    c_diff->Scale(1/c_diff->GetEntries());
	    
	    sc_diff->Scale(100/sc_diff->GetEntries());
	    sc_diff->GetXaxis()->CenterLabels();
	    sc_diff->GetYaxis()->CenterLabels();
	    sc_diff->Write();
	    
	    TCanvas * c1 = new TCanvas ("c1","c1");
	    c1->Divide(2,1);
	    c1->cd(1);
	    s_diff->Draw();
	    c1->cd(2);
	    c_diff->Draw();
	    
	    c1->Write();

		s_diff->Write();
		c_diff->Write();
		


	return 0;
}
