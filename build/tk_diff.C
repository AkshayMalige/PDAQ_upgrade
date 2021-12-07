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
		double s_diff_r;
		double c_diff_r;
};
//Int_t Fee_scan_analysis(const char *  	inputFile, const char* outputFile)


  void normalize(TH1* hist)                                                                                                                                                                               
   {                                                                                                                                                                                                        
      for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)                                                                                                                                                          
      {                                                                                                                                                                                                    
         hist->SetBinContent( j, hist->GetBinContent(j) / hist->GetBinWidth(j) );                                                                                                                          
//         hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );                                                                                                                                  
         hist->SetBinError( j, hist->GetBinError(j) / hist->GetBinWidth(j) );                                                                                                                              
      }                                                                                                                                                                                                    
   } 
   
Int_t tk_diff()
{

	TFile* scan_results = new TFile("beam_marking_track_45915_tk_diff.root", "RECREATE");
	//TFile* scan_results = new TFile("lab_marking_track_44584_tk_diff.root", "RECREATE");

    	ifstream iFile("./ieee/beam_marking_track_45915_pandas.txt");  
    	//ifstream iFile("./lab_marking_track_44584_pandas.txt");  
    	track tk_obj;

    	vector<track> vec_data; 
    	
	    TH1F *h_s_diff = new TH1F("s_diff","s_diff",300,-3,3);
	    TH1F *h_c_diff = new TH1F("c_diff","c_diff",300,-3,3);
	    TH2F *sc_diff = new TH2F("sc_diff","sc_diff",12,-3,3,12,-3,3);
	    
	    TH1F *h_s_slope = new TH1F("h_s_slope","h_s_slope",200,-50,50);
	    TH1F *h_h_slope = new TH1F("h_h_slope","h_h_slope",200,-50,50);
    	    TH1F *h_s_const = new TH1F("h_s_const","h_s_const",400,-100,100);
    	    TH1F *h_h_const = new TH1F("h_h_const","h_h_const",400,-100,100);
 
	  
	    while (!iFile.eof())
	    {
	    	string trigger;
		double s_slope;
		double s_const;
		double h_slope;
		double h_const;
		double s_diff;
		double c_diff;
		double s_diff_r;
		double c_diff_r;
		
		iFile >> trigger >>s_slope >> s_const >> h_slope >> h_const >> s_diff >> c_diff >> s_diff_r >> c_diff_r;
		

		/*tk_obj.trigger = trigger;
		tk_obj.s_slope = s_slope;
		tk_obj.s_const = s_const*0.505;
		tk_obj.h_slope = h_slope;
		tk_obj.h_const = h_const*0.505;
		tk_obj.s_diff = s_diff;
		tk_obj.c_diff = c_diff;
		tk_obj.s_diff_r = s_diff_r;
		tk_obj.c_diff_r = c_diff_r;*/

		//vec_data.push_back(tk_obj);
		
		
	    	h_s_diff->Fill(s_diff_r);
	    	h_c_diff->Fill(c_diff_r);
	    	sc_diff->Fill(s_diff , c_diff);	    
	    	h_s_slope->Fill(s_slope);
	    	h_h_slope->Fill(h_slope);
	    	h_s_const->Fill(s_const*0.505);
	    	h_h_const->Fill(h_const*0.505);

	    }    	    
	    	    
	    	    
	   /* for(int i=0; i< vec_data.size(); i++){
	    	s_diff->Fill(vec_data[i].s_diff_r);
	    	c_diff->Fill(vec_data[i].c_diff_r);
	    	sc_diff->Fill(vec_data[i].s_diff , vec_data[i].c_diff);	    
	    	h_s_slope->Fill(vec_data[i].s_slope);
	    	h_h_slope->Fill(vec_data[i].h_slope);
	    	h_s_const->Fill(vec_data[i].s_const);
	    	h_h_const->Fill(vec_data[i].h_const);
	    	
	    }*/
	    
	    for(int s=0; s<200; s++){
	    	cout<<s<<": "<<h_s_slope->GetBinContent(s)<<"\t";
	    }
	    	h_s_diff->Scale(1/h_s_diff->GetEntries());
	    	h_c_diff->Scale(1/h_c_diff->GetEntries());
	    
	    	sc_diff->Scale(100/sc_diff->GetEntries());
	    	sc_diff->GetXaxis()->CenterLabels();
	    	sc_diff->GetYaxis()->CenterLabels();
	    	sc_diff->Write();
	    
	    	TCanvas * c1 = new TCanvas ("c1","c1");
	    	c1->Divide(2,1);
	    	c1->cd(1);
	    	h_s_diff->Draw();
	    	c1->cd(2);
	    	h_c_diff->Draw();
	    
	    	c1->Write();

		h_s_diff->Write();
		h_c_diff->Write();
		
		
		//normalize(h_s_const);
		//normalize(h_h_const);
		h_s_slope->Write();
		h_h_slope->Write();
		h_s_const->Write();
		h_h_const->Write();
		
		
		/////////////////////////////
		h_s_const->GetXaxis()->SetRangeUser(-2,18);
		h_h_const->GetXaxis()->SetRangeUser(-2,18);
		h_s_const->Scale(1/h_s_const->GetEntries());
		h_h_const->Scale(1/h_h_const->GetEntries());
		h_s_const->SetLineColor(kBlue);
		h_h_const->SetLineColor(kRed);
		h_s_const->SetLineWidth(3);
		h_h_const->SetLineWidth(3);
		
		
		
		h_s_slope->GetXaxis()->SetRangeUser(-4,4);
		h_h_slope->GetXaxis()->SetRangeUser(-4,4);
		h_s_slope->Scale(1/h_s_slope->GetEntries());
		h_h_slope->Scale(1/h_h_slope->GetEntries());
		h_s_slope->SetLineColor(kBlue);
		h_h_slope->SetLineColor(kRed);
		h_s_slope->SetLineWidth(3);
		h_h_slope->SetLineWidth(3);
		
		TCanvas *c2 = new TCanvas("c2","c2",600,700);
   		TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   		pad1->SetBottomMargin(0);
   		pad1->Draw();
   		pad1->cd();
   		h_s_slope->DrawCopy();
   		h_h_slope->Draw("same");
   		c2->cd();
   		TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
   		pad2->SetTopMargin(0);
   		pad2->Draw();
   		pad2->cd();
   		/*h_s_slope->Sumw2();
   		h_s_slope->SetStats(0);
   		h_s_slope->Divide(h_h_slope);
   		h_s_slope->SetMarkerStyle(21);
   		h_s_slope->Draw("ep");*/
   		
   		h_s_diff->SetMarkerStyle(21);
   		h_s_diff->Draw("ep");
   		c2->cd(); 

		TCanvas *c3 = new TCanvas("c3","c3",600,700);
   		TPad *pad3 = new TPad("pad3","pad3",0,0.3,1,1);
   		pad3->SetBottomMargin(0);
   		pad3->Draw();
   		pad3->cd();
   		h_s_const->DrawCopy();
   		h_h_const->Draw("same");
   		c3->cd();
   		TPad *pad4 = new TPad("pad4","pad4",0,0,1,0.3);
   		pad4->SetTopMargin(0);
   		pad4->Draw();
   		pad4->cd();
   		/*h_s_const->Sumw2();
   		h_s_const->SetStats(0);
   		h_s_const->Divide(h_h_const);
   		h_s_const->SetMarkerStyle(21);
   		h_s_const->Draw("ep");*/
   		
   		h_c_diff->SetMarkerStyle(21);
   		h_c_diff->Draw("ep");
   		c3->cd();  
   		
   		c2->Write();
   		c3->Write();  		

   		


	return 0;
}



