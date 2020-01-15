#include <TH2.h>
#include <TH1.h>
#include <iostream>
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include <TF1.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <TGraph.h>
#include "TObject.h"
#include <vector>
#include "TObject.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMultiGraph.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

using namespace std;


int count_line(TString talysName) {
 ifstream f1;
   char c;
   int numchars, numlines;

   f1.open(talysName);

   numchars = 0;
   numlines = 0;
   f1.get(c);
   while (f1) {
     while (f1 && c != '\n') {
       numchars = numchars + 1;
       f1.get(c);
     }
     numlines = numlines + 1;
     f1.get(c);
   }
   cout << "The file has " << numlines << " lines and " 
     << numchars << " characters" << endl;
f1.close();
   return numlines;
}


Bool_t cross_talk()
{

int HV = 5;
int TP = 3;
int TH = 2;    
    //const char* filename[x];
    TString talysName = "crstlk.txt";  //a has all the file names you must read
    fstream talysInput(talysName,ios::in);
    if(!talysInput.is_open())
	{
	cout << "Could not open"<< talysName << endl;
	return kFALSE;
	}



    const Int_t max = count_line(talysName);  // max is the total number of files in the list a.txt

    TString file_name[max];

    for(Int_t i=0; i<max; i++)
    {
        talysInput >> file_name[i];
    }

    talysInput.close();
  
	double index[HV] ;
	double arr_ct_all[TP*TH][HV];
	double arr_ct_more800[TP*TH][HV];
	double arr_ct_less800[TP*TH][HV];

int setting =0;
int hv_cntr =0;

TCanvas * C1 =new TCanvas ( "C1","C1" );
TCanvas * C2 =new TCanvas ( "C2","C2" );
TCanvas * C3 =new TCanvas ( "C3","C3" );

    for (int k=0; k<max; k++)
    {

	   
    	//cout<<file_name[k]<<endl;

    	TString f1 = file_name[k];
	
        TFile* file = new TFile ( f1, "READ" );

        
        if (file->IsOpen()) 
        {
           // cout << Form("File nr. %d opened successfully",k)<< endl;

        } 
        else 
        {
            cout << Form("Error opening file nr. %d",k) << endl;
            exit(-1);
        }
        file->cd();
        TH1F* Cross_TlkCase1;
	TH1F* Cross_TlkCase2;
	TH1F* h_CrossCase;
        Cross_TlkCase1 = ( TH1F* )file->Get ( "Cross_TlkCase1" ); //change here the name of your histogram
	Cross_TlkCase2 = ( TH1F* )file->Get ( "Cross_TlkCase2" );
	h_CrossCase = ( TH1F* )file->Get ( "h_CrossCase" );

        //cout<<h_CrossCase->GetBinContent(3)<<endl;
	
	
	//index[max] = k;
	arr_ct_all[setting][hv_cntr] = ( h_CrossCase->GetBinContent(3) * 100 ) / 8;
	arr_ct_more800[setting][hv_cntr] = ( Cross_TlkCase2->GetBinContent(3) * 100 ) / 8;
	arr_ct_less800[setting][hv_cntr] = ( Cross_TlkCase1->GetBinContent(3) * 100 ) / 8;

	//cout <<index[max] <<"\t"<<arr_ct_all[k] << "\t"<<arr_ct_more800[k] << "\t" << arr_ct_less800[k]<<endl;

	hv_cntr++;
	if (hv_cntr == 5)
	{
		setting++;
		hv_cntr =0;
	}

	file->Close();
        
    }


TString filename = "CT.root";
TFile * hfile = TFile::Open(filename,"RECREATE");
	
double hv_index[5]={1650,1700,1750,1800,1850};

TGraph *gTp15Th6 = new TGraph(HV, hv_index, arr_ct_all[0]);
TGraph *gTp20Th6 = new TGraph(HV, hv_index, arr_ct_all[1]);
TGraph *gTp35Th6 = new TGraph(HV, hv_index, arr_ct_all[2]);
TGraph *gTp15Th20 = new TGraph(HV, hv_index, arr_ct_all[3]);
TGraph *gTp20Th20 = new TGraph(HV, hv_index, arr_ct_all[4]);
TGraph *gTp35Th20 = new TGraph(HV, hv_index, arr_ct_all[5]);

TGraph *gTp15Th6_more = new TGraph(HV, hv_index, arr_ct_more800[0]);
TGraph *gTp20Th6_more = new TGraph(HV, hv_index, arr_ct_more800[1]);
TGraph *gTp35Th6_more = new TGraph(HV, hv_index, arr_ct_more800[2]);
TGraph *gTp15Th20_more = new TGraph(HV, hv_index, arr_ct_more800[3]);
TGraph *gTp20Th20_more = new TGraph(HV, hv_index, arr_ct_more800[4]);
TGraph *gTp35Th20_more = new TGraph(HV, hv_index, arr_ct_more800[5]);

TGraph *gTp15Th6_less = new TGraph(HV, hv_index, arr_ct_less800[0]);
TGraph *gTp20Th6_less = new TGraph(HV, hv_index, arr_ct_less800[1]);
TGraph *gTp35Th6_less = new TGraph(HV, hv_index, arr_ct_less800[2]);
TGraph *gTp15Th20_less = new TGraph(HV, hv_index, arr_ct_less800[3]);
TGraph *gTp20Th20_less = new TGraph(HV, hv_index, arr_ct_less800[4]);
TGraph *gTp35Th20_less = new TGraph(HV, hv_index, arr_ct_less800[5]);

TMultiGraph *mg1 = new TMultiGraph();
TMultiGraph *mg2 = new TMultiGraph();
TMultiGraph *mg3 = new TMultiGraph();

C1->cd();

gTp15Th6->SetMarkerStyle(8);
gTp15Th6->SetMarkerSize(1);
gTp15Th6->SetMarkerColor(kOrange+1);
gTp15Th6->SetLineColor(kOrange+1);
gTp15Th6->SetLineWidth(2);

gTp20Th6->SetMarkerStyle(3);
gTp20Th6->SetMarkerSize(1.5);
gTp20Th6->SetMarkerColor(kRed-7);
gTp20Th6->SetLineColor(kRed-7);
gTp20Th6->SetLineWidth(2); 

gTp35Th6->SetMarkerStyle(22);
gTp35Th6->SetMarkerSize(1.5);
gTp35Th6->SetMarkerColor(kBlue-5);
gTp35Th6->SetLineColor(kBlue-5);
gTp35Th6->SetLineWidth(2); 

gTp15Th20->SetMarkerStyle(8);
gTp15Th20->SetMarkerSize(1);
gTp15Th20->SetMarkerColor(kCyan-5);
gTp15Th20->SetLineColor(kCyan-5);
gTp15Th20->SetLineWidth(2);

gTp20Th20->SetMarkerStyle(3);
gTp20Th20->SetMarkerSize(1.5);
gTp20Th20->SetMarkerColor(kTeal-7);
gTp20Th20->SetLineColor(kTeal-7);
gTp20Th20->SetLineWidth(2); 

gTp35Th20->SetMarkerStyle(22);
gTp35Th20->SetMarkerSize(1.5);
gTp35Th20->SetMarkerColor(kGray+2);
gTp35Th20->SetLineColor(kGray+2);
gTp35Th20->SetLineWidth(2);

mg1->Add(gTp15Th6,"LP");
mg1->Add(gTp20Th6,"LP");
mg1->Add(gTp35Th6,"LP");
mg1->Add(gTp15Th20,"LP");
mg1->Add(gTp20Th20,"LP");
mg1->Add(gTp35Th20,"LP");

   
mg1->Draw("AL");     

TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);

//leg->SetHeader("Input High Voltage");
leg->SetFillColor(1);
leg->AddEntry(gTp15Th6,"Tp15Th6","lep");
leg->AddEntry(gTp20Th6,"Tp20Th6","lep");
leg->AddEntry(gTp35Th6,"Tp35Th6","lep");
leg->AddEntry(gTp15Th20,"Tp15Th20","lep");
leg->AddEntry(gTp20Th20,"Tp20Th20","lep");
leg->AddEntry(gTp35Th20,"Tp35Th20","lep");
leg->SetFillStyle(0);
leg->Draw();

mg1->GetXaxis()->SetTitle("High Voltage [v]");
mg1->GetYaxis()->SetTitle("Cross_talk_all");
mg1->SetTitle("Cross_Talk_all %");

//mg1->Write();
C1->Write();

C2->cd();

gTp15Th6_more->SetMarkerStyle(8);
gTp15Th6_more->SetMarkerSize(1);
gTp15Th6_more->SetMarkerColor(kOrange+1);
gTp15Th6_more->SetLineColor(kOrange+1);
gTp15Th6_more->SetLineWidth(2);

gTp20Th6_more->SetMarkerStyle(3);
gTp20Th6_more->SetMarkerSize(1.5);
gTp20Th6_more->SetMarkerColor(kRed-7);
gTp20Th6_more->SetLineColor(kRed-7);
gTp20Th6_more->SetLineWidth(2); 

gTp35Th6_more->SetMarkerStyle(22);
gTp35Th6_more->SetMarkerSize(1.5);
gTp35Th6_more->SetMarkerColor(kBlue-5);
gTp35Th6_more->SetLineColor(kBlue-5);
gTp35Th6_more->SetLineWidth(2); 

gTp15Th20_more->SetMarkerStyle(8);
gTp15Th20_more->SetMarkerSize(1);
gTp15Th20_more->SetMarkerColor(kCyan-5);
gTp15Th20_more->SetLineColor(kCyan-5);
gTp15Th20_more->SetLineWidth(2);

gTp20Th20_more->SetMarkerStyle(3);
gTp20Th20_more->SetMarkerSize(1.5);
gTp20Th20_more->SetMarkerColor(kTeal-7);
gTp20Th20_more->SetLineColor(kTeal-7);
gTp20Th20_more->SetLineWidth(2); 

gTp35Th20_more->SetMarkerStyle(22);
gTp35Th20_more->SetMarkerSize(1.5);
gTp35Th20_more->SetMarkerColor(kGray+2);
gTp35Th20_more->SetLineColor(kGray+2);
gTp35Th20_more->SetLineWidth(2);

mg2->Add(gTp15Th6_more,"LP");
mg2->Add(gTp20Th6_more,"LP");
mg2->Add(gTp35Th6_more,"LP");
mg2->Add(gTp15Th20_more,"LP");
mg2->Add(gTp20Th20_more,"LP");
mg2->Add(gTp35Th20_more,"LP");

mg2->Draw("AL"); 
leg->Draw();

mg2->GetXaxis()->SetTitle("High Voltage [v]");
mg2->GetYaxis()->SetTitle("Cross_talk_more %");
mg2->SetTitle("Cross_Talk_more");

C2->Write();
//mg2->Write();

C3->cd();

gTp15Th6_less->SetMarkerStyle(8);
gTp15Th6_less->SetMarkerSize(1);
gTp15Th6_less->SetMarkerColor(kOrange+1);
gTp15Th6_less->SetLineColor(kOrange+1);
gTp15Th6_less->SetLineWidth(2);

gTp20Th6_less->SetMarkerStyle(3);
gTp20Th6_less->SetMarkerSize(1.5);
gTp20Th6_less->SetMarkerColor(kRed-7);
gTp20Th6_less->SetLineColor(kRed-7);
gTp20Th6_less->SetLineWidth(2); 

gTp35Th6_less->SetMarkerStyle(22);
gTp35Th6_less->SetMarkerSize(1.5);
gTp35Th6_less->SetMarkerColor(kBlue-5);
gTp35Th6_less->SetLineColor(kBlue-5);
gTp35Th6_less->SetLineWidth(2); 

gTp15Th20_less->SetMarkerStyle(8);
gTp15Th20_less->SetMarkerSize(1);
gTp15Th20_less->SetMarkerColor(kCyan-5);
gTp15Th20_less->SetLineColor(kCyan-5);
gTp15Th20_less->SetLineWidth(2);

gTp20Th20_less->SetMarkerStyle(3);
gTp20Th20_less->SetMarkerSize(1.5);
gTp20Th20_less->SetMarkerColor(kTeal-7);
gTp20Th20_less->SetLineColor(kTeal-7);
gTp20Th20_less->SetLineWidth(2); 

gTp35Th20_less->SetMarkerStyle(22);
gTp35Th20_less->SetMarkerSize(1.5);
gTp35Th20_less->SetMarkerColor(kGray+2);
gTp35Th20_less->SetLineColor(kGray+2);
gTp35Th20_less->SetLineWidth(2);

mg3->Add(gTp15Th6_less,"LP");
mg3->Add(gTp20Th6_less,"LP");
mg3->Add(gTp35Th6_less,"LP");
mg3->Add(gTp15Th20_less,"LP");
mg3->Add(gTp20Th20_less,"LP");
mg3->Add(gTp35Th20_less,"LP");

mg3->Draw("AL"); 
leg->Draw();

mg3->GetXaxis()->SetTitle("High Voltage [v]");
mg3->GetYaxis()->SetTitle("Cross_talk_less %");
mg3->SetTitle("Cross_Talk_less");

//mg3->Write();

C3->Write();

}
