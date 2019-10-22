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
#include "TMultiGraph.h"

using namespace std;

const Int_t max2 =5;


Bool_t resolution()
{
	
TCanvas *c = new TCanvas("c","c");
	TString talysName = "res.txt";
	
	fstream talysInput(talysName,ios::in);

if(!talysInput.is_open())
{
cout << "Could not open"<< talysName << endl;
return kFALSE;
}

TString filename = "restn.root";
TFile * hfile = TFile::Open(filename,"RECREATE");

double tp15th6[5];
double tp20th6[5];
double tp35th6[5];

double tp15th20[5];
double tp20th20[5];
double tp35th20[5];
double hv[5]={1650,1700,1750,1800,1850};

//const Int_t max1 =14;


//for(Int_t i=0; i<1; i++)
//{
//header.ReadLine(talysInput);
//}

for(Int_t i=0; i<max2; i++)
{
	talysInput >> tp15th6[i] >> tp20th6[i] >> tp35th6[i] >> tp15th20[i] >> tp20th20[i] >> tp35th20[i];

}

talysInput.close();
int index =0; 



TGraph *gTp15Th6 = new TGraph(max2, hv, tp15th6);
TGraph *gTp20Th6 = new TGraph(max2, hv, tp20th6);
TGraph *gTp35Th6 = new TGraph(max2, hv, tp35th6);
TGraph *gTp15Th20 = new TGraph(max2, hv, tp15th20);
TGraph *gTp20Th20 = new TGraph(max2, hv, tp20th20);
TGraph *gTp35Th20 = new TGraph(max2, hv, tp35th20);

TMultiGraph *mg1 = new TMultiGraph();
TMultiGraph *mg2 = new TMultiGraph();


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

// //gTalys3->Draw("ALSame");

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

// c->SetLogx();
// c->SetLogy();

mg1->GetXaxis()->SetTitle("High Voltage [v]");
mg1->GetYaxis()->SetTitle("sigma [um]");
mg1->SetTitle("Resolution");
     
//      gPad->Modified();
 mg1->Write();


return kTRUE;
}
