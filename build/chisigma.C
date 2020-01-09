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



Bool_t chisigma()
{
	
//TCanvas *c = new TCanvas("c","c");

TString talysName = "a.txt";
int no_of_tracks = count_line(talysName);
fstream talysInput(talysName,ios::in);
int no_of_sigs = 16;


if(!talysInput.is_open())
{
cout << "Could not open"<< talysName << endl;
return kFALSE;
}

TString filename = "chi.root";
TFile * hfile = TFile::Open(filename,"RECREATE");


double chi_sig[no_of_sigs][no_of_tracks];

for(int a=0; a<no_of_sigs; a++)
{
	for(int b=0; b<no_of_tracks; b++)
	{
		chi_sig[a][b] =0;
	}
}

TH1* h_chi[no_of_sigs];

for (int a=0; a<no_of_sigs; a++)
{
	h_chi[a] = new TH1F ( Form("h_chi%i",(a*10)+150), Form("h_chi%i;chi^{2}",(a*10)+150), 35, 0, 35 );
}

for(Int_t i=0; i<no_of_tracks; i++)
{
	talysInput >> chi_sig[0][i] >> chi_sig[1][i] >> chi_sig[2][i] >> chi_sig[3][i] >> chi_sig[4][i] >> chi_sig[5][i] >> chi_sig[6][i] >> chi_sig[7][i] >> chi_sig[8][i] >> chi_sig[9][i] >> chi_sig[10][i] >> chi_sig[11][i] >> chi_sig[12][i] >> chi_sig[13][i] >> chi_sig[14][i] >> chi_sig[15][i];

}
cout<<"txt file read"<<endl;
talysInput.close();

for (int a=0; a<no_of_sigs; a++)
{
	for(int b=0; b<no_of_tracks; b++)
	{
		if (chi_sig[a][b] ==0 )continue;
		h_chi[a]->Fill(chi_sig[a][b]);
	}

}

TF1* chi2_6dofF = new TF1("chi2_6dofF","ROOT::Math::chisquared_pdf(x,6,0)",0,35);

double sig_arr[no_of_sigs];
double chi_arr[no_of_sigs];

for (int a=0; a<no_of_sigs; a++)
{	
	h_chi[a]->Scale(1/h_chi[a]->Integral(0,35));  
    
    	double chiResult = h_chi[a]->Chisquare(chi2_6dofF);
    	cout<<"Sigma :"<<(a*10)+150<<"\tChi result :"<<chiResult<<endl;
	sig_arr[a]=(a*10)+150;
	chi_arr[a]=chiResult;
	h_chi[a]->Write();

}
chi2_6dofF->Write();

TGraph* g_chi_sig = new TGraph ( no_of_sigs, sig_arr, chi_arr );
g_chi_sig->Draw("LP");
g_chi_sig->Write();
return kTRUE;
}
