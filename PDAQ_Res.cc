#include "PDAQ_Spl_Res.h"
#include "string.h"
#include "TStyle.h"
#include <bits/stdc++.h>

using namespace std;

string convertToString (char *a , int size)
{
    string s(a);
    
    return s;
}

Bool_t PDAQ_Res ( char* intree, char* outtree, int maxEvents )

{
 
    printf("Tree Opened\n");

    TFile inFile ( intree );
    TTree* tree = ( TTree* ) inFile.Get ( "PDAQ_tree" );
    if ( !tree ) {
        std::cerr << "Tree PDAQ_tree was not found\n";
        std::exit ( 1 );
    }
    tree->Print();

    TFile* driftfile1 = new TFile ( outtree, "RECREATE" );   
    TF1* f1 = new TF1 ( "f1", "pol1" );
    
    std::vector<double>* vec_Drifttime = 0;
    std::vector<double>* vec_x = 0;
    std::vector<double>* vec_y = 0;
    std::vector<double>* vec_z = 0;
    std::vector<double>* vec_layer = 0;
    std::vector<double>* vec_straw = 0;
    std::vector<double>* vec_plane = 0;
    std::vector<double>* vec_tot = 0;
    
    tree->SetBranchAddress ( "vec_Drifttime", &vec_Drifttime );
    tree->SetBranchAddress ( "vec_x", &vec_x );
    tree->SetBranchAddress ( "vec_y", &vec_y );
    tree->SetBranchAddress ( "vec_z", &vec_z );
    tree->SetBranchAddress ( "vec_layer", &vec_layer );
    tree->SetBranchAddress ( "vec_straw", &vec_straw );
    tree->SetBranchAddress ("vec_plane",&vec_plane);
    tree->SetBranchAddress ("vec_tot",&vec_tot);
    
    TGraph* gDR;
    gDR = ( TGraph* ) inFile.Get ( "PDAQ_DR" );
    
    Int_t iev = ( Int_t ) tree->GetEntries();
    printf("\n Number of entries in the tree : %d \n",iev);
    
    
    
    return kTRUE;
}

int main ( int argc, char** argv )
{

    if ( argc >= 3 )

        if ( !argv[3] ) 
        {
            printf ( "\n\nNote : One Million Events will be processed. To change "
                     "add the number of events to be processed after the ouput "
                     "file name.\n" );
            sleep ( 2 );
            PDAQ_Res ( argv[1], argv[2], 100000000 );
        } 
        else 
        {
            PDAQ_Res ( argv[1], argv[2], atoi ( argv[3] ) );
        }

    else 
    {
        return 1;
    }

    return 0;
}
