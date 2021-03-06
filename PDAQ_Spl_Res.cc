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

Bool_t PDAQ_Spl_Res ( char* intree, char* outtree, int maxEvents )

//Bool_t PDAQ_Spl_Res ( void )
{

    cout << "Opened" << endl;


    TFile inFile ( intree );
    TTree* tree = ( TTree* ) inFile.Get ( "PDAQ_tree" );
    if ( !tree ) {
        std::cerr << "Tree PDAQ_tree was not found\n";
        std::exit ( 1 );
    }

//     TFile* file = TFile::Open ( "d15.root", "READ" );
//     TTree* tree = 0;
//     file->GetObject ( "PDAQ_tree", tree );

    //     TFile * file = TFile::Open("PDAQ_Stt_Tracks.root", "READ");
    //     TTree* tree = 0;
    //     file->GetObject("PDAQ_EMC_STT_cluster_analysis", tree);
//     if ( !tree ) {
//         std::cerr << "Tree doesn't exists" << std::endl;
//         return 1;
//     }
    tree->Print();

    TF1* f1 = new TF1 ( "f1", "pol1" );
//     TF1* f2 = new TF1 ( "f2", "pol1" );
//     TF1* f3 = new TF1 ( "f3", "pol1" );
//     TF1* f4 = new TF1 ( "f4", "pol1" );

    TFile* driftfile1 = new TFile ( outtree, "RECREATE" );
    // TFile* driftfile2 = new TFile("Y_Coordinates.root", "RECREATE");
    Double_t theta_X = 0;
    
    TH1F* hx = new TH1F ( "hx", "Spatial_Resolution_X", 500, -5, 5 );
    TH1F* hdx = new TH1F ( "hdx", "res_non_corrected", 100, -10, 10 );
    TH1F* hy = new TH1F ( "hy", "Spatial_Resolution_X", 500, -1, 4 );
    TH1F* DR = new TH1F ( "dr", "drift_radius", 350, 0, 350 );
    TH1F* h_theeta_X = new TH1F ( "h_theeta_X", "h_theeta_X", 400, -2, 2 );
    
    TH1F* hx_slope = new TH1F ( "hx_slope", "Slope_X", 1000, -0.5, 0.5 );
    TH1F* hx_const = new TH1F ( "hx_const", "Constant_X", 1000, 0, 50 );

    TH1F* Z_value = new TH1F ( "Z_value", "dummy", 5, 0, 5 );
    
    TH2F* h_XfvsZ = new TH2F ( "h_XfvsZ", "h_XfvsZ;X of track fit X_{f} [mm]; Z [mm]", 500, 0, 500,7,-10,600 );

    
    TH1F* h_dt[2];
    TH1F* h_tot[2];
    TH1F* h_thet_plane[2];
    TH1F* h_thet_res[2];
    TH1F* h_orientation[2];
    TH1F* h_thet_straw_x[2];

    
    
    TH1F* h_dx[18];
    TH1F* h_res_plane[18];
    TH2F* h_Dr_vs_dr[18];
    TH1F* h_plane_str_no[18];

    
    for (int mb=0; mb<2; mb++){
        h_dt[mb] = new TH1F (Form("h_dt%i",mb+1),Form("h_dt%i;Drift Time [ns]",mb+1),500, 0, 500);
        h_tot[mb] = new TH1F (Form("h_tot%i",mb+1),Form("h_tot%i;TOT [ns]",mb+1),500, 0, 500);
        h_thet_plane[mb] = new TH1F (Form("h_thet_plane%i",mb+1),Form("h_thet_plane%i;Plane No",mb+1),18, 0,18);
        h_thet_res[mb] = new TH1F (Form("h_thet_res%i",mb+1),Form("h_thet_res%i;dr [mm]",mb+1),500, -5, 5);
        h_orientation[mb] = new TH1F (Form("h_orientation%i",mb+1),Form("h_orientation%i;L/R",mb+1),3, 0, 3);
        h_thet_straw_x[mb] = new TH1F (Form("h_thet_straw_x%i",mb+1),Form("h_thet_straw_x%i;dr [mm]",mb+1),100, 0, 50);
        
            
    }
    
    for ( int mc = 0; mc < 18; mc++ ) {
        //h_dx[mc] = new TH2F ( Form("h_dx_Plane_%d",mc+1), Form("h_dx_Plane_%d",mc+1), 1000, -1, 1,100,0,10 );
        h_res_plane[mc] = new TH1F ( Form("h_res_plane%d",mc+1), Form("h_res_plane%d;dr [mm]",mc+1), 500, -5, 5 );
        h_dx[mc] = new TH1F ( Form("h_dx_Plane_%d",mc+1), Form("h_dx_Plane_%d;dx [mm]",mc+1),100,-10,10 );
        h_Dr_vs_dr[mc] = new TH2F ( Form("h_Dr_vs_dr_%d",mc+1), Form("h_Dr_vs_dr_%d;Drift radius [mm];dr [mm]",mc+1),60, 0,6,100,-1,1 );
        
        h_plane_str_no[mc] = new TH1F (Form("h_plane_str_no%i",mc+1),Form("h_plane_str_no%i;Straw No.",mc+1),60, 0, 60);

    }
    
    TCanvas * Ct; //Canvas for ToT
    Ct=new TCanvas ( "Ct","Ct" );
    
    TCanvas * xf_z; //Canvas for ToT
    xf_z=new TCanvas ( "xf_z","xf_z" );

    std::vector<double>* vec_Drifttime = 0;
    std::vector<double>* vec_x = 0;
    std::vector<double>* vec_y = 0;
    std::vector<double>* vec_z = 0;
    std::vector<double>* vec_layer = 0;
    std::vector<double>* vec_straw = 0;
    std::vector<double>* vec_plane = 0;
    std::vector<double>* vec_tot = 0;

    //     std::vector<double>* vec_module;
    //     std::vector<double>* vec_fee;
    //     std::vector<double>* vec_fee_ch;
    //     std::vector<double>* vec_tdc_ch;
    std::vector<SttHit*> vec_dr_xaxis[2];
    std::vector<SttHit*> vec_dr_yaxis[2];
    std::vector<double> vec_driftradius;

    std::vector<double> vec_emcX;
    std::vector<double> vec_emcY;

    std::vector<double> vec_Xtheta;
    std::vector<double> vec_Ytheta;

    std::vector<double> vec_Xsr;
    std::vector<double> vec_Ysr;
    
    

    Double_t a1[230];
    Double_t b1[230];

    Double_t zaxis[16];

    Double_t A1[200];
    Double_t A2[200];

    Double_t B1[200];
    Double_t B2[200];

    Double_t xperfect3[16];
    Double_t yperfect3[16];

    Int_t array_length = 8;
    Int_t pair_length = 2;
    Int_t total_pairs = ( array_length / pair_length );
    Int_t no_of_pairs = ( total_pairs * total_pairs );

    tree->SetBranchAddress ( "vec_Drifttime", &vec_Drifttime );
    tree->SetBranchAddress ( "vec_x", &vec_x );
    tree->SetBranchAddress ( "vec_y", &vec_y );
    tree->SetBranchAddress ( "vec_z", &vec_z );
    tree->SetBranchAddress ( "vec_layer", &vec_layer );
    tree->SetBranchAddress ( "vec_straw", &vec_straw );
    tree->SetBranchAddress ("vec_plane",&vec_plane);
    tree->SetBranchAddress ("vec_tot",&vec_tot);


    TGraph* gDR = new TGraph ( 220, a1, b1 );
    gDR = ( TGraph* ) inFile.Get ( "PDAQ_DR" );
   // gDR = ( TGraph* ) file->Get ( "PDAQ_DR" );

    // Double_t strawY[8];
    // Double_t strawZ[8];

    Double_t yperfectmean[16];

    SttHit* new_y1 = 0;
    SttHit* new_y2 = 0;
    int x_min = -10;
    int x_max = 90;
    int x_bins = ( x_max - x_min ) * 5;
    int z_min = -10;
    int z_max = 80;
    int z_bins = ( z_max - z_min ) * 5;
    Int_t iev = ( Int_t ) tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl << endl;

    vector<SttHit*> vec_hits;
    vector<SttHit*> vec_hits_new;
    vec_driftradius.clear();
    vec_emcX.clear();
    vec_emcY.clear();
    // loop over all the vectors in the tree.

    TCanvas* c1 = new TCanvas ( "c1","c1" );
    TH2I* h = new TH2I ( "name", "h;x,y [cm];z [cm]", x_bins, x_min,x_max, z_bins, z_min, z_max );
    
    //c1->cd();
    h->Draw();
    for ( Int_t i = 0; i < iev; i++ ) {
        Int_t count = 0;
        Int_t count1 = 0;
        
        if ( i == maxEvents ) {
            break;
        }
        if ( i % 1000 == 0 ) {
            cout << i << endl;
        }

        tree->GetEntry ( i );
        cout << endl;
        cout << "Track No. " << i << endl;
        Int_t oiv = vec_Drifttime->size();
        // cout << "vecsize  :" << oiv << endl;

        if ( vec_Drifttime->size() > 16 ) {
            cout << "SIZE GREATER THAN 16" << endl;
        }

        for ( int w = 0; w < 2; w++ ) {
            vec_dr_xaxis[w].clear();
            vec_dr_yaxis[w].clear();
            for ( int v = 0; v < vec_dr_xaxis[0].size(); v++ ) {
                delete ( vec_dr_xaxis[w][v] );
                delete ( vec_dr_yaxis[w][v] );
            }

            // vector<SttHit*> vec_hits;
            vec_hits.clear();
            vec_hits_new.clear();



            char buff[200];
            sprintf ( buff, "can_%d", i );
        }
        //TCanvas* c1 = new TCanvas ( buff );

        // sprintf ( buff, "h_%d", i );


        /*        TH2I* h = new TH2I ( buff, "h;x,y [cm];z [cm]", x_bins, x_min,x_max, z_bins, z_min, z_max );
                //c1->cd();
                //h->Draw()*/;
        // c1->Range(-10,-10,60,60);
        //Loop over the vector elements

        for ( int n = 0; n < oiv; n++ ) {
            SttHit* a = new SttHit();
            a->drifttime = vec_Drifttime->at ( n );
            a->x = vec_x->at ( n );
            a->y = vec_y->at ( n );
            a->z = vec_z->at ( n );
            a->layer = vec_layer->at ( n );
            a->tot = vec_tot->at(n);
            //a->fee = vec_fee->at(n);
            //a->fee_channel = vec_fee_ch->at(n);
            //a->channel = vec_tdc_ch->at(n);
            a->plane = vec_plane->at(n);
            a->straw = vec_straw->at ( n );

            vec_hits.push_back ( a );

            float uuu = 0.0;
//             bool is_y = ( a->layer % 2 == 0 );
//             if ( is_y ) {
//                 uuu = a->y;
//             } else {
//                 uuu = a->x;
//             }
            if ( a->layer==1 || a->layer==4 ||a->layer== 5||a->layer==8 ) {
                uuu = a->x;
            } else {
                uuu= a->y;
            }



            //TEllipse* straws = new TEllipse ( uuu, a->z, 0.505, 0.505 );
            //TBox* box = new TBox(25,60,43,78);
            //straws->SetFillColor ( 37 );

            //cout<<"Check :"<<a->drifttime<<"\t"<<a->layer<<"\t"<<a->straw<<"\t"<<a->x<<"\t"<<a->y<<"\t"<<a->z<<endl;
            ///////

//             cout << "\t" << vec_hits.size() << "\t" << vec_hits[n]->drifttime
//             << "\t" << vec_hits[n]->x << "\t" << vec_hits[n]->y << "\t" <<
//             vec_hits[n]->z << "\t" << endl;
            //c1->cd();
            //straws->Draw ( "same" );
            //h->Fill ( uuu, a->z );
            //box->Draw("same");

        }

        // Filter to get only the sets having unique Z values
        Int_t oc_count = 0;

        if ( vec_hits.size() > 0 ) {
            for ( Int_t u = 0; u < vec_hits.size() - 1; u++ ) {
                for ( Int_t uu = u + 1; uu < vec_hits.size(); uu++ ) {
                    if ( fabs ( vec_hits[u]->z - vec_hits[uu]->z ) < 0.5 ) {
                        oc_count++;
                    }
                    // cout<<vec_hits[u]->z<<"\t"<<vec_hits[uu]->z<<endl;
                }
            }
            // cout<<oc_count<<endl;

            TGraph* xfit;
            TGraph* xfit1;
            TGraph* yfit;
            TGraph* yfit1;

            Double_t strawX[8];
            Double_t strawY[8];
            Double_t strawZ[8];
            std::vector<double> vec_strawX;
            std::vector<double> vec_strawY;
            std::vector<double> vec_strawZx;
            std::vector<double> vec_strawZy;

            Z_value->Fill ( oc_count );

            // Take only the events having unique Z coordinate

            if ( oc_count >= 0 ) {
                // cout << "oc_count :";
                // cout<<oc_count<<"\t"<<vec_hits.size()<<endl;
                vec_hits_new = vec_hits;
                // cout<<vec_hits_new.size()<<endl;;
                //Z_value->Fill(1);
            }

            else {

                //Z_value->Fill(2);

                continue;
            }

            ///////////////////////////////////////////////////////////////////

            // Draw an ellipse on the straw position
            vec_strawX.clear();
            vec_strawY.clear();
            vec_strawZx.clear();
            vec_strawZy.clear();
            for ( Int_t j = 0; j < vec_hits_new.size(); j++ ) {
                //
                //                 //cout<<"THIS IS Z
                //                 :"<<vec_hits_new[j]->z<<endl;
                //                 strawX[j] = vec_hits_new[j]->x;
                //                 strawY[j] = vec_hits_new[j]->y;
                //                 strawZ[j] = vec_hits_new[j]->z;
                //
                //                 float uuu = 0.0;
                //                 bool is_y = (vec_hits_new[j]->layer % 2 ==
                //                 0);
                //                 if (is_y)
                //                     uuu = vec_hits_new[j]->y;
                //                 else
                //                     uuu = vec_hits_new[j]->x;
                //                 straws2 = new TEllipse(uuu,
                //                 vec_hits_new[j]->z, 0.505, 0.505);
                //                 if (is_y)
                //                     straws2->SetFillColor(46);
                //                 else
                //                     straws2->SetFillColor(38);
                //                 straws2->SetFillStyle(1001);
                //                 straws2->SetLineColor(1);
                //                 straws2->SetLineWidth(1);
                //                 straws2->Draw("same");
                //
                //                 // Get the drift radius for the corresponding
                //                 drift time for the hits having x-coordinates

                if ( vec_hits_new[j]->layer == 1 || vec_hits_new[j]->layer == 4 ||vec_hits_new[j]->layer == 5 || vec_hits_new[j]->layer == 8 ) {

                    Double_t dr1 = gDR->Eval ( vec_hits_new[j]->drifttime );
                    SttHit* new_x1 = new SttHit();
                    SttHit* new_x2 = new SttHit();

                    vec_driftradius.push_back ( dr1 );

                    new_x1->x = vec_hits_new[j]->x + dr1;
                    new_x1->y = vec_hits_new[j]->y;
                    new_x1->z = vec_hits_new[j]->z;
                    new_x1->layer = vec_hits_new[j]->layer;
                    new_x1->plane = vec_hits_new[j]->plane;
                    new_x1->straw = vec_hits_new[j]->straw;
                    new_x1->tot = vec_hits_new[j]->tot;
                    new_x1->drifttime = vec_hits_new[j]->drifttime;
                    //                     vec_hits_new[j]->module;
                    //                     new_x1->fee = vec_hits_new[j]->fee;
                    //                     new_x1->fee_channel =
                    //                     vec_hits_new[j]->fee_channel;

                    vec_strawX.push_back ( vec_hits_new[j]->x );
                    vec_strawZx.push_back ( vec_hits_new[j]->z );

                    new_x2->x = vec_hits_new[j]->x - dr1;
                    new_x2->y = vec_hits_new[j]->y;
                    new_x2->z = vec_hits_new[j]->z;
                    new_x2->layer = vec_hits_new[j]->layer;
                    new_x2->plane = vec_hits_new[j]->plane;
                    new_x2->straw = vec_hits_new[j]->straw;
                    new_x2->tot = vec_hits_new[j]->tot;
                    new_x2->drifttime = vec_hits_new[j]->drifttime;
                    //                     new_x2->module =
                    //                     vec_hits_new[j]->module;
                    //                     new_x2->fee = vec_hits_new[j]->fee;
                    //                     new_x2->fee_channel =
                    //                     vec_hits_new[j]->fee_channel;
                    //
                    vec_dr_xaxis[0].push_back ( new_x1 );
                    vec_dr_xaxis[1].push_back ( new_x2 );
                    // zaxis[j]= vec_hits[j]->z;

                    // cout<< "Xaxis :
                    // "<<vec_hits[j]->layer<<"\t"<<vec_hits[j]->x<<"\t"<<vec_hits[j]->y<<"\t"<<vec_hits[j]->z<<"\t"<<vec_hits[j]->x<<"\t"<<dr1<<endl;
                    count++;
                }

                // Get the drift radius for the corresponding drift time for the
                // hits having y-coordinates

                if ( vec_hits_new[j]->layer == 2 || vec_hits_new[j]->layer == 3 || vec_hits_new[j]->layer == 6 || vec_hits_new[j]->layer == 7 ) {

                    Double_t dr2 = gDR->Eval ( vec_hits_new[j]->drifttime );

                    SttHit* new_y1 = new SttHit();
                    SttHit* new_y2 = new SttHit();

                    new_y1->y = vec_hits_new[j]->y + dr2;
                    new_y1->z = vec_hits_new[j]->z;
                    new_y1->x = vec_hits_new[j]->x;
                    new_y1->layer = vec_hits_new[j]->layer;
                    new_y1->plane = vec_hits_new[j]->plane;
                    new_y1->straw = vec_hits_new[j]->straw;
                    new_y1->tot = vec_hits_new[j]->tot;
                    new_y1->drifttime = vec_hits_new[j]->drifttime;
                    //                     new_y1->module =
                    //                     vec_hits_new[j]->module;
                    //                     new_y1->fee = vec_hits_new[j]->fee;
                    //                     new_y1->fee_channel =
                    //                     vec_hits_new[j]->fee_channel;
                    //
//                     if (vec_hits_new[j]->y >0)
//                     {
//                     strawY[j] = vec_hits_new[j]->y;
//                     }
                    // cout<< "CHECK ::  "<<vec_hits_new[j]->y<<endl;
                    vec_strawY.push_back ( vec_hits_new[j]->y );
                    vec_strawZy.push_back ( vec_hits_new[j]->z );

                    new_y2->y = vec_hits_new[j]->y - dr2;
                    new_y2->z = vec_hits_new[j]->z;
                    new_y2->x = vec_hits_new[j]->x;
                    new_y2->layer = vec_hits_new[j]->layer;
                    new_y2->plane = vec_hits_new[j]->plane;
                    new_y2->straw = vec_hits_new[j]->straw;
                    new_y2->tot = vec_hits_new[j]->tot;
                    new_y2->drifttime = vec_hits_new[j]->drifttime;
                    //                     new_y2->module =
                    //                     vec_hits_new[j]->module;
                    //                     new_y2->fee = vec_hits_new[j]->fee;
                    //                     new_y2->fee_channel =
                    //                     vec_hits_new[j]->fee_channel;

                    // SttHit* new_y1=vec_hits[j]->y+dr2;
                    // SttHit* new_y2=vec_hits[j]->y-dr2;
                    vec_dr_yaxis[0].push_back ( new_y1 );
                    vec_dr_yaxis[1].push_back ( new_y2 );

                    // yperfectmean[j] = dr2;
                    // zaxis[j]= vec_hits[j]->z;
                    // cout <<"\t"<< "Yaxis
                    // :"<<vec_hits[j]->y<<vec_hits[j]->layer<<endl;
                    //<<"\t"<<vec_hits[j]->layer<<"\t"<<dr2<<endl;

                    count1++;
                    // delete new_y1;
                    // delete new_y2;
                }
            }

//             for (Int_t jx=0; jx<4; jx++)
//             {
//                  cout <<"Y here  :"<<" "<<strawY[jx]<<endl;
//
//             }

            /////////////////////////////////////////FIND THE BEST
            /// COMBINATION//////////////////////////////////////////////////////

            // COMBINATIONS to get least ChiSquar

            Int_t e_index = 0;
            std::vector<Double_t> chiX_array;
            Double_t chi_value;

            vector<vector<SttHit*>> myCombination;
            myCombination.clear();

            Double_t dSlope0 = 0;
            Double_t dConst0 = 0;
            Double_t X_perpX = 0;
            Double_t X_perpY = 0;
            Double_t X_short = 0;
            Double_t centerTotrack = 0;//from the center of the straw as a perpendicular
            Double_t theeta =0;
            Double_t dradius = 0;
            Double_t rotatedX =0;
            Double_t rotatedY =0;
            Double_t staticRes =0;

            Double_t residue =0;

            int no_comb_x = pow ( 2, count );
            chiX_array.resize ( no_comb_x );

            for ( Int_t co = 0; co < no_comb_x; co++ ) {
//                 cout<<no_comb_x<<endl;
//                 cout<<"{";

                std::vector<SttHit*> vt;
                vt.clear();

                for ( Int_t coo = 0; coo < count; coo++ ) {

                    int bit_idx = ( 1 << coo );
                    int comb_idx = ( co & bit_idx ) ? 1 : 0;

                    A1[coo] = vec_dr_xaxis[comb_idx][coo]->x;
                    A2[coo] = vec_dr_xaxis[comb_idx][coo]->z;

//                     cout<<"   STRAW_LOC@@@@@@@@@@@@@@@@&&&&&&&&&:"<<vec_dr_xaxis[comb_idx][coo]->x<<"   " ;

                    vt.push_back ( vec_dr_xaxis[comb_idx][coo] );
                }

                myCombination.push_back ( vt );

//                 cout<<"}";
//                 cout << "\n";
//
//                 cout <<"size1 A1 :"<< ( sizeof ( A1 ) ) /8<<endl;
//                 cout <<"size2 A2:"<< ( sizeof ( A2 ) ) /8<<endl;
                // cout <<"size3 xperfect2:"<<(sizeof(xperfect2))/8<<endl;

//                   for (Int_t cooo=0; cooo< myCombination;cooo++)
//                {
//                    for( Int_t cro=0; cro < vt.size(); cro++)
//                    {
//                        cout << "Hello"<< " ";
//                    }
//                }
                // cout<<"\n"<<endl;

                TGraph* chiX = new TGraph ( vt.size(), A1, A2 );
                chiX->Fit ( f1, "q" );
                chi_value = f1->GetChisquare();
                

                // cout<<"\nX chi value"<<chi_value<<endl;

                chiX_array[co] = chi_value/vt.size();
                delete chiX;
            }

            // Get the index of the combination having the least chisquare.
            Double_t smallest = chiX_array[0];

            Int_t chi_index = 0;
            for ( Int_t ci = 0; ci < no_comb_x; ci++ ) {
                //cout << "chiSquare  :" << chiX_array[ci] << endl;

                if ( smallest > chiX_array[ci] ) {
                    smallest = chiX_array[ci];
                    chi_index = ci;
                }
            }

            Double_t Ex[50];
            SttHit* hit[50];

//             cout << "\n\n"
//                  << "BEST COMBINATION  :";
//             for ( Int_t gdd = 0; gdd < count; gdd++ ) {
//                 //cout << "\t " << myCombination.at ( chi_index ).at ( gdd ) << " ";
//                 Ex[gdd] = myCombination.at ( chi_index ).at ( gdd );
//             }
            //cout << "\n\n\n" << endl;

            //cout << "\n SMALLEST ChiSquar  :" << smallest << "(" << chi_index<< ")" << endl;

            for ( Int_t coo = 0; coo < count; coo++ ) {
                int bit_idx = ( 1 << coo );
                int comb_idx = ( chi_index & bit_idx ) ? 1 : 0;
                
                Ex[coo] = myCombination.at ( chi_index ).at ( coo )->x;
                A1[coo] = vec_dr_xaxis[comb_idx][coo]->x;
                A2[coo] = vec_dr_xaxis[comb_idx][coo]->z;
                hit[coo] = myCombination.at ( chi_index ).at ( coo );
                
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            //             for (Int_t cdf = 0; cdf < vec_hits_new.size(); cdf++)
            //             {
            //                 cout << "  XX  :" << strawX[cdf] << "  ZZ  :" <<
            //                 strawZ[cdf] << endl;
            //             }

            // Plot X vs Z coordinate and fit the points

            if ( vec_dr_xaxis[0].size() != 0 ) {

               /* xfit = new TGraph ( count, Ex, A2 );
                xfit1 = new TGraph ( vec_hits_new.size(), strawX, strawZ );

                xfit->GetXaxis()->SetLimits ( x_min, x_max );
                xfit->GetYaxis()->SetLimits ( z_min, z_max ); */
               
                
                xfit = new TGraph ( count,A2, Ex );
                xfit1 = new TGraph ( vec_hits_new.size(), strawZ, strawX );

                xfit->GetYaxis()->SetLimits ( x_min, x_max );
                xfit->GetXaxis()->SetLimits ( z_min, z_max );

                xfit->GetXaxis()->SetTitle ( "X Axix [cm]" );
                xfit->GetYaxis()->SetTitle ( "Z Axix [cm]" );
                xfit->SetTitle ( "X-Coordinates" );

                xfit->SetMarkerStyle ( 7 );
                xfit->SetMarkerSize ( 1 );
                xfit->SetMarkerColor ( kRed + 2 );
                xfit->SetLineColor ( 0 );
                xfit->SetLineWidth ( 2 );

                xfit1->SetLineColor ( 0 );

                xfit1->SetMarkerStyle ( 24 );
                xfit1->SetMarkerSize ( 2 );
                xfit1->SetMarkerColor ( kRed + 2 );
                // xfit1->SetLineColor(kRed-10);
                // xfit1->SetLineWidth(2);

                xfit->Fit ( f1, "q" );

                f1 = xfit->GetFunction ( "f1" );
                Double_t trck_constant = f1->GetParameter ( 0 ); //constant - c
                Double_t trck_slope = f1->GetParameter ( 1 ); //slope    - m
                theta_X = ((atan ( trck_slope)) * 180) / 3.1415926 ;
                hx_slope->Fill(trck_slope);
                hx_const->Fill(trck_constant);
                // cout<< " THETA X  : "<< theta_X <<endl;
                // emcX_arry[i] = (70 - p0) / p1;
            //    vec_emcX.push_back ( ( 70 - p0 ) / p1 );
              //  cout<<"slope : "<<trck_slope<<"\tconstant : "<<trck_constant<<"\tTheta :"<<theta_X<<endl;
               // if (trck_slope > 80 || trck_slope < -80)
                h_theeta_X->Fill(theta_X);
              //  {
                    for ( Int_t k = 0; k < count; k++ ) {


                        //********************************** Calculate residual ******************************************
                        
                        centerTotrack= fabs((trck_slope * vec_strawZx[k]) - vec_strawX[k] + trck_constant) / ( sqrt ( 1 + ( trck_slope * trck_slope ) ) ) ; //((m * Px) - Py + C)/(sqrt(1+(m*m))) distance b/w a point and line. reference in the downloaded pdf.
                        
                        dradius = fabs ( vec_strawX[k] - A1[k] );
                        
                        residue = (centerTotrack-dradius) * 10; // Cm to mm
                        
                        //********************************** Calculate residual ******************************************
                       // cout<<centerTotrack<<"\t"<<dradius<<"\t"<<residue<<endl;
                   //    cout<<"layer : "<<hit[k]->layer<<"  plane :  "<<hit[k]->plane<<"  straw : "<<hit[k]->straw<<"  tot : "<<hit[k]->tot<<"  dt  : "<<hit[k]->drifttime<<"  c_t : "<<centerTotrack<<"  dr : "<<dradius<<" res : "<<residue<<endl;
//cout<<hit[k]->straw<<endl;
                        
                       // xperfect3[k] = ( A2[k] - trck_constant ) / trck_slope; //x coo of the fit track
                       
                       //xperfect3[k] = ( A1[k] - trck_constant ) / trck_slope;
                       xperfect3[k] = (trck_slope*A2[k])+trck_constant;
                        
                       // cout<<"Z tk: "<<A1[k]<<"\t Cnst tk : "<<trck_constant<<"\t X tk : "<<xperfect3[k]<<"\t Slp tk : "<< trck_slope <<"\t"<<(xperfect3[k]-trck_constant)/trck_slope<<endl; 
                       
                    //   cout<<"Lay :"<<hit[k]->layer<<" Strw :"<<hit[k]->straw<<" Hit X :"<<hit[k]->x<<" Hit Z :"<<hit[k]->z<<" A1 :"<<A1[k]<<" A2 :"<<A2[k]<<" Tk X :"<<xperfect3[k]<<" DT :"<<hit[k]->drifttime<<endl;
 
                       // dSlope0 = - ( 1 / trck_slope ); //gives the slope of the perp line
                        
                       // dConst0 = vec_strawZx[k] + ( vec_strawX[k] / trck_slope ); //get the constant c in y=mx+c
                        
                       // X_perpX = ( dConst0 - trck_constant ) / ( trck_slope - dSlope0 );// at the point where the perpendicular meets the track, the equation of lines becomes the same hence x = (c2 - c1)/(m1-m2)
                        
                      //  X_perpY = ( dSlope0 * X_perpX ) + dConst0; // y coordinate at the point where the perpendicular meets the track line.

                        hx->Fill ( residue );
                        hdx->Fill((xperfect3[k]-vec_strawX[k])*10);
                        
                        int plane_dx = ( myCombination.at ( chi_index ).at ( k )->layer-1 ) *2 +myCombination.at ( chi_index ).at ( k )->plane;
                        h_res_plane[plane_dx]->Fill( residue );                    
                        h_dx[plane_dx]->Fill((xperfect3[k]-vec_strawX[k])*10); //cm to mm
                        h_Dr_vs_dr[plane_dx]->Fill(dradius*10,residue);
                        h_XfvsZ->Fill(xperfect3[k]*10,vec_strawZx[k]*10);
                        h_plane_str_no[plane_dx]->Fill(hit[k]->straw);
                        
                        if(theta_X < 0.3 && theta_X > -0.3 ){
                            h_dt[0]->Fill(hit[k]->drifttime);
                            h_tot[0]->Fill(hit[k]->tot);                            
                            h_thet_plane[0]->Fill(plane_dx);
                            h_thet_res[0]->Fill(residue);
                            h_thet_straw_x[0]->Fill(vec_strawX[k]);
                            if (xperfect3[k] > vec_strawX[k]){
                                h_orientation[0]->Fill(0);
                                //cout<<"\tGOOD_Right :- "<<"DR : "<<dradius<<"\tstrawX "<< vec_strawX[k] <<"\tXe:"<< A1[k] <<"\tXf :"<<xperfect3[k]<<"\tXf-Xe :"<<(xperfect3[k]-A1[k])*10<<endl;
                            }
                            else {
                                h_orientation[0]->Fill(1);
                                //cout<<"\tGOOD_Left :- "<<"DR : "<<dradius<<"\tstrawX "<< vec_strawX[k] <<"\tXe:"<< A1[k] <<"\tXf :"<<xperfect3[k]<<"\tXf-Xe :"<<(xperfect3[k]-A1[k])*10<<endl;
                            }
                        }
                        else {
                            h_dt[1]->Fill(hit[k]->drifttime);
                            h_tot[1]->Fill(hit[k]->tot); 
                            
                            h_thet_plane[1]->Fill(plane_dx);
                            h_thet_res[1]->Fill(residue);
                            h_thet_straw_x[1]->Fill(vec_strawX[k]);
                            if (xperfect3[k] > vec_strawX[k]){
                                h_orientation[1]->Fill(0);
                                //cout<<"\tBAD_Right :- "<<"DR : "<<dradius<<"\tstrawX "<< vec_strawX[k] <<"\tXe:"<< A1[k] <<"\tXf :"<<xperfect3[k]<<"\tXf-Xe :"<<(xperfect3[k]-A1[k])*10<<endl;
                            }
                            else {
                                h_orientation[1]->Fill(1);
                               // cout<<"\tBAD_Left :- "<<"DR : "<<dradius<<"\tstrawX "<< vec_strawX[k] <<"\tXe:"<< A1[k] <<"\tXf :"<<xperfect3[k]<<"\tXf-Xe :"<<(xperfect3[k]-A1[k])*10<<endl;
                            }
                            
                        }

                    /* 
                        if (xperfect3[k] > vec_strawX[k]){
                            cout<<"\tRight :- "<<"DR : "<<dradius<<"\tstrawX "<< vec_strawX[k] <<"\tXe:"<< A1[k] <<"\tXf :"<<xperfect3[k]<<"\tXf-Xe :"<<(xperfect3[k]-A1[k])*10<<endl;
                        }
                        else {
                            cout<<"\tLeft :- "<<"DR : "<<dradius<<"\tstrawX "<< vec_strawX[k] <<"\tXe:"<< A1[k] <<"\tXf :"<<xperfect3[k]<<"\tXf-Xe :"<<(xperfect3[k]-A1[k])*10<<endl;
                        }*/

                    }
                //}
                //cout << "\n\n" << endl;
            }

            ///////////////////////////////////////////////////////////////////////////////////
            // COMBINATIONS to get least Y
            // ChiSquar2############################################

//             Int_t e_index1 = 0;
//             std::vector<Double_t> chiY_array;
//             Double_t chi_value1;
//
//             Double_t dSlope = 0;
//             Double_t dConst = 0;
//             Double_t Y_perpX = 0;
//             Double_t Y_perpY = 0;
//             Double_t Y_short = 0;
//
//             vector<vector<Double_t>> myCombination1;
//             myCombination1.clear();
//
//             int no_comb_y = pow(2, count1);
//             chiY_array.resize(no_comb_y);
//             for (Int_t co = 0; co < no_comb_y; co++)
//             {
//
//                 vector<Double_t> vt1;
//                 vt1.clear();
//
//                 for (Int_t coo = 0; coo < count1; coo++)
//                 {
//                     int bit_idx = (1 << coo);
//                     int comb_idx = (co & bit_idx) ? 1 : 0;
//
//                     B1[coo] = vec_dr_yaxis[comb_idx][coo]->y;
//                     B2[coo] = vec_dr_yaxis[comb_idx][coo]->z;
//                     vt1.push_back(vec_dr_yaxis[comb_idx][coo]->y);
//                 }
//
//                 myCombination1.push_back(vt1);
//
//                 // cout<<"}";
//                 // cout << "\n";
//
//                 // cout << "size1 B1 :" << (sizeof(B1)) / 8 << endl;
//                 // cout << "size2 B2:" << (sizeof(B2)) / 8 << endl;
//                 //                cout << "size3 yperfect2:" <<
//                 //                (sizeof(yperfect2)) / 8 << endl;
//
//                 //                 for (Int_t coooq = 0; coooq < total_pairs;
//                 //                 coooq++) {
//                 //                     cout << "\t\t" << B1[coooq] << "\t";
//                 //                 }
//                 // cout << "\n" << endl;
//
//                 TGraph* chiY = new TGraph(vt1.size(), B1, B2);
//                 chiY->Fit(f4, "q");
//                 chi_value1 = f4->GetChisquare();
//                 //	cout<<"\n\n@#$@#$@#$@#$@#$"<<chi_value<<endl;
//
//                 chiY_array[co] = chi_value1;
//
//                 delete chiY;
//             }
//
//             Double_t smallest1 = chiY_array[0];
//
//             Int_t chi_index1 = 0;
//             for (Int_t ciq = 0; ciq < no_comb_y; ciq++)
//             {
//                 // cout << "chiSquareY  :" << chiY_array[ciq] << endl;
//
//                 if (smallest1 > chiY_array[ciq]) {
//                     smallest1 = chiY_array[ciq];
//                     chi_index1 = ciq;
//                 }
//             }
//
//             Double_t Why[50];
//
//             //             cout << "\n"
//             //                  << "BEST COMBINATION 2 :";
//             for (Int_t gddq = 0; gddq < count1; gddq++)
//             {
//                 // cout << "\t " << myCombination1.at(chi_index1).at(gddq) << "
//                 // ";
//                 Why[gddq] = myCombination1.at(chi_index1).at(gddq);
//             }
//             // cout << "\n\n\n" << endl;
//
//             // cout << "\n SMALLEST ChiSquar 2 :" << smallest1 << "(" <<
//             // chi_index1 << ")" << endl;
//
//             for (Int_t coo = 0; coo < count1; coo++)
//             {
//                 int bit_idx = (1 << coo);
//                 int comb_idx = (chi_index & bit_idx) ? 1 : 0;
//
//                 B1[coo] = vec_dr_yaxis[comb_idx][coo]->y;
//                 B2[coo] = vec_dr_yaxis[comb_idx][coo]->z;
//             }
//
//             if (vec_dr_yaxis[0].size() != 0) {
//
//                 yfit = new TGraph(count1, Why, B2);
//                 yfit1 = new TGraph(vec_hits_new.size(), strawY, strawZ);
//
//                 //				yfit->GetXaxis()->SetLimits(40,100);
//
//                 // TGraph* yfit2 = new
//                 // TGraph(vt1.size(),yperfectmean,yperfect2);
//
//                 yfit->GetXaxis()->SetTitle("Y Axix [cm]");
//                 yfit->GetYaxis()->SetTitle("Z Axix [cm]");
//                 yfit->SetTitle("Y-Coordinates");
//
//                 yfit->SetMarkerStyle(7);
//                 yfit->SetMarkerSize(2);
//                 yfit->SetMarkerColor(kBlue + 2);
//                 yfit->SetLineColor(0);
//                 yfit->SetLineWidth(2);
//
//                 yfit1->SetMarkerStyle(24);
//                 yfit1->SetMarkerSize(2);
//                 yfit1->SetMarkerColor(kBlue + 2);
//                 yfit1->SetLineColor(0);
//
//                 yfit->Fit(f2, "q");
//
//                 f2 = yfit->GetFunction("f2");
//                 Double_t pp0 = f2->GetParameter(0);
//                 Double_t pp1 = f2->GetParameter(1);
//                 // emcY_arry[i] = (70 - pp0) / pp1;
//                 vec_emcY.push_back((70 - pp0) / pp1);
//
//                 // cout << " EMC   :   Y "<< ((70 - pp0) / pp1)<<endl;
//
//                 // yfit->SetPoint(1, 10,10);
//
//                 for (Int_t kk = 0; kk < count1; kk++)
//                 {
//
//                     yperfect3[kk] = (B2[kk] - pp0) / pp1;
//                     // cout << "SPACE Y  :" << fabs(B1[kk] - yperfect3[kk]) <<
//                     // endl;
//                     // hy->Fill(fabs(B1[kk] - yperfect3[kk]));
//
//                     // cout << "Double Check Y :"<< vec_strawY[kk]<<endl;
//
//                     dSlope = -(1 / pp1);
//                     dConst = vec_strawZy[kk] + (vec_strawY[kk] / pp1);
//                     Y_perpX = (dConst - pp0) / (pp1 + (1 / pp1));
//                     Y_perpY = (dSlope * Y_perpX) + dConst;
//                     Y_short = (fabs(sqrt(((vec_strawY[kk] - Y_perpX) *
//                                           (vec_strawY[kk] - Y_perpX)) +
//                                          ((vec_strawZy[kk] - Y_perpY) *
//                                           (vec_strawZy[kk] - Y_perpY))))) -
//                               (fabs(vec_strawY[kk] - B1[kk]));
//
//                     hy->Fill(Y_short);
//
//                     // cout << "Y Real Check  :  " << vec_strawY[kk] << "\t" <<
//                     // vec_strawZy[kk] <<"\t" << Y_perpX <<"\t" <<Y_perpY<<"\t
//                     // PERP : "<<Y_short << endl;
//                 }
//             }
            // }
            // }/*
            //  cout << "\nEMC XYZ   :     "<< emcX << " - "<< emcY<<endl;
            //
            //             //mg->Add(yfit);
            //             //mg->Add(yfit2);
            //             //xfit1->Draw("AP");
            //             //mg->Add(yfit);
            //
            //             mg->Add(xfit);
            // 		mg->Add(xfit1);
            // 		mg->Write("X");
            // 		mg1->Add(yfit);
            // 		mg1->Add(yfit1);
            // 		mg1->Write("Y");
            //
            //             xfit->Draw("same,P");
            //             xfit->Write("X");
            //             yfit->Draw("same,P");
            //             yfit->Write("Y");
            //
            //       //  emcX_arry[i] = emcX;
            //        // emcY_arry[i] = emcY;
            //
            //
        }

        for ( int s = 0; s < vec_hits.size(); s++ ) {
            delete ( vec_hits[s] );
        }
        //
        //       // cout << "\n\nEMC COORDINATES : "<< emcX << "  -  "<< emcY<<
        //       endl;
        //
        //
        //
//         c1->Write();
//         c1->Close();
    }

    Double_t emcX_arry[vec_emcX.size()];
    //Double_t emcY_arry[vec_emcX.size()];
    cout << "Total hits   :" << vec_driftradius.size() << endl;
    // cout<<"vec_emcX size :"<<vec_emcX.size()<<endl;
//       for (Int_t ev = 0; ev < vec_emcX.size(); ev++)
//       {
//           emcX_arry[ev] = vec_emcX[ev];
//           emcY_arry[ev] = vec_emcY[ev];
//       }
    // cout << "\n\n SIZES   :   " << vec_emcX.size()<< "  -  "<<
    // vec_emcX.size()<<endl;
    // cout << "\n\n\n SIZES   :   " << emcX_arry.length<< "  -  "<<
    // emcY_arry.length<<endl;
    //      for ( Int_t abc = 0; abc < vec_emcX.size(); abc++)
    //      {
    //          cout << "\n\n\nEMC COORDINATES : "<< emcX_arry[abc] << "  -  "<<
    //      emcY_arry[abc]<< endl;
    //      }
//       TGraph* emc = new TGraph(vec_emcX.size(), emcX_arry, emcY_arry);
//       emc->SetMarkerStyle(3);
//       emc->SetMarkerSize(1);
//       emc->SetMarkerColor(2);
//       emc->SetLineColor(0);
    // emc->SetLineWidth(2);
    
    TF1* fh = new TF1("fh", "gaus",  -2.5, 2.5);
    hx->Fit(fh,"R");

    int a_size = sizeof(outtree)/sizeof(char);
    TString s_a = convertToString(outtree,a_size)+".png";
    gStyle->SetOptFit(1101);
    hx->GetXaxis()->SetTitle("Residuals dr[mm]");
    hx->Draw();   
    c1->SaveAs(s_a,"png");
    hx->Write();
    
    hdx->Fit(fh,"R");
    hdx->GetXaxis()->SetTitle("dx [mm]");
    hdx->Draw();
    hdx->GetXaxis()->SetTitle("Distance to the track from the straw center dx [mm]");
    hdx->Write();
   
    double plane_res_array[17];
    double plane_indx_array[17];
    for (int dx=0; dx<17; dx++){
        h_dx[dx]->GetXaxis()->SetTitle("Distance to the track from the straw center dx [mm]");
        //h_dx[dx]->GetYaxis()->SetTitle("Shortest distance to the track from the straw center [mm]");
        h_dx[dx]->Write(); 
        
        h_res_plane[dx]->GetXaxis()->SetTitle("Residual dr [mm]");
        h_res_plane[dx]->Fit(fh,"RQ");
        gStyle->SetOptFit(1101);
        h_res_plane[dx]->Draw();
        h_res_plane[dx]->Write();
        plane_indx_array[dx]=dx;
        plane_res_array[dx]= h_res_plane[dx]->GetXaxis()->GetBinCenter(h_res_plane[dx]->GetMaximumBin());
        
        h_Dr_vs_dr[dx]->Write();
        h_plane_str_no[dx]->Write();        
    }
    
                       /*     hx->Fill ( residue );  //dr
                        hdx->Fill((xperfect3[k]-A1[k])*10); //dx
                        h_res_plane[plane_dx]->Fill( residue );          //dr          
                        h_dx[plane_dx]->Fill((xperfect3[k]-A1[k])*10); //cm to mm //dx
                        h_Dr_vs_dr[plane_dx]->Fill(dradius*10,residue); //dr  */
    
    
    TGraph* gPlaneRes = new TGraph ( 17, plane_indx_array, plane_res_array );
    gPlaneRes->SetName ( "Plane_Residuals" );
    gPlaneRes->Draw("AP");
    gPlaneRes->GetXaxis()->SetTitle("Plane no");
    gPlaneRes->GetXaxis()->SetTitle("Maximum of the Residuals [mm]");

    gPlaneRes->GetXaxis()->SetRangeUser(-2,18);
    gPlaneRes->GetYaxis()->SetRangeUser(-1,1);
    gPlaneRes->SetMarkerColor(kOrange+1);
    gPlaneRes->SetMarkerStyle(15);
    gPlaneRes->SetMarkerSize(4);
    
    gPlaneRes->Write();
    
    hx_slope->Write();
    hx_const->Write();

   
    //h->Write();
    //hy->Write();
    h_theeta_X->GetXaxis()->SetTitle("Theta [deg]");
    h_theeta_X->Write();

    //c1->Write();
    //c1->Close();
    // cout <<"THETA SIZE : "<<vec_Xsr.size()<<"  "<<vec_Xtheta.size()<<endl;

    Double_t arr_Xt[vec_Xtheta.size()];
    Double_t arr_Xsr[vec_Xtheta.size()];

    for ( int tt = 0; tt < vec_Xtheta.size(); tt++ ) {
        arr_Xt[tt] = vec_Xtheta[tt];
        arr_Xsr[tt] = vec_Xsr[tt];
    }
    TGraph* theta_SR = new TGraph ( vec_Xtheta.size(), arr_Xt, arr_Xsr );

    theta_SR->Write ( "P" );
//       emc->Draw("P");
//       emc->Write();
    Z_value->Write();
    
    for (int i =0; i<2; i++){
        h_dt[i]->Write();
        h_tot[i]->Write();
        h_thet_plane[i]->Write();
        h_thet_res[i]->Write();
        h_orientation[i]->Write();
        h_thet_straw_x[i]->Write();

    }
    
    TH1F * projh2X;

    Int_t  z_bins_arr[4]={0,2,4,6};
    const EColor colours[] = {kBlue,kRed,kGreen,kBlack};
    Ct->cd();
    
    TLegend* leg = new TLegend ( 0.5,0.7,0.9,0.9 );
    leg->SetHeader ( "ProjX" );
    leg->SetFillColor ( 1 );
    
    for (int b=0; b<4; b++)
    { 
	    projh2X = (TH1F*) h_XfvsZ->ProjectionX(Form("Layer_%i",b),z_bins_arr[b],z_bins_arr[b]+1);
	    projh2X->Draw("same");
        projh2X->SetLineWidth(2);
	    projh2X->SetLineColor(colours[b]);
	    leg->AddEntry (projh2X,Form("Layer_%d",b+1),"lep" );
       
    }
    leg->SetFillStyle ( 0 );
    leg->Draw();
    Ct->Write();
    
    xf_z->cd();
    h_XfvsZ->Draw("colz");

    double ln_x =17.4;
    for (int l=0; l<32; l++)
    {
        TLine *line = new TLine(ln_x*10,0,ln_x*10,600);
        line->Draw("same");
        line->SetLineColor(kRed);
        ln_x = ln_x+0.505;
    }
    xf_z->Write();
    
    
    driftfile1->Close();

    // emc->Write();
    // DR->Draw();
    // cout<<"count in Layer1+layer3 : "<< count<<"\t"<<"count in Layer2+layer4
    // : "<<"\t"<<count1<<"\t"<<"TOTAL : "<<count+count1<<endl;
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
            PDAQ_Spl_Res ( argv[1], argv[2], 100000000 );
        } 
        else 
        {
            PDAQ_Spl_Res ( argv[1], argv[2], atoi ( argv[3] ) );
        }

    else 
    {
        return 1;
    }

    return 0;
}



