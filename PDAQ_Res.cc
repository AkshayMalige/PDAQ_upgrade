#include "PDAQ_Res.h"
#include "string.h"
#include "TStyle.h"
#include <bits/stdc++.h>

using namespace std;

string convertToString (char *a , int size)
{
    string s(a);
    
    return s;
}

std::vector<SttHit*> filter_drift_position (std::vector<SttHit*> vec_hits , TGraph* gDR)
{
    std::vector<SttHit*> vec_dphits;
    int layer_no[4] = {1,4,5,8};
    std::vector<SttHit*>vec_lay;

    for(int j=0; j<4; j++)
    {
        vec_lay.clear();
        for(int i=0; i<vec_hits.size(); i++)
        {
            if (vec_hits[i]->layer == layer_no[j])
            {
                vec_lay.push_back(vec_hits[i]);
               // cout<<j <<"\t"<<layer_no[j]<<"\t"<<vec_hits[i]->plane<<endl;
            }
        }
       // cout<<"size "<<vec_lay.size()<<endl;
        if (vec_lay[0]->x < vec_lay[1]->x)
        {
            vec_lay[0]->driftposition = vec_lay[0]->x + gDR->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->driftposition = vec_lay[1]->x - gDR->Eval ( vec_lay[1]->drifttime );
            vec_lay[0]->DriftRadius = gDR->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->DriftRadius = gDR->Eval ( vec_lay[1]->drifttime );
            
        }
        
        else 
        {
            vec_lay[0]->driftposition = vec_lay[0]->x - gDR->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->driftposition = vec_lay[1]->x + gDR->Eval ( vec_lay[1]->drifttime );
            vec_lay[0]->DriftRadius = gDR->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->DriftRadius = gDR->Eval ( vec_lay[1]->drifttime );
        }
        vec_dphits.push_back(vec_lay[0]);
        vec_dphits.push_back(vec_lay[1]);
    }
    return vec_dphits;    
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
    
  // TFile* gcrr_file = new TFile ( "er.root", "READ" ); 
  // TGraph* gDR = ( TGraph* ) gcrr_file->Get ( "Dt_Correction1" );
    TGraph* gDR = ( TGraph* ) inFile.Get ( "PDAQ_DR" );
    gDR->SetName("gDR");
    
    driftfile1->cd();
    
    TH1F* hx = new TH1F ( "hx", "Spatial_Resolution_X", 500, -5, 5 );
    TH1F* hdx = new TH1F ( "hdx", "res_non_corrected", 500, -5, 5 );
    TH1F* h_dr = new TH1F ( "h_dr", "drift_radius", 350, 0, 350 );
    TH1F* hx_theeta = new TH1F ( "h_theeta_X", "h_theeta_X", 400, -2, 2 );    
    TH1F* hx_slope = new TH1F ( "hx_slope", "Slope_X", 1000, -0.25, 0.25 );
    TH1F* hx_const = new TH1F ( "hx_const", "Constant_X", 1000, 0, 50 );
    TH2F* h_XfvsZ = new TH2F ( "h_XfvsZ", "h_XfvsZ;X of track fit X_{f} [mm]; Z [mm]", 200, 150, 350,7,-10,600 );
   // TH2F* h_Dr_vs_dr = new TH2F ("h_Dr_vs_dr", "h_Dr_vs_dr;Drift radius [mm];dr [mm]",60,0,6,100,-1,1 );
    TH2F* Dt_vs_dr = new TH2F ("Dt_vs_dr", "Dt_vs_dr;Drift time [ns];dr [mm]",200,0,200,100,-1,1 );
    
    
    TH1F* h_dx[18];
    TH1F* h_res_plane[18];
    TH2F* h_Dr_vs_dr[18];
    TH1F* h_plane_str_no[18];
    TH2F* h_dt_vs_dr[18];
    //int plane_no[8] = {0,1,6,7,8,9,14,15};
    
    for ( int mc = 0; mc < 18; mc++ ) 
    {
        h_res_plane[mc] = new TH1F ( Form("h_res_plane%i",mc+1), Form("h_res_plane%i;dr [mm]",mc+1), 500, -5, 5 );
        h_dx[mc] = new TH1F ( Form("h_dx_Plane_%i",mc+1), Form("h_dx_Plane_%i;dx [mm]",mc+1),500,-5,5 );
        h_Dr_vs_dr[mc] = new TH2F ( Form("h_Dr_vs_dr_%i",mc+1), Form("h_Dr_vs_dr_%i;Drift radius [mm];dr [mm]",mc+1),60, 0,6,100,-1,1 );
        
        h_plane_str_no[mc] = new TH1F (Form("h_plane_str_no%i",mc+1),Form("h_plane_str_no%i;Straw No.",mc+1),60, 0, 60);
        
        h_dt_vs_dr[mc] = new TH2F ( Form("h_dt_vs_dr%i",mc+1), Form("h_dt_vs_dr%i;Drift time [ns];dr [mm]",mc+1),200, 0,200,100,-1,1 );

    }
    
    TCanvas * Ct; //Canvas for ToT
    Ct=new TCanvas ( "Ct","Ct" );
    
    TCanvas * xf_z; //Canvas for ToT
    xf_z=new TCanvas ( "xf_z","xf_z" );
    
    TCanvas* c1 = new TCanvas ( "c1","c1" );

    
    int iev = ( int ) tree->GetEntries();
    printf("\n Number of entries in the tree : %d \n",iev);
    
    vector<SttHit*> vec_hits;
    std::vector<std::vector<SttHit*>> vec_layerhits;
    
    for ( int i = 0; i < iev; i++ ) 
    {
        int count = 0;
        
        if ( i == maxEvents ) 
        {
            break;
        }
        if ( i % 1000 == 0 ) 
        {
            cout << i << endl;
        }
        
        tree->GetEntry ( i );
        printf("\nTrack No. %i \n",i);
        
        int oiv = vec_Drifttime->size();
        
        vec_hits.clear();
        vec_layerhits.clear();

        
        for ( int n = 0; n < oiv; n++ ) 
        {
            
            SttHit* a = new SttHit();
            a->drifttime = vec_Drifttime->at ( n );
            a->x = vec_x->at ( n );
            a->y = vec_y->at ( n );
            a->z = vec_z->at ( n );
            a->layer = vec_layer->at ( n );
            a->tot = vec_tot->at(n);
            a->plane = (( vec_layer->at ( n ) -1 ) *2 ) + vec_plane->at ( n );
            a->straw = vec_straw->at ( n );

            vec_hits.push_back ( a );
        }
        
        if (vec_hits.size()<1) continue;
        
        int plane_mult[17];
        
        for (int a =0; a<17; a++){
            plane_mult[a]=0;
        }
        
        for (int u=0; u<vec_hits.size(); u++)
        {

            plane_mult[vec_hits[u]->plane]++;             
            
        }
        
        int plane_mult_counter =0;
        
        for (int aa=0; aa<17; aa++)
        {
            if (plane_mult[aa]>0)
            {
                plane_mult_counter++;
            }
        }
        
        bool clean_track = (plane_mult_counter == 8) ? true : false;
                
        if (clean_track == false) continue;
        
        std::vector<SttHit*> vec_dphits;
        vec_dphits.clear();
        
        vec_dphits = filter_drift_position(vec_hits, gDR);
        //filter_drift_position(vec_hits, gDR);
        double hit_x[vec_dphits.size()];
        double hit_z[vec_dphits.size()];
        int index_p =0;
        
        for (int k=0; k<vec_dphits.size(); k++)
        {
            if (vec_dphits[k]->layer==5)
            {
                vec_dphits[k]->driftposition = vec_dphits[k]->driftposition;
            }
        }
        
        for (int a=0; a<vec_dphits.size(); a++)
        {
           // cout<<vec_dphits[a]->layer<<"\t"<<vec_dphits[a]->plane<<"\t"<<vec_dphits[a]->x<<"\t"<<vec_dphits[a]->DriftRadius<<"\t"<<vec_dphits[a]->driftposition<<endl;
        //    if (vec_dphits[a]->layer !=8)
        //    {
                hit_x[index_p]=vec_dphits[a]->driftposition;
                hit_z[index_p]=vec_dphits[a]->z;
                //cout<<a<<"\t"<<vec_dphits[a]->plane<<"\t"<<vec_dphits.size()<<endl;
                index_p++;
         //   }

        }
    
        TGraph* track = new TGraph ( vec_dphits.size(),hit_z, hit_x );
        
        track->Fit ( f1, "q" );

        f1 = track->GetFunction ( "f1" );
        double trck_constant = f1->GetParameter ( 0 ); //constant - c
        double trck_slope = f1->GetParameter ( 1 ); //slope    - m
        double theta_X = ((atan ( trck_slope)) * 180) / 3.1415926 ;
        hx_slope->Fill(trck_slope);
        hx_const->Fill(trck_constant);
        hx_theeta->Fill(theta_X);
        
        double track_x[vec_dphits.size()];
        
        for ( int k = 0; k < vec_dphits.size(); k++ ) 
        {                        
            double centerTotrack = fabs((trck_slope * vec_dphits[k]->z) - vec_dphits[k]->x + trck_constant) / ( sqrt ( 1 + ( trck_slope * trck_slope ) ) ) ; //((m * Px) - Py + C)/(sqrt(1+(m*m))) distance b/w a point and line. reference in the downloaded pdf.
                        
            double dradius = fabs ( vec_dphits[k]->x - vec_dphits[k]->driftposition );                        
            double residue = (centerTotrack-dradius) * 10; // Cm to mm                        
            track_x[k] = (trck_slope * vec_dphits[k]->z) + trck_constant;
          //  hdx->Fill((track_x[k] - vec_dphits[k]->x)*10);
            hx->Fill ( residue );
            
           // if (vec_dphits[k]->driftposition > vec_dphits[k]->x)
           // {
                hdx->Fill((track_x[k] - vec_dphits[k]->driftposition)*10);
              //  cout<<"Right : straw x:"<<vec_dphits[k]->x<<"\tXf :"<<track_x[k]<<"\tXe :"<<vec_dphits[k]->driftposition<<"\tdx :"<<(track_x[k] - vec_dphits[k]->driftposition)*10<<endl;
                
           // }
           /* else
            {
                hdx->Fill((vec_dphits[k]->driftposition - track_x[k] )*10); 
              //  cout<<"Left : straw x:"<<vec_dphits[k]->x<<"\tXf :"<<track_x[k]<<"\tXe :"<<vec_dphits[k]->driftposition<<"\tdx :"<<(track_x[k] - vec_dphits[k]->driftposition)*10<<endl;
            }*/
            h_dr->Fill(vec_dphits[k]->DriftRadius * 10);
            h_XfvsZ->Fill(track_x[k]*10,vec_dphits[k]->z*10); 
            Dt_vs_dr->Fill(vec_dphits[k]->drifttime,residue);
            
            h_res_plane[vec_dphits[k]->plane]->Fill( residue );  
         //   if (vec_dphits[k]->driftposition > vec_dphits[k]->x)
           // {
                h_dx[vec_dphits[k]->plane]->Fill((track_x[k] - vec_dphits[k]->driftposition)*10);
           /* }
            else
            {
                h_dx[vec_dphits[k]->plane]->Fill((vec_dphits[k]->driftposition - track_x[k] )*10);  
            }*/
         //   h_dx[vec_dphits[k]->plane]->Fill((track_x[k] - vec_dphits[k]->driftposition)*10); 
            h_Dr_vs_dr[vec_dphits[k]->plane]->Fill(dradius*10,residue);
            h_plane_str_no[vec_dphits[k]->plane]->Fill(vec_dphits[k]->straw);
            //cout<<vec_dphits[k]->plane<<"\t"<<vec_dphits[k]->drifttime<<"\t"<<residue<<endl;
            h_dt_vs_dr[vec_dphits[k]->plane]->Fill(vec_dphits[k]->drifttime,residue);
        }
                        
    }
    
    TH1F * projhDT;
    TF1* ft = new TF1("fh", "gaus");
    double dt_time[200];
    double dt_correted[200];
    double crr_mean =0;
    cout<<"DT200 "<<gDR->Eval(200)<<endl;
    for(int t=0; t<200; t++)
    {
        projhDT = (TH1F*) Dt_vs_dr->ProjectionY(Form("timeBin%i",t),t,t+1);
        projhDT->Fit(ft, "q");
        dt_time[t]= t;
        double crr_mean = (t<30 || t>170) ? 0 : ft->GetParameter(1);
        dt_correted[t]= gDR->Eval(t) + (crr_mean/10);   
        cout<<t<<"\t"<<crr_mean<<"\t"<<dt_correted[t]<<endl;
    }
    delete ft;
    
    TGraph* gDriftRadius = new TGraph ( 200, dt_time, dt_correted );
    gDriftRadius->SetName ( "Dt_Correction" );
    gDriftRadius->Write();
    
    TF1* fh = new TF1("fh", "gaus",  -2.5, 2.5);
    hx->Fit(fh,"R");
    int a_size = sizeof(outtree)/sizeof(char);
    TString s_a = convertToString(outtree,a_size)+".png";
    gStyle->SetOptFit(1101);
    hx->GetXaxis()->SetTitle("Residuals dr[mm]");
    hx->Draw();   
    //c1->SaveAs(s_a,"png");
    hx->Write();
    
    hdx->Fit(fh,"R");
    hdx->GetXaxis()->SetTitle("dx [mm]");
    hdx->Draw();
    hdx->GetXaxis()->SetTitle("X_{f} - X_{e} [mm]");
    hdx->Write();
    hx_slope->Write();
    hx_const->Write();
    hx_theeta->GetXaxis()->SetTitle("#Theta [deg]");
    hx_theeta->Write();
    h_dr->GetXaxis()->SetTitle("Drift radius [mm]");
    h_dr->Write();
    Dt_vs_dr->Write();
    gDR->Write();
    
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
    
    double plane_res_array[18];
    double plane_indx_array[18];
    for (int dx=0; dx<18; dx++){
        h_dx[dx]->GetXaxis()->SetTitle("X_{f} - X_{e} [mm]");
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
        h_dt_vs_dr[dx]->Write();
    }
    
    TGraph* gPlaneRes = new TGraph ( 8, plane_indx_array, plane_res_array );
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

    
    driftfile1->Close();
    
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
