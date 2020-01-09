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

std::vector<SttHit*> filter_drift_position (std::vector<SttHit*> vec_hits , TGraph* gDR[8])
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
            vec_lay[0]->driftposition = vec_lay[0]->x + gDR[vec_lay[0]->layer-1]->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->driftposition = vec_lay[1]->x - gDR[vec_lay[1]->layer-1]->Eval ( vec_lay[1]->drifttime );
            vec_lay[0]->DriftRadius = gDR[vec_lay[0]->layer-1]->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->DriftRadius = gDR[vec_lay[1]->layer-1]->Eval ( vec_lay[1]->drifttime );
            
        }
        
        else 
        {
            vec_lay[0]->driftposition = vec_lay[0]->x - gDR[vec_lay[0]->layer-1]->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->driftposition = vec_lay[1]->x + gDR[vec_lay[1]->layer-1]->Eval ( vec_lay[1]->drifttime );
            vec_lay[0]->DriftRadius = gDR[vec_lay[0]->layer-1]->Eval ( vec_lay[0]->drifttime );
            vec_lay[1]->DriftRadius = gDR[vec_lay[1]->layer-1]->Eval ( vec_lay[1]->drifttime );
        }
        vec_dphits.push_back(vec_lay[0]);
        vec_dphits.push_back(vec_lay[1]);
    }
    return vec_dphits;    
}

bool track_recon (int iev , histograms* h, TGraph* gDR[8] , int maxEvents,vec * v,TTree* tree, bool write_file)
{
    ofstream myfile;
    myfile.open ("a.txt");
    std::vector<SttHit*> vec_hits;
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
            cout <<"Track No "<< i << endl;
        }
        
        tree->GetEntry ( i );
       // printf("\nTrack No. %i \n",i);
        
        int oiv = v->vec_Drifttime->size();
        
        vec_hits.clear();
        vec_layerhits.clear();

        
        for ( int n = 0; n < oiv; n++ ) 
        {
            
            SttHit* a = new SttHit();
            a->drifttime = v->vec_Drifttime->at ( n );
            a->x = v->vec_x->at ( n );
            a->y = v->vec_y->at ( n );
            a->z = v->vec_z->at ( n );
            a->layer = v->vec_layer->at ( n );
            a->tot = v->vec_tot->at(n);
            a->plane = (( v->vec_layer->at ( n ) -1 ) *2 ) + v->vec_plane->at ( n );
            a->straw = v->vec_straw->at ( n );

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
        
        double hit_x[vec_dphits.size()];
        double hit_z[vec_dphits.size()];
        int index_p =0;
        TF1* f1 = new TF1 ( "f1", "pol1" );
        
        for (int a=0; a<vec_dphits.size(); a++)
        {
            hit_x[index_p]=vec_dphits[a]->driftposition;
            hit_z[index_p]=vec_dphits[a]->z;
            index_p++;
        }
    
        TGraph* track = new TGraph ( vec_dphits.size(),hit_z, hit_x );
        
        track->Fit ( f1, "q" );

        f1 = track->GetFunction ( "f1" );
        double trck_constant = f1->GetParameter ( 0 ); //constant - c
        double trck_slope = f1->GetParameter ( 1 ); //slope    - m
        double theta_X = ((atan ( trck_slope)) * 180) / 3.1415926 ;
        h->hx_slope->Fill(trck_slope);
        h->hx_const->Fill(trck_constant);
        h->hx_theeta->Fill(theta_X);
        
        double track_x[vec_dphits.size()];
       // if(f1->GetChisquare()>0.02) continue;
        h->hx_chi->Fill(f1->GetChisquare());
        double calc_chi[16];
        
        for(int a=0; a<16; a++)
        {
            calc_chi[a]=0;
        }
        
        for ( int k = 0; k < vec_dphits.size(); k++ ) 
        {                        
            double centerTotrack = fabs((trck_slope * vec_dphits[k]->z) - vec_dphits[k]->x + trck_constant) / ( sqrt ( 1 + ( trck_slope * trck_slope ) ) ) ; //((m * Px) - Py + C)/(sqrt(1+(m*m))) distance b/w a point and line. reference in the downloaded pdf.
            double dradius = fabs ( vec_dphits[k]->x - vec_dphits[k]->driftposition );                        
            double residue = (centerTotrack-dradius) * 10; // Cm to mm                        
            track_x[k] = (trck_slope * vec_dphits[k]->z) + trck_constant;
            h->hx->Fill ( residue );
            h->hdx->Fill((track_x[k] - vec_dphits[k]->driftposition)*10);  
            h->h_dr->Fill(vec_dphits[k]->DriftRadius * 10);
            h->h_XfvsZ->Fill(track_x[k]*10,vec_dphits[k]->z*10); 
            h->Dt_vs_dr->Fill(vec_dphits[k]->drifttime,residue);            
            h->h_res_plane[vec_dphits[k]->plane]->Fill( residue );  
            h->h_dx[vec_dphits[k]->plane]->Fill((track_x[k] - vec_dphits[k]->driftposition)*10);
            h->h_Dr_vs_dr[vec_dphits[k]->plane]->Fill(dradius*10,residue);
            h->h_plane_str_no[vec_dphits[k]->plane]->Fill(vec_dphits[k]->straw);
            h->h_dt_vs_dr[vec_dphits[k]->plane]->Fill(vec_dphits[k]->drifttime,residue);
            if(residue < 0.6 && residue > -0.6){ h->h_Dt_vs_dr_Lay[vec_dphits[k]->layer-1]->Fill(vec_dphits[k]->drifttime,residue);          
            }
            
            double ini_sig = 0.015;
            for(int a=0; a<16; a++)
            {
                calc_chi[a]+=pow((centerTotrack-dradius),2)/(pow((ini_sig),2));
            //  calc_chi[a]+=pow((centerTotrack-dradius),2)/pow(0.022,2);
            //  cout<<"\t"<<ini_sig<<"\t"<<ini_sig<<endl;
                ini_sig=ini_sig+0.001;
            
            }
            ini_sig = 0.015;

           
        }
        for (int a=0; a<16; a++)
        {
            myfile<<calc_chi[a]<<"\t";
        }
        myfile<<endl;

        
        delete f1;
        delete track;
        h->h_cal_chi->Fill(calc_chi[6]);
    }
    myfile.close();    
    return true;
}

bool reset_hist(histograms* h,TH2F* h_Dt_vs_dr_Lay_copy[8])
{
    h->hx->Reset("ICESM");
    h->hdx->Reset("ICESM");
    h->h_dr->Reset("ICESM");
    h->hx_theeta->Reset("ICESM");
    h->hx_slope->Reset("ICESM");
    h->hx_const->Reset("ICESM");
    h->h_XfvsZ->Reset("ICESM");
    h->Dt_vs_dr->Reset("ICESM");
    h->hx_chi->Reset("ICESM");
    h->h_cal_chi->Reset("ICESM");
    
    for (int hd =0; hd<18; hd++)
    {
        h->h_res_plane[hd]->Reset("ICESM");
        h->h_Dr_vs_dr[hd]->Reset("ICESM");
        h->h_plane_str_no[hd]->Reset("ICESM");
        h->h_dt_vs_dr[hd]->Reset("ICESM");
        h->h_dx[hd]->Reset("ICESM");
    }
    
    for(int a=0; a<8; a++)
    {
        h_Dt_vs_dr_Lay_copy[a]= h->h_Dt_vs_dr_Lay[a];
        h_Dt_vs_dr_Lay_copy[a]->Write();
        h->h_Dt_vs_dr_Lay[a]->Reset("ICESM");
    }
    
    return true;
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
    histograms* h = new histograms();
    vec * v = new vec();         
    tree->SetBranchAddress ( "vec_Drifttime", &v->vec_Drifttime );
    tree->SetBranchAddress ( "vec_x", &v->vec_x );
    tree->SetBranchAddress ( "vec_y", &v->vec_y );
    tree->SetBranchAddress ( "vec_z", &v->vec_z );
    tree->SetBranchAddress ( "vec_layer", &v->vec_layer );
    tree->SetBranchAddress ( "vec_straw", &v->vec_straw );
    tree->SetBranchAddress ("vec_plane",&v->vec_plane);
    tree->SetBranchAddress ("vec_tot",&v->vec_tot);

//     TGraph* gDR = ( TGraph* ) inFile.Get ( "PDAQ_DR" );
//     gDR->SetName("gDR");
    TGraph* gDR[8]; 
    for(int a=0; a<8; a++)
    {
        gDR[a] = ( TGraph* ) inFile.Get ( Form("PDAQ_DR%i",a+1) );
        gDR[a]->SetName(Form("gDR%i",a+1));
        
        h->h_Dt_vs_dr_Lay[a] = new TH2F(Form("h_Dt_vs_dr_Lay%i",a+1), Form("h_Dt_vs_dr_Lay%i;Drift time [ns];dr [mm]",a+1),200, 0,200,100,-1,1);
    }

    
    driftfile1->cd(); 
    
    for ( int mc = 0; mc < 18; mc++ ) 
    {
        h->h_res_plane[mc] = new TH1F ( Form("h_res_plane%i",mc+1), Form("h_res_plane%i;dr [mm]",mc+1), 500, -5, 5 );
        h->h_dx[mc] = new TH1F ( Form("h_dx_Plane_%i",mc+1), Form("h_dx_Plane_%i;dx [mm]",mc+1),500,-5,5 );
        h->h_Dr_vs_dr[mc] = new TH2F ( Form("h_Dr_vs_dr_%i",mc+1), Form("h_Dr_vs_dr_%i;Drift radius [mm];dr [mm]",mc+1),60, 0,6,100,-1,1 );
        
        h->h_plane_str_no[mc] = new TH1F (Form("h_plane_str_no%i",mc+1),Form("h_plane_str_no%i;Straw No.",mc+1),60, 0, 60);
        
        h->h_dt_vs_dr[mc] = new TH2F ( Form("h_dt_vs_dr%i",mc+1), Form("h_dt_vs_dr%i;Drift time [ns];dr [mm]",mc+1),200, 0,200,100,-1,1 );

    }
    
    
    TCanvas * Ct; //Canvas for ToT
    Ct=new TCanvas ( "Ct","Ct" );
    
    TCanvas * xf_z; //Canvas for ToT
    xf_z=new TCanvas ( "xf_z","xf_z" );
    
    TCanvas* c1 = new TCanvas ( "c1","c1" );

    
    int iev = ( int ) tree->GetEntries();
    printf("\n Number of entries in the tree : %d \n",iev);

    int no_of_it =7;
    
    TGraph* Dt_vs_drCrr[no_of_it];
     TH2F* h_Dt_vs_dr_Lay_copy[8];
//     
//     for(int a=0; a<no_of_it; a++)
//     {
        for(int b=0; b<8; b++)
        {
            h_Dt_vs_dr_Lay_copy[b] = new TH2F ( Form("h_Dt_vs_dr_Lay_copy%i",b), Form("h_Dt_vs_dr_Lay_copy%i;Drift time [ns];dr [mm]",b),200, 0,200,100,-1,1 );
        }
//     }

    for (int it=0; it<no_of_it; it++)
    {
        bool write_file = (it == no_of_it -1) ? true : false;
        track_recon (iev ,  h,  gDR , maxEvents, v, tree, write_file );    
        
        for (int s=0; s<8; s++)
        {
            TH1F * projhDT;
            TF1* ft = new TF1("fh", "gaus");
            double dt_time[200];
            double dt_correted[200];
            double correction[200];
            double crr_mean =0;
            for(int t=0; t<200; t++)
            {
                projhDT = (TH1F*) h->h_Dt_vs_dr_Lay[s]->ProjectionY(Form("timeBin%i",t),t,t+1);
              //  projhDT = (TH1F*) h->Dt_vs_dr->ProjectionY(Form("timeBin%i",t),t,t+1);
                projhDT->Fit(ft, "q");
                dt_time[t]= t;
              //  double crr_mean = (t<30 || t>160) ? 0 : ft->GetParameter(1);
                double crr_mean = (t<30 || t>160) ? 0 : projhDT->GetMean();
              //  cout<<ft->GetParameter(1)<<"\t"<<projhDT->GetMean()<<endl;

                dt_correted[t]= gDR[s]->Eval(t) + (crr_mean/10);
                correction[t] = crr_mean;
                
            }
            delete ft;                
                    
            gDR[s] = new TGraph ( 200, dt_time, dt_correted );
            gDR[s]->SetName(Form("gDR%i",s+1));
            
            Dt_vs_drCrr[it] = new TGraph ( 200, dt_time, correction );
            Dt_vs_drCrr[it]->SetName(Form("Dt_vs_drCrr%i",it));
            
        }
        if (it != no_of_it-1) reset_hist(h,h_Dt_vs_dr_Lay_copy);

    }
    
    
    for(int a=0; a< no_of_it; a++)
    {
      Dt_vs_drCrr[a]->Write();  
    }
    
    for(int a=0; a<8; a++)
    {
        gDR[a]->Write();   
        h->h_Dt_vs_dr_Lay[a]->Write();
    }
    
    
    TF1* chi2_6dofF = new TF1("chi2_6dofF","ROOT::Math::chisquared_pdf(x,6,0)",0,50);
    
    h->h_cal_chi->Scale(1/h->h_cal_chi->Integral(0,50));  
    
    double chiResult = h->h_cal_chi->Chisquare(chi2_6dofF);
    cout<<"Chi fit "<<chiResult<<endl;
    
    //h->h_cal_chi->Fit(chi2_6dofF,"R");
    
    h->h_cal_chi->Write();    
    chi2_6dofF->Write();
    
    TF1* fh = new TF1("fh", "gaus",  -2.5, 2.5);
    h->hx->Fit(fh,"R");
    int a_size = sizeof(outtree)/sizeof(char);
    TString s_a = convertToString(outtree,a_size)+".png";
    gStyle->SetOptFit(1101);
    h->hx->GetXaxis()->SetTitle("Residuals dr[mm]");
    h->hx->Draw();   
    //c1->SaveAs(s_a,"png");
    h->hx->Write();
    
    
    h->hdx->Fit(fh,"R");
    h->hdx->GetXaxis()->SetTitle("dx [mm]");
    h->hdx->Draw();
    h->hdx->GetXaxis()->SetTitle("X_{f} - X_{e} [mm]");
    h->hdx->Write();
    h->hx_slope->Write();
    h->hx_const->Write();
    h->hx_theeta->GetXaxis()->SetTitle("#Theta [deg]");
    h->hx_theeta->Write();
    h->h_dr->GetXaxis()->SetTitle("Drift radius [mm]");
    h->h_dr->Write();
    h->Dt_vs_dr->Write();
    h->hx_chi->Write();
    
   // 

    //gDR->Write();
    
    TH1F * projh2X;
    Int_t  z_bins_arr[4]={0,2,4,6};
    const EColor colours[] = {kBlue,kRed,kGreen,kBlack};
    Ct->cd();
    
    TLegend* leg = new TLegend ( 0.5,0.7,0.9,0.9 );
    leg->SetHeader ( "ProjX" );
    leg->SetFillColor ( 1 );
    
    for (int b=0; b<4; b++)
    { 
	    projh2X = (TH1F*) h->h_XfvsZ->ProjectionX(Form("Layer_%i",b),z_bins_arr[b],z_bins_arr[b]+1);
	    projh2X->Draw("same");
        projh2X->SetLineWidth(2);
	    projh2X->SetLineColor(colours[b]);
	    leg->AddEntry (projh2X,Form("Layer_%d",b+1),"lep" );
       
    }
    
    leg->SetFillStyle ( 0 );
    leg->Draw();
    Ct->Write();
    
    xf_z->cd();
    h->h_XfvsZ->Draw("colz");

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
        h->h_dx[dx]->GetXaxis()->SetTitle("X_{f} - X_{e} [mm]");
        h->h_dx[dx]->Write();
        
        h->h_res_plane[dx]->GetXaxis()->SetTitle("Residual dr [mm]");
        h->h_res_plane[dx]->Fit(fh,"RQ");
        gStyle->SetOptFit(1101);
        h->h_res_plane[dx]->Draw();
        h->h_res_plane[dx]->Write();
        plane_indx_array[dx]=dx;
        plane_res_array[dx]= h->h_res_plane[dx]->GetXaxis()->GetBinCenter(h->h_res_plane[dx]->GetMaximumBin());
        
        h->h_Dr_vs_dr[dx]->Write();
        h->h_plane_str_no[dx]->Write();   
        h->h_dt_vs_dr[dx]->Write();
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
