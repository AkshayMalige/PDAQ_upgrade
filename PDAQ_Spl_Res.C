#include "PDAQ_Spl_Res.h"
using namespace std;


Bool_t PDAQ_Spl_Res(void)
{

 cout<<"Opened"<<endl; 

 
    TFile * file = TFile::Open("Drift_Radius_test.root", "READ");
    TTree* tree = 0;
    file->GetObject("DR_Tree",tree);
    
//     TFile * file = TFile::Open("PDAQ_Stt_Tracks.root", "READ");
//     TTree* tree = 0;
//     file->GetObject("PDAQ_EMC_STT_cluster_analysis", tree);
    if (!tree) {
        std::cerr << "Tree doesn't exists" << std::endl;
        return 1;
    }
    tree->Print();
    
    
    
    TF1* f1 = new TF1("f1", "pol1");
    TF1* f2 = new TF1("f2", "pol1");
    TF1* f3 = new TF1("f3", "pol1");
    TF1* f4 = new TF1("f4", "pol1");

    TFile* driftfile1 = new TFile("XYZ_Coordinates.root", "RECREATE");
    //TFile* driftfile2 = new TFile("Y_Coordinates.root", "RECREATE");
    Double_t theta_X =0;

    TH1F* hx = new TH1F("hx", "dummy_resolutionX", 500, -1, 4);
    TH1F* hy = new TH1F("hy", "dummy_resolutionY", 500, -1, 4);
    TH1F* DR = new TH1F("dr", "drift_radius", 350, 0, 350);
    TH2F *h_theeta_X =new TH2F("h_theeta_X","h_theeta_X",50,0,5,20,0,2);

    TH1F* Z_value = new TH1F("Z_value", "dummy", 5, 0, 5);

    std::vector<double>* Vec_o_test =0;
    std::vector<double>* vec_x =0;
    std::vector<double>* vec_y=0;
    std::vector<double>* vec_z=0;
    std::vector<double>* vec_layer =0;
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

    Double_t A1[100];
    Double_t A2[100];

    Double_t B1[100];
    Double_t B2[100];



    Double_t xperfect3[16];
    Double_t yperfect3[16];

    Int_t array_length = 8;
    Int_t pair_length = 2;
    Int_t total_pairs = (array_length / pair_length);
    Int_t no_of_pairs = (total_pairs * total_pairs);

    tree->SetBranchAddress("Vec_o_test", &Vec_o_test);
    tree->SetBranchAddress("vec_x", &vec_x);
    tree->SetBranchAddress("vec_y", &vec_y);
    tree->SetBranchAddress("vec_z", &vec_z);
    tree->SetBranchAddress("vec_layer", &vec_layer);
//     tree->SetBranchAddress("vec_module", &vec_module);
//     tree->SetBranchAddress("vec_fee", &vec_fee);
//     tree->SetBranchAddress("vec_fee_ch", &vec_fee_ch);
//     tree->SetBranchAddress("vec_tdc_ch", &vec_tdc_ch);
// 
    TGraph* gDR = new TGraph(220, a1, b1);
    gDR = (TGraph*)file->Get("PDAQ_DR");

    //Double_t strawY[8];
    //Double_t strawZ[8];

    Double_t yperfectmean[16];

    SttHit* new_y1 = 0;
    SttHit* new_y2 = 0;

    Int_t iev = (Int_t)tree->GetEntries();
    cout << "number of entries in tree:" << iev << endl
         << endl;

    vector<SttHit*> vec_hits;
    vector<SttHit*> vec_hits_new;
    vec_driftradius.clear();
    vec_emcX.clear();
    vec_emcY.clear();
    //loop over all the vectors in the tree.

    for (Int_t i = 0; i < iev; i++) {
	Int_t count = 0;
    	Int_t count1 = 0;

        tree->GetEntry(i);
        cout << endl;
        cout << "entry no. " << i << endl;
        Int_t oiv = Vec_o_test->size();
        cout << "vecsize  :" << oiv << endl;

        if (Vec_o_test->size() > 8) {
            cout << "SIZE GREATER THAN 8" << endl;
        }

        //vector<SttHit*> vec_hits;
        vec_hits.clear();
        vec_hits_new.clear();
        vec_dr_xaxis[0].clear();
        vec_dr_xaxis[1].clear();
        vec_dr_yaxis[0].clear();
        vec_dr_yaxis[1].clear();

        Int_t oc_count = 0;

//         char buff[200];
//         sprintf(buff, "can_%d", i);
//         TCanvas* c1 = new TCanvas(buff);
//         sprintf(buff, "h_%d", i);
        int x_min = -10;
        int x_max = 90;
        int x_bins = (x_max - x_min) * 5;
        int z_min = -10;
        int z_max = 80;
        int z_bins = (z_max - z_min) * 5;
// 
//         TH2I* h = new TH2I(buff, "h;x,y [cm];z [cm]", x_bins, x_min, x_max, z_bins, z_min, z_max);
//         c1->cd();
//         h->Draw();
//         // 		c1->Range(-10,-10,60,60);
// 
// 
//         //Loop over the vector elements
// 
// 
// 
        for (int n = 0; n < oiv; n++) {
            SttHit* a = new SttHit();

            a->drifttime = Vec_o_test->at(n);
            a->x = vec_x->at(n);
            a->y = vec_y->at(n);
            a->z = vec_z->at(n);
            a->layer = vec_layer->at(n);
//             a->module = vec_module->at(n);
//             a->fee = vec_fee->at(n);
//             a->fee_channel = vec_fee_ch->at(n);
//             a->channel = vec_tdc_ch->at(n);
 
            vec_hits.push_back(a);

            float uuu = 0.0;
            bool is_y = (a->layer % 2 == 0);
            if (is_y)
                uuu = a->y;
            else
                uuu = a->x;
// 
//             straws = new TEllipse(uuu, vec_z->at(n), 0.505, 0.505);
//           //  box = new TBox(25,60,43,78);
// 		 //	box->SetFillColor(37);
// 
//             //cout<<"see :"<<a->drifttime<<"\t"<<a->layer<<"\t"<<a->module<<"\t"<<a->fee<<"\t"<<a->fee_channel<<"\t"<<a->x<<"\t"<<a->y<<"\t"<<a->z<<endl;
//             /////////
// 
            cout << "\t" << vec_hits.size() << "\t" << vec_hits[n]->drifttime << "\t" << vec_hits[n]->x << "\t" << vec_hits[n]->y << "\t" << vec_hits[n]->z << "\t" << endl;
//             straws->Draw("same");
//           //  box->Draw("same");
// 
         }

        //Filter to get only the sets having unique Z values

        if (vec_hits.size() > 0) {
            for (Int_t u = 0; u < vec_hits.size() - 1; u++) {
                for (Int_t uu = u + 1; uu < vec_hits.size(); uu++) {
                    if (fabs(vec_hits[u]->z - vec_hits[uu]->z) < 0.5)
                        oc_count++;

                    //cout<<vec_hits[u]->z<<"\t"<<vec_hits[uu]->z<<endl;
                }
            }
            //cout<<oc_count<<endl;

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

            Z_value->Fill(3);

            //Take only the events having unique Z coordinate

            if (oc_count >= 0) {
                //cout << "oc_count :";
                //cout<<oc_count<<"\t"<<vec_hits.size()<<endl;
                vec_hits_new = vec_hits;
                //cout<<vec_hits_new.size()<<endl;;
                Z_value->Fill(1);
            }

            else {

                Z_value->Fill(2);

                continue;
            }

            ///////////////////////////////////////////////////////////////////

            // Draw an ellipse on the straw position
            vec_strawX.clear();
            vec_strawY.clear();
            vec_strawZx.clear();
            vec_strawZy.clear();
             for (Int_t j = 0; j < vec_hits_new.size(); j++) {
// 
//                 //cout<<"THIS IS Z :"<<vec_hits_new[j]->z<<endl;
//                /* strawX[j] = vec_hits_new[j]->x;
//                 strawY[j] = vec_hits_new[j]->y;
//                 strawZ[j] = vec_hits_new[j]->z;*/
// 
//                 float uuu = 0.0;
//                 bool is_y = (vec_hits_new[j]->layer % 2 == 0);
//                 if (is_y)
//                     uuu = vec_hits_new[j]->y;
//                 else
//                     uuu = vec_hits_new[j]->x;
//                 straws2 = new TEllipse(uuu, vec_hits_new[j]->z, 0.505, 0.505);
//                 if (is_y)
//                     straws2->SetFillColor(46);
//                 else
//                     straws2->SetFillColor(38);
//                 straws2->SetFillStyle(1001);
//                 straws2->SetLineColor(1);
//                 straws2->SetLineWidth(1);
//                 straws2->Draw("same");
// 
//                 // Get the drift radius for the corresponding drift time for the hits having x-coordinates
 
                if (vec_hits_new[j]->layer == 1 || vec_hits_new[j]->layer == 3) {

                    Double_t dr1 = gDR->Eval(vec_hits_new[j]->drifttime);
                    SttHit* new_x1 = new SttHit();
                    SttHit* new_x2 = new SttHit();

                    vec_driftradius.push_back(dr1);

                    new_x1->x = vec_hits_new[j]->x + dr1;
                    new_x1->y = vec_hits_new[j]->y;
                    new_x1->z = vec_hits_new[j]->z;
                    new_x1->layer = vec_hits_new[j]->layer;
//                     new_x1->module = vec_hits_new[j]->module;
//                     new_x1->fee = vec_hits_new[j]->fee;
//                     new_x1->fee_channel = vec_hits_new[j]->fee_channel;

                    vec_strawX.push_back(vec_hits_new[j]->x);
                    vec_strawZx.push_back(vec_hits_new[j]->z);

                    new_x2->x = vec_hits_new[j]->x - dr1;
                    new_x2->y = vec_hits_new[j]->y;
                    new_x2->z = vec_hits_new[j]->z;
                    new_x2->layer = vec_hits_new[j]->layer;
//                     new_x2->module = vec_hits_new[j]->module;
//                     new_x2->fee = vec_hits_new[j]->fee;
//                     new_x2->fee_channel = vec_hits_new[j]->fee_channel;
// 
                    vec_dr_xaxis[0].push_back(new_x1);
                    vec_dr_xaxis[1].push_back(new_x2);
                    //zaxis[j]= vec_hits[j]->z;

                    //cout<< "Xaxis :  "<<vec_hits[j]->layer<<"\t"<<vec_hits[j]->x<<"\t"<<vec_hits[j]->y<<"\t"<<vec_hits[j]->z<<"\t"<<vec_hits[j]->x<<"\t"<<dr1<<endl;
                    count++;
                 }

                // Get the drift radius for the corresponding drift time for the hits having y-coordinates

                if (vec_hits_new[j]->layer == 2 || vec_hits_new[j]->layer == 4) {
                   
                    Double_t dr2 = gDR->Eval(vec_hits_new[j]->drifttime);

                    SttHit* new_y1 = new SttHit();
                    SttHit* new_y2 = new SttHit();

                    new_y1->y = vec_hits_new[j]->y + dr2;
                    new_y1->z = vec_hits_new[j]->z;
                    new_y1->x = vec_hits_new[j]->x;
                    new_y1->layer = vec_hits_new[j]->layer;
//                     new_y1->module = vec_hits_new[j]->module;
//                     new_y1->fee = vec_hits_new[j]->fee;
//                     new_y1->fee_channel = vec_hits_new[j]->fee_channel;
// 
                    /*if (vec_hits_new[j]->y >0)
                    {
                    strawY[j] = vec_hits_new[j]->y;
                    }*/
                    //cout<< "CHECK ::  "<<vec_hits_new[j]->y<<endl;
                    vec_strawY.push_back(vec_hits_new[j]->y);
                    vec_strawZy.push_back(vec_hits_new[j]->z);

                    new_y2->y = vec_hits_new[j]->y - dr2;
                    new_y2->z = vec_hits_new[j]->z;
                    new_y2->x = vec_hits_new[j]->x;
                    new_y2->layer = vec_hits_new[j]->layer;
//                     new_y2->module = vec_hits_new[j]->module;
//                     new_y2->fee = vec_hits_new[j]->fee;
//                     new_y2->fee_channel = vec_hits_new[j]->fee_channel;

                    //SttHit* new_y1=vec_hits[j]->y+dr2;
                    //SttHit* new_y2=vec_hits[j]->y-dr2;
                    vec_dr_yaxis[0].push_back(new_y1);
                    vec_dr_yaxis[1].push_back(new_y2);

                    yperfectmean[j] = dr2;
                    //zaxis[j]= vec_hits[j]->z;
                    //cout <<"\t"<< "Yaxis :"<<vec_hits[j]->y<<vec_hits[j]->layer<<endl;
                    //<<"\t"<<vec_hits[j]->layer<<"\t"<<dr2<<endl;

                    count1++;
                }
 
             }

            /*for (Int_t jx=0; jx<4; jx++)
            {
            	 cout <<"Y here  :"<<" "<<strawY[jx]<<endl;

            }*/


            /////////////////////////////////////////FIND THE BEST COMBINATION//////////////////////////////////////////////////////

            //COMBINATIONS to get least ChiSquar

            Int_t e_index = 0;
            std::vector<Double_t> chiX_array;
            Double_t chi_value;

            vector<vector<Double_t> > myCombination;
            myCombination.clear();


            Double_t dSlope0 =0; 
            Double_t dConst0 =0;
            Double_t X_perpX = 0;
            Double_t X_perpY = 0;
            Double_t X_short = 0;
            int no_comb_x = pow(2, count);
            chiX_array.resize(no_comb_x);

	    for (Int_t co = 0; co < no_comb_x; co++) {
/*                cout<<no_comb_x<<endl;
                cout<<"{"*/;

                vector<Double_t> vt;
                vt.clear();

		for (Int_t coo = 0; coo < count; coo++) {
		  
                	int bit_idx = (1 << coo);
                	int comb_idx = (co & bit_idx) ? 1 : 0;

                    A1[coo] = vec_dr_xaxis[comb_idx][coo]->x;
                    A2[coo] = vec_dr_xaxis[comb_idx][coo]->z;

                    //cout<<"   STRAW_LOC  @@@@@@@@@@@@@@@@&&&&&&&&&:"<<vec_dr_xaxis[comb_idx][coo]->x<<"  ";

                    vt.push_back(vec_dr_xaxis[comb_idx][coo]->x);
                }

                myCombination.push_back(vt);

                //cout<<"}";
                cout << "\n";

                //cout <<"size1 A1 :"<<(sizeof(A1))/8<<endl;
                //cout <<"size2 A2:"<<(sizeof(A2))/8<<endl;
                //cout <<"size3 xperfect2:"<<(sizeof(xperfect2))/8<<endl;

             /*   for (Int_t cooo=0; cooo< myCombination;cooo++)
			{
				for( Int_t cro=0; cro < vt.size(); cro++)
				{
					cout << "Hello"<< " ";
				}
			}*/
			cout<<"\n"<<endl;
			
                TGraph* chiX = new TGraph(vt.size(), A1, A2);
                chiX->Fit(f3);
                chi_value = f3->GetChisquare();
                	cout<<"\n\nX chi value"<<chi_value<<endl;

                chiX_array[co] = chi_value;
            }

            //Get the index of the combination having the least chisquare.
            Double_t smallest = chiX_array[0];

            Int_t chi_index = 0;
            for (Int_t ci = 0; ci < no_comb_x; ci++) {
                cout << "chiSquare  :" << chiX_array[ci] << endl;

                if (smallest > chiX_array[ci]) {
                    smallest = chiX_array[ci];
                    chi_index = ci;
                }
            }

            Double_t Ex[50];

            cout << "\n\n"
                 << "BEST COMBINATION  :";
            for (Int_t gdd = 0; gdd < count; gdd++) {
                cout << "\t " << myCombination.at(chi_index).at(gdd) << " ";
                Ex[gdd]=myCombination.at(chi_index).at(gdd);
            }
            cout << "\n\n\n" << endl;

            cout << "\n SMALLEST ChiSquar  :" << smallest << "(" << chi_index << ")" << endl;

            for (Int_t coo = 0; coo < count; coo++) {
             	int bit_idx = (1 << coo);
               	int comb_idx = (chi_index & bit_idx) ? 1 : 0;

                A1[coo] = vec_dr_xaxis[comb_idx][coo]->x;
                A2[coo] = vec_dr_xaxis[comb_idx][coo]->z;
            }
             //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            for (Int_t cdf = 0; cdf < vec_hits_new.size(); cdf++) {
                cout << "  XX  :" << strawX[cdf] << "  ZZ  :" << strawZ[cdf] << endl;
            }


            //Plot X vs Z coordinate and fit the points

            if (vec_dr_xaxis[0].size() != 0) {

                xfit = new TGraph(count, Ex, A2);
                xfit1 = new TGraph(vec_hits_new.size(), strawX, strawZ);

                xfit->GetXaxis()->SetLimits(x_min, x_max);
                xfit->GetYaxis()->SetLimits(z_min, z_max);

                xfit->GetXaxis()->SetTitle("X Axix [cm]");
                xfit->GetYaxis()->SetTitle("Z Axix [cm]");
                xfit->SetTitle("X-Coordinates");

                xfit->SetMarkerStyle(7);
                xfit->SetMarkerSize(1);
                xfit->SetMarkerColor(kRed + 2);
                xfit->SetLineColor(0);
                xfit->SetLineWidth(2);

                xfit1->SetLineColor(0);

                xfit1->SetMarkerStyle(24);
                xfit1->SetMarkerSize(2);
                xfit1->SetMarkerColor(kRed + 2);
                //xfit1->SetLineColor(kRed-10);
                //xfit1->SetLineWidth(2);

                xfit->Fit(f1);

                f1 = xfit->GetFunction("f1");
                Double_t p0 = f1->GetParameter(0);
                Double_t p1 = f1->GetParameter(1);
                theta_X = atan (p1);
                cout<< " THETA X  : "<< theta_X <<endl;
        		//emcX_arry[i] = (70 - p0) / p1;
        		vec_emcX.push_back((70 - p0) / p1);

                for (Int_t k = 0; k < count; k++) {

                    xperfect3[k] = (A2[k] - p0) / p1;

                    cout << "SPACE X  :" << fabs(A1[k] - xperfect3[k]) << endl;

                    dSlope0 = -(1/p1);
                    dConst0 = vec_strawZx[k]+(vec_strawX[k]/p1);
                    X_perpX = (dConst0 - p0) / (p1 + (1/p1));
                    X_perpY = (dSlope0 * X_perpX) + dConst0;
                    X_short = (fabs(sqrt(((vec_strawX[k] - X_perpX) * (vec_strawX[k] - X_perpX)) + ((vec_strawZx[k] - X_perpY) * (vec_strawZx[k] - X_perpY)) ))) - (fabs(vec_strawX[k] - A1[k]));

                    hx -> Fill(X_short);
                    h_theeta_X->Fill(theta_X,X_short);

                    vec_Xtheta.push_back(theta_X);
                    vec_Xsr.push_back(X_short);


                }
                cout << "\n\n" << endl;
            }

            ///////////////////////////////////////////////////////////////////////////////////
            //COMBINATIONS to get least Y ChiSquar2############################################

            Int_t e_index1 = 0;
            std::vector<Double_t> chiY_array;
            Double_t chi_value1;

            Double_t dSlope =0; 
            Double_t dConst =0;
            Double_t Y_perpX = 0;
            Double_t Y_perpY = 0;
            Double_t Y_short = 0;

            vector<vector<Double_t> > myCombination1;
            myCombination1.clear();

            int no_comb_y = pow(2, count1);
            chiY_array.resize(no_comb_y);
            for (Int_t co = 0; co < no_comb_y; co++) {
                
                vector<Double_t> vt1;
                vt1.clear();

                for (Int_t coo = 0; coo < count1; coo++) {
                	int bit_idx = (1 << coo);
                	int comb_idx = (co & bit_idx) ? 1 : 0;

                    B1[coo] = vec_dr_yaxis[comb_idx][coo]->y;
                    B2[coo] = vec_dr_yaxis[comb_idx][coo]->z;
                    vt1.push_back(vec_dr_yaxis[comb_idx][coo]->y);
                }

                myCombination1.push_back(vt1);

                //cout<<"}";
                cout << "\n";

                cout << "size1 B1 :" << (sizeof(B1)) / 8 << endl;
                cout << "size2 B2:" << (sizeof(B2)) / 8 << endl;
//                cout << "size3 yperfect2:" << (sizeof(yperfect2)) / 8 << endl;

                for (Int_t coooq = 0; coooq < total_pairs; coooq++) {
                    cout << "\t\t" << B1[coooq] << "\t";
                }
                cout << "\n" << endl;

                TGraph* chiY = new TGraph(vt1.size(), B1, B2);
                chiY->Fit(f4);
                chi_value1 = f4->GetChisquare();
                //	cout<<"\n\n@#$@#$@#$@#$@#$"<<chi_value<<endl;

                chiY_array[co] = chi_value1;
            }

            Double_t smallest1 = chiY_array[0];

			Int_t chi_index1 = 0;
            for (Int_t ciq = 0; ciq < no_comb_y; ciq++) {
                cout << "chiSquareY  :" << chiY_array[ciq] << endl;

                if (smallest1 > chiY_array[ciq]) {
                    smallest1 = chiY_array[ciq];
                    chi_index1 = ciq;
                }
            }

            Double_t Why[50];

            cout << "\n\n"
                 << "BEST COMBINATION 2 :";
            for (Int_t gddq = 0; gddq < count1; gddq++) {
                cout << "\t " << myCombination1.at(chi_index1).at(gddq) << " ";
                Why[gddq]=myCombination1.at(chi_index1).at(gddq);
            }
            cout << "\n\n\n" << endl;

            cout << "\n SMALLEST ChiSquar 2 :" << smallest1 << "(" << chi_index1 << ")" << endl;

           for (Int_t coo = 0; coo < count1; coo++) {
             	int bit_idx = (1 << coo);
               	int comb_idx = (chi_index & bit_idx) ? 1 : 0;

                B1[coo] = vec_dr_yaxis[comb_idx][coo]->y;
                B2[coo] = vec_dr_yaxis[comb_idx][coo]->z;
            }

            if (vec_dr_yaxis[0].size() != 0) {

                yfit = new TGraph(count1, Why, B2);
                yfit1 = new TGraph(vec_hits_new.size(), strawY, strawZ);

                //				yfit->GetXaxis()->SetLimits(40,100);

                //TGraph* yfit2 = new TGraph(vt1.size(),yperfectmean,yperfect2);

                yfit->GetXaxis()->SetTitle("Y Axix [cm]");
                yfit->GetYaxis()->SetTitle("Z Axix [cm]");
                yfit->SetTitle("Y-Coordinates");

                yfit->SetMarkerStyle(7);
                yfit->SetMarkerSize(2);
                yfit->SetMarkerColor(kBlue + 2);
                yfit->SetLineColor(0);
                yfit->SetLineWidth(2);

                yfit1->SetMarkerStyle(24);
                yfit1->SetMarkerSize(2);
                yfit1->SetMarkerColor(kBlue + 2);
                yfit1->SetLineColor(0);

                yfit->Fit(f2);

                f2 = yfit->GetFunction("f2");
                Double_t pp0 = f2->GetParameter(0);
                Double_t pp1 = f2->GetParameter(1);
        		//emcY_arry[i] = (70 - pp0) / pp1;
        	    vec_emcY.push_back((70 - pp0) / pp1);
	
        		//cout << " EMC   :   Y "<< ((70 - pp0) / pp1)<<endl;

                //yfit->SetPoint(1, 10,10);

                for (Int_t kk = 0; kk < count1; kk++) {

                    yperfect3[kk] = (B2[kk] - pp0) / pp1;
                    cout << "SPACE Y  :" << fabs(B1[kk] - yperfect3[kk]) << endl;
                   // hy->Fill(fabs(B1[kk] - yperfect3[kk]));

                   // cout << "Double Check Y :"<< vec_strawY[kk]<<endl;



                    dSlope = -(1/pp1);
                    dConst = vec_strawZy[kk]+(vec_strawY[kk]/pp1);
                    Y_perpX = (dConst - pp0) / (pp1 + (1/pp1));
                    Y_perpY = (dSlope * Y_perpX) + dConst;
                    Y_short = (fabs(sqrt(((vec_strawY[kk] - Y_perpX) * (vec_strawY[kk] - Y_perpX)) + ((vec_strawZy[kk] - Y_perpY) * (vec_strawZy[kk] - Y_perpY)) ))) - (fabs(vec_strawY[kk] - B1[kk]));

                     hy -> Fill(Y_short);

                    cout << "Y Real Check  :  " << vec_strawY[kk] << "\t" << vec_strawZy[kk] <<"\t" << Y_perpX <<"\t" <<Y_perpY<<"\t PERP : "<<Y_short << endl;
                }
            }

          //  cout << "\nEMC XYZ   :     "<< emcX << " - "<< emcY<<endl;
// 
//             //mg->Add(yfit);
//             //mg->Add(yfit2);
//             //xfit1->Draw("AP");
//             //mg->Add(yfit);
// 
//             /*mg->Add(xfit);
// 		mg->Add(xfit1);
// 		mg->Write("X");
// 		mg1->Add(yfit);
// 		mg1->Add(yfit1);
// 		mg1->Write("Y");*/
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
// 
//       // cout << "\n\nEMC COORDINATES : "<< emcX << "  -  "<< emcY<< endl;
//  
// 
// 
//         c1->Write();
//         c1->Close();



    }
    Double_t emcX_arry[500];
    Double_t emcY_arry[500];
    cout << "all driftradius   :" << vec_driftradius.size() << endl;
    
    for ( Int_t ev=0; ev< vec_emcX.size(); ev++)
    {
emcX_arry[ev] = vec_emcX[ev];
emcY_arry[ev] = vec_emcY[ev];
    }
         cout << "\n\n\n SIZES   :   " << vec_emcX.size()<< "  -  "<< vec_emcX.size()<<endl;
     	// cout << "\n\n\n SIZES   :   " << emcX_arry.length<< "  -  "<< emcY_arry.length<<endl;
           /*  for ( Int_t abc = 0; abc < vec_emcX.size(); abc++)
		    {
		    	cout << "\n\n\nEMC COORDINATES : "<< emcX_arry[abc] << "  -  "<< emcY_arry[abc]<< endl;
		    }  */ 
    TGraph* emc = new TGraph(vec_emcX.size(), emcX_arry, emcY_arry);
	emc->SetMarkerStyle(3);
	emc->SetMarkerSize(1);
	emc->SetMarkerColor(2);
	emc->SetLineColor(0);
	//emc->SetLineWidth(2);
    hx->Write();
    hy->Write();
    h_theeta_X->Write();

cout <<"THETA SIZE : "<<vec_Xsr.size()<<"  "<<vec_Xtheta.size()<<endl;

Double_t arr_Xt[500];
Double_t arr_Xsr[500];

for (int tt =0; tt< 500; tt++)
{
arr_Xt[tt]=vec_Xtheta[tt];
arr_Xsr[tt]=vec_Xsr[tt];
 }
   TGraph* theta_SR = new TGraph(500, arr_Xt, arr_Xsr);

    theta_SR->Write("P");
    Z_value->Write();
    emc->Draw("P");
    emc->Write();
    //emc->Write();
    //DR->Draw();
    //cout<<"count in Layer1+layer3 : "<< count<<"\t"<<"count in Layer2+layer4 : "<<"\t"<<count1<<"\t"<<"TOTAL : "<<count+count1<<endl;
    return kTRUE;
 
}


int main() {

	return PDAQ_Spl_Res();
}