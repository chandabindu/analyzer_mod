#include "TFile.h"
#include "TH1F.h"
#include "TChain.h"
#include "TH2F.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <math.h>
#include <utility>
#include <map>
#include <TCanvas.h>
#include <iterator>
const Int_t xch = 256;
const Int_t ych = 512;
vector < double> xcentroid;
vector <double> xsigma;
vector < double> ycentroid;
vector <double> ysigma;
void Sorting()
{
	cout<<" Start time "<<time(NULL)<<endl;
	int xapv,yapv,ij;
	Double_t xcentroid1;
	Double_t ycentroid1;
	Double_t xcharge;
	Double_t ycharge;
	//Pedestal information
	TString infile = "SelectPulse.txt";
	ifstream fh_infile;
	fh_infile.open(infile);
	if(!fh_infile)
	{
		cout<<"Cannot open pedestal file : "<<infile<<endl;
		return;
	}
	Double_t dumy1,dumy2;
	ij=0;
	while(fh_infile && ij<768)
	{
		fh_infile>>dumy1>>dumy2;
		if(ij<256)
		{
			xcentroid.push_back(dumy1);
			xsigma.push_back(dumy2);
			//cout<<"xij "<<ij<<" dumy1 "<<dumy1<<"  dumy2 "<<dumy2<<endl;

		}
		else
		{
			ycentroid.push_back(dumy1);
			ysigma.push_back(dumy2);
			//cout<<"yij "<<ij<<" dumy1 "<<dumy1<<"  dumy2 "<<dumy2<<endl;
		}
		ij++;
	
	}
	/*
	for(int pqr = 0; pqr<256;pqr++)
	{
		cout<<"xpqr "<<pqr<<"  xcentroid "<<xcentroid[pqr]<<" xsigma "<<xsigma[pqr]<<endl;
	}
	for(int pqr = 0; pqr<512;pqr++)
	{
		cout<<"ypqr "<<pqr<<"  ycentroid "<<ycentroid[pqr]<<" ysigma "<<ysigma[pqr]<<endl;
	}
	*/

	fh_infile.close();
	int NX = xcentroid.size();
	int NY = ycentroid.size();
    	xapv = NX/128;
	yapv = NY/128;
    	cout<<"NX "<<NX<<"  xapv  "<<xapv<<" NY "<<NY<<" yapv "<<yapv<<endl;
    	Int_t argc=746;
    	cout<<"Enter run number"<<endl;
    	cin>>argc;
    	TChain *T = new TChain("T");
    	T->Add(Form("~/rootfiles/test_%d.root",argc));
   
	Double_t x_strip, y_strip;//# of strips in X/Y
	Double_t xadc[6][xch];
    	Double_t yadc[6][ych];//6 Samples for x and Y
    	Double_t xstripID[xch];  
    	Double_t ystripID[ych];
  	Double_t dummy = 0;  
    	Double_t xcmn[2] = {-1.e9}; 
    	Double_t ycmn[4] = {-1.e9}; 
     	// Double_t integral[6][128] = {-1e10};
	vector<double>  xintegral;
	vector<double>  yintegral;
	
    	TH1D* h_xcmn[2];
    	TH1D* h_ycmn[4];
    	for(int pk=0;pk<xapv;pk++)
    	{
    		h_xcmn[pk] = new TH1D(Form("h_xcmn_%d",pk),Form("CMN for xapv %d",pk),300,-300,300);
    	}
    	for(int pk=0;pk<yapv;pk++)
    	{
    		h_ycmn[pk] = new TH1D(Form("h_ycmn_%d",pk),Form("CMN for yapv %d",pk),300,-300,300);
    	}

	TH2D *hits = new TH2D("XY hits","XY hits",512,-0.5,511.5,256,-0.5,255.5);
	TH2D *h_charge = new TH2D("charge","YX charge",100,-0.5,4999.5,100,-0.5,4999.5);
	
    	T->SetBranchAddress("sbs.gems.x1.nch", &x_strip);
	T->SetBranchAddress("sbs.gems.y1.nch", &y_strip);
    	T->SetBranchAddress("sbs.gems.x1.strip", xstripID);
    	T->SetBranchAddress("sbs.gems.y1.strip", ystripID);
    	for(int ij=0;ij<6;ij++)
	{
      		T->SetBranchAddress(Form("sbs.gems.x1.adc%d",ij),&xadc[ij]);
      		T->SetBranchAddress(Form("sbs.gems.y1.adc%d",ij),&yadc[ij]);
  	}
	cout<<"Enter event number "<<endl;
        int nevent = 2;
	cin>>nevent;
	int entries = T->GetEntries();
    	//for(int ij=0;ij<nevent;ij++)
	for(int ij=0;ij<entries;ij++)
	{
    		T->GetEntry(ij);
    		//cout<<"x_strip "<<x_strip<<"  y_strip "<<y_strip<<endl;
		for(int pq=0;pq<x_strip;pq++)
    		{
    			dummy=0;
    			for(int ik=0;ik<6;ik++)
    			{
				//cout<<"pq "<<pq<<" ik  "<<ik<<" xadc "<<xadc[ik][pq]<<endl;
    				dummy += xadc[ik][pq];
    			}
    			dummy = dummy/6.;
			xintegral.push_back(dummy);
			//cout<<"pq "<<pq<<" xdummy "<<dummy<<endl;
	     		// hx[pq]->Fill(dummy);
		}

    		vector <double> xadc_temp;
		for(int pk=0;pk<xch;pk++)
    		{
    			xadc_temp.push_back(xintegral[pk]-xcentroid[pk]);
			//cout<<"pk "<<pk<<"X Integral "<<xintegral[pk]<<" xcentroid "<<xcentroid[pk]<<"  xadc_temp "<<xadc_temp[pk]<<endl;	        
		}

    		for(int pk=0;pk<xapv;pk++)
    		{
    			xcmn[pk]=0;
    			vector<double> xadc_temp_sort;
    			//xadc_temp_sort.insert(xadc_temp_sort.end(),&xadc_temp[128*pk],&xadc_temp[128*(pk+1)]);
    			xadc_temp_sort.insert(xadc_temp_sort.begin(),&xadc_temp[128*pk],&xadc_temp[128*(pk+1)]);
    		/*	int xyz=0;
			for(std::vector<float>::iterator itv = xadc_temp_sort.begin();itv!=xadc_temp_sort.end();++itv)
			{
				cout<<"pk "<<pk<<" xyz "<<xyz<<" xadc "<<*itv<<endl;
				xyz++;
			}*/
			sort(xadc_temp_sort.begin(),xadc_temp_sort.end());
    			/*for(int ip=0;ip<128;ip++)
			 {
				 cout<<"ip  "<<ip<<" adc_temp_sort "<<adc_temp_sort[ip]<<endl;
				 
	     		 }*/
			
    			for(int ik=28;ik<100;ik++)
    			{
    				xcmn[pk] +=xadc_temp_sort[ik];
    			}
    			xcmn[pk] = xcmn[pk]/72.;
    			h_xcmn[pk]->Fill(xcmn[pk]);
    			//cout<<"x   apv  "<<pk<<"  cmn "<<xcmn[pk]<<endl;
			xadc_temp_sort.clear();
	    	}

    		for(int pq=0;pq<y_strip;pq++)
    		{
    			dummy = 0;
    			for(int ik = 0;  ik<6; ik++)
    			{
				//cout<<"pq "<<pq<<" ik  "<<ik<<" yadc "<<yadc[ik][pq]<<endl;
    				dummy += yadc[ik][pq];
    			}
    			dummy = dummy/6;
    			yintegral.push_back(dummy);
    			//hy[pq]->Fill(dummy);
	        }

    		vector <double> yadc_temp;
		for(int pk=0;pk<ych;pk++)
    		{
    			yadc_temp.push_back(yintegral[pk]-ycentroid[pk]);
    			//cout<<"pk "<<pk<<"Y Integral "<<yintegral[pk]<<" ycentroid "<<ycentroid[pk]<<"  yadc_temp "<<yadc_temp[pk]<<endl;
	        }
    		for(int pk=0;pk<yapv;pk++)
    		{
    			ycmn[pk]=0;
    			vector<double> yadc_temp_sort;
    			yadc_temp_sort.insert(yadc_temp_sort.begin(),&yadc_temp[128*pk],&yadc_temp[128*(pk+1)]);
    			sort(yadc_temp_sort.begin(),yadc_temp_sort.end());
    			/*for(int ip=0;ip<128;ip++)
			 {
				 cout<<"ip  "<<ip<<" adc_temp_sort "<<yadc_temp_sort[ip]<<endl;
				 
	     		 }*/
			
    			for(int ik=28;ik<100;ik++)
    			{
    				ycmn[pk] +=yadc_temp_sort[ik];
    			}
    			ycmn[pk] = ycmn[pk]/72.;
    			h_ycmn[pk]->Fill(ycmn[pk]);
    			//cout<<"y apv  "<<pk<<"  cmn "<<cmn[pk]<<endl;
			yadc_temp_sort.clear();
	    	}

    		map<int, double> hitx;
		map<int, TH1F*> xhitHisto;
		map<int, TCanvas*> C_xhits;
		int nxhits=0;
    		for(int pk=0;pk<xapv;pk++)
    		{
    			double corrch = 0;
    			for(int ik=0;ik<128;ik++)
    			{
    				corrch = xadc_temp[128*pk+ik]-xcmn[pk];
				//cout<<"pk "<<pk<<" ik "<<ik<<" xadc "<<xadc_temp[128*pk+ik]<<" xcmn "<<xcmn[pk]<<endl;
    				//if(-1*corrch>5*xsigma[128*pk+ik])
    				if(xintegral[128*pk+ik]<(xcentroid[128*pk+ik]-5*xsigma[128*pk+ik]))
    				//if(corrch>0)
    				{
					if(ij==nevent-1){
					xhitHisto[nxhits] = new TH1F(Form("hx_%d_%d",ij,(Int_t)xstripID[128*pk+ik]),Form("XEvent %d strip %d",ij,(Int_t)xstripID[128*pk+ik]),6,0,6); 
					for(int mn = 0;mn<6;mn++)
					{
						xhitHisto[nxhits]->SetBinContent(mn+1,xcentroid[128*pk+ik]-xadc[mn][128*pk+ik]+xcmn[pk]);
						//xhitHisto[nxhits]->SetBinContent(mn+1,xadc[mn][128*pk+ik]);
					}
					C_xhits[nxhits] = new TCanvas(Form("Cx_%d_%d",ij,(Int_t)xstripID[128*pk+ik]),Form("Canvas XEvent %d strip %d",ij,(Int_t)xstripID[128*pk+ik]),500,500); 
					xhitHisto[nxhits]->Draw("");
					nxhits++;}
					hitx.insert(pair <int, double> ((Int_t)xstripID[128*pk+ik],-1*corrch));
					//cout<<"Event "<<ij<<" X corrch "<<corrch<<" sigma  "<<xsigma[128*pk+ik]<<"  strip no "<<xstripID[128*pk+ik]<<endl; 
    				}
    			}
    		}
		map<int, double>::iterator itm;
		xcentroid1=0;
		ycentroid1=0;
		xcharge = 0;
		ycharge = 0;
		for(itm= hitx.begin(); itm!=hitx.end();++itm)
		{
			xcharge +=itm->second; 
			xcentroid1 += (itm->first)*(itm->second);
    			//cout<<"x strip id "<<itm->first<<" charge "<<itm->second<<endl;
    		}
		xcentroid1 = xcentroid1/xcharge;

    		map<int, double> hity;
		map<int, TH1F*>yhitHisto;
                map<int, TCanvas*> C_yhits;
                int nyhits=0;
    		for(int pk=0;pk<yapv;pk++)
    		{
    			double corrch = 0;
    			for(int ik=0;ik<128;ik++)
    			{
    				corrch = yadc_temp[128*pk+ik]-ycmn[pk];
    				//if(-1*corrch>5*ysigma[128*pk+ik])
    				if(yintegral[128*pk+ik]<(ycentroid[128*pk+ik]-5*ysigma[128*pk+ik]))
    				//if(corrch>0)
    				{
					if(ij==nevent-1){
					yhitHisto[nyhits] = new TH1F(Form("hy_%d_%d",ij,(Int_t)ystripID[128*pk+ik]),Form("YEvent %d strip %d",ij,(Int_t)ystripID[128*pk+ik]),6,0,6); 
					for(int mn = 0;mn<6;mn++)
					{
						yhitHisto[nyhits]->SetBinContent(mn+1,ycentroid[128*pk+ik]-yadc[mn][128*pk+ik]+ycmn[pk]);
						//yhitHisto[nxhits]->SetBinContent(mn+1,yadc[mn][128*pk+ik]);
					}
					C_yhits[nyhits] = new TCanvas(Form("Cy_%d_%d",ij,(Int_t)ystripID[128*pk+ik]),Form("Canvas YEvent %d strip %d",ij,(Int_t)ystripID[128*pk+ik]),500,500);
					yhitHisto[nyhits]->Draw("");
					nyhits++;}

					hity.insert(pair <int, double> ((Int_t)ystripID[128*pk+ik],-1*corrch));
					//cout<<"Event "<<ij<<" y corrch "<<corrch<<" sigma  "<<ysigma[128*pk+ik]<<"  strip no "<<ystripID[128*pk+ik]<<endl; 
    				}
    			}
    		}
		for(itm= hity.begin(); itm!=hity.end();++itm)
		{
			ycharge += itm->second;
			ycentroid1 += (itm->first)*(itm->second);
    			//cout<<"y strip id "<<itm->first<<" charge "<<itm->second<<endl;
    		}
		ycentroid1 = ycentroid1/ycharge;
		//if(xcentroid1>0 && ycentroid1>0)
		{
		//	cout<<"Entries  "<<ij<<"  xcentroid "<<xcentroid1<<"  ycentroid "<<ycentroid1<<" xcharge "<<xcharge<<" ycharge "<<ycharge<<endl;
		        hits->Fill(ycentroid1,xcentroid1);
			h_charge->Fill(xcharge,ycharge);
		}



    		xadc_temp.clear();
    		xintegral.clear();
    		yadc_temp.clear();
    		yintegral.clear();
    		hitx.clear();
    		hity.clear();
		//cout<<endl;
    	}

	TCanvas *C_cmn = new TCanvas("CMN","CMN",900,600);
    	C_cmn->Divide(2,3);
    	for(int ip=0;ip<6;ip++)
    	{
    		C_cmn->cd(ip+1);
    		if(ip<2)
			h_xcmn[ip]->Draw();
		else
			h_ycmn[ip-2]->Draw();

    	}

	TCanvas *C_hits = new TCanvas("Hits","Hits",800,600);
	hits->Draw("COLZ");
	TCanvas *C_charge = new TCanvas("charge","charge",800,600);
	h_charge->Draw("COLZ");
}
