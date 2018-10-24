#include "TFile.h"
#include "TH1F.h"
#include <assert.h>
#include <utility>
#include "TF1.h"
#include "TChain.h"
#include "TMath.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

#include <TCanvas.h>

Double_t fitf(Double_t *v, Double_t *par)
{
  	Double_t arg = 0;
  	if (par[2] != 0) arg = (v[0] - par[1])/par[2];
	Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
   	return fitval;
}


void PedCal()
//void PlotGEMData()
{
    Int_t argc=-10;
    cout<<"Enter run number"<<endl;
    cin>>argc;
    TChain *T = new TChain("T");
    T->Add(Form("~/rootfiles/test_%d.root",argc));

    ofstream fh_out;
    TString outfile="SelectPulse.txt";
    fh_out.open(outfile);
    if(!fh_out)
    {
	    cout<<"Cannot open output file " <<outfile<<endl;
	    return;
    }
  //  cout<<"hi Chandan 1"<<endl;
    const int xch = 256;
    const int ych = 512;
    Double_t x_strip, y_strip;//# of strips in X/Y

    Double_t xadc[6][xch];
    Double_t yadc[6][ych];//6 Samples for x and Y
    Double_t xstripID[xch];  
    Double_t ystripID[ych];
    Double_t dummy = 0;  
    

    T->SetBranchAddress("sbs.gems.x1.nch", &x_strip);
    T->SetBranchAddress("sbs.gems.y1.nch", &y_strip);
    T->SetBranchAddress("sbs.gems.x1.strip", xstripID);
    T->SetBranchAddress("sbs.gems.y1.strip", ystripID);
    for(int ij=0;ij<6;ij++)
      {
  	  T->SetBranchAddress(Form("sbs.gems.x1.adc%d",ij),&xadc[ij]);
  	  T->SetBranchAddress(Form("sbs.gems.y1.adc%d",ij),&yadc[ij]);
      }
    TH1F *hx[xch],*hy[ych];
    Double_t xcentroid[xch]={-10000};
    Double_t ycentroid[ych]={-10000};
    Double_t xsigma[xch] = {-10};
    Double_t ysigma[ych] = {-10};
    Double_t xamp[xch] = {-100};
    Double_t yamp[ych] = {-100};


    for(int ij=0;ij<xch;ij++)
    {
	 hx[ij] = new TH1F(Form("hx_%d",ij+1),Form("Pedestal x_%d",ij+1),1500,0,1500);
    }	 
    for(int ij=0;ij<ych;ij++)
    {
	 hy[ij] = new TH1F(Form("hy_%d",ij+1),Form("Pedestal y_%d",ij+1),1500,0,1500);
    }	 
    int entries = T->GetEntries();
    

    for(int ij=0;ij<entries;ij++)
    //for(int ij=0;ij<10;ij++)
    //if(event<entries)
    {
	    T->GetEntry(ij);
	    //cout<<"x_strip "<<x_strip<<"  y_strip "<<y_strip<<endl;
	    for(int pq=0;pq<x_strip;pq++)
    	    {
    		    dummy=0;
		    for(int ik=0;ik<6;ik++)
		    {
			    dummy += xadc[ik][pq];
		    }
		    dummy = dummy/6.;
		    hx[pq]->Fill(dummy);

	    }
	    for(int pq=0;pq<y_strip;pq++)
    	    {
		    dummy = 0;
		    for(int ik = 0;  ik<6; ik++)
		    {
			    dummy += yadc[ik][pq];
		    }
		    dummy = dummy/6;
		    hy[pq]->Fill(dummy);
	    }
    }
    
    TCanvas *plotx[xch],*ploty[ych];

    for(int ij=0; ij<xch;ij++)
    {
	    hx[ij]->GetXaxis()->SetRangeUser(1,1500);
	    xcentroid[ij]=hx[ij]->GetMaximumBin()-1;
	    xamp[ij]=hx[ij]->GetBinContent(xcentroid[ij]);
	    xsigma[ij]=20;
	    TF1 *func = new TF1("fit",fitf,xcentroid[ij]-100,xcentroid[ij]+100,3);
	    /*if(ij%25==0)
	    cout<<"  xch  "<<ij<<" xcentroid  "<<xcentroid[ij]<<
	    "  xamp "<<xamp[ij]<<"  xsigma  "<<xsigma[ij]<<endl;*/
	    func->SetParameters(xamp[ij],xcentroid[ij],xsigma[ij]);
	    func->SetParNames("Constant","Mean_value","Sigma");
	    hx[ij]->Fit("fit","RQ0","",xcentroid[ij]-100,xcentroid[ij]+100);
	    fh_out<<func->GetParameter(1)<<"\t"<<fabs(func->GetParameter(2))<<endl;
	    if(ij%25==0){
	    cout<<"After fitting  xcentroid "<<func->GetParameter(1)<<"  xamplitude  "<<func->GetParameter(0)<<" xsigma  "<< func->GetParameter(2)<<endl;
	    plotx[ij]= new TCanvas(Form("xplot_%d",ij+1),Form("Ped x strip %d",ij+1),800,800);
	    plotx[ij]->SetLogy();
	    hx[ij]->Draw();
	    func->Draw("same");
	    }
    }
    for(int ij=0; ij<ych;ij++)
    {
	    hy[ij]->GetXaxis()->SetRangeUser(1,1500);
	    ycentroid[ij]=hy[ij]->GetMaximumBin()-1;
	    yamp[ij]=hy[ij]->GetBinContent(ycentroid[ij]);
	    ysigma[ij]=20;
	    TF1 *func1 = new TF1("fit1",fitf,ycentroid[ij]-100,ycentroid[ij]+100,3);
	    /*if(ij%50==0)
	    cout<<"  ych  "<<ij<<" ycentroid  "<<ycentroid[ij]<<
	    "  yamp "<<yamp[ij]<<"  ysigma  "<<ysigma[ij]<<endl;*/
	    func1->SetParameters(yamp[ij],ycentroid[ij],ysigma[ij]);
	    func1->SetParNames("Constant","Mean_value","Sigma");
	    hy[ij]->Fit("fit1","RQ0","",ycentroid[ij]-100,ycentroid[ij]+100);
	    //cout<<"After fitting  centroid "<<func1->GetParameter(1)<<"  yamplitude  "<<func1->GetParameter(0)<<" ysigma  "<< func1->GetParameter(2)<<endl;
	    fh_out<<func1->GetParameter(1)<<"\t"<<fabs(func1->GetParameter(2))<<endl;
	    if(ij%25==0){
	    cout<<"After fitting  centroid "<<func1->GetParameter(1)<<"  yamplitude  "<<func1->GetParameter(0)<<" ysigma  "<< func1->GetParameter(2)<<endl;
	    ploty[ij]= new TCanvas(Form("yplot_%d",ij+1),Form("Ped y strip %d",ij+1),800,800);
	    ploty[ij]->SetLogy();
	    hy[ij]->Draw();
	    func1->Draw("same");
	    }
    }
	fh_out.close();

    /*
    TCanvas *plot1 = new TCanvas("plot1","plot1",800,800);
    plot1->Divide(2,4);
    plot1->cd(1)->SetLogy();
    hx[25]->Draw();
    plot1->cd(2)->SetLogy();
    hx[50]->Draw();
    plot1->cd(3)->SetLogy();
    hx[75]->Draw();
    plot1->cd(4)->SetLogy();
    hx[100]->Draw();
    plot1->cd(5)->SetLogy();
    hx[125]->Draw();
    plot1->cd(6)->SetLogy();
    hx[150]->Draw();
    plot1->cd(7)->SetLogy();
    hx[175]->Draw();
    plot1->cd(8)->SetLogy();
    hx[200]->Draw();
		
    TCanvas *plot2 = new TCanvas("plot2","plot2",800,800);
    plot2->Divide(2,4);
    plot2->cd(1)->SetLogy();
    hy[50]->Draw();
    plot2->cd(2)->SetLogy();
    hy[100]->Draw();
    plot2->cd(3)->SetLogy();
    hy[150]->Draw();
    plot2->cd(4)->SetLogy();
    hy[200]->Draw();
    plot2->cd(5)->SetLogy();
    hy[250]->Draw();
    plot2->cd(6)->SetLogy();
    hy[300]->Draw();
    plot2->cd(7)->SetLogy();
    hy[350]->Draw();
    plot2->cd(8)->SetLogy();
    hy[400]->Draw();
    */
	
}
