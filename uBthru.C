/**
 * Load an MC sample and populate kinematic histograms.
 *
 * \param name A unique string identifier
 * \param filename GENIE GST file to load
 * \returns A Sample with histograms
 */
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cmath>


bool FidVol(double x, double y, double z){
	double radius   = 168.1;  //cm
	double y_min    = -200.; //cm
	double y_max    = 200.;  //cm
	double z_center = 168.1; //cm
	double y_offset = 14.46; //cm
	if(y+y_offset > y_min && y+y_offset < y_max && radius > std::sqrt((z - z_center)*(z - z_center) + x*x)){
		return true;
	}
	else {
		return false;
	}
}

/** The main function. */
void uBthru() {
//  gROOT->SetBatch(true);
//  gStyle->SetOptStat(0);

//  TCanvas c;

  TFile *f = new TFile("/pnfs/annie/persistent/users/jminock/dirt/overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt.root", "READ");
//  TFile *f = new TFile("gntp.control.gst.root", "READ");
//  TFile *f = new TFile("gntp.v5.gst.root", "READ");
//  TFile *f = new TFile("gntp.v6.gst.root", "READ");

  TDirectory* dir = (TDirectoryFile*)f->Get("nuselection");

  TTree* nu = (TTree*)dir->Get("NeutrinoSelectionFilter");
  TTree* subrun = (TTree*)dir->Get("SubRun");

  vector<int>* mc_pdg = new vector<int>();
  vector<float>* mc_E = new vector<float>();
  vector<float>* mc_vx = new vector<float>();
  vector<float>* mc_vy = new vector<float>();
  vector<float>* mc_vz = new vector<float>();
  vector<float>* mc_px = new vector<float>();
  vector<float>* mc_py = new vector<float>();
  vector<float>* mc_pz = new vector<float>();

  float pot;


  nu->SetBranchAddress("mc_E", &mc_E);
  nu->SetBranchAddress("mc_pdg", &mc_pdg);
  nu->SetBranchAddress("mc_vx", &mc_vx);
  nu->SetBranchAddress("mc_vy", &mc_vy);
  nu->SetBranchAddress("mc_vz", &mc_vz);
  nu->SetBranchAddress("mc_px", &mc_px);
  nu->SetBranchAddress("mc_py", &mc_py);
  nu->SetBranchAddress("mc_pz", &mc_pz);
  subrun->SetBranchAddress("pot", &pot);

//  TCanvas *c1 = new TCanvas("TruevRecoTrackComp","Simple Reco Muon E 1 Track",900,600);
//  c1->SetGrid();

  TH1D* hmcvtxx = new TH1D("hmcvtxx","MC World Vertex X", 100, -840, 1200);
  TH1D* hmcvtxy = new TH1D("hmcvtxy","MC World Vertex Y", 100, -800, 1600);
  TH1D* hmcvtxz = new TH1D("hmcvtxz","MC World Vertex Z", 100, -1800, 1300);

  double muon_m = 105.7;
  int count_intank = 0;
  int count_cc = 0;
  int count_ccintank = 0;
  int count_FMV = 0;
  float xmin = 9999;
  float xmax = -9999;
  float ymin = 9999;
  float ymax = -9999;
  float zmin = 9999;
  float zmax = -9999;

  // Loop over MC events
  long nEntries = nu->GetEntries();
  for (long i=0; i<nEntries; i++) {
    nu->GetEntry(i);

    for(int j = 0; j < mc_pdg->size(); j++){
      if(mc_pdg->at(j) == 13){
        if(mc_vx->at(j) < xmin) xmin = mc_vx->at(j);
        if(mc_vx->at(j) > xmax) xmax = mc_vx->at(j);

        if(mc_vy->at(j) < ymin) ymin = mc_vy->at(j);
        if(mc_vy->at(j) > ymax) ymax = mc_vy->at(j);

        if(mc_vz->at(j) < zmin) zmin = mc_vz->at(j);
        if(mc_vz->at(j) > zmax) zmax = mc_vz->at(j);

//        hmcvtxx->Fill(mc_vx->at(j));
//        hmcvtxy->Fill(mc_vy->at(j));
//        hmcvtxz->Fill(mc_vz->at(j));

        bool hitsFMV = false;
        if(mc_vz->at(j) <= 0.){
          float time = std::abs(mc_vz->at(j))/(mc_pz->at(j)/muon_m);
          float FMVx = mc_vx->at(j) + (mc_px->at(j)/muon_m)*time;
          float FMVy = mc_vy->at(j) + (mc_py->at(j)/muon_m)*time;
          if(FMVx < 230. && FMVx > -90. && FMVy < 200. && FMVy > -200.) hitsFMV = true;
        }
        if(hitsFMV) count_FMV++;
      }
    }
  }

  long nEntries2 = subrun->GetEntries();
  float pot_total = 0;
  for(long i=0; i<nEntries2; i++){
    subrun->GetEntry(i);
    pot_total += pot;
  }

  std::cout << "Total:     " << nEntries << std::endl;
  std::cout << "POT:       " << pot_total << std::endl;
/*  std::cout << "InTank:    " << count_intank << std::endl;
  std::cout << "Total CC:  " << count_cc << std::endl;
  std::cout << "CC InTank: " << count_ccintank << std::endl;
*/  std::cout << "Hit FMV:   " << count_FMV << std::endl;
  std::cout << "X min: " << xmin << ", X max: " << xmax << std::endl;
  std::cout << "Y min: " << ymin << ", Y max: " << ymax << std::endl;
  std::cout << "Z min: " << zmin << ", Z max: " << zmax << std::endl;

/*  hmcvtxz->Draw();
  hmcvtxx->GetXaxis()->SetTitle("Vtx X [cm]");
  hmcvtxy->GetXaxis()->SetTitle("Vtx Y [cm]");
  hmcvtxz->GetXaxis()->SetTitle("Vtx Z [cm]");


  auto l1 = new TLine(0,0,0,800);
  l1->SetLineColor(kRed);
  l1->SetLineWidth(2);
  l1->Draw();

  auto l2 = new TLine(3.36,0,3.36,800);
  l2->SetLineColor(kRed);
  l2->SetLineWidth(2);
  l2->Draw();

  c1->Update();
*/
  nu->ResetBranchAddresses();
  subrun->ResetBranchAddresses();

}
