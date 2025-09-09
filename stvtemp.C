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
        double radius   = 100.;  //cm
        double y_min    = -100.; //cm
        double y_max    = 100.;  //cm
        double z_center = 168.1; //cm
        double y_offset = 14.46; //cm
        if(y+y_offset > y_min && y+y_offset < y_max && radius > std::sqrt((z - z_center)*(z - z_center) + x*x)){
                return true;
        }
        else {
                return false;
        }
}

void stvtemp(){
        //Open file and trees
        TFile *f = new TFile("PhaseIITree_1mil_temp_nuwro_ntuple.root","update");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tTrig = (TTree*)f->Get("phaseIITriggerTree");

        //insert variables here
        int rnTrig;
        int rnMRD;
        int evNTrig;
        int evNMRD;
        int trigNum;
        int trigword, hasTank, hasMRD, tankMRDCoinc, noveto;
        int nrings;
        double nuE, nuvtxx, nuvtxy, nuvtxz, nupx, nupy, nupz;
        double fslpx, fslpy, fslpz, mcfslE;
        double fslvtxx, fslvtxy, fslvtxz;
        int isCCINC, isQEL, isRES, isDIS, isCOH, isMEC, mcfslpdg;
        double mcvtxx, mcvtxy, mcvtxz, mcdirx, mcdiry, mcdirz, mcangle, mcmuonE, mctanktracklength, mcmrdtracklength;
        int hasPi0, hasPiP, hasPiM, hasPiPC, hasPiMC, hasKP, hasKM, hasKPC, hasKMC;
        int numMRDTracks;
        int mcentersmrd, mcexitsmrd, mcpenetratesmrd;

	int isCC;

        bool mcfv;
	double mcp, mcct;
	bool mc_no_mesons, no_followers;
	bool recofv, recoMRDInc, recoMRD0pi, reco0pi;
	double recop, recopc, recoPE;

        int simpleflag, simplefv;
        double simpleenergy, simplecostheta, simplept, simplemrdenergy, simplemrdtrack, simpletanktrack;
        double simplevtxx, simplevtxy, simplevtxz, simplestopvtxx, simplestopvtxy, simplestopvtxz;
        double simplemrdstartx, simplemrdstarty, simplemrdstartz, simplemrdstopx, simplemrdstopy, simplemrdstopz;
        double rcsr, rcmr;
	double Qij, PE;
        double recovtxx, recovtxy, recovtxz;

//      vector<double>* MRDTrackLengthTrig = new vector<double>();
//      vector<double>* MRDTrackLengthMRD = new vector<double>();
        //TBranch *bwgts = 0;
        //Set branch addresses
//      tTrig->SetBranchAddress("weight_All_UBGenie",&All_weight);
//      tTrig->SetBranchAddress("weight_TunedCentralValue_UBGenie",&TCV_weight);
//      TBranch *TCV = tTrig->Branch("weight_TunedCentralValue_UBGenie",&TCV_weight);

        tTrig->SetBranchAddress("IsCCINC",&isCCINC);


//      tTrig->SetBranchAddress("truePenetratesMRD",&mcpenetratesmrd);
        TBranch *TFV    = tTrig->Branch("trueCC",&isCC);


//      tMRD->SetBranchAddress("runNumber",&rnMRD);
//      tMRD->SetBranchAddress("eventNumber",&evNMRD);
//      tMRD->SetBranchAddress("MRDTrackLength",&MRDTrackLengthMRD);


        double muon_m = 105.66;
        Long64_t nentriesTrig = tTrig->GetEntries();
        std::cout << "TriggerTree: " << nentriesTrig << std::endl;
        //fill histograms
        for (Long64_t i = 0; i < nentriesTrig; i++) {
//      for (Long64_t i = 0; i < 2; i++) {
                tTrig->GetEntry(i);
                if(i%10000 == 0) std::cout << i << std::endl;

		isCC = isCCINC;

                TFV->Fill();
	}


        tTrig->Write("",TObject::kOverwrite);
//      tTrig->ResetBranchAddresses();
        delete f;
}
