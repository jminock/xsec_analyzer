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
#include <algorithm>

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

void stvPrepData(){
        //insert variables here
        int rnTrig;
        int rnMRD;
        int evNTrig;
        int evNMRD;
        int trigword, hasTank, hasMRD, tankMRDCoinc, noveto;
        int numMRDTracks;
        int mcentersmrd, mcexitsmrd, mcpenetratesmrd;

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

        vector<double>* MRDTrackAngle = new vector<double>();
        vector<double>* MRDTrackAngleError = new vector<double>();
        vector<double>* MRDPenetrationDepth = new vector<double>();
        vector<double>* MRDTrackLength = new vector<double>();
        vector<double>* MRDEntryPointRadius = new vector<double>();
        vector<double>* MRDEnergyLoss = new vector<double>();
        vector<double>* MRDEnergyLossError = new vector<double>();
        vector<bool>* MRDSide = new vector<bool>();
        vector<bool>* MRDStop = new vector<bool>();
        vector<bool>* MRDThrough = new vector<bool>();
//      vector<double>* MRDTrackLengthTrig = new vector<double>();
//      vector<double>* MRDTrackLengthMRD = new vector<double>();
 

        //Open file and trees

        //Open file and trees
//        string file_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_2percent_data_ntuple.root";
        string file_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_0.0.0.root";

	//check if files exist
	if(gSystem->AccessPathName(file_path.c_str())){
		std::cout << "WARNING: " << file_path << " does not exist. Skipping." << std::endl;
	}
 
        TFile *f = new TFile(file_path.c_str(),"update");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tTrig = (TTree*)f->Get("phaseIITriggerTree");

        //Set branch addresses
//      tTrig->SetBranchAddress("weight_All_UBGenie",&All_weight);
//      tTrig->SetBranchAddress("weight_TunedCentralValue_UBGenie",&TCV_weight);
//      TBranch *TCV = tTrig->Branch("weight_TunedCentralValue_UBGenie",&TCV_weight);

        tTrig->SetBranchAddress("runNumber",&rnTrig);
        tTrig->SetBranchAddress("eventNumber",&evNTrig);
        tTrig->SetBranchAddress("trigword",&trigword);
        tTrig->SetBranchAddress("HasTank",&hasTank);
        tTrig->SetBranchAddress("HasMRD",&hasMRD);
        tTrig->SetBranchAddress("TankMRDCoinc",&tankMRDCoinc);
        tTrig->SetBranchAddress("NoVeto",&noveto);
        tTrig->SetBranchAddress("numMRDTracks",&numMRDTracks);
        tTrig->SetBranchAddress("MRDTrackAngle",&MRDTrackAngle);
        tTrig->SetBranchAddress("MRDTrackAngleError",&MRDTrackAngleError);
        tTrig->SetBranchAddress("MRDPenetrationDepth",&MRDPenetrationDepth);
        tTrig->SetBranchAddress("MRDTrackLength",&MRDTrackLength);
        tTrig->SetBranchAddress("MRDEntryPointRadius",&MRDEntryPointRadius);
        tTrig->SetBranchAddress("MRDEnergyLoss",&MRDEnergyLoss);
        tTrig->SetBranchAddress("MRDEnergyLossError",&MRDEnergyLossError);
        tTrig->SetBranchAddress("MRDSide",&MRDSide);
        tTrig->SetBranchAddress("MRDStop",&MRDStop);
        tTrig->SetBranchAddress("MRDThrough",&MRDThrough);

//      tTrig->SetBranchAddress("XSecWeights",&xsecweights);
//      tTrig->SetBranchAddress("FluxWeights",&fluxweights);

        //Simple Reco
        tTrig->SetBranchAddress("simpleRecoFlag",&simpleflag);
        tTrig->SetBranchAddress("simpleRecoEnergy",&simpleenergy);
        tTrig->SetBranchAddress("simpleRecoVtxX",&simplevtxx);
        tTrig->SetBranchAddress("simpleRecoVtxY",&simplevtxy);
        tTrig->SetBranchAddress("simpleRecoVtxZ",&simplevtxz);
        tTrig->SetBranchAddress("simpleRecoCosTheta",&simplecostheta);
        tTrig->SetBranchAddress("simpleRecoPt",&simplept);
        tTrig->SetBranchAddress("simpleRecoFV",&simplefv);
        tTrig->SetBranchAddress("simpleRecoMrdEnergyLoss",&simplemrdenergy);
        tTrig->SetBranchAddress("simpleRecoTrackLengthInMRD",&simplemrdtrack);
        tTrig->SetBranchAddress("simpleRecoTrackLengthInTank",&simpletanktrack);

        //Ring Counting Reco
//        tTrig->SetBranchAddress("RCSRPred",&rcsr);
//        tTrig->SetBranchAddress("RCMRPred",&rcmr);
	tTrig->SetBranchAddress("Qij",&Qij);
	tTrig->SetBranchAddress("promptMuonTotalPE",&PE);

        //true Muon info from WCSim
//      tTrig->SetBranchAddress("trueEntersMRD",&mcentersmrd);
//      tTrig->SetBranchAddress("trueExitsMRD",&mcexitsmrd);
//      tTrig->SetBranchAddress("truePenetratesMRD",&mcpenetratesmrd);
	TBranch *RFV    = tTrig->Branch("recoFV",&recofv);
	TBranch *RPE    = tTrig->Branch("recoPE",&recoPE);
	TBranch *TCT    = tTrig->Branch("trueCosTheta",&mcct);
	TBranch *RMuP   = tTrig->Branch("simpleRecoMomentum",&recop);
	TBranch *RMuPC  = tTrig->Branch("simpleRecoMomentumCor",&recopc);
	TBranch *RinMRDInc = tTrig->Branch("recoInc_contained_in_MRD",&recoMRDInc);
	TBranch *RinMRD0pi = tTrig->Branch("reco0pi_contained_in_MRD",&recoMRD0pi);
	TBranch *R0pi   = tTrig->Branch("reco_0pi",&reco0pi);

        double muon_m = 105.7;
        Long64_t nentriesTrig = tTrig->GetEntries();
        std::cout << "TriggerTree: " << nentriesTrig << std::endl;
        //fill histograms
        for (Long64_t i = 0; i < nentriesTrig; i++) {
                tTrig->GetEntry(i);
//                if(i%1000 == 0) std::cout << i << std::endl;
                bool simplerecoFV = FidVol(simplevtxx*100.,simplevtxy*100.,simplevtxz*100.);

		recofv = simplerecoFV;
		recoPE = ((PE > 500) && (PE < 3000));
		double recototE = simpleenergy + muon_m;
		recop  = std::sqrt(recototE*recototE - muon_m*muon_m);
		recopc = std::sqrt(recototE*recototE - muon_m*muon_m)*0.82 + 160.;
		recoMRD0pi = numMRDTracks == 1 ? MRDStop->at(0) : false;
		recoMRDInc = numMRDTracks > 0 ? MRDStop->at(0) : false;
		reco0pi = ((PE > 200.*Qij*Qij) && (PE < 2000.*std::cbrt(4.5-Qij)+1500.));

		RFV->Fill();
		RPE->Fill();
		RMuP->Fill();
		RMuPC->Fill();
		RinMRDInc->Fill();
		RinMRD0pi->Fill();
		R0pi->Fill();
	}


        tTrig->Write("",TObject::kOverwrite);
//      tTrig->ResetBranchAddresses();
        delete f;
}
