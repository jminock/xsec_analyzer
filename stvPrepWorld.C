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
#include <fstream>
#include <sstream>
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

void stvPrepWorld(){
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
        int isCC, isQEL, isRES, isDIS, isCOH, isMEC, mcfslpdg;
        double mcvtxx, mcvtxy, mcvtxz, mcdirx, mcdiry, mcdirz, mcangle, mcmuonE, mctanktracklength, mcmrdtracklength;
        int hasPi0, hasPiP, hasPiM, hasPiPC, hasPiMC, hasKP, hasKM, hasKPC, hasKMC;
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
	double mrd_eff; //weight for MRD Efficiency correction

        vector<double>* MRDTrackAngle = new vector<double>();
        vector<double>* MRDTrackAngleError = new vector<double>();
        vector<double>* MRDPenetrationDepth = new vector<double>();
        vector<double>* MRDTrackLength = new vector<double>();
        vector<double>* MRDEntryPointRadius = new vector<double>();
        vector<double>* MRDEnergyLoss = new vector<double>();
        vector<double>* MRDEnergyLossError = new vector<double>();
        vector<double>* MRDTrackStartX = new vector<double>();
        vector<double>* MRDTrackStartY = new vector<double>();
        vector<double>* MRDTrackStartZ = new vector<double>();
        vector<double>* MRDTrackStopX = new vector<double>();
        vector<double>* MRDTrackStopY = new vector<double>();
        vector<double>* MRDTrackStopZ = new vector<double>();
        vector<bool>* MRDSide = new vector<bool>();
        vector<bool>* MRDStop = new vector<bool>();
        vector<bool>* MRDThrough = new vector<bool>();
	vector<int>* mcFolPPDG = new vector<int>();
//      vector<double>* MRDTrackLengthTrig = new vector<double>();
//      vector<double>* MRDTrackLengthMRD = new vector<double>();
 

        //Open file and trees

        //Open file and trees
        string file_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_400k_test_ntuple.root";
	string mrd_cal_file = "/exp/annie/app/users/jminock/ANNIE_Aux_files/.txt"

	//check if files exist
	if(gSystem->AccessPathName(file_path.c_str())){
		std::cout << "WARNING: " << file_path << " does not exist. Skipping." << std::endl;
	}
	if(gSystem->AccessPathName(mrd_cal_file.c_str())){
		std::cout << "WARNING: " << mrd_cal_file << " does not exist. Stopping." << std::endl;
		return false;
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
        tTrig->SetBranchAddress("MRDTrackStartX",&MRDTrackStartX);
        tTrig->SetBranchAddress("MRDTrackStartY",&MRDTrackStartY);
        tTrig->SetBranchAddress("MRDTrackStartZ",&MRDTrackStartZ);
        tTrig->SetBranchAddress("MRDSide",&MRDSide);
        tTrig->SetBranchAddress("MRDStop",&MRDStop);
        tTrig->SetBranchAddress("MRDThrough",&MRDThrough);

        tTrig->SetBranchAddress("trueCC",&isCC);
        tTrig->SetBranchAddress("trueQEL",&isQEL);
        tTrig->SetBranchAddress("trueRES",&isRES);
        tTrig->SetBranchAddress("trueDIS",&isDIS);
        tTrig->SetBranchAddress("trueCOH",&isCOH);
        tTrig->SetBranchAddress("trueMEC",&isMEC);
        tTrig->SetBranchAddress("trueNuIntxVtx_X",&nuvtxx);
        tTrig->SetBranchAddress("trueNuIntxVtx_Y",&nuvtxy);
        tTrig->SetBranchAddress("trueNuIntxVtx_Z",&nuvtxz);
        tTrig->SetBranchAddress("trueNeutrinoEnergy",&nuE);
        tTrig->SetBranchAddress("trueNeutrinoMomentum_X",&nupx);
        tTrig->SetBranchAddress("trueNeutrinoMomentum_Y",&nupy);
        tTrig->SetBranchAddress("trueNeutrinoMomentum_Z",&nupz);
        tTrig->SetBranchAddress("truePi0",&hasPi0);
        tTrig->SetBranchAddress("truePiPlus",&hasPiP);
        tTrig->SetBranchAddress("truePiMinus",&hasPiM);
        tTrig->SetBranchAddress("truePiPlusCher",&hasPiPC);
        tTrig->SetBranchAddress("truePiMinusCher",&hasPiMC);
        tTrig->SetBranchAddress("trueKPlus",&hasKP);
        tTrig->SetBranchAddress("trueKMinus",&hasKM);
        tTrig->SetBranchAddress("trueKPlusCher",&hasKPC);
        tTrig->SetBranchAddress("trueKMinusCher",&hasKMC);
        tTrig->SetBranchAddress("trueFSLEnergy",&mcfslE);
        tTrig->SetBranchAddress("trueFSLPdg",&mcfslpdg);
	tTrig->SetBranchAddress("trueFollowerParentPDG",&mcFolPPDG);
        tTrig->SetBranchAddress("triggerNumber",&trigNum);

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
        tTrig->SetBranchAddress("trueMuonEnergy",&mcmuonE);
        tTrig->SetBranchAddress("trueVtxX",&mcvtxx);
        tTrig->SetBranchAddress("trueVtxY",&mcvtxy);
        tTrig->SetBranchAddress("trueVtxZ",&mcvtxz);
        tTrig->SetBranchAddress("trueDirX",&mcdirx);
        tTrig->SetBranchAddress("trueDirY",&mcdiry);
        tTrig->SetBranchAddress("trueDirZ",&mcdirz);
        tTrig->SetBranchAddress("trueTrackLengthInWater",&mctanktracklength);
        tTrig->SetBranchAddress("trueTrackLengthInMRD",&mcmrdtracklength);
        tTrig->SetBranchAddress("trueAngle",&mcangle);
//      tTrig->SetBranchAddress("trueEntersMRD",&mcentersmrd);
//      tTrig->SetBranchAddress("trueExitsMRD",&mcexitsmrd);
//      tTrig->SetBranchAddress("truePenetratesMRD",&mcpenetratesmrd);
        TBranch *TFV    = tTrig->Branch("trueFV",&mcfv);
	TBranch *RFV    = tTrig->Branch("recoFV",&recofv);
	TBranch *RPE    = tTrig->Branch("recoPE",&recoPE);
	TBranch *TMuP   = tTrig->Branch("trueMuonMomentum",&mcp);
	TBranch *TCT    = tTrig->Branch("trueCosTheta",&mcct);
	TBranch *RMuP   = tTrig->Branch("simpleRecoMomentum",&recop);
	TBranch *RMuPC  = tTrig->Branch("simpleRecoMomentumCor",&recopc);
	TBranch *TnoPi  = tTrig->Branch("true_no_mesons",&mc_no_mesons);
	TBranch *TnoF   = tTrig->Branch("true_no_followers",&no_followers);
	TBranch *RinMRDInc = tTrig->Branch("recoInc_contained_in_MRD",&recoMRDInc);
	TBranch *RinMRD0pi = tTrig->Branch("reco0pi_contained_in_MRD",&recoMRD0pi);
	TBranch *R0pi   = tTrig->Branch("reco_0pi",&reco0pi);
	TBranch *MRDEff = tTrig->Branch("MRDEff",&mrd_eff);

	//read in calibration file
	std::vector<double> bins;
	std::vector<double> factor;

        double muon_m = 105.7;
        Long64_t nentriesTrig = tTrig->GetEntries();
        std::cout << "TriggerTree: " << nentriesTrig << std::endl;
        //fill histograms
        for (Long64_t i = 0; i < nentriesTrig; i++) {
                tTrig->GetEntry(i);
//                if(i%1000 == 0) std::cout << i << std::endl;
                bool inFV = FidVol(nuvtxx,nuvtxy,nuvtxz);
                bool simplerecoFV = FidVol(simplevtxx*100.,simplevtxy*100.,simplevtxz*100.);
                bool hasPi = ((hasPi0) || (hasPiP) || (hasPiM));
                bool hasVisPi = ((hasPi0) || (hasPiPC) || (hasPiMC));
                bool hasNonMuon = ((hasPi) || (hasKP) || (hasKM));
                bool hasVisNonMuon = ((hasVisPi) || (hasKPC) || (hasKMC));

                mcfv = inFV;
		recofv = simplerecoFV;
		recoPE = ((PE > 500) && (PE < 3000));
		mcp = std::sqrt(mcmuonE*mcmuonE - muon_m*muon_m);
		mcct = std::cos(mcangle*M_PI/180.);
		double recototE = simpleenergy + muon_m;
		recop  = std::sqrt(recototE*recototE - muon_m*muon_m);
		recopc = std::sqrt(recototE*recototE - muon_m*muon_m)*0.82 + 160.;
		mc_no_mesons = !(hasVisNonMuon);
		no_followers = (std::find(mcFolPPDG->begin(), mcFolPPDG->end(), 211) != mcFolPPDG->end() || std::find(mcFolPPDG->begin(), mcFolPPDG->end(), -211) != mcFolPPDG->end()) ? false : true;
		recoMRD0pi = numMRDTracks == 1 ? MRDStop->at(0) : false;
		recoMRDInc = numMRDTracks > 0 ? MRDStop->at(0) : false;
		reco0pi = ((PE > 200.*Qij*Qij) && (PE < 2000.*std::cbrt(4.5-Qij)+1500.));


                TFV->Fill();
		RFV->Fill();
		RPE->Fill();
		TMuP->Fill();
		TCT->Fill();
		RMuP->Fill();
		RMuPC->Fill();
		TnoPi->Fill();
		TnoF->Fill();
		RinMRDInc->Fill();
		RinMRD0pi->Fill();
		R0pi->Fill();
		MRDEff->Fill();
	}


        tTrig->Write("",TObject::kOverwrite);
//      tTrig->ResetBranchAddresses();
        delete f;
}
