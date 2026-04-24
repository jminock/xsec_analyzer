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

//Script that strips all weights and unnessecary variables out of ANNIE MC files
//Used for DV CV files
void skinny_legend_got_me_quakin(){

	int runs = 5000;//5000;
	int subruns = 1;
	//Loop through runs
        //Open file and trees
        string file_path = "/pnfs/annie/persistent/users/jminock/v1_3_4_world_stv_DV_ntuples/PhaseIITree_CV_stv_ntuple.root";
        string out_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_FD_stv_ntuple.root";

	//check if files exist
	if(gSystem->AccessPathName(file_path.c_str())){
		std::cout << "WARNING: " << file_path << " does not exist. Skipping." << std::endl;
	}
        //Open file and trees
        TFile *f = new TFile(file_path.c_str(),"read");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *T = (TTree*)f->Get("phaseIITriggerTree");

	T->SetBranchStatus("*", 0);
	//Reco
	T->SetBranchStatus("trigword", 1);
	T->SetBranchStatus("HasTank", 1);
	T->SetBranchStatus("HasMRD", 1);
	T->SetBranchStatus("TankMRDCoinc", 1);
	T->SetBranchStatus("NoVeto", 1);
	T->SetBranchStatus("simpleRecoFlag", 1);
	T->SetBranchStatus("reco0pi_contained_in_MRD", 1);
	T->SetBranchStatus("recoInc_contained_in_MRD", 1);
	T->SetBranchStatus("recoFV", 1);
	T->SetBranchStatus("recoPE", 1);
	T->SetBranchStatus("reco_0pi", 1);
	T->SetBranchStatus("simpleRecoCosTheta", 1);
	T->SetBranchStatus("simpleRecoMomentumCor", 1);
        T->SetBranchStatus("simpleRecoEnergy", 1);
	//Misc
	T->SetBranchStatus("MRDEff", 1);
	T->SetBranchStatus("DirtMu", 1);

        TFile *fnew = new TFile(out_path.c_str(),"recreate");
	auto Tnew = T->CloneTree();

	fnew->Write();
	delete f;
	delete fnew;
}
