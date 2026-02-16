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

//Script that strips all weights (except TunedCentralValue) out of ANNIE MC files
//Used for non numuMC files that do not need weights
void stripWeights(){

	int runs = 4000;//4000;
	int subruns = 1;
	//Loop through runs
	for(int rn = 3500; rn < runs; rn++){
		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
        //Open file and trees
        string file_path = "/pnfs/annie/persistent/users/jminock/v1_3_3_world_stv_ntuples/PhaseIITree_0." + std::to_string(rn) + ".0.root";
        string out_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_0." + std::to_string(rn) + ".0.root";

	//check if files exist
	if(gSystem->AccessPathName(file_path.c_str())){
		std::cout << "WARNING: " << file_path << " does not exist. Skipping." << std::endl;
		continue;
	}
        //Open file and trees
        TFile *f = new TFile(file_path.c_str(),"read");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *T = (TTree*)f->Get("phaseIITriggerTree");

	T->SetBranchStatus("weight_All0_UBGenie", 0);
	T->SetBranchStatus("weight_All1_UBGenie", 0);
	T->SetBranchStatus("weight_All2_UBGenie", 0);
	T->SetBranchStatus("weight_All3_UBGenie", 0);
	T->SetBranchStatus("weight_All4_UBGenie", 0);
	T->SetBranchStatus("weight_All5_UBGenie", 0);
//	T->SetBranchStatus("weight_All_UBGenie", 0);
//	T->SetBranchStatus("weight_AxFFCCQEshape_UBGenie", 0);
//	T->SetBranchStatus("weight_DecayAngMEC_UBGenie", 0);
//	T->SetBranchStatus("weight_NormCCCOH_UBGenie", 0);
//	T->SetBranchStatus("weight_NormNCCOH_UBGenie", 0);
//	T->SetBranchStatus("weight_RPA_CCQE_UBGenie", 0);
//	T->SetBranchStatus("weight_RootinoFix_UBGenie", 0);
//	T->SetBranchStatus("weight_ThetaDelta2NRad_UBGenie", 0);
//	T->SetBranchStatus("weight_Theta_Delta2Npi_UBGenie", 0);
//	T->SetBranchStatus("weight_VecFFCCQEshape_UBGenie", 0);
//	T->SetBranchStatus("weight_XSecShape_CCMEC_UBGenie", 0);
//	T->SetBranchStatus("weight_flux_all", 0);
	T->SetBranchStatus("weight_horncurrent_FluxUnisim", 0);
	T->SetBranchStatus("weight_expskin_FluxUnisim", 0);
	T->SetBranchStatus("weight_pioninexsec_FluxUnisim", 0);
	T->SetBranchStatus("weight_piontotxsec_FluxUnisim", 0);
	T->SetBranchStatus("weight_pionqexsec_FluxUnisim", 0);
	T->SetBranchStatus("weight_nucleoninexsec_FluxUnisim", 0);
	T->SetBranchStatus("weight_nucleontotxsec_FluxUnisim", 0);
	T->SetBranchStatus("weight_nucleonqexsec_FluxUnisim", 0);
	T->SetBranchStatus("weight_piplus_PrimaryHadronSWCentralSplineVariation", 0);
	T->SetBranchStatus("weight_piminus_PrimaryHadronSWCentralSplineVariation", 0);
	T->SetBranchStatus("weight_kminus_PrimaryHadronNormalization", 0);
	T->SetBranchStatus("weight_kplus_PrimaryHadronFeynmanScaling", 0);
	T->SetBranchStatus("weight_kzero_PrimaryHadronSanfordWang", 0);

        TFile *fnew = new TFile(out_path.c_str(),"recreate");
	auto Tnew = T->CloneTree();

	fnew->Write();
	}
}
