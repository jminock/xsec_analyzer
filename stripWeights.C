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
        //Open file and trees
        TFile *f = new TFile("PhaseIITree_40k_ntuple.root","read");
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
	T->SetBranchStatus("weight_All_UBGenie", 0);
	T->SetBranchStatus("weight_AxFFCCQEshape_UBGenie", 0);
	T->SetBranchStatus("weight_DecayAngMEC_UBGenie", 0);
	T->SetBranchStatus("weight_NormCCCOH_UBGenie", 0);
	T->SetBranchStatus("weight_Norm_NCCOH_UBGenie", 0);
	T->SetBranchStatus("weight_RPA_CCQE_UBGenie", 0);
	T->SetBranchStatus("weight_RootinoFix_UBGenie", 0);
	T->SetBranchStatus("weight_ThetaDelta2NRad_UBGenie", 0);
	T->SetBranchStatus("weight_Theta_Delta2Npi_UBGenie", 0);
	T->SetBranchStatus("weight_VecFFCCQEshape_UBGenie", 0);
	T->SetBranchStatus("weight_XSecShape_CCMEC_UBGenie", 0);
	T->SetBranchStatus("weight_flux_all", 0);
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

        TFile *fnew = new TFile("PhaseIITree_40k_DVCV_ntuple.root","recreate");
	auto Tnew = T->CloneTree();

	fnew->Write();
}
