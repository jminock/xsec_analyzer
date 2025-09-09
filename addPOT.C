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

void addPOT(){
        //insert variables here
        double pot;

        //Open file and trees
        string file_path = "/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_1mil_stv_ntuple_1.root";

	//check if files exist
	if(gSystem->AccessPathName(file_path.c_str())){
		std::cout << "WARNING: " << file_path << " does not exist. Skipping." << std::endl;
	}
 
        TFile *f = new TFile(file_path.c_str(),"update");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tTrig = (TTree*)f->Get("phaseIITriggerTree");

//      tTrig->SetBranchAddress("truePenetratesMRD",&mcpenetratesmrd);
        TBranch *tPOT    = tTrig->Branch("true_pot",&pot);

        double muon_m = 105.7;
        Long64_t nentriesTrig = tTrig->GetEntries();
        std::cout << "TriggerTree: " << nentriesTrig << std::endl;
        //fill histograms
        for (Long64_t i = 0; i < nentriesTrig; i++) {
                tTrig->GetEntry(i);
                if(i%1000 == 0) std::cout << i << std::endl;

		pot = 3.706115e20;

                tPOT->Fill();
	}

        tTrig->Write("",TObject::kOverwrite);
//      tTrig->ResetBranchAddresses();
        delete f;
}
