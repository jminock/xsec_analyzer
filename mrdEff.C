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

void breakCSV(string line, vector<double> &tokens){
	//change line into stream
	std::stringstream line_in(line);
	std::string temp_token;
	double token;
	while (line_in.good()){
		std::getline(line_in, temp_token, ','); //break up by comma
		std::stringstream to_token(temp_token);
		to_token >> token;          //turn into double
		tokens.push_back(token);    //save to array
	}
}

int findBin(double Y, int iter, vector<double> const &bins){
	if(iter == bins.size()) return 0;
	if(Y >= bins[iter]) return findBin(Y, iter+1, bins);
	else return iter;
}
void mrdEff(){
        //insert variables here
	double mrd_eff; //weight for MRD Efficiency correction

        vector<double>* MRDTrackStartX = new vector<double>();
        vector<double>* MRDTrackStartY = new vector<double>();
        vector<double>* MRDTrackStartZ = new vector<double>();
        vector<bool>* MRDStop = new vector<bool>();
        
        //Open file and trees

        //Open file and trees
        string file_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_tank_ntuple.root";
	string mrd_cal_file = "/exp/annie/app/users/jminock/ANNIE_AuxFiles/MRDEffCal.txt";

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
        tTrig->SetBranchAddress("MRDTrackStartX",&MRDTrackStartX);
        tTrig->SetBranchAddress("MRDTrackStartY",&MRDTrackStartY);
        tTrig->SetBranchAddress("MRDTrackStartZ",&MRDTrackStartZ);
        tTrig->SetBranchAddress("MRDStop",&MRDStop);

	TBranch *MRDEff = tTrig->Branch("MRDEff",&mrd_eff);

	//read in calibration file
	std::vector<double> bins;
	std::vector<double> factor;
	string bins_str = "";
	string factor_str = "";

	std::ifstream cal_file(mrd_cal_file.c_str(), ios::in);
	cal_file >> bins_str;
	cal_file >> factor_str;
	cal_file.close();

	breakCSV(bins_str, bins);
	breakCSV(factor_str, factor);

        double muon_m = 105.7;
        Long64_t nentriesTrig = tTrig->GetEntries();
        std::cout << "TriggerTree: " << nentriesTrig << std::endl;
        //fill histograms
        for(Long64_t i = 0; i < nentriesTrig; i++) {
                tTrig->GetEntry(i);
                if(i%1000 == 0) std::cout << i << std::endl;
                //Establish MRD Y variable to use if there are multiple tracks
		int ntracks = MRDTrackStartY->size();
		if(ntracks <= 0){ //assign 0 for no MRD tracks
			mrd_eff = 0.0;
			MRDEff->Fill();
			continue;
		}
		double MRD_Y = -9999;
		for(int j = 0; j < ntracks; j++){
			if(MRDStop->at(j)) MRD_Y = MRDTrackStartY->at(j); //confirm its stopping track
		}
		//Establish which bin the event falls into
		int bin = findBin(MRD_Y, 0, bins);
		//Assign weight based on bin
		mrd_eff = factor[bin];
		MRDEff->Fill();
	}


        tTrig->Write("",TObject::kOverwrite);
//      tTrig->ResetBranchAddresses();
        delete f;
}
