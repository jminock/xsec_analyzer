/*Author: P. Englezos <p.englezos@physics.rutgers.edu>
 *
 * Usage: ./chi_square_cc0pi files.txt
 *
 */

//#include "chi_square_cc0pi_christian.hh"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"
#include "EventCategory.hh"
#include "TreeUtils.hh"
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include <cassert>
#include <set>
#include <vector>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TPad.h>
#include "TColor.h"
#include "TInterpreter.h"
#include <algorithm>
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
//#include "PlotUtils.hh"
#include "TMatrixT.h"
#include "WienerSVDUnfolder.hh"
//#include "FiducialVolume.hh"
#include <iomanip>
#include "TMatrixD.h"
//#include "MatrixUtils.hh"
//#include "SliceBinning.hh"
//#include "SliceHistogram.hh"
#include "TLatex.h"
//#include "HistUtils.hh"
#include "RooStats/RooStatsUtils.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"


#include "includes/PlotUtils.hh"
#include "includes/AnnieGeometryTools.hh"

void chi_square_all_gens(std::string infiles);
void multiply_1d_hist_by_matrix(TMatrixD *mat, TH1 *hist);
bool FidVol(double x, double y, double z);
void Test_plot();


void multiply_1d_hist_by_matrix(TMatrixD *mat, TH1 *hist)
{
   int num_bins = mat->GetNcols();
    TMatrixD hist_mat(num_bins, 1);
    for (int r = 0; r < num_bins; ++r)
    {
        hist_mat(r, 0) = hist->GetBinContent(r + 1);
    }
   
    TMatrixD hist_mat_transformed(*mat, TMatrixD::EMatrixCreatorsOp2::kMult, hist_mat);

    for (int r = 0; r < num_bins; ++r)
    {
        double val = hist_mat_transformed(r, 0);
        hist->SetBinContent(r + 1, val);
    }
}

bool FidVol(double x, double y, double z){
  double radius   = 100.;  //cm
  double y_min    = -100.; //cm
  double y_max    = 100.;  //cm
  double z_center = 168.1;  //cm
  if(y > y_min && y < y_max && radius > std::sqrt((z - z_center)*(z - z_center) + x*x)){
    return true;
  }
  else {
    return false;
  }
}



///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
void chi_square_all_gens(std::string infiles) {
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);

   //  std::vector< double > nbins = {-1.,-0.775,-0.675,-0.575,-0.475,-0.4,-0.325,-0.25,-0.175,-0.1,-0.025,0.025,0.1,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.};
   
   //  std::vector<double> nbins = {0.8,0.95,1.0};
     std::vector<double> nbins = {600.,810.,910.,1005.,1100.,1200};
   
     //-------Most up-to-date root file for Genie Closure test----------//
     //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_costheta_closure_runs_1-3_April_17th.root", "systcalc.conf" ); 
   
     //-------Most up-to-date root file for NuWro fake data test (without AltCV at xsec uncertainties)----------//
     //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/without_nuwro_uncertainty/univmake_costheta_no_nuwro_uncertainty_runs_1-3_April_17th.root", "systcalc.conf" ); 
   
      //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_Philip_nuwro_Runs_1_to_3_momentum_all_uncertainties_scaled_as_usual_ATTEMPT_1.root", "systcalc.conf" );
       //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_Philip_genie_Runs_1_to_3_momentum_ATTEMPT_1.root", "systcalc.conf" );
       //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_Philip_nuwro_fake_data_study_Runs_1_to_3_momentum_ATTEMPT_1.root", "systcalc.conf" );
   
   //-------Most up-to-date root file for NuWro fake data test (including AltCV at xsec uncertainties)----------//
      //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_nuwro_fake_data_runs_1-3_momentum_old.root","systcalc.conf" );
      //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_nuwro_fake_data_runs_1-3_angle_old.root","systcalc.conf" );
   
      //after the flipped tracks rejection:
      //  auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_nuwro_fake_data_runs_1-3_momentum.root", "systcalc.conf" );
      //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_nuwro_fake_data_runs_1-3_costheta_2.root", "systcalc.conf" );
   
      //------------Genie with all uncertainties------------//
     //auto* mcc9 = new MCC9SystematicsCalculator("/uboone/data/users/englezos/jointxsec/ru_github/stv-analysis-new/stv-analysis-new/univmake_genie_with_uncertainties_runs1-3_muon_angle.root","systcalc.conf" );
      //auto* mcc9 = new MCC9SystematicsCalculator("./output_files/univmake_corrected-muon-angle_xsec_stat_nuwro_uncertainties.root", "systcalc.conf" );    

 

 //\\\\\\\\\\\\\\\\\\\\\\\ -> Up-to-Date//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  //NuWro with full uncertainties
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_Philip_nuwro_Runs_1_to_3_momentum_all_uncertainties_scaled_as_usual_ATTEMPT_1.root", "systcalc.conf" );
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_Philip_GITHUB_nuwro_with_all_uncertainties_Runs_1_to_3_April_26th_first_bin_up_to_0175_pmu.root", "systcalc.conf" );
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_Philip_nuwro_Runs_1_to_3_costheta_all_uncertainties_scaled_as_usual_April_26th_ATTEMPT_1.root", "systcalc.conf" );
  

  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_nuwro_ALL_uncertainties_pmu_May_10th_TEST_on_May_13th.root", "systcalc.conf" );
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_May_10th_all_uncertainties_nuwro_pmu_upto02.root", "systcalc.conf" );      

   //  auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_nuwro_ALL_uncertainties_costheta_May_10th_TEST.root", "systcalc.conf" );

  //Genie closure test
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_Philip_genie_Runs_1_to_3_momentum_ATTEMPT_1.root", "systcalc.conf" );
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_Philip_GITHUB_upto0175_pmu_April_26th_corrections.root", "systcalc.conf" );
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_Philip_genie_Runs_1_to_3_costheta_April_26th_ATTEMPT_1.root", "systcalc.conf" );

  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_genie_closure_May_10th_pmu.root", "systcalc.conf" );
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_May_10th_genie_closure_pmu_upto02.root", "systcalc.conf" );

  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/genie_closure_test/univmake_genie_May_10th_costheta.root", "systcalc.conf" );

  //NuWro fake data test
 //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_Philip_nuwro_fake_data_study_Runs_1_to_3_momentum_ATTEMPT_1.root", "systcalc.conf" );
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_Philip-GITHUB-_nuwro_Runs_1_to_3_momentum_all_uncertainties_scaled_as_usual_ATTEMPT_1.root","systcalc.conf");
  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/with_extra_uncertainties/univmake_Philip-GITHUB-_nuwro_Runs_1_to_3_momentum_all_uncertainties_scaled_as_usual_ATTEMPT_1_up_to_02.root","systcalc.conf");
   //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_Philip-GITHUB-_nuwro_Runs_1_to_3_momentum_all_uncertainties_scaled_as_usual_ATTEMPT_1_up_to_0175.root","systcalc.conf");
  // auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_Philip_nuwro_fake_data_study_Runs_1_to_3_costheta_April_26th_ATTEMPT_1.root", "systcalc.conf" );
   
   //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_nuwro_xsec_stats_May_10th_pmu.root","systcalc.conf");
   //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_May_10th_nuwro_fake_data_pmu_upto02.root","systcalc.conf"); 

  //auto* mcc9 = new MCC9SystematicsCalculator("./univmake_verified/nuwro_fake_data/only_xsec_and_mc_uncertainties/univmake_nuwro_xsec_stats_May_10th_costheta.root", "systcalc.conf" );   

  //ANNIE
  auto* mcc9 = new MCC9SystematicsCalculator(
  //    "/exp/annie/data/users/jminock/stv-analysis/stv-univmake-output-20k-pmu.root",
    "output.root",
    "systcalc.conf" );

  const auto &syst = *mcc9;

  TH1D* genie_cv_truth_vals = new TH1D("genie_cv_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  TH1D* fake_data_truth_vals = new TH1D("fake_data_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  TH1D* unfolded_events_vals = new TH1D("unfolded_events_vals", ";cos#theta_{#mu}; Scaled Events", nbins.size() - 1, nbins.data());

  std::cout<<"Now examining real data."<<std::endl;

  double total_pot = mcc9->total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  double num_Ar = num_O_targets_in_FV(); //5.25e28;
  double conv_factor = (num_Ar * integ_flux)/1e38;

  TH1D* genie_cv_truth = mcc9->cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();
  TH1D* unfolded_events = dynamic_cast< TH1D* >(genie_cv_truth->Clone("unfolded_events") );
  unfolded_events->Reset(); 

  const auto& fake_data_univ = mcc9->fake_data_universe();
  TH1D* fake_data_truth = fake_data_univ->hist_true_.get(); 

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto smearcept_ptr = syst.get_cv_smearceptance_matrix();
  const TMatrixD& smearcept = *smearcept_ptr;
  auto true_signal = syst.get_cv_true_signal();
  auto meas = syst.get_measured_events();
  const auto& data_signal = meas.reco_signal_;
  const auto& data_covmat = meas.cov_matrix_;
  std::cout << "N columns: " << data_covmat->GetNcols() << std::endl;

  auto inv_data_covmat = invert_matrix( *data_covmat );
  TDecompChol chol( *inv_data_covmat );
  TMatrixD Q( chol.GetU() );
  TMatrixD R = Q * smearcept;
  TMatrixD R_tr( TMatrixD::EMatrixCreatorsOp1::kTransposed, R );

  TH2D *inv_data_covmat_hist = new TH2D(smearcept);
  TCanvas* c4 = new TCanvas("c4");
  c4->cd();
  gStyle->SetPaintTextFormat("4.2f");
  inv_data_covmat_hist->SetMarkerColor(kRed);
  inv_data_covmat_hist->SetMarkerSize(0.8);
  inv_data_covmat_hist->GetXaxis()->SetTitle("True");
  inv_data_covmat_hist->GetYaxis()->SetTitle("Reco");
  inv_data_covmat_hist->Draw("colz text");
  c4->SaveAs("inv_data_covmat.pdf"); 
  //  c4->Draw(); 
  //  c4->Update(); 
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
    std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv ) );
    auto result = unfolder->unfold( *mcc9 );
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //TH2D *h_error = new TH2D(*result.err_prop_matrix_);

  TH2D *h_error = new TH2D(*result.unfolding_matrix_); //AAAA
  TCanvas* c40 = new TCanvas("c40");
  c40->cd();
  gStyle->SetPaintTextFormat("4.2f");
  h_error->SetMarkerColor(kRed);
  h_error->SetMarkerSize(0.8);
  h_error->GetXaxis()->SetTitle("True");
  h_error->GetYaxis()->SetTitle("Reco");
  h_error->Draw("colz text");
  c40->SaveAs("unfolding_mat.pdf");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout << "num_true_bins: " << num_true_bins << std::endl;
  std::cout << "nbins size: " << nbins.size() << std::endl;
  for ( int t = 0; t < num_true_bins; ++t ) { 
    double evts = 0.;
    double error = 0.;
    if ( t < nbins.size() - 1) {
          
        ///////////////////////////////////////////////////
      std::cout<<"background subtracted data [reco_signal]: "<<data_signal->operator()( t, 0 )<<std::endl;
      std::cout<<"CV_true_signal: "<<true_signal->operator()( t, 0 )<<std::endl;
      std::cout<<"Scaled CV_true_signal to BNB POT: "<<true_signal->operator()( t, 0 )*0.259<<std::endl;//WHERE IS THAT FACTOR COMING FROM???
         //TH2D *h = new TH2D(*result.err_prop_matrix_);
         /////////////////////////////////////////////////         

      evts = result.unfolded_signal_->operator()( t, 0 );
      error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );

      std::cout << "evts: " << evts << std::endl;
      unfolded_events->SetBinContent( t + 1, evts );
      unfolded_events->SetBinError( t + 1, error );
      unfolded_events_vals->SetBinContent( t + 1, evts/conv_factor/(nbins[t+1] - nbins[t]) );
      unfolded_events_vals->SetBinError( t + 1, error/conv_factor/(nbins[t+1] - nbins[t]) );
      std::cout<<"ERROR: "<<std::pow(error/conv_factor/(nbins[t+1] - nbins[t]),2)<<std::endl;        
 
      //std::cout<<"unfolded signal events: "<<evts/*<<" pre-scaled error: "<<error*/<<std::endl;
      //std::cout<<"Scale: "<<conv_factor/(nbins[t+1] - nbins[t])<<" SCALED events: "<<evts/conv_factor/(nbins[t+1] - nbins[t])<<" Scaled error:  "<<error/conv_factor/(nbins[t+1] - nbins[t])<<"\n"<<std::endl;

    }
  }

  TMatrixD *A_C = result.add_smear_matrix_.get();
  //TH2D* ac_matrix_bins = new TH2D(*A_C);  

  auto true_covmat_ptr = result.cov_matrix_.get();
  const TMatrixD& true_covmat = *true_covmat_ptr; 
  TH2D* ac_matrix_bins = new TH2D(true_covmat);
  
  TH2D* ac_matrix_hist = new TH2D("ac_matrix_hist", "; cos#theta^{true}_{#mu}; cos#theta^{true}_{#mu}", nbins.size() - 1, nbins.data(), nbins.size() - 1, nbins.data());
  for ( int iter_row = 0; iter_row < nbins.size() - 1; ++iter_row ) {
      for ( int iter_col = 0; iter_col < nbins.size() - 1; ++iter_col ) {
          ac_matrix_hist->SetBinContent(iter_row + 1, iter_col + 1, ac_matrix_bins->GetBinContent(iter_row + 1, iter_col + 1)/conv_factor/conv_factor/(nbins[ iter_row + 1 ] -  nbins[iter_row])/(nbins[iter_col + 1] -  nbins[iter_col]) );
      }
  } 

   ///////////////
  TH2D* total_correlation_matrix_hist = new TH2D( "total_correlation_matrix_hist", "; p^{true}_{#mu}; p^{true}_{#mu}", nbins.size() - 1, nbins.data(), nbins.size() - 1, nbins.data());
   std::cout << "Size: " << nbins.size() << std::endl;

  for ( int iter_row = 0; iter_row < nbins.size() - 1; ++iter_row ) {
      for ( int iter_col = 0; iter_col < nbins.size() - 1; ++iter_col ) {
         total_correlation_matrix_hist->SetBinContent(iter_row + 1, iter_col + 1, ac_matrix_hist->GetBinContent(iter_row + 1, iter_col + 1)/std::sqrt(ac_matrix_hist->GetBinContent(iter_row + 1, iter_row + 1))/std::sqrt(ac_matrix_hist->GetBinContent(iter_col + 1, iter_col + 1)));
      }
  }

   ///////////////

  TCanvas* c3 = new TCanvas;
  gStyle->SetPalette(kViridis);
  //ac_matrix_hist->Draw("colz");
  //ac_matrix_hist->GetZaxis()->SetRangeUser(-10000.,715000.);
  gStyle->SetPaintTextFormat("4.2f");
  ac_matrix_hist->SetMarkerColor(kRed);
  ac_matrix_hist->SetMarkerSize(0.2);
  ac_matrix_hist->Draw("colz text");
  c3->SaveAs("ac_matrix_hist.pdf"); 

  TCanvas* c7 = new TCanvas;
  gStyle->SetPalette(kViridis);
  total_correlation_matrix_hist->Draw("colz");
  total_correlation_matrix_hist->GetZaxis()->SetRangeUser(-1.,1.);
  c7->SaveAs("total_correlation_matrix_hist.pdf");


  multiply_1d_hist_by_matrix(A_C, genie_cv_truth);
  multiply_1d_hist_by_matrix(A_C, fake_data_truth);

  for ( int t = 0; t < nbins.size() - 1; ++t ) { 
      fake_data_truth_vals->SetBinContent( t + 1, fake_data_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
      genie_cv_truth_vals->SetBinContent( t + 1, genie_cv_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
  }

 // for ( int t = 0; t < result.cov_matrix_->GetNcols() ; ++t ) {
 //     for ( int u = 0; u < result.cov_matrix_->GetNrows() ; ++u ) {
  std::cout << "number of columns: " << result.cov_matrix_->GetNcols() << std::endl;
  for ( int t = 0; t < nbins.size() - 1; ++t ) { 
     for ( int u = 0; u < nbins.size() - 1 ; ++u ) {
         result.cov_matrix_->operator()(u,t) = result.cov_matrix_->operator()(u,t)/conv_factor/conv_factor/(nbins[u+1] -  nbins[u])/(nbins[t+1] -  nbins[t]);
      }
  }

  auto inv_cov_mat = invert_matrix(*result.cov_matrix_, 1e-4 );
   //  auto inv_cov_mat = invert_matrix(*result.cov_matrix_);

  std::vector<std::string> filenames;
  std::ifstream myfile(infiles);
  std::copy(std::istream_iterator<std::string>(myfile),
            std::istream_iterator<std::string>(),
            std::back_inserter(filenames));
  std::cout << "File Count: " << filenames.size() << std::endl;

  TCanvas* c2 = new TCanvas("c2");
  TLegend *leg=new TLegend(0.6,0.5,0.85,0.85);  //0.88
  //TLegend *leg=new TLegend(0.65,0.65,0.88,0.88);

  std::cout<<"Now examining generator true data"<<std::endl; 
 
  for (int i = 0; i < filenames.size(); i++) {

  std::cout << "File: " << filenames.at(i) << std::endl;
  TFile* file = TFile::Open((filenames.at(i)).c_str());
  TTree *tree = (TTree*)file->Get("phaseIITriggerTree");
	//insert variables here
  double fSF = 4.27095e-40;
  int rnTrig;
  int rnMRD;
  int evNTrig;
  int evNMRD;
  int trigNum;
  int trigword, hasTank, hasMRD, tankMRDCoinc, noveto;
  double nuE, nuvtxx, nuvtxy, nuvtxz, nupx, nupy, nupz;
  double fslpx, fslpy, fslpz, mcfslE;
  double fslvtxx, fslvtxy, fslvtxz;
  int isCC, isQEL, isRES, isDIS, isCOH, isMEC, mcfslpdg;
  double mcvtxx, mcvtxy, mcvtxz, mcdirx, mcdiry, mcdirz, mcangle, mcmuonE, mctanktracklength, mcmrdtracklength;
  int hasPi0, hasPiP, hasPiM, hasPiPC, hasPiMC, hasKP, hasKM, hasKPC, hasKMC;
  int numMRDTracks;
  int mcentersmrd, mcexitsmrd, mcpenetratesmrd;
  int simpleflag, simplefv;
  double simpleenergy, simplecostheta, simplept, simplemrdenergy, simplemrdtrack, simpletanktrack;
  double simplevtxx, simplevtxy, simplevtxz, simplestopvtxx, simplestopvtxy, simplestopvtxz;
  double simplemrdstartx, simplemrdstarty, simplemrdstartz, simplemrdstopx, simplemrdstopy, simplemrdstopz;
  std::vector<double>* MRDTrackAngle = new std::vector<double>();
  std::vector<double>* MRDTrackAngleError = new std::vector<double>();
  std::vector<double>* MRDPenetrationDepth = new std::vector<double>();
  std::vector<double>* MRDTrackLength = new std::vector<double>();
  std::vector<double>* MRDEntryPointRadius = new std::vector<double>();
  std::vector<double>* MRDEnergyLoss = new std::vector<double>();
  std::vector<double>* MRDEnergyLossError = new std::vector<double>();
  std::vector<double>* MRDTrackStartX = new std::vector<double>();
  std::vector<double>* MRDTrackStartY = new std::vector<double>();
  std::vector<double>* MRDTrackStartZ = new std::vector<double>();
  std::vector<double>* MRDTrackStopX = new std::vector<double>();
  std::vector<double>* MRDTrackStopY = new std::vector<double>();
  std::vector<double>* MRDTrackStopZ = new std::vector<double>();
  std::vector<bool>* MRDSide = new std::vector<bool>();
  std::vector<bool>* MRDStop = new std::vector<bool>();
  std::vector<bool>* MRDThrough = new std::vector<bool>();
 //  std::vector<double>* MRDTrackLengthTrig = new std::vector<double>();
 //  std::vector<double>* MRDTrackLengthMRD = new std::vector<double>();
   std::vector<double>* All_weight = new std::vector<double>();
 //  map<string,vector<double>>* xsecweights = new map<string,vector<double>>();
 //  map<string,vector<double>>* fluxweights = new map<string,vector<double>>();
   //TBranch *bwgts = 0;
   //Set branch addresses
 //  tree->SetBranchAddress("weight_All_UBGenie",&All_weight);
   tree->SetBranchAddress("runNumber",&rnTrig);
   tree->SetBranchAddress("eventNumber",&evNTrig);
   tree->SetBranchAddress("trigword",&trigword);
   tree->SetBranchAddress("HasTank",&hasTank);
   tree->SetBranchAddress("HasMRD",&hasMRD);
   tree->SetBranchAddress("TankMRDCoinc",&tankMRDCoinc);
   tree->SetBranchAddress("NoVeto",&noveto);
   tree->SetBranchAddress("numMRDTracks",&numMRDTracks);
  tree->SetBranchAddress("MRDTrackAngle",&MRDTrackAngle);
  tree->SetBranchAddress("MRDTrackAngleError",&MRDTrackAngleError);
  tree->SetBranchAddress("MRDPenetrationDepth",&MRDPenetrationDepth);
  tree->SetBranchAddress("MRDTrackLength",&MRDTrackLength);
  tree->SetBranchAddress("MRDEntryPointRadius",&MRDEntryPointRadius);
  tree->SetBranchAddress("MRDEnergyLoss",&MRDEnergyLoss);
  tree->SetBranchAddress("MRDEnergyLossError",&MRDEnergyLossError);
  tree->SetBranchAddress("MRDTrackStartX",&MRDTrackStartX);
  tree->SetBranchAddress("MRDTrackStartY",&MRDTrackStartY);
  tree->SetBranchAddress("MRDTrackStartZ",&MRDTrackStartZ);
  tree->SetBranchAddress("MRDTrackStopX",&MRDTrackStopX);
  tree->SetBranchAddress("MRDTrackStopY",&MRDTrackStopY);
  tree->SetBranchAddress("MRDTrackStopZ",&MRDTrackStopZ);
  tree->SetBranchAddress("MRDSide",&MRDSide);
  tree->SetBranchAddress("MRDStop",&MRDStop);
  tree->SetBranchAddress("MRDThrough",&MRDThrough);

  tree->SetBranchAddress("trueCC",&isCC);
  tree->SetBranchAddress("trueQEL",&isQEL);
  tree->SetBranchAddress("trueRES",&isRES);
  tree->SetBranchAddress("trueDIS",&isDIS);
  tree->SetBranchAddress("trueCOH",&isCOH);
  tree->SetBranchAddress("trueMEC",&isMEC);
  tree->SetBranchAddress("trueNuIntxVtx_X",&nuvtxx);
  tree->SetBranchAddress("trueNuIntxVtx_Y",&nuvtxy);
  tree->SetBranchAddress("trueNuIntxVtx_Z",&nuvtxz);
  tree->SetBranchAddress("trueNeutrinoEnergy",&nuE);
  tree->SetBranchAddress("trueNeutrinoMomentum_X",&nupx);
  tree->SetBranchAddress("trueNeutrinoMomentum_Y",&nupy);
  tree->SetBranchAddress("trueNeutrinoMomentum_Z",&nupz);
  tree->SetBranchAddress("truePi0",&hasPi0);
  tree->SetBranchAddress("truePiPlus",&hasPiP);
  tree->SetBranchAddress("truePiMinus",&hasPiM);
  tree->SetBranchAddress("truePiPlusCher",&hasPiPC);
  tree->SetBranchAddress("truePiMinusCher",&hasPiMC);
  tree->SetBranchAddress("trueKPlus",&hasKP);
  tree->SetBranchAddress("trueKMinus",&hasKM);
  tree->SetBranchAddress("trueKPlusCher",&hasKPC);
  tree->SetBranchAddress("trueKMinusCher",&hasKMC);
  tree->SetBranchAddress("trueFSLMomentum_X",&fslpx);
  tree->SetBranchAddress("trueFSLMomentum_Y",&fslpy);
  tree->SetBranchAddress("trueFSLMomentum_Z",&fslpz);
  tree->SetBranchAddress("trueFSLVtx_X",&fslvtxx);
  tree->SetBranchAddress("trueFSLVtx_Y",&fslvtxy);
  tree->SetBranchAddress("trueFSLVtx_Z",&fslvtxz);
  tree->SetBranchAddress("trueFSLEnergy",&mcfslE);
  tree->SetBranchAddress("trueFSLPdg",&mcfslpdg);
  tree->SetBranchAddress("triggerNumber",&trigNum);
  //  tree->SetBranchAddress("XSecWeights",&xsecweights);
  //  tree->SetBranchAddress("FluxWeights",&fluxweights);

  //Simple Reco
  tree->SetBranchAddress("simpleRecoFlag",&simpleflag);
  tree->SetBranchAddress("simpleRecoEnergy",&simpleenergy);
  tree->SetBranchAddress("simpleRecoVtxX",&simplevtxx);
  tree->SetBranchAddress("simpleRecoVtxY",&simplevtxy);
  tree->SetBranchAddress("simpleRecoVtxZ",&simplevtxz);
  tree->SetBranchAddress("simpleRecoStopVtxX",&simplestopvtxx);
  tree->SetBranchAddress("simpleRecoStopVtxY",&simplestopvtxy);
  tree->SetBranchAddress("simpleRecoStopVtxZ",&simplestopvtxz);
  tree->SetBranchAddress("simpleRecoCosTheta",&simplecostheta);
  tree->SetBranchAddress("simpleRecoPt",&simplept);
  tree->SetBranchAddress("simpleRecoFV",&simplefv);
  tree->SetBranchAddress("simpleRecoMrdEnergyLoss",&simplemrdenergy);
  tree->SetBranchAddress("simpleRecoTrackLengthInMRD",&simplemrdtrack);
  tree->SetBranchAddress("simpleRecoTrackLengthInTank",&simpletanktrack);
  tree->SetBranchAddress("simpleRecoMRDStartX",&simplemrdstartx);
  tree->SetBranchAddress("simpleRecoMRDStartY",&simplemrdstarty);
  tree->SetBranchAddress("simpleRecoMRDStartZ",&simplemrdstartz);
  tree->SetBranchAddress("simpleRecoMRDStopX",&simplemrdstopx);
  tree->SetBranchAddress("simpleRecoMRDStopY",&simplemrdstopy);
  tree->SetBranchAddress("simpleRecoMRDStopZ",&simplemrdstopz);

  //true Muon info from WCSim
  tree->SetBranchAddress("trueMuonEnergy",&mcmuonE);
  tree->SetBranchAddress("trueVtxX",&mcvtxx);
  tree->SetBranchAddress("trueVtxY",&mcvtxy);
  tree->SetBranchAddress("trueVtxZ",&mcvtxz);
  tree->SetBranchAddress("trueDirX",&mcdirx);
  tree->SetBranchAddress("trueDirY",&mcdiry);
  tree->SetBranchAddress("trueDirZ",&mcdirz);
  tree->SetBranchAddress("trueTrackLengthInWater",&mctanktracklength);
  tree->SetBranchAddress("trueTrackLengthInMRD",&mcmrdtracklength);
  tree->SetBranchAddress("trueAngle",&mcangle);
  tree->SetBranchAddress("trueEntersMRD",&mcentersmrd);
  tree->SetBranchAddress("trueExitsMRD",&mcexitsmrd);
  tree->SetBranchAddress("truePenetratesMRD",&mcpenetratesmrd);

  TH1D* h_muons_gen = new TH1D("h_muons_gen", "; p_{#mu}; Number of Events", nbins.size() - 1, nbins.data());
  double muon_m = 105.66;

  for (int i_tree=0; i_tree<tree->GetEntries(); i_tree++) {//Loop over the entries.
    tree->GetEntry(i_tree);
    bool inFV = FidVol(nuvtxx,nuvtxy,nuvtxz);
    bool recoFV = FidVol(simplevtxx*100.,simplevtxy*100.,simplevtxz*100.);
    double fslp = std::sqrt(fslpx*fslpx + fslpy*fslpy + fslpz*fslpz);
    bool exceedsCher = (160. < mcmuonE);
    bool ismuon = (mcfslpdg == 13);
    double mcmuonp = std::sqrt(mcmuonE*mcmuonE - muon_m*muon_m);
    double mcmuonpx = mcmuonp*mcdirx;
    double mcmuonpy = mcmuonp*mcdiry;
    double mcmuonpz = mcmuonp*mcdirz;
    double MRDtotE = 0;
    for( int j = 0; j < numMRDTracks; j++ ){
      MRDtotE += MRDEnergyLoss->at(j);
    }

    //Reco cuts
 //    if((numMRDTracks == 1) && recoFV && (MRDStop->at(0)) && (simpleenergy > 600.) && (simpleenergy < 1200.) && (simplecostheta > 0.8)){
     if(inFV && isCC && ismuon && (std::cos(mcangle * M_PI / 180.) > 0.8) && (mcmuonE > 600.) && (mcmuonE < 1200.)){
       h_muons_gen->Fill(simpleenergy);
     }
   }
 
 //    h_muons_gen->FillRandom("gaus",10000);
   multiply_1d_hist_by_matrix(A_C, h_muons_gen);  
   for ( int t = 0; t < nbins.size() - 1; ++t ) { 
     h_muons_gen->SetBinContent( t + 1, h_muons_gen->GetBinContent(t + 1)*16.*fSF*1e38/(nbins[t+1] - nbins[t]));
   }
 
   h_muons_gen->SetStats(0);
   h_muons_gen->SetLineWidth(3);
 
   if (i==0){
     //h_muons_gen->GetYaxis()->SetTitle("d#sigma/dp_{#mu} [10^{-38} cm^{2}/(GeV/c)/Ar]");
     //h_muons_gen->GetYaxis()->SetTitle("#frac{d#sigma}{dcos#theta_{#mu}} [ 10^{-38} #frac{cm^{2}}{Ar} ]");
     h_muons_gen->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{MeV/c O} ]");        
 
 //  h_muons_gen->SetAxisRange(0., 30.,"Y");
     h_muons_gen->SetLineColor( kOrange + 2);
     h_muons_gen->Draw(); //DRAW
   }
   else{
     if (i==1) h_muons_gen->SetLineColor(kCyan - 3); 
     if (i==2) h_muons_gen->SetLineColor(kRed - 4);
     if (i==3) h_muons_gen->SetLineColor(kGreen + 3);
     //if (i==3) h_muons_gen->Draw();
     h_muons_gen->Draw("same");
   } 
   if (i==0) leg->AddEntry(h_muons_gen,"Genie 3.0.6","l");
   if (i==1) leg->AddEntry(h_muons_gen,"NEUT 5.0.4","l");
   if (i==2) leg->AddEntry(h_muons_gen,"Genie 2.12.10","l");
   if (i==3) leg->AddEntry(h_muons_gen,"NuWro 19.02","l");
   
   double chi_square = 0;
   for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
        for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
            chi_square += (unfolded_events_vals->GetBinContent(k+1)-h_muons_gen->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1)-h_muons_gen->GetBinContent(j+1));
            //std::cout<<chi_square<<std::endl;
        }
    }
     
   double dof = inv_cov_mat->GetNrows();
   double p_value = TMath::Prob(chi_square, dof);
   //double sigma = TMath::Sqrt( TMath::ChisquareQuantile( 1-p_value, dof ) );
  double sigma = RooStats::PValueToSignificance(p_value);
  std::cout<<"chi_square is: "<<chi_square<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  leg->AddEntry((TObject*)0, TString::Format("#chi^{2} / d.o.f = %g / %g", chi_square, dof), "");
  //leg->AddEntry((TObject*)0, TString::Format("p = %g || #sigma = %g ", p_value, sigma), ""); 
  leg->AddEntry((TObject*)0, TString::Format("p = %.2f", p_value), "");
 }

  double chi_square_cv = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square_cv += (unfolded_events_vals->GetBinContent(k+1) - genie_cv_truth_vals->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1) - genie_cv_truth_vals->GetBinContent(j+1));
       }
   }
  
  leg->AddEntry(genie_cv_truth_vals,"MicroBooNE Tune","l");
  double dof = inv_cov_mat->GetNrows();
  double p_value = TMath::Prob(chi_square_cv, dof);
  std::cout<<"chi_square is: "<<chi_square_cv<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  double sigma = RooStats::PValueToSignificance(p_value);
  leg->AddEntry(genie_cv_truth_vals, TString::Format("#chi^{2} / d.o.f = %.2f / %g", chi_square_cv, dof, p_value), "");
  leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %.2f",  p_value), "");
  //leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %g || #sigma = %g",  p_value, sigma), "");

  double chi_square_fake = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square_fake += (unfolded_events_vals->GetBinContent(k+1) - fake_data_truth_vals->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1) - fake_data_truth_vals->GetBinContent(j+1));
        }
   }

  leg->AddEntry(fake_data_truth_vals,"Truth (NuWro)","l");
  p_value = TMath::Prob(chi_square_fake, dof);
  sigma = RooStats::PValueToSignificance(p_value);
  std::cout<<"chi_square is: "<<chi_square_fake<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  leg->AddEntry(fake_data_truth_vals, TString::Format("#chi^{2} / d.o.f = %.2f / %g ", chi_square_fake, dof), "");
  //leg->AddEntry(fake_data_truth_vals, TString::Format("p = %g || #sigma = %g", p_value, sigma), "");
  leg->AddEntry(fake_data_truth_vals, TString::Format("p = %.2f", p_value), "");

  unfolded_events_vals->SetStats(0);
  unfolded_events_vals->SetLineWidth(3);
  unfolded_events_vals->SetLineColor(kBlack);
  //unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{GeV/c Ar} ]"); 
  //unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dcos#theta_{#mu}} [ 10^{-38} #frac{cm^{2}}{Ar} ]"); 
  unfolded_events_vals->GetYaxis()->SetTitle("Number of events"); 
  //unfolded_events_vals->GetYaxis()->SetRangeUser(4000.,116000.);
  unfolded_events_vals->Draw("e same"); //DRAW
  //unfolded_events_vals->Draw("e");
  fake_data_truth_vals->SetLineColor( kBlue );
  fake_data_truth_vals->SetLineWidth( 3 );
  fake_data_truth_vals->SetLineStyle( 2 );
  fake_data_truth_vals->Draw( "hist same" ); //DRAW
  genie_cv_truth_vals->SetLineColor( kMagenta - 3 );
  genie_cv_truth_vals->SetLineWidth( 3 );
  genie_cv_truth_vals->SetLineStyle( 2 );
  genie_cv_truth_vals->Draw( "hist same" ); //DRAW
  leg->AddEntry(unfolded_events_vals,"Fake BNB data","l");
  leg->Draw();
  c2->SaveAs("overlay_chi_sqr_CC0pi.pdf");  
}
///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
void Test_plot() {
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);

    char text_title_pdf1[2024];
    char text_title_pdf2[2024];
    char text_title_pdf3[2024];
    char text_title_pdf4[2024];
    char RootName[1024];
  sprintf(text_title_pdf1, "Make_Plots_Annie_unimake_test_3_setTriggerto1_changedPOT.pdf(","" );

  sprintf(text_title_pdf2, "Make_Plots_Annie_unimake_test_3_setTriggerto1_changedPOT.pdf","" );
  sprintf(text_title_pdf3, "Make_Plots_Annie_unimake_test_3_setTriggerto1_changedPOT.pdf)","" );
  sprintf(text_title_pdf4, "Make_Plots_Annie_unimake_test_3_setTriggerto1_changedPOT","" );

std::string text_title_pdf2_string(text_title_pdf2);

    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties.txt" );
   TCanvas* c4 = new TCanvas("c4");
   c4-> Print(text_title_pdf1);
   //  std::vector< double > nbins = {-1.,-0.775,-0.675,-0.575,-0.475,-0.4,-0.325,-0.25,-0.175,-0.1,-0.025,0.025,0.1,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.};
   
   //  std::vector<double> nbins = {0.8,0.95,1.0};
     std::vector<double> nbins = {600.,810.,910.,1005.,1100.,1200};
   
   
   
   
  //ANNIE
  auto* mcc9 = new MCC9SystematicsCalculator(
  //    "/exp/annie/data/users/jminock/stv-analysis/stv-univmake-output-20k-pmu.root",
    "output.root",
    "systcalc.conf" );

  const auto &syst = *mcc9;

  TH1D* genie_cv_truth_vals = new TH1D("genie_cv_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  TH1D* fake_data_truth_vals = new TH1D("fake_data_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  TH1D* unfolded_events_vals = new TH1D("unfolded_events_vals", ";cos#theta_{#mu}; Scaled Events", nbins.size() - 1, nbins.data());
  TH1D* fake_data_reco_vals = new TH1D("fake_data_reco_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());


  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  bool NormArea = false;
   bool Setgrid = true; 
   bool BinWidthNorm = false; 

  Draw_HIST(
   reco_bnb_hist,
  "reco_bnb_hist",
   reco_ext_hist,
   "reco_ext_hist",
   "",
   "KE_{#mu} ",
   "Events",
    NormArea, 
    Setgrid,
    BinWidthNorm,
   -99,
   c4,
   text_title_pdf2_string);
   
   TH1D* reco_bnb_hist_clone = (TH1D*)reco_bnb_hist->Clone("reco_bnb_hist_clone");
   TH1D* reco_ext_hist_clone = (TH1D*)reco_ext_hist->Clone("reco_ext_hist_clone");
   
   reco_bnb_hist_clone->Divide(reco_ext_hist_clone);
   reco_ext_hist_clone->Divide(reco_ext_hist_clone);
   
     Draw_HIST(
   reco_bnb_hist_clone,
  "reco_bnb_hist/reco_ext_hist",
   reco_ext_hist_clone,
   "reco_ext_hist/reco_ext_hist",
   "",
   "KE_{#mu} ",
   "#frac{reco_bnb_hist}{reco_ext_hist}",
    NormArea, 
    Setgrid,
    BinWidthNorm,
   -99,
   c4,
   text_title_pdf2_string);
   

  //syst.data_hists_.at( NFT::kOnBNB ).get()->Add(syst.data_hists_.at( NFT::kExtBNB ).get(),-1);



  std::cout<<"Now examining real data."<<std::endl;

  double total_pot = mcc9->total_bnb_data_pot_;
  std::cout << "total_pot: " << total_pot << std::endl;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  std::cout << "integ_flux: " << integ_flux << std::endl;
  double num_Ar = num_O_targets_in_FV(); //5.25e28;
  std::cout << "num_Ar: " << num_Ar << std::endl;
  double conv_factor = (num_Ar * integ_flux)/1e38;
  std::cout << "conv_factor: " << conv_factor << std::endl;

  TH1D* genie_cv_truth = mcc9->cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();
  TH1D* unfolded_events = dynamic_cast< TH1D* >(genie_cv_truth->Clone("unfolded_events") );
  unfolded_events->Reset(); 
   
   std::cout<<" Total POT ="<< total_pot << std::endl;
  
  const auto& fake_data_univ = mcc9->fake_data_universe();
  TH1D* fake_data_truth = fake_data_univ->hist_true_.get(); 
  TH1D* fake_data_reco = fake_data_univ->hist_reco_.get(); 
  std::cout << "AAAA" << std::endl;

  Draw_STACK_HIST(
   reco_bnb_hist,
  "reco_bnb_hist",
   reco_ext_hist,
   "reco_ext_hist",
   fake_data_reco,
   "fake_data_reco",
   "",
   "E_{#mu} ",
   "Events",
    NormArea, 
    Setgrid,
    BinWidthNorm,
   -99,
   c4,
   text_title_pdf2_string);
 

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto smearcept_ptr = syst.get_cv_smearceptance_matrix();
  const TMatrixD& smearcept = *smearcept_ptr;
  auto true_signal = syst.get_cv_true_signal();
  auto meas = syst.get_measured_events();
  const auto& data_signal = meas.reco_signal_;
  const auto& data_covmat = meas.cov_matrix_;
  std::cout << "N columns: " << data_covmat->GetNcols() << std::endl;

  auto inv_data_covmat = invert_matrix( *data_covmat );
  TDecompChol chol( *inv_data_covmat );
  TMatrixD Q( chol.GetU() );
  TMatrixD R = Q * smearcept;
  TMatrixD R_tr( TMatrixD::EMatrixCreatorsOp1::kTransposed, R );

  TH2D *inv_data_covmat_hist = new TH2D(smearcept);

  c4->cd();
  gStyle->SetPaintTextFormat("4.2f");
  inv_data_covmat_hist->SetMarkerColor(kRed);
  inv_data_covmat_hist->SetMarkerSize(0.8);
  inv_data_covmat_hist->GetXaxis()->SetTitle("True");
  inv_data_covmat_hist->GetYaxis()->SetTitle("Reco");
  inv_data_covmat_hist->Draw("colz text");
  c4->SaveAs(text_title_pdf2); 
  //  c4->Draw(); 
  //  c4->Update(); 
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
    std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv ) );
    
    auto result = unfolder->unfold( *mcc9 );
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //TH2D *h_error = new TH2D(*result.err_prop_matrix_);

  TH2D *h_error = new TH2D(*result.unfolding_matrix_); //AAAA
  TCanvas* c40 = new TCanvas("c40");
  c40->cd();
  gStyle->SetPaintTextFormat("4.2f");
  h_error->SetMarkerColor(kRed);
  h_error->SetMarkerSize(0.8);
  h_error->GetXaxis()->SetTitle("True");
  h_error->GetYaxis()->SetTitle("Reco");
  h_error->Draw("colz text");
  c40->SaveAs(text_title_pdf2);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout << "num_true_bins: " << num_true_bins << std::endl;
  std::cout << "nbins size: " << nbins.size() << std::endl;
  for ( int t = 0; t < num_true_bins; ++t ) { 
    double evts = 0.;
    double error = 0.;
    if ( t < nbins.size() - 1) {
          
        ///////////////////////////////////////////////////
      std::cout<<"background subtracted data [reco_signal]: "<<data_signal->operator()( t, 0 )<<std::endl;
      std::cout<<"CV_true_signal: "<<true_signal->operator()( t, 0 )<<std::endl;
      std::cout<<"Scaled CV_true_signal to BNB POT: "<<true_signal->operator()( t, 0 )*0.259<<std::endl;//WHERE IS THAT FACTOR COMING FROM???
         //TH2D *h = new TH2D(*result.err_prop_matrix_);
         /////////////////////////////////////////////////         

      evts = result.unfolded_signal_->operator()( t, 0 );
      error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );

      std::cout << "evts: " << evts << std::endl;
      unfolded_events->SetBinContent( t + 1, evts );
      unfolded_events->SetBinError( t + 1, error );
      unfolded_events_vals->SetBinContent( t + 1, evts/conv_factor/(nbins[t+1] - nbins[t]) );
      unfolded_events_vals->SetBinError( t + 1, error/conv_factor/(nbins[t+1] - nbins[t]) );
      std::cout<<"ERROR: "<<std::pow(error/conv_factor/(nbins[t+1] - nbins[t]),2)<<std::endl;        
 
      //std::cout<<"unfolded signal events: "<<evts/*<<" pre-scaled error: "<<error*/<<std::endl;
      //std::cout<<"Scale: "<<conv_factor/(nbins[t+1] - nbins[t])<<" SCALED events: "<<evts/conv_factor/(nbins[t+1] - nbins[t])<<" Scaled error:  "<<error/conv_factor/(nbins[t+1] - nbins[t])<<"\n"<<std::endl;

    }
  }

  TMatrixD *A_C = result.add_smear_matrix_.get();
  //TH2D* ac_matrix_bins = new TH2D(*A_C);  

  auto true_covmat_ptr = result.cov_matrix_.get();
  const TMatrixD& true_covmat = *true_covmat_ptr; 
  TH2D* ac_matrix_bins = new TH2D(true_covmat);
  
  TH2D* ac_matrix_hist = new TH2D("ac_matrix_hist", "; cos#theta^{true}_{#mu}; cos#theta^{true}_{#mu}", nbins.size() - 1, nbins.data(), nbins.size() - 1, nbins.data());
  for ( int iter_row = 0; iter_row < nbins.size() - 1; ++iter_row ) {
      for ( int iter_col = 0; iter_col < nbins.size() - 1; ++iter_col ) {
          ac_matrix_hist->SetBinContent(iter_row + 1, iter_col + 1, ac_matrix_bins->GetBinContent(iter_row + 1, iter_col + 1)/conv_factor/conv_factor/(nbins[ iter_row + 1 ] -  nbins[iter_row])/(nbins[iter_col + 1] -  nbins[iter_col]) );
      }
  } 

   ///////////////
  TH2D* total_correlation_matrix_hist = new TH2D( "total_correlation_matrix_hist", "; p^{true}_{#mu}; p^{true}_{#mu}", nbins.size() - 1, nbins.data(), nbins.size() - 1, nbins.data());
   std::cout << "Size: " << nbins.size() << std::endl;

  for ( int iter_row = 0; iter_row < nbins.size() - 1; ++iter_row ) {
      for ( int iter_col = 0; iter_col < nbins.size() - 1; ++iter_col ) {
         total_correlation_matrix_hist->SetBinContent(iter_row + 1, iter_col + 1, ac_matrix_hist->GetBinContent(iter_row + 1, iter_col + 1)/std::sqrt(ac_matrix_hist->GetBinContent(iter_row + 1, iter_row + 1))/std::sqrt(ac_matrix_hist->GetBinContent(iter_col + 1, iter_col + 1)));
      }
  }

   ///////////////

  TCanvas* c3 = new TCanvas;
  gStyle->SetPalette(kViridis);
  //ac_matrix_hist->Draw("colz");
  //ac_matrix_hist->GetZaxis()->SetRangeUser(-10000.,715000.);
  gStyle->SetPaintTextFormat("4.2f");
  ac_matrix_hist->SetMarkerColor(kRed);
  ac_matrix_hist->SetMarkerSize(0.2);
  ac_matrix_hist->Draw("colz text");
  c3->SaveAs(text_title_pdf2); 

  TCanvas* c7 = new TCanvas;
  gStyle->SetPalette(kViridis);
  total_correlation_matrix_hist->Draw("colz");
  total_correlation_matrix_hist->GetZaxis()->SetRangeUser(-1.,1.);
  c7->SaveAs(text_title_pdf2);


  multiply_1d_hist_by_matrix(A_C, genie_cv_truth);
  multiply_1d_hist_by_matrix(A_C, fake_data_truth);

  for ( int t = 0; t < nbins.size() - 1; ++t ) { 
      fake_data_truth_vals->SetBinContent( t + 1, fake_data_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
      genie_cv_truth_vals->SetBinContent( t + 1, genie_cv_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
  }

 // for ( int t = 0; t < result.cov_matrix_->GetNcols() ; ++t ) {
 //     for ( int u = 0; u < result.cov_matrix_->GetNrows() ; ++u ) {
  std::cout << "number of columns: " << result.cov_matrix_->GetNcols() << std::endl;
  for ( int t = 0; t < nbins.size() - 1; ++t ) { 
     for ( int u = 0; u < nbins.size() - 1 ; ++u ) {
         result.cov_matrix_->operator()(u,t) = result.cov_matrix_->operator()(u,t)/conv_factor/conv_factor/(nbins[u+1] -  nbins[u])/(nbins[t+1] -  nbins[t]);
      }
  }

  auto inv_cov_mat = invert_matrix(*result.cov_matrix_, 1e-4 );
   //  auto inv_cov_mat = invert_matrix(*result.cov_matrix_);


  TCanvas* c2 = new TCanvas("c2");
  TLegend *leg=new TLegend(0.6,0.5,0.85,0.85);  //0.88
  //TLegend *leg=new TLegend(0.65,0.65,0.88,0.88);

  std::cout<<"Now examining generator true data"<<std::endl; 


  double chi_square_cv = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square_cv += (unfolded_events_vals->GetBinContent(k+1) - genie_cv_truth_vals->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1) - genie_cv_truth_vals->GetBinContent(j+1));
       }
   }
  
  leg->AddEntry(genie_cv_truth_vals,"MicroBooNE Tune","l");
  double dof = inv_cov_mat->GetNrows();
  double p_value = TMath::Prob(chi_square_cv, dof);
  std::cout<<"chi_square is: "<<chi_square_cv<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  double sigma = RooStats::PValueToSignificance(p_value);
  leg->AddEntry(genie_cv_truth_vals, TString::Format("#chi^{2} / d.o.f = %.2f / %g", chi_square_cv, dof, p_value), "");
  leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %.2f",  p_value), "");
  //leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %g || #sigma = %g",  p_value, sigma), "");

  double chi_square_fake = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square_fake += (unfolded_events_vals->GetBinContent(k+1) - fake_data_truth_vals->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1) - fake_data_truth_vals->GetBinContent(j+1));
        }
   }

  leg->AddEntry(fake_data_truth_vals,"Truth (NuWro)","l");
  p_value = TMath::Prob(chi_square_fake, dof);
  sigma = RooStats::PValueToSignificance(p_value);
  std::cout<<"chi_square is: "<<chi_square_fake<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  leg->AddEntry(fake_data_truth_vals, TString::Format("#chi^{2} / d.o.f = %.2f / %g ", chi_square_fake, dof), "");
  //leg->AddEntry(fake_data_truth_vals, TString::Format("p = %g || #sigma = %g", p_value, sigma), "");
  leg->AddEntry(fake_data_truth_vals, TString::Format("p = %.2f", p_value), "");

  unfolded_events_vals->SetStats(0);
  unfolded_events_vals->SetLineWidth(3);
  unfolded_events_vals->SetLineColor(kBlack);
  //unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{GeV/c Ar} ]"); 
  //unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dcos#theta_{#mu}} [ 10^{-38} #frac{cm^{2}}{Ar} ]"); 
  auto max = genie_cv_truth_vals->GetMaximum()*1.8; 
  unfolded_events_vals->GetYaxis()->SetTitle("Number of events"); 
  //unfolded_events_vals->GetYaxis()->SetRangeUser(4000.,116000.);
  unfolded_events_vals->SetMaximum(max);
  unfolded_events_vals->Draw("e"); //DRAW
  
  
     // Calculate the errors
    TH1D *h_upper = (TH1D*)unfolded_events_vals->Clone("h_upper");
    TH1D *h_lower = (TH1D*)unfolded_events_vals->Clone("h_lower");

    for (int i = 1; i <= unfolded_events_vals->GetNbinsX(); ++i) {
        double value = unfolded_events_vals->GetBinContent(i);
        double error = unfolded_events_vals->GetBinError(i);
        
        h_upper->SetBinContent(i, value + error);
        h_lower->SetBinContent(i, value - error);
    }
  
    h_upper->SetLineColorAlpha(kRed, 0.55);
    h_upper->SetLineStyle(9); // Dashed line for upper bound
    h_upper->Draw("HIST SAME"); // Draw the upper bound histogram
    h_upper->SetLineWidth( 2 );
    h_lower->SetLineColorAlpha(kRed, 0.55);
    h_lower->SetLineStyle(9); // Dashed line for lower bound
    h_lower->Draw("HIST SAME"); // Draw the lower bound histogram
    h_lower->SetLineWidth( 2 );
     leg->AddEntry(h_upper, "Error Bounds", "l");
    //leg->AddEntry(h_lower, "Lower Bound", "l");
  
  //unfolded_events_vals->Draw("e");
  fake_data_truth_vals->SetLineColor( kBlue );
  fake_data_truth_vals->SetLineWidth( 3 );
  fake_data_truth_vals->SetLineStyle( 2 );
  fake_data_truth_vals->Draw( "hist same" ); //DRAW
  genie_cv_truth_vals->SetLineColor( kMagenta - 3 );
  genie_cv_truth_vals->SetLineWidth( 3 );
  genie_cv_truth_vals->SetLineStyle( 2 );
  genie_cv_truth_vals->Draw( "hist same" ); //DRAW
  
  leg->AddEntry(unfolded_events_vals,"Fake BNB data","l");
leg->Draw();
  c2->SaveAs(text_title_pdf2);  



  Draw_HIST(
   fake_data_truth_vals,
  "fake_data_truth_vals",
   unfolded_events_vals,
   "unfolded_events_vals",
   "",
   "KE_{#mu} ",
   "Events",
    NormArea, 
    Setgrid,
    BinWidthNorm,
   -99,
   c2,
   text_title_pdf2 );
   
   
   
   
      TH1D* unfolded_events_vals_clone = (TH1D*)unfolded_events_vals->Clone("unfolded_events_vals_clone");
      TH1D* fake_data_truth_vals_clone = (TH1D*)fake_data_truth_vals->Clone("fake_data_truth_vals_clone");


fake_data_truth_vals_clone->Divide(unfolded_events_vals_clone);
unfolded_events_vals_clone->Divide(unfolded_events_vals_clone);
  
  Draw_HIST(
   fake_data_truth_vals_clone,
  "fake_data_truth_vals/unfolded_events_vals",
   unfolded_events_vals_clone,
   "unfolded_events_vals/unfolded_events_vals",
   "",
   "KE_{#mu} ",
   "Ratio",
    NormArea, 
    Setgrid,
    BinWidthNorm,
   -99,
   c2,
   text_title_pdf2 );




c40->SaveAs(text_title_pdf3);


}



int main(int argc, char* argv[]) {
   //chi_square_all_gens(argv[1]);
   Test_plot();
   return 0;
}
