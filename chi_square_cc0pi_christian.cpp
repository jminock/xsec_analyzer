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
//#include "DAgostiniUnfolder.hh"
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
  sprintf(text_title_pdf1, "Make_Plots_Annie_unimake_0pi_full.pdf(","" );
  sprintf(text_title_pdf2, "Make_Plots_Annie_unimake_0pi_full.pdf","" );
  sprintf(text_title_pdf3, "Make_Plots_Annie_unimake_0pi_full.pdf)","" );
  sprintf(text_title_pdf4, "Make_Plots_Annie_unimake_0pi_full","" );

std::string text_title_pdf2_string(text_title_pdf2);

    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties.txt" );
   TCanvas* c4 = new TCanvas("c4");
   c4-> Print(text_title_pdf1);
   //  std::vector< double > nbins = {-1.,-0.775,-0.675,-0.575,-0.475,-0.4,-0.325,-0.25,-0.175,-0.1,-0.025,0.025,0.1,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.};
   
   //  std::vector<double> nbins = {0.8,0.95,1.0};
     std::vector<double> nbins = {600.,740.,860.,960.,1080.,1200};
//     std::vector<double> nbins = {600.,710.,860.,960.,1080.,1200.};
//     std::vector<double> nbins = {600.,700.,800.,900.,1000.,1100.,1200.};
//     std::vector<double> nbins = {500.,680.,770.,860.,950.,1040.,1160.,1400.};


   
   
   
  //ANNIE
  auto* mcc9 = new MCC9SystematicsCalculator(
  //    "/exp/annie/data/users/jminock/stv-analysis/stv-univmake-output-20k-pmu.root",
    "output_0pi_full.root",
    "systcalc_data.conf" );

  const auto &syst = *mcc9;

  TH1D* genie_cv_truth_vals = new TH1D("genie_cv_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
//  TH1D* fake_data_truth_vals = new TH1D("fake_data_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  TH1D* unfolded_events_vals = new TH1D("unfolded_events_vals", ";p_{#mu}; Scaled Events", nbins.size() - 1, nbins.data());
//  TH1D* fake_data_reco_vals = new TH1D("fake_data_reco_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());


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
  
/*  const auto& fake_data_univ = mcc9->fake_data_universe();
  TH1D* fake_data_truth = fake_data_univ->hist_true_.get(); 
  TH1D* fake_data_reco = fake_data_univ->hist_reco_.get(); 

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
 */

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
//    std::unique_ptr< Unfolder > unfolder (new DAgostiniUnfolder( DAgostiniUnfolder::ConvergenceCriterion::FigureOfMerit, 0.025 ) );

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
      std::cout<<"Scaled CV_true_signal to BNB POT: "<<true_signal->operator()( t, 0 )*(1.68/2.7)<<std::endl;//WHERE IS THAT FACTOR COMING FROM???
         //TH2D *h = new TH2D(*result.err_prop_matrix_);
         /////////////////////////////////////////////////         

      evts = result.unfolded_signal_->operator()( t, 0 );
      error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );
//      error = std::max(0., result.cov_matrix_->operator()( t, t ));

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
  TH2D* ac_matrix_bins = new TH2D(*A_C);  

  auto true_covmat_ptr = result.cov_matrix_.get();
  const TMatrixD& true_covmat = *true_covmat_ptr; 
  //TH2D* ac_matrix_bins = new TH2D(true_covmat);
  
  TH2D* ac_matrix_hist = new TH2D("ac_matrix_hist", "A_{C} Matrix; p^{true}_{#mu}; p^{true}_{#mu}", nbins.size() - 1, nbins.data(), nbins.size() - 1, nbins.data());
  for ( int iter_row = 0; iter_row < nbins.size() - 1; ++iter_row ) {
      for ( int iter_col = 0; iter_col < nbins.size() - 1; ++iter_col ) {
//          ac_matrix_hist->SetBinContent(iter_row + 1, iter_col + 1, ac_matrix_bins->GetBinContent(iter_row + 1, iter_col + 1)/conv_factor/conv_factor/(nbins[ iter_row + 1 ] -  nbins[iter_row])/(nbins[iter_col + 1] -  nbins[iter_col]) );
//          ac_matrix_hist->SetBinContent(iter_row + 1, iter_col + 1, ac_matrix_bins->GetBinContent(iter_row + 1, iter_col + 1)/(nbins[ iter_row + 1 ] -  nbins[iter_row])/(nbins[iter_col + 1] -  nbins[iter_col]) );
          ac_matrix_hist->SetBinContent(iter_row + 1, iter_col + 1, ac_matrix_bins->GetBinContent(iter_row + 1, iter_col + 1) );
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
  ac_matrix_hist->SetMarkerSize(0.8);
  ac_matrix_hist->Draw("colz text");
  c3->SaveAs(text_title_pdf2); 

  TCanvas* c7 = new TCanvas;
  gStyle->SetPalette(kViridis);
  total_correlation_matrix_hist->Draw("colz");
  total_correlation_matrix_hist->GetZaxis()->SetRangeUser(-1.,1.);
  c7->SaveAs(text_title_pdf2);


  multiply_1d_hist_by_matrix(A_C, genie_cv_truth);
//  multiply_1d_hist_by_matrix(A_C, fake_data_truth);

  for ( int t = 0; t < nbins.size() - 1; ++t ) { 
//      fake_data_truth_vals->SetBinContent( t + 1, fake_data_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
      genie_cv_truth_vals->SetBinContent( t + 1, genie_cv_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
  }

 // for ( int t = 0; t < result.cov_matrix_->GetNcols() ; ++t ) {
 //     for ( int u = 0; u < result.cov_matrix_->GetNrows() ; ++u ) {
  std::cout << "number of columns: " << result.cov_matrix_->GetNcols() << std::endl;
  for ( int t = 0; t < nbins.size() - 1; ++t ) { 
     for ( int u = 0; u < nbins.size() - 1 ; ++u ) {
         result.cov_matrix_->operator()(u,t) = result.cov_matrix_->operator()(u,t)/conv_factor/conv_factor/(nbins[u+1] -  nbins[u])/(nbins[t+1] -  nbins[t]);
         std::cout << result.cov_matrix_->operator()(u,t) << " ";
      }
     std::cout << std::endl;
  }

  auto inv_cov_mat = invert_matrix(*result.cov_matrix_, 1e-4 );
//     auto inv_cov_mat = invert_matrix(*result.cov_matrix_);


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
  
//  leg->AddEntry(genie_cv_truth_vals,"MicroBooNE Tune","l");
  leg->AddEntry(genie_cv_truth_vals,"GENIE CV Truth","l");
  double dof = inv_cov_mat->GetNrows();
  double p_value = TMath::Prob(chi_square_cv, dof);
  std::cout<<"chi_square is: "<<chi_square_cv<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  double sigma = RooStats::PValueToSignificance(p_value);
  leg->AddEntry(genie_cv_truth_vals, TString::Format("#chi^{2} / d.o.f = %.2f / %g", chi_square_cv, dof, p_value), "");
  leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %.2f",  p_value), "");
  //leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %g || #sigma = %g",  p_value, sigma), "");

/*  double chi_square_fake = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square_fake += (unfolded_events_vals->GetBinContent(k+1) - fake_data_truth_vals->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1) - fake_data_truth_vals->GetBinContent(j+1));
        }
   }

//  leg->AddEntry(fake_data_truth_vals,"Truth (NuWro)","l");
  leg->AddEntry(fake_data_truth_vals,"Truth (Fake Data)","l");
  p_value = TMath::Prob(chi_square_fake, dof);
  sigma = RooStats::PValueToSignificance(p_value);
  std::cout<<"chi_square is: "<<chi_square_fake<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  leg->AddEntry(fake_data_truth_vals, TString::Format("#chi^{2} / d.o.f = %.2f / %g ", chi_square_fake, dof), "");
  //leg->AddEntry(fake_data_truth_vals, TString::Format("p = %g || #sigma = %g", p_value, sigma), "");
  leg->AddEntry(fake_data_truth_vals, TString::Format("p = %.2f", p_value), "");
*/
  unfolded_events_vals->SetStats(0);
  unfolded_events_vals->SetLineWidth(3);
  unfolded_events_vals->SetLineColor(kBlack);
  //unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{GeV/c Ar} ]"); 
  //unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dcos#theta_{#mu}} [ 10^{-38} #frac{cm^{2}}{Ar} ]"); 
  auto max = genie_cv_truth_vals->GetMaximum()*1.8; 
//  unfolded_events_vals->GetYaxis()->SetTitle("Number of events"); 
  unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{MeV/c O} ]"); 
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
//  fake_data_truth_vals->SetLineColor( kBlue );
//  fake_data_truth_vals->SetLineWidth( 3 );
//  fake_data_truth_vals->SetLineStyle( 2 );
//  fake_data_truth_vals->Draw( "hist same" ); //DRAW
  genie_cv_truth_vals->SetLineColor( kMagenta - 3 );
  genie_cv_truth_vals->SetLineWidth( 3 );
  genie_cv_truth_vals->SetLineStyle( 2 );
  genie_cv_truth_vals->Draw( "hist same" ); //DRAW
  
//  leg->AddEntry(unfolded_events_vals,"Fake BNB data","l");
  leg->AddEntry(unfolded_events_vals,"Unfolded BNB data","l");
leg->Draw();
  c2->SetLeftMargin(0.2);
  c2->SaveAs(text_title_pdf2);  


/*
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
  */ 
   
   
   
      TH1D* unfolded_events_vals_clone = (TH1D*)unfolded_events_vals->Clone("unfolded_events_vals_clone");
//      TH1D* fake_data_truth_vals_clone = (TH1D*)fake_data_truth_vals->Clone("fake_data_truth_vals_clone");


//fake_data_truth_vals_clone->Divide(unfolded_events_vals_clone);
unfolded_events_vals_clone->Divide(unfolded_events_vals_clone);
/*  
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
*/



c40->SaveAs(text_title_pdf3);


}



int main(int argc, char* argv[]) {
   //chi_square_all_gens(argv[1]);
   Test_plot();
   return 0;
}
