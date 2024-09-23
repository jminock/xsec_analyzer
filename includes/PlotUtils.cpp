///////////////////////
// Centerized Place for Plotting function 
// Christian Nguyen 12/18/2023
//
////////////////////
#include "PlotUtils.hh"

bool DeBug = true;
// Helper function that produces the standard MicroBooNE plot legend title
// with the BNB POT displayed
std::string get_legend_title( double bnb_pot ) {
  // Print the data POT exposure to a precision of 3 decimal digits, then
  // separate out the base-ten exponent for the legend header
  std::stringstream temp_ss;
  temp_ss << std::scientific << std::setprecision(3) << bnb_pot;

  std::string pot_str = temp_ss.str();
  size_t e_pos = pot_str.find( 'e' );

  // Digits not including 'e' and the base-ten exponent
  std::string pot_digits_str = pot_str.substr( 0, e_pos );
  // pot_str now contains just the base-ten exponent
  pot_str.erase( 0, e_pos + 1u );
  // If there's a leading '+' in the exponent, erase that too
  if ( pot_str.front() == '+' ) pot_str.erase( 0, 1u );

  std::string legend_title = " MicroBooNE " + pot_digits_str + " #times 10^{"
    + pot_str + "} POT, Preliminary";

  return legend_title;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
// Return a unique string on every call, so we can name root hists
// without it clobbering them
TString uniq()
{
  static int i=0;
  return TString::Format("uniq%d", i++);
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

void dump_1d_histogram( const std::string& hist_col_prefix,
  const TH1D& hist,
  std::map< std::string, std::vector<std::string> >& pgf_plots_hist_table,
  bool include_yerror, bool include_x_coords  )
{
  // TODO: consider adding a check for pre-existing duplicate columns
  // (to avoid accidentally overwriting information)

  const std::string bin_col_name = "bin";
  if ( include_x_coords ) {
    pgf_plots_hist_table[ bin_col_name ] = std::vector< std::string >();
  }

  std::string y_col_name = hist_col_prefix;
  pgf_plots_hist_table[ y_col_name ] = std::vector< std::string >();

  std::string yerror_col_name;
  if ( include_yerror ) {
    yerror_col_name = hist_col_prefix + "_error";
    pgf_plots_hist_table[ yerror_col_name ] = std::vector< std::string >();
  }

  // Include the ordinary bins and also the overflow bin. The latter
  // (presumably empty) will help certain pgfplots plotting styles look good
  int overflow_bin_idx = hist.GetNbinsX() + 1;
  for ( int b = 1; b <= overflow_bin_idx; ++b ) {
    double y = hist.GetBinContent( b );

    pgf_plots_hist_table.at( y_col_name )
      .push_back( std::to_string(y) );

    if ( include_yerror ) {
      double yerror = hist.GetBinError( b );
      pgf_plots_hist_table.at( yerror_col_name )
        .push_back( std::to_string(yerror) );
    }

    if ( include_x_coords ) {
      pgf_plots_hist_table.at( bin_col_name )
        .push_back( std::to_string(b) );

      double low_edge = hist.GetBinLowEdge( b );
      double half_width = 0.5 * hist.GetBinWidth( b );

      std::string x_col_name = hist_col_prefix;
      if ( !hist_col_prefix.empty() ) x_col_name += '_';
      x_col_name += 'x';
      auto end = pgf_plots_hist_table.end();
      auto iter = pgf_plots_hist_table.find( x_col_name );
      if ( iter == end ) {
        pgf_plots_hist_table[ x_col_name ]
          = std::vector< std::string > { std::to_string(low_edge) };
      }
      else iter->second.push_back( std::to_string(low_edge) );

      std::string x_hw_col_name = x_col_name + "_halfwidth";
      auto end2 = pgf_plots_hist_table.end();
      auto iter2 = pgf_plots_hist_table.find( x_hw_col_name );
      if ( iter2 == end2 ) {
        pgf_plots_hist_table[ x_hw_col_name ]
          = std::vector<std::string> { std::to_string(half_width) };
      }
      else iter2->second.push_back( std::to_string(half_width) );
    } // include x coordinates
  } // bin number
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void dump_tgraph( const std::string& col_prefix, const TGraph& graph,
  std::map< std::string, std::vector<std::string> >& pgf_plots_table )
{
  // TODO: consider adding a check for pre-existing duplicate columns
  // (to avoid accidentally overwriting information)

  std::string x_col_name = col_prefix + "_x";
  pgf_plots_table[ x_col_name ] = std::vector< std::string >();

  std::string y_col_name = col_prefix + "_y";
  pgf_plots_table[ y_col_name ] = std::vector< std::string >();

  int num_points = graph.GetN();
  for ( int p = 0; p < num_points; ++p ) {
    double x = graph.GetX()[ p ];
    double y = graph.GetY()[ p ];

    pgf_plots_table.at( x_col_name ).push_back( std::to_string(x) );
    pgf_plots_table.at( y_col_name ).push_back( std::to_string(y) );
  } // point index
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void write_pgfplots_file( const std::string& out_filename,
  std::map< std::string, std::vector<std::string> >& pgfplots_table )
{
  std::ofstream out_file( out_filename );

  // Column headings
  for ( const auto& pair : pgfplots_table ) {
    const std::string& col_name = pair.first;
    out_file << "  " << col_name;
  }
  out_file << '\n';

  // Data rows (all columns are assumed to have the same number of rows)
  size_t num_rows = pgfplots_table.cbegin()->second.size();
  for ( size_t r = 0u; r < num_rows; ++r ) {
    for ( const auto& pair : pgfplots_table ) {
      out_file << "  " << pair.second.at( r );
    }
    out_file << '\n';
  }
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void ScaleHistogramsInStack(THStack* stack, double scaleFactor, char *option ) {
    TList* list = stack->GetHists();
    
    for (int i = 0; i < list->GetSize(); i++) {
        TH1* hist = static_cast<TH1*>(list->At(i));
        if (hist) {
            hist->Scale(scaleFactor, option);
        }
    }
}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void ScaleHistogramsInStack(THStack* stack, double scaleFactor ) {
    TList* list = stack->GetHists();
    
    for (int i = 0; i < list->GetSize(); i++) {
        TH1* hist = static_cast<TH1*>(list->At(i));
        if (hist) {
            hist->Scale(scaleFactor);
        }
    }
}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void colNormalize(TH2& hist, bool includeFlows)
{
  int first_bin = includeFlows ? 0 : 1;
  int last_bin = includeFlows ? hist.GetNbinsX()+1 : hist.GetNbinsX();
  //Int_t nbins = includeFlows ? hist.GetNbinsX()+2 : hist.GetNbinsX();
  //const auto nXBins = nbins;
  //const auto nYBins = nbins;
  for(int col = first_bin; col <= last_bin; ++col)
  {
    const double denom = hist.Integral(col, col, 0, last_bin);
    for(int row = first_bin; row <= last_bin; ++row)
    {
    
      double input =   hist.GetBinContent(col, row)/denom;
      bool isnan_bool = std::isnan(input);
      if(isnan_bool==false || input <  .001)  hist.SetBinContent(col, row,input);
      else hist.SetBinContent(col, row, 0);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void rowNormalize(TH2& hist, bool includeFlows)
{

 int first_bin = includeFlows ? 0 : 1;
  int last_bin = includeFlows ? hist.GetNbinsX()+1 : hist.GetNbinsX();

  for(int row = first_bin; row <= last_bin; ++row)
  {
    const double denom = hist.Integral(0, last_bin, row, row);
    for(int col = first_bin; col <= last_bin; ++col)
    {
      double input =   hist.GetBinContent(col, row)/denom;
      bool isnan_bool = std::isnan(input);
      if(isnan_bool==false || input <  .001)  hist.SetBinContent(col, row, input);
      else hist.SetBinContent(col, row, 0);

    }
  }
}//////////////////
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void Draw_heatMap(
  TH2D *h_migration_input, 
  const char* xaxislabel,
  const char* yaxislabel,
  const char* Title,
  const char* pdf,
  int rownormtype,
  TCanvas *can,
  bool includeFlows)
{
  TH2D *h_migration = (TH2D*)h_migration_input->Clone("h_migration");

  int first_bin = includeFlows ? 0 : 1;
  int last_bin_x = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
  int last_bin_y = includeFlows ? h_migration->GetNbinsY()+1 : h_migration->GetNbinsY();
  Int_t nbins_x = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();
  Int_t nbins_y = includeFlows ? h_migration->GetNbinsY()+2 : h_migration->GetNbinsY();

  h_migration->GetXaxis()->SetTitle(xaxislabel);
  h_migration->GetXaxis()->CenterTitle();
  h_migration->GetYaxis()->SetTitle(yaxislabel);
  h_migration->GetYaxis()->CenterTitle();
  h_migration->SetTitle(Title);  
  TMatrixD m_migration(nbins_x, nbins_y);
  TH2D tmp(*h_migration);
  tmp.Reset();

  for (int y = first_bin; y <= last_bin_y; ++y){
    for (int x = first_bin; x <= last_bin_x; ++x){

      double NumberperBin =  h_migration->GetBinContent(x,y);
      bool isnan_bool = std::isnan(NumberperBin);

      if(NumberperBin > 0.0 && isnan_bool == false)
      {
        if (includeFlows)
        {
          m_migration[x][y] = NumberperBin; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = NumberperBin; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, NumberperBin);
      }

      else if (isnan_bool == true || NumberperBin < 0.0)
      {
        if (includeFlows){
          m_migration[x][y] = 0.0; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = 0.0; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, 0.0);

      }

    }
  }


  if(rownormtype==0)  colNormalize(tmp,includeFlows);
  else if(rownormtype==1)  rowNormalize(tmp,includeFlows);
  else{std::cout<<"Normtype none"<<std::endl;}



  gStyle->SetPaintTextFormat("2.2f");
  gStyle->SetOptStat(0);
  tmp.SetMarkerColor(kRed);
  tmp.SetMarkerSize(0.52);
  tmp.SetMinimum(.01);
  tmp.DrawCopy("colz text");
  can->Print(pdf);
  
}
////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void Draw_heatMap_notext(
  TH2D *h_migration_input, 
  const char* xaxislabel,
  const char* yaxislabel,
  const char* Title,
  const char* pdf,
  int rownormtype,
  TCanvas *can,
  bool includeFlows)
{
  TH2D *h_migration = (TH2D*)h_migration_input->Clone("h_migration");

  int first_bin = includeFlows ? 0 : 1;
  int last_bin_x = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
  int last_bin_y = includeFlows ? h_migration->GetNbinsY()+1 : h_migration->GetNbinsY();
  Int_t nbins_x = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();
  Int_t nbins_y = includeFlows ? h_migration->GetNbinsY()+2 : h_migration->GetNbinsY();

  h_migration->GetXaxis()->SetTitle(xaxislabel);
  h_migration->GetXaxis()->CenterTitle();
  h_migration->GetYaxis()->SetTitle(yaxislabel);
  h_migration->GetYaxis()->CenterTitle();
  h_migration->SetTitle(Title);  
  TMatrixD m_migration(nbins_x, nbins_y);
  TH2D tmp(*h_migration);
  tmp.Reset();

  for (int y = first_bin; y <= last_bin_y; ++y){
    for (int x = first_bin; x <= last_bin_x; ++x){

      double NumberperBin =  h_migration->GetBinContent(x,y);
      bool isnan_bool = std::isnan(NumberperBin);

      if(NumberperBin > 0.0 && isnan_bool == false)
      {
        if (includeFlows)
        {
          m_migration[x][y] = NumberperBin; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = NumberperBin; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, NumberperBin);
      }

      else if (isnan_bool == true || NumberperBin < 0.0)
      {
        if (includeFlows){
          m_migration[x][y] = 0.0; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = 0.0; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, 0.0);

      }

    }
  }


  if(rownormtype==0)  colNormalize(tmp,includeFlows);
  else if(rownormtype==1)  rowNormalize(tmp,includeFlows);
  else{std::cout<<"Normtype none"<<std::endl;}



  //gStyle->SetPaintTextFormat("2.2f");
  gStyle->SetOptStat(0);
  tmp.SetMarkerColor(kRed);
  tmp.SetMarkerSize(0.52);
  tmp.SetMinimum(.01);
  tmp.DrawCopy("colz");
  can->Print(pdf);
  
}
////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
/*
void Draw_heatMap_4DMigration(
  TH2D *h_migration_input, 
  const char* xaxislabel,
  const char* yaxislabel,
  const char* Title,
  const char* pdf,
  int rownormtype,
  TCanvas *can,
  bool includeFlows,
  bool IsInclusive,
  std::map<double,std::vector<double>> Bin_edges, 
  double X1_window, double Y1_window, double X2_window, double Y2_window)
{
  TH2D *h_migration = (TH2D*)h_migration_input->Clone("h_migration");


  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();



  int first_bin = includeFlows ? 0 : 1;
  int last_bin_x = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
  int last_bin_y = includeFlows ? h_migration->GetNbinsY()+1 : h_migration->GetNbinsY();
  Int_t nbins_x = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();
  Int_t nbins_y = includeFlows ? h_migration->GetNbinsY()+2 : h_migration->GetNbinsY();

  h_migration->GetXaxis()->SetTitle(xaxislabel);
  h_migration->GetXaxis()->CenterTitle();
  h_migration->GetYaxis()->SetTitle(yaxislabel);
  h_migration->GetYaxis()->CenterTitle();
  h_migration->SetTitle(Title);  
  TMatrixD m_migration(nbins_x, nbins_y);
  TH2D tmp(*h_migration);
  tmp.Reset();

  for (int y = first_bin; y <= last_bin_y; ++y){
    for (int x = first_bin; x <= last_bin_x; ++x){

      double NumberperBin =  h_migration->GetBinContent(x,y);
      bool isnan_bool = std::isnan(NumberperBin);

      if(NumberperBin > 0.0 && isnan_bool == false)
      {
        if (includeFlows)
        {
          m_migration[x][y] = NumberperBin; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = NumberperBin; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, NumberperBin);
      }

      else if (isnan_bool == true || NumberperBin < 0.0)
      {
        if (includeFlows){
          m_migration[x][y] = 0.0; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = 0.0; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, 0.0);

      }

    }
  }


  if(rownormtype==0)  colNormalize(tmp,includeFlows);
  else if(rownormtype==1)  rowNormalize(tmp,includeFlows);
  else{std::cout<<"Normtype none"<<std::endl;}



  gStyle->SetPaintTextFormat("2.2f");
  gStyle->SetOptStat(0);
  tmp.SetMarkerColor(kRed);
  tmp.SetMarkerSize(0.52);
  tmp.SetMinimum(.01);
  tmp.DrawCopy("colz text");
  
   TPad *null=new TPad("null", "null", X1_window, Y1_window, X2_window, Y2_window);
   null->SetFillStyle(0);
   null->SetFrameFillStyle(0);
   null->Draw();
   null->cd();

    // Shrink the TPad by 30%
    //Double_t widthShrinkFactor = 0.7;
    //Double_t heightShrinkFactor = 0.7;
    //null->SetPad(0.7, 0.5, 0.7 + widthShrinkFactor * (1 - 0.7), 0.5 + heightShrinkFactor * (0.75 - 0.5));



    // Draw something on the transparent pad
    // For example, draw a line
    null->cd();
    gPad->Update();
    std::map<int , BinMap> TH1Poly_binMap;
    
    if(IsInclusive==true){
        UBTH2Poly* originalPlot = Make2DHist_inclusive_UB(Bin_edges, TH1Poly_binMap, "originalPlot");
        originalPlot->SetTitle("");
        originalPlot->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
        originalPlot->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
        originalPlot->Draw("text");
        DrawBinningInfo(TH1Poly_binMap);
        DrawBinningNum(TH1Poly_binMap);
        }
    else {
        UBTH2Poly* originalPlot = Make2DHist_UB(Bin_edges, TH1Poly_binMap, "originalPlot");
        originalPlot->SetTitle("");
        originalPlot->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
        originalPlot->GetYaxis()->SetTitle("P_{#mu} [GeV/c]");
        originalPlot->Draw("text");
        DrawBinningInfo(TH1Poly_binMap);
        DrawBinningNum(TH1Poly_binMap);
    }
    
    
    pad1->cd();
   can->Print(pdf);
   can->Closed();
}
*/
////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////


void Draw_heatMap_BinCheck(
  TH2D *h_migration_input, 
  const char* xaxislabel,
  const char* yaxislabel,
  const char* Title,
  const char* pdf,
  TCanvas *can,
  bool includeFlows)
{
  TH2D *h_migration = (TH2D*)h_migration_input->Clone("h_migration");

  int first_bin = includeFlows ? 0 : 1;
  int last_bin_x = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
  int last_bin_y = includeFlows ? h_migration->GetNbinsY()+1 : h_migration->GetNbinsY();
  Int_t nbins_x = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();
  Int_t nbins_y = includeFlows ? h_migration->GetNbinsY()+2 : h_migration->GetNbinsY();

  h_migration->GetXaxis()->SetTitle(xaxislabel);
  h_migration->GetXaxis()->CenterTitle();
  h_migration->GetYaxis()->SetTitle(yaxislabel);
  h_migration->GetYaxis()->CenterTitle();
  h_migration->SetTitle(Title);  
  TMatrixD m_migration(nbins_x, nbins_y);
  TH2D tmp(*h_migration);
  tmp.Reset();


  for (int y = first_bin; y <= last_bin_y; ++y){
    for (int x = first_bin; x <= last_bin_x; ++x){

      double NumberperBin =  h_migration->GetBinContent(x,y);
      bool isnan_bool = std::isnan(NumberperBin);

      if(NumberperBin > 0.0 && isnan_bool == false)
      {
        if (includeFlows)
        {
          m_migration[x][y] = NumberperBin; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = NumberperBin; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, 0);
      }

      else if (isnan_bool == true || NumberperBin < 0.0)
      {
        if (includeFlows){
          m_migration[x][y] = 0.0; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = 0.0; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, 0.0);

      }

    }
  }
//superimpose lines at the xbins positions TLine l; c1->Update(); Double_t ymin = c1->GetUymin(); Double_t ymax = c1->GetUymax(); l.SetLineStyle(2); for (Int_t bin=1;nxbins-1;bin++) { l.DrawLine(xbins[bin],ymin,xbins[bin],ymax); } 

  gStyle->SetPaintTextFormat("2.0f");
  gStyle->SetOptStat(0);
  //tmp.SetMarkerColor(kRed);
  //tmp.SetMarkerSize(0.7);
  tmp.DrawCopy(" ");
  
  
   for (Int_t i = 1; i <= h_migration->GetNbinsX(); ++i) {
        for (Int_t j = 1; j <= h_migration->GetNbinsY(); ++j) {
            Int_t globalBin = h_migration->GetBin(i, j);
            Double_t binContent = h_migration->GetBinContent(globalBin);
Double_t     binNumb = -10 + i*10 + j ; 

Double_t xLowEdge = h_migration->GetXaxis()->GetBinLowEdge(i);
    Double_t xUpEdge = h_migration->GetXaxis()->GetBinUpEdge(i);
    Double_t yLowEdge = h_migration->GetYaxis()->GetBinLowEdge(j);
    Double_t yUpEdge = h_migration->GetYaxis()->GetBinUpEdge(j);

    // Print the bin edges
    std::cout << "Bin Edges for Global Bin " << globalBin << ":\n";
    std::cout << "X-axis: [" << xLowEdge << ", " << xUpEdge << "]\n";
    std::cout << "Y-axis: [" << yLowEdge << ", " << yUpEdge << "]\n";


            // Offset values for demonstration
            Double_t offset1 = 0.025;
            Double_t offset2 = -0.025;

            Double_t offset2y = -0.01;
            // Add first number as text with offset
            TText* text1 = new TText(h_migration->GetXaxis()->GetBinCenter(i) + offset1,
                                     h_migration->GetYaxis()->GetBinCenter(j),
                                     Form("%2.0f",binContent));
            text1->SetTextSize(0.0125);
            text1->SetTextAlign(22);  // Centered
            text1->SetTextColor(kRed);
            text1->Draw();

            // Add second number as text with offset
            TText* text2 = new TText(h_migration->GetXaxis()->GetBinCenter(i) + offset2,
                                     h_migration->GetYaxis()->GetBinCenter(j),
                                     Form("%1.0f", binNumb));  // Example: Double the bin content
            text2->SetTextSize(0.0125);
            text2->SetTextAlign(22);  // Centered
            text2->SetTextColor(kBlue);
            text2->Draw();
            
                 TText* textrange = new TText(h_migration->GetXaxis()->GetBinCenter(i),
                                     h_migration->GetYaxis()->GetBinCenter(j)+ offset2,
                                     Form("(%1.2f,%1.2f) X (%1.2f,%1.2f) ",xLowEdge,xUpEdge,yLowEdge,yUpEdge));
            textrange->SetTextSize(0.008);
            textrange->SetTextAlign(22);  // Centered
            textrange->SetTextColor(22);
            textrange->Draw();
             
             
             
        }
    }
  
  
  TLine l; can->Update();
  l.SetLineStyle(2);
    for (Int_t i = 1; i <=  h_migration->GetNbinsX(); ++i) {
        for (Int_t j = 1; j <=  h_migration->GetNbinsY(); ++j) {
            Double_t xLowEdge =  h_migration->GetXaxis()->GetBinLowEdge(i);
            Double_t xUpEdge =  h_migration->GetXaxis()->GetBinUpEdge(i);
            Double_t yLowEdge =  h_migration->GetYaxis()->GetBinLowEdge(j);
            Double_t yUpEdge =  h_migration->GetYaxis()->GetBinUpEdge(j);

            // Create TLine for the bin edges
            TLine* lineXLow = new TLine(xLowEdge, yLowEdge, xLowEdge, yUpEdge);
            TLine* lineXUp = new TLine(xUpEdge, yLowEdge, xUpEdge, yUpEdge);
            TLine* lineYLow = new TLine(xLowEdge, yLowEdge, xUpEdge, yLowEdge);
            TLine* lineYUp = new TLine(xLowEdge, yUpEdge, xUpEdge, yUpEdge);

            // Set line properties (color, style, etc.) if needed
            // lineXLow->SetLineColor(kRed);
            // lineXUp->SetLineColor(kBlue);
            // lineYLow->SetLineColor(kGreen);
            // lineYUp->SetLineColor(kMagenta);

            // Draw TLines on the canvas
            lineXLow->Draw();
            lineXUp->Draw();
            lineYLow->Draw();
            lineYUp->Draw();
        }
    }

  
  
  can->Print(pdf);
  
}//////////////////
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////


/*
void DrawMagration_heatMap_LabelBinNumber_new(
  MnvH2D *h_mig,  const char* xaxislabel,
  const char* yaxislabel, const char* Title,
  const char* pdf, TCanvas *can, MnvPlotter *plotter,
  MnvH2D *h_binning, int bini, int binningtype,
  int rownormtype, double text_size )

{

  TH2D *h_migration = (TH2D*)h_mig->Clone("h_migration");

  //RooUnfoldResponse UnfoldTool;

  int first_bin =  1;
  int last_bin_x =  h_migration->GetNbinsX();
  int last_bin_y =  h_migration->GetNbinsY();

  Int_t nbins_x =  h_migration->GetNbinsX();
  Int_t nbins_y =  h_migration->GetNbinsY();


  h_migration->GetXaxis()->SetTitle(xaxislabel);
  h_migration->GetXaxis()->CenterTitle();
  h_migration->GetYaxis()->SetTitle(yaxislabel);
  h_migration->GetYaxis()->CenterTitle();
  h_migration->SetMinimum(0.001);

  TMatrixD m_migration(nbins_x, nbins_y);
  TMatrixD m_migration2(nbins_x, nbins_y);
  TH2D tmp(*h_migration);
  TH2D tmp2(*h_migration);
  tmp.Reset();
  tmp2.Reset();

  for (int y = first_bin; y <= last_bin_y; ++y){
    for (int x = first_bin; x <= last_bin_x; ++x){

      double bcx = ((TAxis*)h_migration->GetXaxis())->GetBinCenter(x);
      double bcy = ((TAxis*)h_migration->GetYaxis())->GetBinCenter(y);
      //int global_bin =   UnfoldTool.FindBin(h_migration, bcx,  bcy);
      double NumberperBin =  h_migration->GetBinContent(x,y);
      int global_bin = h_migration->GetBin(x, y);
      tmp.SetBinContent(x, y, global_bin);
      bool isnan_bool = std::isnan(NumberperBin);
      if(isnan_bool == false)tmp2.SetBinContent(x, y, NumberperBin);
      else {tmp2.SetBinContent(x, y, 0.0);}

    }
  }



  tmp.SetMinimum(0.001);
  gStyle->SetTickLength(0.0,"x");
  gStyle->SetTickLength(0.0,"y");
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetPalette(kBird);
  tmp.SetMarkerSize(.42);
  tmp.SetMarkerColor(kBlack);
  tmp.GetXaxis()->SetTickLength(0.);
  tmp.GetYaxis()->SetTickLength(0.);
  tmp.SetBarOffset(0.05);

  tmp2.SetMarkerSize(.42);
  tmp2.SetMarkerColor(kBlack);
  tmp2.GetXaxis()->SetTickLength(0.);
  tmp2.GetYaxis()->SetTickLength(0.);
  tmp2.SetBarOffset(-0.05);
  tmp2.GetZaxis()->SetLabelSize(0.01);

  if(rownormtype==2)  colNormalize(tmp2);
  else if(rownormtype==3)  rowNormalize(tmp2);
  else{std::cout<<"Normtype none"<<std::endl;}

  tmp2.DrawCopy("colz text");
  //tmp.DrawCopy("text Same");
  //drawBinRange(h_binning, binningtype, bini, xaxislabel, ".2f", true);
  drawBinRange(h_binning, binningtype, bini, xaxislabel, text_size, ".2f", true);
  //h_mig->Draw("COLZ");
  //plotter->AddHistoTitle(Title, .04);
  //can->Print(pdf);
}//end of function
*/

void Draw_heatMap(
  TH2Poly *h_migration_input, 
  const char* xaxislabel,
  const char* yaxislabel,
  const char* Title,
  const char* pdf,
  int rownormtype,
  TCanvas *can,
  bool includeFlows)
{
 TH2Poly *h_migration = (TH2Poly*)h_migration_input->Clone("h_migration");

  int first_bin = includeFlows ? 0 : 1;
  int last_bin_x = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
  int last_bin_y = includeFlows ? h_migration->GetNbinsY()+1 : h_migration->GetNbinsY();
  Int_t nbins_x = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();
  Int_t nbins_y = includeFlows ? h_migration->GetNbinsY()+2 : h_migration->GetNbinsY();

  h_migration->GetXaxis()->SetTitle(xaxislabel);
  h_migration->GetXaxis()->CenterTitle();
  h_migration->GetYaxis()->SetTitle(yaxislabel);
  h_migration->GetYaxis()->CenterTitle();
  h_migration->SetTitle(Title);  
  //TMatrixD m_migration(nbins_x, nbins_y);
 // TH2D tmp(*h_migration);
  //tmp.Reset();
/*
  for (int y = first_bin; y <= last_bin_y; ++y){
    for (int x = first_bin; x <= last_bin_x; ++x){

      double NumberperBin =  h_migration->GetBinContent(x,y);
      bool isnan_bool = std::isnan(NumberperBin);

      if(NumberperBin > 0.0 && isnan_bool == false)
      {
        if (includeFlows)
        {
          m_migration[x][y] = NumberperBin; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = NumberperBin; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, NumberperBin);
      }

      else if (isnan_bool == true || NumberperBin < 0.0)
      {
        if (includeFlows){
          m_migration[x][y] = 0.0; //yeah that's right  y/x
        }
        else
        {
          m_migration[x-1][y-1] = 0.0; //yeah that's right  y/x
        }

        tmp.SetBinContent(x, y, 0.0);

      }

    }
  }
*/

  //if(rownormtype==0)  colNormalize(tmp,includeFlows);
  //else if(rownormtype==1)  rowNormalize(tmp,includeFlows);
  //else{std::cout<<"Normtype none"<<std::endl;}



  gStyle->SetPaintTextFormat("2.2f");
  gStyle->SetOptStat(0);
  //tmp.SetMarkerColor(kRed);
  //tmp.SetMarkerSize(0.7);
  //tmp.DrawCopy("colz text");
  
  h_migration->Draw("");
  can->Print(pdf);
  
}//////////////////
//////////////////////////////
void Draw_HIST_Resolution(
  TH1D *hist_1_input,
   char *histotitle,
   std::string xaxislabel,
   std::string yaxislabel,
   std::string pdf_name,
   bool NormArea, 
   bool Setgrid,
   bool BinWidthNorm,
   double Ymax,
   TCanvas *cE)
{

  char textplace[1024];
  TH1D *hist_1 = (TH1D *)hist_1_input->Clone("");
  

  cE->SetGrid();

  //histHelium->GetXaxis()->SetTitleSize(0.035);
  hist_1->GetXaxis()->SetTitle(Form("%s",xaxislabel.c_str()));
  hist_1->GetYaxis()->SetTitle(Form("%s",yaxislabel.c_str()));
  hist_1->SetTitle(histotitle);
  hist_1->SetTitleOffset(1.2);

   if(NormArea==true){
   hist_1->Scale((1.0 / (hist_1->Integral())));
   }


  if(BinWidthNorm== true){
    hist_1->Scale(1.0,"width");
  }

  if(Ymax != -99){
    hist_1->SetMaximum(Ymax);
  }  
   else{
        double max1 = hist_1->GetMaximum();
         hist_1->SetMaximum(max1* 1.25);
       }
       
  int nbins_1 = hist_1->GetNbinsX();
  double Mean = hist_1->GetMean(1);
  double stdDev = hist_1->GetStdDev(1);
  
  hist_1->GetXaxis()->SetTitleSize(0.03);
  hist_1->SetLineColor(kAzure - 1);
  hist_1->Draw("hist e");
  TLatex* text = new TLatex;
  text->SetNDC();
  text->SetTextSize(0.03);
  text->SetTextColor(kRed);
  sprintf(textplace, "Mean = %.2f",Mean );
  text->DrawLatex(0.15, 0.85, textplace);
  sprintf(textplace, "stdDev = %.2f",stdDev );
  text->DrawLatex(0.15, 0.8, textplace);
  
  std::string plotname = Form("%s",pdf_name.c_str());
 cE-> Print(plotname.c_str());

}//////////////////
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void DrawStackMCandData(
    TH1* h_data_input, 
    TH1* h_MC_Total_input,
    THStack* Stack_input,
    char *Yaxis_title,
    char *Xaxis_title,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax,
    TCanvas *Canvas)
{

 TH1D* h_Data =(TH1D*)h_data_input->Clone("h_Data");
 TH1D* h_Total_MC =(TH1D*)h_MC_Total_input->Clone("h_Total_MC");
 THStack* Stack = (THStack*)Stack_input->Clone("Stack");

  if(DoBinWidthNorm==true){
    h_Data->Scale(1.0,"width");
    h_Total_MC->Scale(1.0,"width");
    ScaleHistogramsInStack(Stack, 1.0, "width" );
    //double ymaxnew = std::max( h_Data->GetMaximum(),
    //h_Total_MC->GetMaximum() ) * 1.5;
    //h_Data->GetYaxis()->SetRangeUser( 0., ymaxnew );
  }

  h_Data->SetLineColor( kBlack );
  h_Data->SetLineWidth( 3 );
  h_Data->SetMarkerStyle( kFullCircle );
  h_Data->SetMarkerSize( 0.8 );
  h_Data->SetStats( false );
  double ymax = std::max( h_Data->GetMaximum(), h_Total_MC->GetMaximum())* 1.7;
  if(YMax==99){YMax = ymax; }
 // slice_mc_plus_ext->hist_->GetMaximum() ) * 1.4;
  //if(DoBinWidt//hNorm==false) 
  h_Data->SetMaximum(YMax);  
  h_Data->GetYaxis()->SetRangeUser( 0., YMax );
  h_Data->GetYaxis()->SetTitle( Yaxis_title );
  h_Data->GetXaxis()->SetTitle( Xaxis_title );
  h_Data->GetYaxis()->SetLabelSize(.03);
  h_Data->GetXaxis()->SetLabelSize(.025);
  h_Data->GetXaxis()->SetTitleSize(0.035);
  h_Data->GetXaxis()->CenterTitle(kFALSE);
  h_Data->SetTitleOffset (1.01,"Y");
  h_Data->SetTitleOffset (.95,"X");
  h_Data->Draw( "e" );
  Stack->Draw( "hist same" );
  h_Total_MC->SetLineWidth( 1 );
  h_Total_MC->Draw( "same hist e" );
  h_Data->Draw( "same e" );

}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void DrawStackMCandData_WithMCBand(
    TH1* h_data_input, 
    TH1* h_MC_Total_input,
    THStack* Stack_input,
    char *Yaxis_title,
    char *Xaxis_title,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax,
    TCanvas *Canvas)
{

 TH1D* h_Data =(TH1D*)h_data_input->Clone("h_Data");
 TH1D* h_Total_MC =(TH1D*)h_MC_Total_input->Clone("h_Total_MC");
 THStack* Stack = (THStack*)Stack_input->Clone("Stack");

  if(DoBinWidthNorm==true){
    h_Data->Scale(1.0,"width");
    h_Total_MC->Scale(1.0,"width");
    ScaleHistogramsInStack(Stack, 1.0, "width" );
    //double ymaxnew = std::max( h_Data->GetMaximum(),
    //h_Total_MC->GetMaximum() ) * 1.5;
    //h_Data->GetYaxis()->SetRangeUser( 0., ymaxnew );
  }

  h_Data->SetLineColor( kBlack );
  h_Data->SetLineWidth( 3 );
  h_Data->SetMarkerStyle( kFullCircle );
  h_Data->SetMarkerSize( 0.8 );
  h_Data->SetStats( false );
  double ymax = std::max( h_Data->GetMaximum(), h_Total_MC->GetMaximum())* 1.7;
  if(YMax==99){YMax = ymax; }
 // slice_mc_plus_ext->hist_->GetMaximum() ) * 1.4;
  //if(DoBinWidt//hNorm==false) 
  h_Data->SetMaximum(YMax);  
  h_Data->GetYaxis()->SetRangeUser( 0., YMax );
  h_Data->GetYaxis()->SetTitle( Yaxis_title );
  h_Data->GetXaxis()->SetTitle( Xaxis_title );
  h_Data->GetYaxis()->SetLabelSize(.02);
  h_Data->GetXaxis()->SetLabelSize(.025);
  h_Data->GetXaxis()->SetTitleSize(0.035);
  h_Data->GetXaxis()->CenterTitle(kFALSE);
  h_Data->SetTitleOffset (1.01,"Y");
  h_Data->SetTitleOffset (.95,"X");
  h_Data->Draw( "e" );
  Stack->Draw( "hist same" );
  h_Total_MC->SetLineWidth( 0 );
  //h_Total_MC->Draw( "same hist e" );
    h_Total_MC->SetFillColorAlpha(kRed, 0.35);
    h_Total_MC->SetFillStyle(3001);
    h_Total_MC->SetMarkerStyle(0);
    h_Total_MC->DrawCopy("SAME E2");
    h_Total_MC->Draw("SAME AXIS");

  
  
  h_Data->Draw( "same e" );

}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////

void DrawStackMCandData(
    TH1* h_data_input, 
    TH1* h_MC_Total_input,
    THStack* Stack_input,
    bool DoBinWidthNorm, 
    double YMax)
{

 TH1D* h_Data =(TH1D*)h_data_input->Clone("h_Data");
 TH1D* h_Total_MC =(TH1D*)h_MC_Total_input->Clone("h_Total_MC");
 //TH1D* h_extBNB =(TH1D*)BG_beamOFF_input->Clone("h_extBNB");
 THStack* Stack = (THStack*)Stack_input->Clone("Stack");

  if(DoBinWidthNorm==true){
    h_Data->Scale(1.0,"width");
    h_Total_MC->Scale(1.0,"width");
    //h_extBNB->Scale(1.0,"width");
    ScaleHistogramsInStack(Stack, 1.0, "width" );
  }

  h_Data->SetLineColor( kBlack );
  h_Data->SetLineWidth( 3 );
  h_Data->SetMarkerStyle( kFullCircle );
  h_Data->SetMarkerSize( 0.8 );
  h_Data->SetStats( false ); 
  h_Data->GetYaxis()->SetRangeUser( 0., YMax );
  h_Data->Draw( "e" );
  Stack->Draw( "hist same" );
  h_Total_MC->SetLineWidth( 1 );
  h_Total_MC->Draw( "same hist e" );
  h_Data->Draw( "same e" );

}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void DrawStackMCandData_withBand(
    TH1* h_data_input, 
    TH1* h_MC_Total_input,
    TH1* BG_beamOFF_input,
    THStack* Stack_input,
    bool DoBinWidthNorm, 
    double YMax)
{

 TH1D* h_Data =(TH1D*)h_data_input->Clone("h_Data");
 TH1D* h_Total_MC =(TH1D*)h_MC_Total_input->Clone("h_Total_MC");
 //TH1D* h_extBNB =(TH1D*)BG_beamOFF_input->Clone("h_extBNB");
 THStack* Stack = (THStack*)Stack_input->Clone("Stack");

  if(DoBinWidthNorm==true){
    h_Data->Scale(1.0,"width");
    h_Total_MC->Scale(1.0,"width");
    //h_extBNB->Scale(1.0,"width");
    ScaleHistogramsInStack(Stack, 1.0, "width" );
  }

  h_Data->SetLineColor( kBlack );
  h_Data->SetLineWidth( 3 );
  //:h_Data->SetMarkerStyle( kFullCircle );
  h_Data->SetMarkerSize( 0.8 );
  
    h_Data->SetMarkerStyle(20);
    h_Data->SetMarkerSize(.5);
    h_Data->SetMarkerColor(1);
    h_Data->SetLineWidth(1.0);
    h_Data->SetLineStyle(1.0);
    h_Data->SetLineColor(1);
    h_Data->DrawCopy("SAME E1 X0");
  
  
  
  
  h_Data->SetStats( false ); 
  h_Data->GetYaxis()->SetRangeUser( 0., YMax );
  h_Data->Draw( "e" );
  Stack->Draw( "hist same" );
  
  h_Total_MC->SetLineWidth( 1 );
  //h_Total_MC->Draw( "same hist e" );
  
  h_Total_MC->SetFillColorAlpha(kRed, 0.45);
  h_Total_MC->SetFillStyle(3001);
  h_Total_MC->SetMarkerStyle(0);

  h_Total_MC->DrawCopy("SAME E2");
  h_Total_MC->Draw("SAME AXIS");
  h_Total_MC->SetFillColor(0);
  h_Total_MC->SetLineColor(2);
  h_Total_MC->SetLineStyle(1);
  h_Total_MC->SetLineWidth(2);
  h_Data->DrawCopy("SAME E1 X0");

}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void DrawStack(
    TH1* h_input, 
    THStack* Stack_input,
    char *Yaxis_title,
    char *Xaxis_title,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax
    )
{

 TH1D* histTotal =(TH1D*)h_input->Clone("histTotal");
 THStack* Stack = (THStack*)Stack_input->Clone("Stack");

  if(DoBinWidthNorm==true){
   histTotal->Scale(1.0,"width");
    //h_extBNB->Scale(1.0,"width");
    ScaleHistogramsInStack(Stack, 1.0, "width" );
  }
  
  histTotal->SetTitle(Title);
  histTotal->GetYaxis()->SetTitle( Yaxis_title );
  histTotal->GetXaxis()->SetTitle( Xaxis_title );
  histTotal->GetYaxis()->SetLabelSize(.024);
  histTotal->GetXaxis()->SetLabelSize(.025);
  histTotal->GetXaxis()->SetTitleSize(0.035);
  histTotal->GetXaxis()->CenterTitle();
  histTotal->SetLineColor( kBlack );
  histTotal->SetLineWidth( 3 );
  histTotal->SetStats( false ); 
  histTotal->GetYaxis()->SetRangeUser( 0., YMax );
  histTotal->Draw( "e" );
  Stack->Draw( "hist same" );
  histTotal->Draw( "e Same" );

}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void DrawStack(
    TH1* h_input, 
    THStack* Stack_input,
    bool DoBinWidthNorm, 
    char *Title,
    double YMax
    )
{

 TH1D* histTotal =(TH1D*)h_input->Clone("histTotal");
 THStack* Stack = (THStack*)Stack_input->Clone("Stack");

  if(DoBinWidthNorm==true){
   histTotal->Scale(1.0,"width");
    //h_extBNB->Scale(1.0,"width");
    ScaleHistogramsInStack(Stack, 1.0, "width" );
  }
  
  histTotal->SetTitle(Title);
  histTotal->GetYaxis()->SetLabelSize(.024);
  histTotal->GetXaxis()->SetLabelSize(.025);
  histTotal->GetXaxis()->SetTitleSize(0.035);
  histTotal->GetXaxis()->CenterTitle();
  histTotal->GetYaxis()->SetRangeUser( 0., YMax );
  histTotal->Draw( "e X0" );
  Stack->Draw( "hist same" );
  histTotal->Draw( "e X0 Same" );

}//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////

void drawBinRange(TH2* h, int axis, int bin, const char* varName, double text_size, const char* numFormatStr, bool left )
{
  double varmin=axis==1 ? h->GetXaxis()->GetBinLowEdge(bin) : h->GetYaxis()->GetBinLowEdge(bin);
  double varmax=axis==1 ? h->GetXaxis()->GetBinUpEdge(bin) :  h->GetYaxis()->GetBinUpEdge(bin);

  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(text_size);
  la->Draw();
}
//////////////////////////////////////////////////////////////////////////////
/// End of Funtion
//////////////////////////////////////////////////////////////////////////////
void drawString(std::string inputString,double text_size,  bool left )
{
 
  TLatex* la=0;
  TString text(inputString);
  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  la->SetTextSize(text_size);
  la->Draw();
}


//======================================================================
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void PlotDataStackedMC2D_ProjY(
	TH2* data,
	TH2* BG_beamOFF_input,
	TH2* Total_MC_input,
	TObjArray stack_input,
	std::vector<int> fillColors,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *yaxislabel, char *zaxislabel_units,
	double Ymax, bool setMaxY, bool doMultipliers,
	std::vector<double> YMultipliers,
	bool do_bin_width_norm,
	double text_size,
	double POT_DATA)
{
	
	TH2D *dataHist = (TH2D*)data->Clone(uniq());
	TH2D *dataHist_BG = (TH2D*)BG_beamOFF_input->Clone(uniq());
	//TH2 *dataHist_clone = data->Clone("dataHist_clone");
	TH2D* h_MC_total = (TH2D*)Total_MC_input->Clone(uniq()); 
	
	
	//TObjArray *stack =  &stack_input;
	TObjArray *stack = (TObjArray*) stack_input.Clone(uniq());
	stack->SetOwner(kTRUE);
	int nVars = stack->GetEntries();
	int nbins_X = data->GetNbinsX();
	int nbins_Y = data->GetNbinsY();
	std::cout << "The NEntries for TObjArray: " << nVars << std::endl;

	TLegend *leg = new TLegend(0.41, 0.12, 0.96, 0.32);
    std::string Leg_title_string = get_legend_title(POT_DATA);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.025);
	leg->SetNColumns(3);
    leg->SetHeader(Leg_title_string.c_str());
	//std::cout<<"nbins_X = " << nbins_X << std::endl;
	//std::cout<<"nbins_Y = " << nbins_Y << std::endl;
	//std::cout<<"nbins_X_mc = " << nbins_X_mc << std::endl;
	//std::cout<<"nbins_Y_mc = " << nbins_Y_mc << std::endl;

	double min_XAxis = data->GetYaxis()->GetXmin();
	double max_XAxis = data->GetYaxis()->GetBinUpEdge(data->GetNbinsY()-1)* 1.09;


	double maxmax = 1.25* GetMaxFromProjectionY(dataHist, do_bin_width_norm);

	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	grid_x = sqrt(nbins_X) + 1;
	grid_y = nbins_X / (grid_x - 1);
	Nbins = dataHist->GetNbinsX();

	if (DeBug == true) std::cout << nbins_X - grid_x * grid_y << std::endl;
	if (DeBug == true) std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");
	gc->SetRightMargin(0.01);
	gc->SetLeftMargin(0.1);
	gc->ResetPads();

	for (int i = 1; i <= Nbins; ++i)
	{
		TObjArray *mcprojyArr = new TObjArray();
		mcprojyArr->SetOwner(kTRUE);

		TH1D *hdatayproj = dataHist->ProjectionY(uniq(),i, i);
		
		TH1D *hdata_BGyproj = dataHist_BG->ProjectionY(uniq(), i, i);
		      hdata_BGyproj->SetFillColor( 28 );
              hdata_BGyproj->SetLineColor( 28 );
              hdata_BGyproj->SetLineWidth( 2 );
              hdata_BGyproj->SetFillStyle( 3005 );
              hdata_BGyproj->SetStats( false );
		
		
		
		TH1D *hMC_yproj = h_MC_total->ProjectionY(uniq(), i, i);
			
		hdatayproj->SetMaximum(maxmax);
		hdata_BGyproj->SetMaximum(maxmax);
		hMC_yproj->SetMaximum(maxmax);
		
		mcprojyArr->Add(hdata_BGyproj);

		
		if (i == 1) {
		leg->AddEntry(hdatayproj, "Data[Stat]", "lpe");
		leg->AddEntry(hMC_yproj, "Total MC[Stat]", "l");
		}


		 for (int iHist = 0; iHist < nVars; ++iHist)
		 {
			TH2D *mcHist_stack = (TH2D*) stack->At(iHist)->Clone("mcHist_stack");
			//int ndf;
			std::string words = mcHist_stack->GetTitle();
			std::cout << "iHist  = " << iHist << " Title = " << words << std::endl;

			TH1D *mcHist_stack_ProjectionY = mcHist_stack->ProjectionY(uniq(), i, i);
			//if (do_bin_width_norm == true) mcHist_stack_ProjectionY->Scale(1., "width");
			if (doMultipliers == true) mcHist_stack_ProjectionY->Scale(YMultipliers.at(i - 1));

			mcHist_stack_ProjectionY->SetTitle(words.c_str());
			mcHist_stack_ProjectionY->SetXTitle("");

			//mcHist_stack_ProjectionY->SetFillColor(fillColors.at(iHist));
			mcHist_stack_ProjectionY->SetLineWidth(0);

			if (i == 1)
			{
				char words_char[words.length() + 1];
				strcpy(words_char, words.c_str());
				leg->AddEntry(mcHist_stack_ProjectionY, words_char, "f");
			}

         mcprojyArr->Add(mcHist_stack_ProjectionY);
		 }
		 
		if (i == 1) {
		
		leg->AddEntry(hdata_BGyproj, "BeamOFF[Stat]", "f");
		}
		
		//leg->AddEntry(mcHist_stack_ProjectionY, words_char, "f");
		

		
		
		gc->cd(i);
		THStack* Projection_Stack = CreateStackFromTObjArray(mcprojyArr);


		if (do_bin_width_norm == true) {
		hdatayproj->Scale(1., "width");
		//hdata_BGyproj->Scale(1., "width");
		hMC_yproj->Scale(1., "width");
		ScaleHistogramsInStack(Projection_Stack, 1.0, "width" );
		}
		
		if (doMultipliers == true) {
		hdatayproj->Scale(YMultipliers.at(i - 1));
		//hdata_BGyproj->Scale(YMultipliers.at(i - 1));
		hMC_yproj->Scale(YMultipliers.at(i - 1));
		ScaleHistogramsInStack(Projection_Stack, YMultipliers.at(i - 1));
		}
		

           DrawStackMCandData(
             hdatayproj, 
             hMC_yproj,
            Projection_Stack,
            /*do_bin_width_norm*/false, 
            maxmax);




		//drawBinRange(dataHist_clone, 1, i, xaxislabel, ".1f", false);
		drawBinRange(dataHist, 1, i, xaxislabel, text_size, ".2f", false);

		if (doMultipliers && YMultipliers.at(i - 1) != 1)
		{
			auto pad = gc->cd(i);
			TLatex *la2 = new TLatex(1 - pad->GetRightMargin() - 0.02,
				1 - pad->GetTopMargin() - 0.08,
				TString::Format("#times %.1f", YMultipliers.at(i - 1)));
			la2->SetTextAlign(33);	// top right
			la2->SetNDC();
			la2->SetTextFont(42);
			la2->SetTextSize(0.03);
			la2->Draw();
		}

		//mcprojyArr->Clear();
	}	/// End of Loop

    std::cout<<"Finished 2D drawing loop "<< std::endl;

	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.02);
	/*
	double MaxY = gc->GetPadMax();
	std::cout<<"Got Maxium Pad"<< std::endl;
	if (setMaxY == false) gc->SetYLimits(0, MaxY *1.35);
	else
	{
		gc->SetYLimits(0, Ymax);
	}
	*/
    gc->SetYLimits(0, Ymax);
    gc->SetXLimits(min_XAxis, max_XAxis);
	gc->SetXTitle(yaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel_units);
	gc->SetHistTexts();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);
	//leg->Clear();
	delete gc;

}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void PlotDataStackedMC2D_ProjX(
	TH2* data,
	TH2* BG_beamOFF_input,
	TH2* Total_MC_input,
	TObjArray stack_input,
	std::vector<int> fillColors,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *yaxislabel, char *zaxislabel_units,
	double Ymax, bool setMaxY, bool doMultipliers,
	std::vector<double> YMultipliers,
	bool do_bin_width_norm,
	double text_size,
	double POT_DATA)
{
	
	TH2D *dataHist = (TH2D*)data->Clone(uniq());
	TH2D *dataHist_BG = (TH2D*)BG_beamOFF_input->Clone(uniq());
	//TH2 *dataHist_clone = data->Clone("dataHist_clone");
	TH2D* h_MC_total = (TH2D*)Total_MC_input->Clone(uniq()); 
	
	
	//TObjArray *stack =  &stack_input;
	TObjArray *stack = (TObjArray*) stack_input.Clone(uniq());
	stack->SetOwner(kTRUE);
	int nVars = stack->GetEntries();
	int nbins_X = data->GetNbinsX();
	int nbins_Y = data->GetNbinsY();
	std::cout << "The NEntries for TObjArray: " << nVars << std::endl;

	TLegend *leg = new TLegend(0.41, 0.12, 0.96, 0.32);
    std::string Leg_title_string = get_legend_title(POT_DATA);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.025);
	leg->SetNColumns(3);
    leg->SetHeader(Leg_title_string.c_str());
	//std::cout<<"nbins_X = " << nbins_X << std::endl;
	//std::cout<<"nbins_Y = " << nbins_Y << std::endl;
	//std::cout<<"nbins_X_mc = " << nbins_X_mc << std::endl;
	//std::cout<<"nbins_Y_mc = " << nbins_Y_mc << std::endl;

	double min_XAxis = data->GetXaxis()->GetXmin();
	double max_XAxis = data->GetXaxis()->GetBinUpEdge(data->GetNbinsX()-1)* 1.09;


	double maxmax = 1.25* GetMaxFromProjectionX(dataHist, do_bin_width_norm);

	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	grid_x = sqrt(nbins_Y) + 1;
	grid_y = nbins_Y / (grid_x - 1);
	Nbins = dataHist->GetNbinsY();

	if (DeBug == true) std::cout << nbins_X - grid_x * grid_y << std::endl;
	if (DeBug == true) std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");
	gc->SetRightMargin(0.01);
	gc->SetLeftMargin(0.1);
	gc->ResetPads();

	for (int i = 1; i <= Nbins; ++i)
	{
		TObjArray *mcprojyArr = new TObjArray();
		mcprojyArr->SetOwner(kTRUE);

		TH1D *hdatayproj = dataHist->ProjectionX(uniq(),i, i);
		
		TH1D *hdata_BGyproj = dataHist_BG->ProjectionX(uniq(), i, i);
		      hdata_BGyproj->SetFillColor( 28 );
              hdata_BGyproj->SetLineColor( 28 );
              hdata_BGyproj->SetLineWidth( 2 );
              hdata_BGyproj->SetFillStyle( 3005 );
              hdata_BGyproj->SetStats( false );
		
		
		
		TH1D *hMC_yproj = h_MC_total->ProjectionX(uniq(), i, i);
			
		hdatayproj->SetMaximum(maxmax);
		hdata_BGyproj->SetMaximum(maxmax);
		hMC_yproj->SetMaximum(maxmax);
		
		mcprojyArr->Add(hdata_BGyproj);

		
		if (i == 1) {
		leg->AddEntry(hdatayproj, "Data[Stat]", "lpe");
		leg->AddEntry(hMC_yproj, "Total MC[Stat]", "l");
		}


		 for (int iHist = 0; iHist < nVars; ++iHist)
		 {
			TH2D *mcHist_stack = (TH2D*) stack->At(iHist)->Clone("mcHist_stack");
			//int ndf;
			std::string words = mcHist_stack->GetTitle();
			std::cout << "iHist  = " << iHist << " Title = " << words << std::endl;

			TH1D *mcHist_stack_ProjectionY = mcHist_stack->ProjectionX(uniq(), i, i);
			//if (do_bin_width_norm == true) mcHist_stack_ProjectionY->Scale(1., "width");
			if (doMultipliers == true) mcHist_stack_ProjectionY->Scale(YMultipliers.at(i - 1));

			mcHist_stack_ProjectionY->SetTitle(words.c_str());
			mcHist_stack_ProjectionY->SetXTitle("");

			//mcHist_stack_ProjectionY->SetFillColor(fillColors.at(iHist));
			mcHist_stack_ProjectionY->SetLineWidth(0);

			if (i == 1)
			{
				char words_char[words.length() + 1];
				strcpy(words_char, words.c_str());
				leg->AddEntry(mcHist_stack_ProjectionY, words_char, "f");
			}

         mcprojyArr->Add(mcHist_stack_ProjectionY);
		 }
		 
		if (i == 1) {
		
		leg->AddEntry(hdata_BGyproj, "BeamOFF[Stat]", "f");
		}
		
		//leg->AddEntry(mcHist_stack_ProjectionY, words_char, "f");
		

		
		
		gc->cd(i);
		THStack* Projection_Stack = CreateStackFromTObjArray(mcprojyArr);


		if (do_bin_width_norm == true) {
		hdatayproj->Scale(1., "width");
		//hdata_BGyproj->Scale(1., "width");
		hMC_yproj->Scale(1., "width");
		ScaleHistogramsInStack(Projection_Stack, 1.0, "width" );
		}
		
		if (doMultipliers == true) {
		hdatayproj->Scale(YMultipliers.at(i - 1));
		//hdata_BGyproj->Scale(YMultipliers.at(i - 1));
		hMC_yproj->Scale(YMultipliers.at(i - 1));
		ScaleHistogramsInStack(Projection_Stack, YMultipliers.at(i - 1));
		}
		

           DrawStackMCandData(
             hdatayproj, 
             hMC_yproj,
            Projection_Stack,
            /*do_bin_width_norm*/false, 
            maxmax);




		//drawBinRange(dataHist_clone, 1, i, xaxislabel, ".1f", false);
		drawBinRange(dataHist, 0, i, yaxislabel, text_size, ".2f", true);

		if (doMultipliers && YMultipliers.at(i - 1) != 1)
		{
			auto pad = gc->cd(i);
			TLatex *la2 = new TLatex(1 - pad->GetRightMargin() - 0.02,
				1 - pad->GetTopMargin() - 0.08,
				TString::Format("#times %.1f", YMultipliers.at(i - 1)));
			la2->SetTextAlign(33);	// top right
			la2->SetNDC();
			la2->SetTextFont(42);
			la2->SetTextSize(0.03);
			la2->Draw();
		}

		//mcprojyArr->Clear();
	}	/// End of Loop

    std::cout<<"Finished 2D drawing loop "<< std::endl;

	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.02);
	/*
	double MaxY = gc->GetPadMax();
	std::cout<<"Got Maxium Pad"<< std::endl;
	if (setMaxY == false) gc->SetYLimits(0, MaxY *1.35);
	else
	{
		gc->SetYLimits(0, Ymax);
	}
	*/
    gc->SetYLimits(0, Ymax);
    gc->SetXLimits(min_XAxis, max_XAxis);
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel_units);
	gc->SetHistTexts();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);
	//leg->Clear();
	delete gc;

}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
double GetMaxFromProjectionY(TH2 *h1,TH2 *h2, bool checkBinwidth){

  TH2D* hist1 = (TH2D*)h1->Clone("hist1");
  TH2D* hist2 = (TH2D*)h2->Clone("hist2");

  int Nbins1 = hist1->GetNbinsX();
  int Nbins2 = hist2->GetNbinsX();

  if(Nbins1 !=Nbins2){std::cout<<"Bins are not matching for GetMaxFromProjectionX "<<std::endl;
    std::cout<< " Nbins1 =  "<< Nbins1 << "  Nbins2 = "<< Nbins2 << std::endl;
    std::cout<< " Title1 =  "<< hist1->GetTitle()<< "  title2 = "<< hist2->GetTitle() << std::endl;

    assert(false);}
  double MaxValue = -1;

  for(int i = 1; i <=Nbins1; ++i) {
    TH1D* hproj_1 =
    hist1->ProjectionY(uniq(), i, i);

    TH1D* hproj_2 =
    hist2->ProjectionY(uniq(), i, i);


    if(checkBinwidth==true){
      hproj_1->Scale(1,"width");
      hproj_2->Scale(1,"width");
    }

    double max1 = hproj_1->GetMaximum();
    double max2 = hproj_2->GetMaximum();
    double Maxpro = (max1 < max2) ? max2 : max1;

    if (Maxpro > MaxValue){MaxValue = Maxpro;}

  }

  return MaxValue;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
double GetMaxFromProjectionX(TH2 *h1,TH2 *h2, bool checkBinwidth){

  TH2D* hist1 = (TH2D*)h1->Clone("hist1");
  TH2D* hist2 = (TH2D*)h2->Clone("hist2");

  int Nbins1 = hist1->GetNbinsX();
  int Nbins2 = hist2->GetNbinsX();

  if(Nbins1 !=Nbins2){std::cout<<"Bins are not matching for GetMaxFromProjectionX "<<std::endl;
    std::cout<< " Nbins1 =  "<< Nbins1 << "  Nbins2 = "<< Nbins2 << std::endl;
    std::cout<< " Title1 =  "<< hist1->GetTitle()<< "  title2 = "<< hist2->GetTitle() << std::endl;

    assert(false);}
  double MaxValue = -1;

  for(int i = 1; i <=Nbins1; ++i) {
    TH1D* hproj_1 =
    hist1->ProjectionX(uniq(), i, i);

    TH1D* hproj_2 =
    hist2->ProjectionX(uniq(), i, i);


    if(checkBinwidth==true){
      hproj_1->Scale(1,"width");
      hproj_2->Scale(1,"width");
    }

    double max1 = hproj_1->GetMaximum();
    double max2 = hproj_2->GetMaximum();
    double Maxpro = (max1 < max2) ? max2 : max1;

    if (Maxpro > MaxValue){MaxValue = Maxpro;}

  }

  return MaxValue;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
double GetMaxFromProjectionY(TH2 * h1, bool checkBinwidth){

  TH2D* hist1 = (TH2D*)h1->Clone("hist1");

  int Nbins1 = hist1->GetNbinsY();
  double MaxValue = -1;

  for(int i = 1; i <=Nbins1; ++i) {
   TH1D* hproj_1 =
    hist1->ProjectionY(uniq(), i, i);

    if(checkBinwidth==true){
      hproj_1->Scale(1,"width");
    }

    double Maxpro = hproj_1->GetMaximum();

    if (Maxpro > MaxValue){MaxValue = Maxpro;}

  }

  return MaxValue;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
double GetMaxFromProjectionX(TH2 * h1, bool checkBinwidth){

  TH2D* hist1 = (TH2D*)h1->Clone("hist1");

  int Nbins1 = hist1->GetNbinsX();
  double MaxValue = -1;

  for(int i = 1; i <=Nbins1; ++i) {
   TH1D* hproj_1 =
    hist1->ProjectionX(uniq(), i, i);

    if(checkBinwidth==true){
      hproj_1->Scale(1,"width");
    }

    double Maxpro = hproj_1->GetMaximum();

    if (Maxpro > MaxValue){MaxValue = Maxpro;}

  }

  return MaxValue;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
THStack* CreateStackFromTObjArray(TObjArray* histArray) {
    // Create a THStack
    THStack* stack = new THStack(uniq(), "Stacked Histograms");

    // Loop over the histograms in the TObjArray and add them to the stack
    for (int i = 0; i < histArray->GetEntries(); ++i) {
        TH1D* hist = dynamic_cast<TH1D*>(histArray->At(i));
        if (hist) {
            stack->Add(hist);
        } else {
            // Handle the case where the object in the TObjArray is not a TH1
            // Print an error message or handle it according to your needs
            std::cerr << "Object at index " << i << " is not a TH1." << std::endl;
        }
    }

    return stack;
}//////End of Function 
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

void AddCutArrow(
       double cut_location,
       double y1,
       double y2,
       double arrow_length,
      bool ArrowHeadGoLeft,
       int arrow_line_width,
       int arrow_line_style,
       int arrow_line_color)
{

    double arrow_tip = arrow_length;
    if (true == ArrowHeadGoLeft) {
        arrow_tip *= -1.0;
    }

    TLine line;
    line.SetLineWidth(arrow_line_width);
    line.SetLineStyle(arrow_line_style);
    line.SetLineColor(arrow_line_color);
    line.DrawLine(cut_location,y1,cut_location,y2);

   std::string arrow_type = "|>";
    double arrow_size = 0.01;
    TArrow arrow;
    arrow.SetLineWidth(arrow_line_width);
    arrow.SetLineStyle(arrow_line_style);
    arrow.SetLineColor(arrow_line_color);
    arrow.DrawArrow(cut_location,y2,cut_location+arrow_tip,y2,arrow_size,arrow_type.c_str());
}/////end of function////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

/*
void DrawStack(
  TH2D* category_hist_input,
  TH1D* Data_input,
  TH1D* Total_MC_input,
  TH1D* extBNB_input,
  SliceHistogram* slice_bnb,
  SliceHistogram* slice_ext,
  SliceHistogram* slice_mc_plus_ext,
  const Slice& slice,
  std::string pdfTitle,
  std::string TotalTitle,
  TCanvas *c1,
  double ymax)
{
  const auto& EventInterp = EventCategoryInterpreter::Instance();
  TLegend* lg1_stacked = new TLegend( 0.35, 0.6, 0.85, 0.85 );
  lg1_stacked->SetNColumns(2);
  lg1_stacked->SetBorderSize(0);
  Double_t defaultTextSize = lg1_stacked->GetTextSize();
  lg1_stacked->SetTextSize(.03); //

  TH2D* category_hist = (TH2D*)category_hist_input->Clone("category_hist");
  TH1D* h_Data =(TH1D*)slice_bnb->hist_.get()->Clone("h_Data");
  lg1_stacked->AddEntry(h_Data, "Data", "pe" );
  TH1D* h_Total_MC =(TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_Total_MC");
  lg1_stacked->AddEntry(h_Total_MC, "Total MC", "le" );
  
  TH1D* h_extBNB =(TH1D*)extBNB_input->Clone("h_extBNB");
  EventInterp.set_ext_histogram_style( h_extBNB );
  
  EventInterp.set_ext_histogram_style( slice_ext->hist_.get() );
  auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );

  THStack* slice_pred_stack = new THStack( "mc+ext", "" );
  slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

  const auto& cat_map = EventInterp.label_map();
  int cat_bin_index = cat_map.size();

  for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
  {
    EventCategory cat = iter->first;
    TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
      cat_bin_index, cat_bin_index );
    temp_mc_hist->SetDirectory( nullptr );
    
    SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
      *temp_mc_hist, slice  );

    EventInterp.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

    slice_pred_stack->Add( temp_slice_mc->hist_.get() );

    std::string lg1_label = EventInterp.label(cat);
    lg1_stacked->AddEntry(temp_slice_mc->hist_.get(),lg1_label.c_str() , "f" );
    //std::string cat_col_prefix = "MC" + std::to_string( cat );

    --cat_bin_index;
  }

 bool makeNormWidth = true; 

  DrawStackMCandData(
    h_Data, 
    h_Total_MC,
    h_extBNB,
    slice_pred_stack,
    "NEvents / Bin Width",
    "",
    makeNormWidth, 
    "title",
    ymax,
    c1);


   lg1_stacked->AddEntry(slice_ext->hist_.get(),"EXT BNB" , "f" );
   lg1_stacked->Draw( "same" );
  //  std::cout<<"printing chi values "<< std::endl;
  TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kRed);
    char textplace[1024];
    sprintf(textplace, "#chi^{2}/ndf = %.2f / %i",chi2_result.chi2_ , chi2_result.num_bins_ );
    text->DrawLatex(0.15, 0.85, textplace);

    sprintf(textplace, "p-value = %.4f",chi2_result.p_value_ );
    text->DrawLatex(0.15, 0.80, textplace);
    //AddHistoTitle(TotalTitle.c_str(), .035);
  ////////////////////////////////////////
  ///// Finished with First Plot
  ////////////////////////////////////////
  char pdf_title[1024];
  
  sprintf(pdf_title, "%s.pdf", pdfTitle.c_str());
  c1 -> Print(pdf_title);

}
*/

void DrawOverlayCanvas(UBTH2Poly* originalPlot, 
char *drawoption,  double X_windowPos, double Y_windowPos,
double X_windowSize, double Y_windowSize) {
    // Create a new canvas
    TCanvas* overlayCanvas = new TCanvas("overlayCanvas", "Overlay Canvas", 800, 600);
    UBTH2Poly* originalPlot_BinNumber = originalPlot->GetCopyWithBinNumbers("originalPlot_BinNumber");
    // Set the position and size of the overlay canvas
    overlayCanvas->SetWindowPosition(X_windowPos, Y_windowPos);
    overlayCanvas->SetWindowSize(X_windowSize, Y_windowSize);

    // Draw the original TH2D on the overlay canvas
  overlayCanvas->cd(); 

    // Create a smaller TH2D for the overlay
 char input[1024]; 
  sprintf(input, "%s SAME",drawoption );
    // Draw the smaller TH2D on top of the original plot
    originalPlot_BinNumber->Draw(input);
}
///////////////////////////////////////////////////
 /*
 void PlotMC2D(
 std::map<Binning2D, TH1D*> InputHistmap,
 std::map<Binning2D , std::string > BinStringMap,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *legendTitle, char *zaxislabel_units,
	double Ymax, bool setMaxY, double Ymin,  bool doMultipliers,
	std::vector<double> YMultipliers,
	bool do_bin_width_norm,
	double text_size,
	double POT_DATA)
{
	auto BinVector = GetProjectBinVector();

	int nbins_X = InputHistmap.size()-1;
	int nbins_Y = InputHistmap.size()-1;
	//std::cout << "The NEntries for TObjArray: " << nVars << std::endl;

	 TLegend *leg = new TLegend(0.41, 0.1, 0.96, 0.25);
    std::string Leg_title_string = get_legend_title(POT_DATA);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.02);
	leg->SetNColumns(3);
  leg->SetHeader("MicroBooNE Simulation");
	

	double min_XAxis = InputHistmap[BinVector.at(0)]->GetXaxis()->GetXmin();
	double max_XAxis = InputHistmap[BinVector.at(0)]->GetXaxis()->GetBinUpEdge(InputHistmap[BinVector.at(0)]->GetNbinsX());
 std::cout<< "max_XAxis = "<< max_XAxis<< std::endl;
 double maxmax = -1; 


 if(setMaxY==true){
  maxmax= Ymax; 
 }
 else {
   double max =  1.15* MaxYofMap(InputHistmap);
 }


	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	grid_x = sqrt(nbins_X + 1);
	grid_y = nbins_X / (grid_x - 1);
	Nbins = InputHistmap.size();

	if (DeBug == true) std::cout << nbins_X - grid_x * grid_y << std::endl;
	if (DeBug == true) std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");

	

	for (int i = 1; i <= Nbins; ++i)
	{
		TH1D *hist = (TH1D *)InputHistmap[BinVector.at(i-1)]->Clone(uniq());
		hist->SetMaximum(maxmax);
		hist->SetMinimum(0.0);
		hist->GetXaxis()->SetRangeUser(min_XAxis,max_XAxis);
		
		if (i == 1) {
		leg->AddEntry(hist, legendTitle, "l");
		}
		
	//////////////////////////////	
		gc->cd(i);
	/////////////////////////////


		if (do_bin_width_norm == true) {
		hist->Scale(1., "width");
		}
		
		if (doMultipliers == true) {
		hist->Scale(YMultipliers.at(i - 1));
		}
		
		hist->Draw("Hist");
   
   drawString(BinStringMap[BinVector.at(i-1)], text_size, false );
   
		if (doMultipliers && YMultipliers.at(i - 1) != 1)
		{
			auto pad = gc->cd(i);
			TLatex *la2 = new TLatex(1 - pad->GetRightMargin() - 0.02,
				1 - pad->GetTopMargin() - 0.08,
				TString::Format("#times %.1f", YMultipliers.at(i - 1)));
			la2->SetTextAlign(33);	// top right
			la2->SetNDC();
			la2->SetTextFont(42);
			la2->SetTextSize(0.03);
			la2->Draw();
		}

		//mcprojyArr->Clear();
	}	/// End of Loop

    std::cout<<"Finished 2D drawing loop "<< std::endl;

	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.015);


	
	gc->SetInterpadSpace(.005);
	//gc->SetRightMargin(0.05);
	//gc->SetLeftMargin(0.08);
  gc->SetYLimits(Ymin, Ymax);
  gc->SetXLimits(min_XAxis,max_XAxis);
  gc->ResetPads();
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel_units);
	gc->SetTitleAlignmentFor6Hist();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);
	//leg->Clear();
	delete gc;

}
*/
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
 /*
 void PlotMC2D_purity_Eff(
 std::map<Binning2D, TH1D*> InputHistmap1,
 std::map<Binning2D, TH1D*> InputHistmap2,
 std::map<Binning2D , std::string > BinStringMap,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *zaxislabel_units,
	double Ymax, bool setMaxY, double Ymin,  bool doMultipliers,
	std::vector<double> YMultipliers,
	bool do_bin_width_norm,
	double text_size)
{
	auto BinVector = GetProjectBinVector();

	int nbins_X = InputHistmap1.size()-1;
	int nbins_Y = InputHistmap1.size()-1;
	//std::cout << "The NEntries for TObjArray: " << nVars << std::endl;

	 TLegend *leg = new TLegend(0.4, 0.08, 0.7, 0.2);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(1);
  leg->SetHeader("MicroBooNE Simulation");
	

	double min_XAxis = InputHistmap1[BinVector.at(0)]->GetXaxis()->GetXmin();
	double max_XAxis = InputHistmap2[BinVector.at(0)]->GetXaxis()->GetBinUpEdge(InputHistmap2[BinVector.at(0)]->GetNbinsX());
 std::cout<< "max_XAxis = "<< max_XAxis<< std::endl;
 double maxmax = -1; 




 if(setMaxY==true){
  maxmax= Ymax; 
 }
 else {
   double max1 =  1.15* MaxYofMap(InputHistmap1);
   double max2 =  1.15* MaxYofMap(InputHistmap2);
   maxmax = 1.0;
   
 }


	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	grid_x = sqrt(nbins_X + 1);
	grid_y = nbins_X / (grid_x - 1);
	Nbins = InputHistmap1.size();

	if (DeBug == true) std::cout << nbins_X - grid_x * grid_y << std::endl;
	if (DeBug == true) std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");

	gc->SetInterpadSpace(.005);
	gc->SetBottomMargin(.00);
  gc->SetTopMargin(.02);
  gc->SetRightMargin(.05);
  gc->ResetPads();


	for (int i = 1; i <= Nbins; ++i)
	{
		TH1D *hist1 = (TH1D *)InputHistmap1[BinVector.at(i-1)]->Clone(uniq());
		TH1D *hist3 = (TH1D *)InputHistmap1[BinVector.at(i-1)]->Clone(uniq());
		TH1D *hist2 = (TH1D *)InputHistmap2[BinVector.at(i-1)]->Clone(uniq());
		hist3->Multiply(hist2);
		hist1->SetMaximum(maxmax);
		hist1->SetMinimum(Ymin);
		hist1->GetXaxis()->SetRangeUser(min_XAxis,max_XAxis);
		hist1->SetLineColor(8);
    hist2->SetLineColor(9); 
    
    hist3->SetLineColor(2);

		hist1->SetTitle(histotitle);
		
		hist1->SetLineWidth(2);
    hist2->SetLineWidth(2);
    hist3->SetLineWidth(2);
		
		if (i == 1) {
		leg->AddEntry(hist1, "Efficiency", "l");
		leg->AddEntry(hist2, "Purity", "l");
		leg->AddEntry(hist3, "Eff*Purity", "l");
		
		}
		
	//////////////////////////////	
		gc->cd(i);
	/////////////////////////////


		if (do_bin_width_norm == true) {
		hist1->Scale(1., "width");
		hist2->Scale(1., "width");
		hist3->Scale(1., "width");
		}
		
		if (doMultipliers == true) {
		hist1->Scale(YMultipliers.at(i - 1));
		hist2->Scale(YMultipliers.at(i - 1));
		hist3->Scale(YMultipliers.at(i - 1));
		}
		
		hist1->Draw("Hist");
		hist2->Draw("Hist SAME");
		hist3->Draw("Hist SAME");
   
   drawString(BinStringMap[BinVector.at(i-1)], text_size, false );
   
		if (doMultipliers && YMultipliers.at(i - 1) != 1)
		{
			auto pad = gc->cd(i);
			TLatex *la2 = new TLatex(1 - pad->GetRightMargin() - 0.02,
				1 - pad->GetTopMargin() - 0.08,
				TString::Format("#times %.1f", YMultipliers.at(i - 1)));
			la2->SetTextAlign(33);	// top right
			la2->SetNDC();
			la2->SetTextFont(42);
			la2->SetTextSize(0.03);
			la2->Draw();
		}

		//mcprojyArr->Clear();
	}	/// End of Loop

    std::cout<<"Finished 2D drawing loop "<< std::endl;

	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.015);


	

	//gc->SetRightMargin(0.05);
	//gc->SetLeftMargin(0.08);
  gc->SetYLimits(Ymin, Ymax);
  gc->SetXLimits(min_XAxis,max_XAxis);
  gc->ResetPads();
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel_units);
	gc->SetTitleAlignmentFor6Hist();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);
	//leg->Clear();
	delete gc;

}
*/
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
/*
void PlotFractionofEvents2D_ProjX(
	PlotUtils::MnvH2D *mc, const TObjArray *stack,
	std::vector<int> fillColors, char *dataname, char *dataname_total,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *yaxislabel, char *zaxislabel_units
)
{
	PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCInclusiveHeliumStyle);
	TH2 *mc_clone = mc->Clone("mcth1");
	PlotUtils::MnvH2D *mcHist_clone = mc->Clone("mcHist_clone");
	double CompleteTotal = mcHist_clone->Integral();
	int nVars = stack->GetEntries();
	int nbins_X = mcHist_clone->GetNbinsX();
	int nbins_Y = mcHist_clone->GetNbinsY();
	std::cout << "The NEntries for TObjArray: " << nVars << std::endl;
	int nbins_X_mc = mcHist_clone->GetNbinsX();
	int nbins_Y_mc = mcHist_clone->GetNbinsY();
	double MinX = -.003; //mcHist_clone->GetXaxis()->GetXmin();
	//if(MinX==0){MinX = .01; }
	double MaxX = mcHist_clone->GetXaxis()->GetBinUpEdge(mcHist_clone->GetNbinsX()-1)*1.09;
	//std::cout<<"nbins_X = " << nbins_X << std::endl;
	//std::cout<<"nbins_Y = " << nbins_Y << std::endl;
	//std::cout<<"nbins_X_mc = " << nbins_X_mc << std::endl;
	//std::cout<<"nbins_Y_mc = " << nbins_Y_mc << std::endl;

	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	//   grid_x = sqrt(nbins_Y)+1;
	// grid_y = nbins_Y/(grid_x-1);
	// Nbins = dataHist->GetNbinsY();

	grid_x = sqrt(nbins_Y) + 1;
	grid_y = nbins_Y / (grid_x - 1);
	Nbins = mcHist_clone->GetNbinsY();

	TLegend * leg;
	if (Nbins == 5 || Nbins == 11)
	{
		leg = new TLegend(0.78, 0.1, 0.98, 0.35);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.015);
		leg->SetNColumns(1);
	}
	else
	{
		leg = new TLegend(0.4, 0.08, 0.95, 0.25);

		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.025);
		leg->SetNColumns(3);
		std::string minervatitle = GetMinervaLegendTitle();
        leg->SetHeader(minervatitle.c_str());
	}

	if (DeBug == true) std::cout << nbins_X - grid_x * grid_y << std::endl;
	if (DeBug == true) std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	std::vector<int> bins;	// = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");
	gc->SetRightMargin(0.06);
	gc->SetLeftMargin(0.06);
	gc->SetBottomMargin(0.004);
	gc->ResetPads();

	for (int i = 1; i <= Nbins; ++i)
	{
		bins.push_back(i);

		TObjArray *mcprojyArr = new TObjArray();

		PlotUtils::MnvH1D *hMCxproj =
			mcHist_clone->ProjectionX(uniq(), i, i)->Clone();
		PlotUtils::MnvH1D *hMCxproj_clone = (PlotUtils::MnvH1D *) hMCxproj->Clone("hMCxproj_clone");

		hMCxproj->ClearAllErrorBands();

		double area = hMCxproj->Integral(1, hMCxproj->GetNbinsX());
		hMCxproj->Scale(1.0 / area);
		hMCxproj_clone->Scale(1.0 / CompleteTotal);

		hMCxproj->SetMaximum(1.25);
		hMCxproj->SetMarkerStyle(20);
		hMCxproj_clone->SetMarkerStyle(21);
		hMCxproj_clone->SetMarkerSize(.5);
		if (i == 1)
		{
			leg->AddEntry(hMCxproj, dataname, "p");

			leg->AddEntry(hMCxproj_clone, dataname_total, "p");
		}

		for (int iHist = 0; iHist < nVars; ++iHist)
		{
			PlotUtils::MnvH2D *mcHist_stack = (PlotUtils::MnvH2D *) stack->At(iHist)->Clone("mcHist_stack");
			//int ndf;
			std::string words = mcHist_stack->GetTitle();
			std::cout << "iHist  = " << iHist << " words = " << words << std::endl;

			PlotUtils::MnvH1D *mcHist_stack_ProjectionY = mcHist_stack->ProjectionX(uniq(), i, i)->Clone();

			mcHist_stack_ProjectionY->SetTitle(words.c_str());
			mcHist_stack_ProjectionY->SetXTitle("");

			mcHist_stack_ProjectionY->SetFillStyle(1001);
			mcHist_stack_ProjectionY->SetFillColor(fillColors.at(iHist));
			mcHist_stack_ProjectionY->SetLineWidth(0);

			if (i == 1)
			{
				char words_char[words.length() + 1];
				strcpy(words_char, words.c_str());
				leg->AddEntry(mcHist_stack_ProjectionY, words_char, "f");
			}

			mcprojyArr->Add(mcHist_stack_ProjectionY);
		}

		gc->cd(i);
		mnvPlotter.axis_label_size = 0.02;
		mnvPlotter.legend_text_size = 0.02;
		mnvPlotter.data_marker_size = 0.4;
		const Int_t mcFillStyle = 1001;
		const Double_t mcScale = 1.0;
		//const std::string& legPos = "N";
		mnvPlotter.axis_maximum = 1.35;

		BinNormalizeTOFractionOF_Events_mvnH1D(*mcprojyArr);

		mnvPlotter.mc_line_width = 0;
		mnvPlotter.DrawDataStackedMC(hMCxproj, mcprojyArr, 1.0, "N", "", -1, -1, 1001);
		hMCxproj_clone->Draw("SAME E1 X0");

		drawBinRange(mc_clone, 2, i, yaxislabel, ".2f", false);

		//if (i == 1) mnvPlotter.WritePreliminary("TL", .015, -.01, -.06, false);
	}	/// End of Loop

	// gc->DrawBinRanges(dataHist, 2, bins, dataHist->GetYaxis()->GetTitle(), binCenter, .02);
	gc->SetTitleAlignmentFor6Hist();    
	gc->SetXLimits(MinX, MaxX);
	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.02);
	gc->SetYLimits(0, 1.245);
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(20);
	gc->SetXTitleSize(20);
	gc->SetYTitle(zaxislabel_units);
	gc->SetHistTexts();
	leg->Draw("SAME");
	mnvPlotter.AddHistoTitle(histotitle, .04);
	gc->Modified();
	gc->Print(pdf_label);

	delete gc;

}
*/
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

/*
 void PlotMC2D_Stack(
 std::map<Binning2D, TH1D*> InputHistmap_DATA,
 std::map<Binning2D, TH1D*> InputHistmap_DATA_BeamOff,
 std::map<Binning2D, TH1D*> InputHistmap_MC,
 std::map<std::pair<Binning2D, EventCategory>, TH1D*> StackMap_input,
 std::map<Binning2D , std::string > BinStringMap,
	char *pdf_label, char *histotitle, char *xaxislabel,
	char *legendTitle, char *zaxislabel_units,
	double Ymax, bool setMaxY, double Ymin,  bool doMultipliers,
	std::vector<double> YMultipliers,
	bool dontDo_bin_width_norm,
	double text_size,
	float POT_DATA, 
float POT_scaler_MC,
float BG_Trigger_scaler)
{

 std::cout<<"inside:PlotMC2D_Stack "<< std::endl;
	auto BinVector = GetProjectBinVector();
  auto& CategoryInterpreter = EventCategoryInterpreter::Instance();
  std::vector<EventCategory>  Category_vector =  CategoryInterpreter.ReturnCategoryVector();
	int nbins_X = InputHistmap_DATA.size()-1;
	int nbins_Y = InputHistmap_DATA.size()-1;
	//std::cout << "The NEntries for TObjArray: " << nVars << std::endl;

	 TLegend *leg = new TLegend(0.41, 0.1, 0.96, 0.25);
    std::string Leg_title_string = get_legend_title(POT_DATA);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.02);
	leg->SetNColumns(3);
  leg->SetHeader(Leg_title_string.c_str());
	

	double min_XAxis = InputHistmap_DATA[BinVector.at(0)]->GetXaxis()->GetXmin();
	double max_XAxis = InputHistmap_DATA[BinVector.at(0)]->GetXaxis()->GetBinUpEdge(InputHistmap_DATA[BinVector.at(0)]->GetNbinsX());
 std::cout<< "max_XAxis = "<< max_XAxis<< std::endl;
 double maxmax = -1; 





	int grid_x;
	int grid_y;
	int Nbins;
	std::string xaxislabel_string(xaxislabel);

	grid_x = sqrt(nbins_X + 1);
	grid_y = nbins_X / (grid_x - 1);
	Nbins = InputHistmap_DATA.size();

	 std::cout << nbins_X - grid_x * grid_y << std::endl;
	 std::cout << "Plotting plotYAxis1D with a grid of " << nbins_X << "\t" << nbins_Y << "\t" << grid_x << "\t" << grid_y << std::endl;

	GridCanvas *gc = new GridCanvas(uniq(), grid_x, grid_y, 800, 550);
	TCanvas cE("c1", "c1");

	std::cout<<"starting Loop "<< std::endl;

	for (int i = 1; i <= Nbins; ++i)
	{
	
	   std::cout<<"Bin index: "<< i << " BinVector.at(i-1) = "<< BinVector.at(i-1) << std::endl; 
	
		TH1D *h_Data_BeamOn = (TH1D *)InputHistmap_DATA[BinVector.at(i-1)]->Clone(uniq());
		std::cout<<" 1"<< std::endl; 
		h_Data_BeamOn->SetTitle(histotitle);
		TH1D *h_Data_BeamOFF = (TH1D *)InputHistmap_DATA_BeamOff[BinVector.at(i-1)]->Clone(uniq());
		std::cout<<"2 "<< std::endl; 
		TH1D *h_MC_Total = (TH1D *)InputHistmap_MC[BinVector.at(i-1)]->Clone(uniq());
		
		 std::cout<<"Made clones "<< std::endl; 
		
		h_MC_Total->Scale(POT_scaler_MC); 
    h_Data_BeamOFF->Scale(BG_Trigger_scaler); 
		
		 THStack* EventCategory_stack = new THStack( "mc+ext", "" );
		
		 
		
		  CategoryInterpreter.set_ext_histogram_style(h_Data_BeamOFF);
      CategoryInterpreter.set_bnb_data_histogram_style( h_Data_BeamOn );
		
		//h_MC_Total->Add(h_Data_BeamOFF);
		
		std::map<EventCategory, TH1D*> h_p_EventCategory_StackMap;
		//std::cout<<" about to do this loop : Category_vector "<< std::endl; 
		for(auto cat: Category_vector){
		  	//std::cout<<" Getting (,) =  ("<< BinVector.at(i-1)<< " , "<< cat<<std::endl; 
		std::string leg_stack_check = CategoryInterpreter.label( cat );
     	//std::cout<<" Getting (,) =  ("<< BinVector.at(i-1)<< " , "<< cat<<std::endl; 
       TH1D* hist = (TH1D*) StackMap_input[{BinVector.at(i-1),cat}]->Clone(uniq());
       //hist->SetTitle(map.second->GetTitle());
       h_p_EventCategory_StackMap.insert(std::pair<EventCategory, TH1D*>(cat,hist)); 
     }
		
	  for(auto &HistInMap : h_p_EventCategory_StackMap){HistInMap.second->Scale(POT_scaler_MC); }

		
		/////////////////
		/// Want Beam Off First 
		////////////////
		//EventCategory_stack->Add( h_Data_BeamOFF ); 
		
	for (auto it = h_p_EventCategory_StackMap.rbegin(); it != h_p_EventCategory_StackMap.rend(); ++it) {
    // 'it' is a reverse iterator pointing to the current element in reverse order
    auto& HistInMap = *it;

    CategoryInterpreter.set_mc_histogram_style( HistInMap.first, HistInMap.second );
    EventCategory_stack->Add( HistInMap.second );
 }
		
		
		 if(i ==1){
		for(auto &HistInMap : h_p_EventCategory_StackMap){
     std::string leg_stack_label = CategoryInterpreter.label( HistInMap.first );
    leg->AddEntry(HistInMap.second, leg_stack_label.c_str() , "f" );
    }
    
  }
		
		
		if (i == 1) {
		leg->AddEntry(h_Data_BeamOn, "Data", "pe" );
    leg->AddEntry(h_MC_Total, "Total MC", "le" );
		leg->AddEntry(h_Data_BeamOFF, "bnbext", "f");
		
		}
		
		
		
	//////////////////////////////	
		gc->cd(i);
	/////////////////////////////


		//if (dontDo_bin_width_norm == true) {
		//h_Data->Scale(1., "width");
		//h_DataBeamOff->Scale(1., "width");
		//h_mc->Scale(1., "width");
		//}
		
		if (doMultipliers == true) {
		h_Data_BeamOn ->Scale(YMultipliers.at(i - 1));
		scaleTHStack(EventCategory_stack, YMultipliers.at(i - 1));
		//EventCategory_stack->Scale(YMultipliers.at(i - 1));
		
		h_MC_Total->Scale(YMultipliers.at(i - 1));
		}
		
		
		
		DrawStackMCandData(
    h_Data_BeamOn, 
    h_MC_Total,
    	EventCategory_stack,
    dontDo_bin_width_norm, 
    Ymax);
   
   
   
   
   drawString(BinStringMap[BinVector.at(i-1)], text_size, false );
		
		if (doMultipliers && YMultipliers.at(i - 1) != 1)
		{
			auto pad = gc->cd(i);
			TLatex *la2 = new TLatex(1 - pad->GetRightMargin() - 0.02,
				1 - pad->GetTopMargin() - 0.05,
				TString::Format("#times %.1f", YMultipliers.at(i - 1)));
			la2->SetTextAlign(33);	// top right
			la2->SetNDC();
			la2->SetTextFont(42);
			la2->SetTextSize(0.03);
			la2->Draw();
		}

 if (i == 1){DrawFakeData();}


		//mcprojyArr->Clear();
	}	/// End of Loop

    std::cout<<"Finished 2D drawing loop "<< std::endl;

	gc->SetYLabel_Size(.015);
	gc->SetXLabel_Size(.015);


	
	gc->SetInterpadSpace(.005);
	//gc->SetRightMargin(0.05);
	//gc->SetLeftMargin(0.08);
  gc->SetYLimits(Ymin, Ymax);
  gc->SetXLimits(min_XAxis,max_XAxis);
  gc->ResetPads();
	gc->SetXTitle(xaxislabel);
	gc->SetYTitleSize(25);
	gc->SetXTitleSize(20);  
	gc->SetYTitle(zaxislabel_units);
	gc->SetTitleAlignmentFor6Hist();
	leg->Draw("SAME");
	gc->Modified();
	gc->Print(pdf_label);
	//leg->Clear();
	delete gc;



}
*/

//////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////
void scaleTHStack(THStack* stack, double scaleFactor) {
    if (!stack) {
        std::cerr << "Error: Invalid THStack pointer." << std::endl;
        return;
    }

    // Get the list of histograms from the THStack
    TList* histograms = stack->GetHists();

    if (!histograms) {
        std::cerr << "Error: No histograms found in THStack." << std::endl;
        return;
    }

    // Iterate through the histograms and scale each one
    TIter next(histograms);
    TObject* obj;

    while ((obj = next())) {
        if (obj->IsA()->InheritsFrom("TH1")) {
            TH1* hist = dynamic_cast<TH1*>(obj);
            hist->Scale(scaleFactor);
        }
    }
}

//////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////

void DrawFakeData(){
    TLatex* textfake = new TLatex;
    textfake->SetNDC();
    textfake->SetTextSize(0.03);
    textfake->SetTextColor(kBlue);
    textfake->DrawLatex(0.12, 0.86, "USING FAKE DATA");

}


void DrawAnnie_MC(){
    TLatex* textfake = new TLatex;
    textfake->SetNDC();
    textfake->SetTextSize(0.03);
    textfake->SetTextColor(kBlue);
    textfake->DrawLatex(0.12, 0.86, "preliminary ANNIE");

}
void DrawTGraph(
  TGraph *g_TGraph1,
  const char* xaxislabel, const char* yaxislabel,
  const char* Title, const char* legend_Title1,
 const char* pdf,
  TCanvas *can,
  bool MakeXaxisLOG, bool MakeYaxisLOG, bool doMax )
  {
     TLegend *legend = new TLegend (.1 , .8, .2 , .9 );

     if(MakeXaxisLOG==true){
       gPad->SetLogx();
     }
     if(MakeYaxisLOG==true){
       gPad->SetLogy();
     }

     can->SetGrid();

    std::string TotalTitle(Title);

    g_TGraph1 -> SetTitle("");
    g_TGraph1 -> GetXaxis() -> SetTitle(xaxislabel);
    g_TGraph1 -> GetYaxis() -> SetTitle(yaxislabel);
    g_TGraph1 -> GetXaxis() -> CenterTitle();
    g_TGraph1 -> GetYaxis() -> CenterTitle();
    g_TGraph1 -> GetXaxis() -> SetTitleSize(0.03);
    g_TGraph1 -> GetYaxis() -> SetTitleSize(0.03);
    g_TGraph1->GetYaxis()->SetLabelSize(0.025);
    g_TGraph1->GetXaxis()->SetLabelSize(0.025);

  
    g_TGraph1 -> SetMarkerColor(4);
    g_TGraph1-> SetMarkerStyle(6);
    
    
    if(doMax==true){
      double SetMaxY = g_TGraph1->GetHistogram()->GetMaximum();
      double SetMinY = g_TGraph1->GetHistogram()->GetMinimum();
      g_TGraph1->SetMaximum(SetMaxY*1.1);
      g_TGraph1->SetMinimum(SetMinY*1.1);
    
    }

    g_TGraph1->SetMarkerStyle(20); // Set marker style
    g_TGraph1->SetMarkerSize(1.0); // Set marker size
    g_TGraph1->SetMarkerColor(kBlue); // Set marker color
    g_TGraph1 -> Draw("AP");
   

   legend -> SetFillColor(0);
   legend -> AddEntry(g_TGraph1, legend_Title1, "p");
   legend -> Draw();


   gPad->Update();
   can->Modified();
   can->Print(pdf);
   can->Closed();

   if(MakeXaxisLOG==true){
     gPad->SetLogx(0);
   }
   if(MakeYaxisLOG==true){
     gPad->SetLogy(0);
   }

   can->SetGrid(0,0);


    return;

}



void DrawTGraph(
  TGraph *g_TGraph1, TGraph *g_TGraph2,
  const char* xaxislabel, const char* yaxislabel,
  const char* Title, const char* legend_Title1,
  const char* legend_Title2, const char* pdf,
  TCanvas *can,
  bool MakeXaxisLOG, bool MakeYaxisLOG, bool doMax )
  {
     TLegend *legend = new TLegend (.4 , .5, .8 , .8 );

     if(MakeXaxisLOG==true){
       gPad->SetLogx();
     }
     if(MakeYaxisLOG==true){
       gPad->SetLogy();
     }

     can->SetGrid();

    std::string TotalTitle(Title);

    g_TGraph2 -> SetTitle("");
    g_TGraph2 -> GetXaxis() -> SetTitle(xaxislabel);
    g_TGraph2 -> GetYaxis() -> SetTitle(yaxislabel);
    g_TGraph2 -> GetXaxis() -> CenterTitle();
    g_TGraph2 -> GetYaxis() -> CenterTitle();
    g_TGraph2 -> GetXaxis() -> SetTitleSize(0.03);
    g_TGraph2 -> GetYaxis() -> SetTitleSize(0.03);
    g_TGraph2->GetYaxis()->SetLabelSize(0.025);
    g_TGraph2->GetXaxis()->SetLabelSize(0.025);

    g_TGraph1 -> SetLineColor(2);
    g_TGraph1 -> SetMarkerColor(2);
    g_TGraph1-> SetMarkerStyle(6);
    g_TGraph1->SetLineStyle(2);
    g_TGraph1->SetLineWidth(2);
    //g_TGraph1->SetMarkerSize(5);
    g_TGraph2 -> SetLineColor(4);
    g_TGraph2 -> SetMarkerColor(4);
    g_TGraph2-> SetMarkerStyle(6);
    g_TGraph2->SetLineStyle(2);
    g_TGraph2->SetLineWidth(2);
    //g_TGraph2->SetMarkerSize(5);

    double SetMaxY, SetMinY;
    if(doMax==true){
      double Amax = g_TGraph1->GetHistogram()->GetMaximum();
      double Bmax = g_TGraph2->GetHistogram()->GetMaximum();
      if(Amax  >= Bmax){SetMaxY = Amax;}
      else{SetMaxY = Bmax;}
      g_TGraph2->SetMaximum(SetMaxY*1.1);

      double Amin = g_TGraph1->GetHistogram()->GetMinimum();
      double Bmin = g_TGraph2->GetHistogram()->GetMinimum();

      if(Amin  <= Bmin){SetMinY = Amin;}
      else{SetMinY = Bmin;}
      g_TGraph2->SetMinimum(SetMinY*.9);
    }

    g_TGraph1->SetFillColor(42);
    g_TGraph2->SetFillColor(38);
    g_TGraph2 -> Draw("alf");
    g_TGraph1 -> Draw("lf");



   legend -> SetTextSize(0.035);
   legend -> SetFillColor(0);
   legend -> AddEntry(g_TGraph1, legend_Title1, "p");
   legend -> AddEntry(g_TGraph2, legend_Title2, "p");
   legend -> Draw();


   gPad->Update();
   can->Modified();
   can->Print(pdf);
   can->Closed();

   if(MakeXaxisLOG==true){
     gPad->SetLogx(0);
   }
   if(MakeYaxisLOG==true){
     gPad->SetLogy(0);
   }

   can->SetGrid(0,0);


    return;

}

/////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////
void DrawTGraph(
  TGraph *g_TGraph1, TGraph *g_TGraph2, TGraph *g_TGraph3,
  const char* xaxislabel, const char* yaxislabel,
  const char* Title, const char* legend_Title1,
  const char* legend_Title2,   
  const char* legend_Title3, const char* pdf,
  TCanvas *can,
  bool MakeXaxisLOG, bool MakeYaxisLOG, bool doMax )
  {
     TLegend *legend = new TLegend (.1 , .75, .3 , .9 );

     if(MakeXaxisLOG==true){
       gPad->SetLogx();
     }
     if(MakeYaxisLOG==true){
       gPad->SetLogy();
     }

     can->SetGrid();

     std::string TotalTitle(Title);

    g_TGraph1 -> SetTitle("");
    g_TGraph1 -> GetXaxis() -> SetTitle(xaxislabel);
    g_TGraph1 -> GetYaxis() -> SetTitle(yaxislabel);
    g_TGraph1 -> GetXaxis() -> CenterTitle();
    g_TGraph1 -> GetYaxis() -> CenterTitle();
    g_TGraph1 -> GetXaxis() -> SetTitleSize(0.03);
    g_TGraph1 -> GetYaxis() -> SetTitleSize(0.03);
    g_TGraph1->GetYaxis()->SetLabelSize(0.025);
    g_TGraph1->GetXaxis()->SetLabelSize(0.025);


    g_TGraph1 -> SetMarkerColor(2);
    g_TGraph1-> SetMarkerStyle(20);
    g_TGraph1->SetMarkerSize(.8); // Set marker size

    //g_TGraph1->SetMarkerSize(5);

    g_TGraph2 -> SetMarkerColor(4);
    g_TGraph2-> SetMarkerStyle(20);
    g_TGraph3->SetMarkerSize(.8); // Set marker size

    
    g_TGraph3->SetMarkerColor(1);
    g_TGraph3 ->SetMarkerStyle(20);
    g_TGraph3->SetMarkerSize(.8); // Set marker size
    //g_TGraph2->SetMarkerSize(5);

    double SetMaxY, SetMinY;
    if(doMax==true){
      double Amax = g_TGraph1->GetHistogram()->GetMaximum();
      double Bmax = g_TGraph2->GetHistogram()->GetMaximum();
      if(Amax  >= Bmax){SetMaxY = Amax;}
      else{SetMaxY = Bmax;}
      g_TGraph2->SetMaximum(SetMaxY*1.35);

      double Amin = g_TGraph1->GetHistogram()->GetMinimum();
      double Bmin = g_TGraph2->GetHistogram()->GetMinimum();

      if(Amin  <= Bmin){SetMinY = Amin;}
      else{SetMinY = Bmin;}
      g_TGraph2->SetMinimum(SetMinY*.9);
    }

    //g_TGraph1->SetFillColor(42);
    //g_TGraph2->SetFillColor(38);
    g_TGraph1->Draw("AP");
    g_TGraph2->Draw("P SAME");
    g_TGraph3->Draw("P SAME");


   legend -> SetTextSize(0.035);
   legend -> SetFillColor(0);
   legend -> AddEntry(g_TGraph1, legend_Title1, "p");
   legend -> AddEntry(g_TGraph2, legend_Title2, "p");
   legend -> AddEntry(g_TGraph3, legend_Title3, "p");
   legend -> Draw();


   gPad->Update();
   can->Modified();
   can->Print(pdf);
   can->Closed();

   if(MakeXaxisLOG==true){
     gPad->SetLogx(0);
   }
   if(MakeYaxisLOG==true){
     gPad->SetLogy(0);
   }

   can->SetGrid(0,0);


    return;

}

/////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////
void Draw_HIST(
  TH1D *hist_1_input,
   char *histotitle,
   std::string xaxislabel,
   std::string yaxislabel,
   bool NormArea, 
   bool Setgrid,
   bool BinWidthNorm,
   double Ymax,
   TCanvas *cE)
{

  char textplace[1024];
  TH1D *hist_1 = (TH1D *)hist_1_input->Clone("");
  

  if (Setgrid==true) cE->SetGrid();

  //histHelium->GetXaxis()->SetTitleSize(0.035);
  hist_1->GetXaxis()->SetTitle(Form("%s",xaxislabel.c_str()));
  hist_1->GetYaxis()->SetTitle(Form("%s",yaxislabel.c_str()));
  hist_1->SetTitle(histotitle);
  hist_1->SetTitleOffset(1.2);

   if(NormArea==true){
   hist_1->Scale((1.0 / (hist_1->Integral())));
   }


  if(BinWidthNorm== true){
    hist_1->Scale(1.0,"width");
  }

  if(Ymax != -99){
    hist_1->SetMaximum(Ymax);
  }  
   else{
        double max1 = hist_1->GetMaximum();
         hist_1->SetMaximum(max1* 1.25);
       }
       
  int nbins_1 = hist_1->GetNbinsX();
  double Mean = hist_1->GetMean(1);
  double stdDev = hist_1->GetStdDev(1);
  
  hist_1->GetXaxis()->SetTitleSize(0.03);
  hist_1->SetLineColor(kAzure - 1);
  hist_1->Draw("hist e");
  

}//////////////////

/////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////
void Draw_HIST(
   TH1D *hist_1_input,
  char *legend_Title1,
   TH1D *hist_2_input,
   char *legend_Title2,
   char *histotitle,
   std::string xaxislabel,
   std::string yaxislabel,
   bool NormArea, 
   bool Setgrid,
   bool BinWidthNorm,
   double Ymax,
   TCanvas *cE,
   std::string pdf_name )
{

  char textplace[1024];
  TH1D *hist_1 = (TH1D *)hist_1_input->Clone("");
  TH1D *hist_2 = (TH1D *)hist_2_input->Clone("");

   TLegend *legend = new TLegend (.1 , .75, .3 , .9 );
  



  if (Setgrid==true) cE->SetGrid();

  //histHelium->GetXaxis()->SetTitleSize(0.035);
  hist_1->GetXaxis()->SetTitle(Form("%s",xaxislabel.c_str()));
  hist_1->GetYaxis()->SetTitle(Form("%s",yaxislabel.c_str()));
  hist_1->SetTitle(histotitle);
  hist_1->SetTitleOffset(1.2);

   if(NormArea==true){
   hist_1->Scale((1.0 / (hist_1->Integral())));
   hist_2->Scale((1.0 / (hist_1->Integral())));
   }




  if(BinWidthNorm== true){
    hist_1->Scale(1.0,"width");
    hist_2->Scale(1.0,"width");
  }

  if(Ymax != -99){
    hist_1->SetMaximum(Ymax);
  }  
   else{
        double max1 = hist_1->GetMaximum();
        double max2 = hist_2->GetMaximum();
        double Max = (max1 > max2 ) ? max1 : max2; 
        
        
         hist_1->SetMaximum(Max* 1.4);
       }
       
       
  int nbins_1 = hist_1->GetNbinsX();
  double Mean = hist_1->GetMean(1);
  double stdDev = hist_1->GetStdDev(1);
  
  hist_1->SetLineColor(kAzure - 1);
  hist_2->SetLineColor(kRed);
  
   double area_1 = hist_1->Integral();
   double area_2 = hist_2->Integral();
   
   
   legend -> AddEntry(hist_1, legend_Title1, "l");
   legend -> AddEntry(hist_2, legend_Title2, "l");
  
  TLatex* text = new TLatex;
  text->SetNDC();
  text->SetTextSize(0.03);
  text->SetTextColor(kRed);

  
  
  
  hist_1->GetXaxis()->SetTitleSize(0.03);
  hist_1->SetMinimum(0.0);
  hist_1->Draw("hist e");
  hist_2->Draw("SAME hist e");
  legend->Draw("SAME");  

  sprintf(textplace, "area hist1= %.2f",area_1 );
  text->DrawLatex(0.7, 0.85, textplace);
  
  sprintf(textplace, "area hist2= %.2f",area_2 );
    text->DrawLatex(0.7, 0.8, textplace);
    
  sprintf(textplace, " hist1/ hist2= %.2f",area_1/area_2 );
    text->DrawLatex(0.7, 0.75, textplace);


    std::string plotname = Form("%s",pdf_name.c_str());
    cE-> Print(plotname.c_str());



}//////////////////
/////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////
void Draw_STACK_HIST(
   TH1D *hist_1_input,
  char *legend_Title1,
   TH1D *hist_2_input,
   char *legend_Title2,
   TH1D *hist_3_input,
   char *legend_Title3,
   char *histotitle,
   std::string xaxislabel,
   std::string yaxislabel,
   bool NormArea, 
   bool Setgrid,
   bool BinWidthNorm,
   double Ymax,
   TCanvas *cE,
   std::string pdf_name )
{

  std::cout << "BBBBBB" << std::endl;
  char textplace[1024];
  TH1D *hist_1_temp = (TH1D *)hist_1_input->Clone("");
  TH1D *hist_2_temp = (TH1D *)hist_2_input->Clone("");
  TH1D *hist_3_temp = (TH1D *)hist_3_input->Clone("");

  THStack* hs = new THStack("hs","Event Rates");
  std::vector<double> nbins = {600.,810.,910.,1005.,1100.,1200};

  std::cout << "CCCCC" << std::endl;
  TH1D *hist_1 = new TH1D("hist_1", ";E_{muon};", nbins.size() - 1, nbins.data());
  TH1D *hist_2 = new TH1D("hist_2", ";E_{muon};", nbins.size() - 1, nbins.data());
  TH1D *hist_3 = new TH1D("hist_3", ";E_{muon};", nbins.size() - 1, nbins.data());

  std::cout << "DDDDD" << std::endl;
  for(int i = 1; i < nbins.size(); i++){
    hist_1->SetBinContent(i, hist_1_temp->GetBinContent(i));
    hist_2->SetBinContent(i, hist_2_temp->GetBinContent(i));
    hist_3->SetBinContent(i, hist_3_temp->GetBinContent(i));
  }

  std::cout << "EEEEE" << std::endl;
  TLegend *legend = new TLegend (.1 , .75, .3 , .9 );
  
  if (Setgrid==true) cE->SetGrid();

  //histHelium->GetXaxis()->SetTitleSize(0.035);
  hist_1->GetXaxis()->SetTitle(Form("%s",xaxislabel.c_str()));
  hist_1->GetYaxis()->SetTitle(Form("%s",yaxislabel.c_str()));
  hist_1->SetTitle(histotitle);
  hist_1->SetTitleOffset(1.2);

  if(NormArea==true){
    hist_1->Scale((1.0 / (hist_1->Integral())));
    hist_2->Scale((1.0 / (hist_1->Integral())));
    hist_3->Scale((1.0 / (hist_1->Integral())));
  }


  if(BinWidthNorm== true){
    hist_1->Scale(1.0,"width");
    hist_2->Scale(1.0,"width");
    hist_3->Scale(1.0,"width");
  }

  std::cout << "FFFF" << std::endl;
  if(Ymax != -99){
    hs->SetMaximum(Ymax);
  }  
  else{
    double max1 = hist_1->GetMaximum();
    double max2 = hist_2->GetMaximum();
    double Max = (max1 > max2 ) ? max1 : max2; 
    
    hs->SetMaximum(Max* 1.4);
  }
       
  std::cout << "GGGG" << std::endl;     
  int nbins_1 = hist_1->GetNbinsX();
  double Mean = hist_1->GetMean(1);
  double stdDev = hist_1->GetStdDev(1);
  
  hist_1->SetLineColor(kAzure - 1);
  hist_2->SetLineColor(kRed);
  hist_3->SetLineColor(kGreen);
  
  double area_1 = hist_1->Integral();
  double area_2 = hist_2->Integral();
   
   std::cout << "HHHHHH" << std::endl;
  legend -> AddEntry(hist_1, legend_Title1, "l");
  legend -> AddEntry(hist_2, legend_Title2, "l");
  legend -> AddEntry(hist_3, legend_Title3, "l");

  TLatex* text = new TLatex;
  text->SetNDC();
  text->SetTextSize(0.03);
  text->SetTextColor(kRed);

  std::cout << "IIIIII" << std::endl;
//  hs->GetXaxis()->SetTitleSize(0.03);
  std::cout << "WWWWWW" << std::endl;
  hs->SetMinimum(0.0);

  std::cout << "JJJJJJ" << std::endl;
  hs->Add(hist_1);
  hs->Add(hist_2);
  hs->Add(hist_3);

  std::cout << "ZZZZZ" << std::endl;
  hs->Draw("nostack e");
//  hist_2->Draw("SAME hist e");
  legend->Draw("SAME");  

  sprintf(textplace, "area hist1= %.2f",area_1 );
  text->DrawLatex(0.7, 0.85, textplace);
  
  sprintf(textplace, "area hist2= %.2f",area_2 );
    text->DrawLatex(0.7, 0.8, textplace);
    
  sprintf(textplace, " hist1/ hist2= %.2f",area_1/area_2 );
    text->DrawLatex(0.7, 0.75, textplace);

  std::string plotname = Form("%s",pdf_name.c_str());
  cE-> Print(plotname.c_str());

}//////////////////
