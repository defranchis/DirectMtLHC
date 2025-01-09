#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TMatrixDSymEigen.h"
#include "TMath.h"
#include "TF2.h"
#include "TStyle.h"
#include "TVectorD.h"
#include "RooMultiVarGaussian.h"
#include <iostream>
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/LikelihoodInterval.h"
//#include "AtlasLabels.C"

using namespace RooStats;

using namespace std;
  

// small function to test if point is inside an ellipse that is centered at x0, y0
// with major and minor radii a and b; angle theta (degrees)
bool inEllipse(double x, double y,
                     double x0, double y0,
                     double a, double b,
                     double theta) // (in degrees)
{
  // shift the center
  x -= x0;
  y -= y0;
  // un-rotate the axes
  theta *= TMath::Pi() / 180.0; // degrees -> radians
  double v = x;
  x = x * std::cos(theta) + y * std::sin(theta);
  y = y * std::cos(theta) - v * std::sin(theta);
  bool inside = (x*x/(a*a) + y*y/(b*b)) <= 1;
  return inside;
}

void Exp2DPlot() {

  const double lhc = 172.515;
  const double lhc_err = 0.329;
  
  const double twoD_atlas_err = 0.468;
  const double twoD_cms_err = 0.406;
  const double rho = 0.153;
  const double twoD_atlas = 172.723;
  const double twoD_cms = 172.367;
  
  const double oneD_atlas_err = 0.482;
  const double oneD_cms_err = 0.415;
  const double oneD_atlas = 172.706;
  const double oneD_cms = 172.522;
  
  const int defaultFont = 42;
  const double text_size = gStyle->GetTextSize();
  const double ymax = 175.3;
  const double xmax = ymax;
  const double ymin = 170.1;
  const double xmin = ymin;
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetTitleOffset(1.8,"Y");
  TH1D::SetDefaultSumw2(true);
  
  RooRealVar atlas("ATLAS","m_{t}^{ATLAS} [GeV]",twoD_atlas,150,200);
  RooRealVar cms("CMS","m_{t}^{CMS} [GeV]",twoD_cms,150,200);

  RooRealVar atlas_mu("ATLAS_mu","ATLAS_mu",twoD_atlas, 150, 200);
  atlas_mu.setConstant(true);
  RooRealVar cms_mu("CMS_mu","CMS_mu",twoD_cms, 150, 200);
  cms_mu.setConstant(true);


  
  TGraphErrors* atlas_cms = new TGraphErrors(1);
  atlas_cms->SetPoint(0, atlas_mu.getVal(), cms_mu.getVal());
  atlas_cms->SetPointError(0, twoD_atlas_err, twoD_cms_err);
  atlas_cms->SetMarkerStyle(29);

  TGraph* g_lhc = new TGraph(1);
  g_lhc->SetPoint(0, lhc, lhc);
  g_lhc->SetMarkerStyle(8);
  TGraphErrors* g_lhc_err = new TGraphErrors(1);
  g_lhc_err->SetPoint(0, lhc, lhc);
  g_lhc_err->SetPointError(0, lhc_err, lhc_err);
  g_lhc_err->SetMarkerStyle(8);
  g_lhc_err->SetLineWidth(3);
  
  
  TMatrixDSym cov(2);
  cov[0][0] = twoD_atlas_err*twoD_atlas_err;
  cov[1][1] = twoD_cms_err*twoD_cms_err;
  cov[1][0] = cov[0][1] = rho*twoD_atlas_err*twoD_cms_err;
  cout << cov[0][0] << " " << cov[0][1] << endl;
  cout << cov[1][0] << " " << cov[1][1] << endl;

  RooMultiVarGaussian *pdf = new RooMultiVarGaussian("pdf","", RooArgList(atlas, cms), RooArgList(atlas_mu, cms_mu), cov);

  TCanvas* c = new TCanvas("c","c",500,500);
  //TH2* pdf_2d = (TH2*) pdf->createHistogram("ATLAS,CMS",200,200) ;
  TH2* pdf_2d = (TH2*)pdf->createHistogram("ATLAS",atlas, RooFit::Binning(100,xmin,xmax), RooFit::YVar(cms,RooFit::Binning(100,ymin,ymax))) ;
  pdf_2d->SetMarkerStyle(21);
  pdf_2d->SetMarkerSize(0.5);
  pdf_2d->SetContour(50);
  pdf_2d->Draw("colz");
  g_lhc->Draw("p");
  c->Print("2dv1.pdf");
  
  TMatrixDSymEigen eigen(cov);
  const TVectorD eVals = eigen.GetEigenValues();
  const TMatrixD eVecs = eigen.GetEigenVectors();
  cout << "Eigen vals: " << eVals[0] << " " << eVals[1] << endl;
  cout << "sqrt(Eigen vals): " << sqrt(eVals[0]) << " " << sqrt(eVals[1]) << endl;
  cout << eVecs.GetNcols () << " " << eVecs.GetNrows() << endl;
  cout << "EV0 : " << eVecs[0][0] << " " << eVecs[0][1] << endl;
  cout << "EV1 : " << eVecs[1][0] << " " << eVecs[1][1] << endl;

  const double angleRad = TMath::ATan(eVecs[0][1] / eVecs[0][0]);
  const double angleDeg = 180.0*angleRad/TMath::Pi();
  std::cout << angleDeg << std::endl;
  
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->cd();
  pdf_2d->Draw("axis");
  pdf_2d->GetYaxis()->SetNdivisions(9, 5, 0);
  pdf_2d->GetXaxis()->SetNdivisions(9, 5, 0);

  // Draw ATLAS and CMS individual combinations, we want bands for this
  TBox* box_atlasI = new TBox(oneD_atlas-oneD_atlas_err, ymin, oneD_atlas+oneD_atlas_err, ymax);
  box_atlasI->SetFillColor(kBlue);
  box_atlasI->SetFillStyle(3005);
  box_atlasI->Draw();
  TLine *line_atlasI = new TLine(oneD_atlas, ymin, oneD_atlas, ymax);
  line_atlasI->SetLineStyle(2);
  line_atlasI->SetLineColor(kBlue);
  line_atlasI->Draw();

  TBox* box_cmsI = new TBox(xmin, oneD_cms-oneD_cms_err, xmax, oneD_cms+oneD_cms_err);
  box_cmsI->SetFillColor(kRed);
  box_cmsI->SetFillStyle(3005);
  box_cmsI->Draw();
  TLine *line_cmsI = new TLine(xmin, oneD_cms, xmax, oneD_cms);
  line_cmsI->SetLineStyle(2);
  line_cmsI->SetLineColor(kRed);
  line_cmsI->Draw();
  
  // elipse for 95%
  TEllipse* el95 = new TEllipse(atlas_mu.getVal(), cms_mu.getVal(), sqrt(5.991*eVals[0]), sqrt(5.991*eVals[1]), 0.0, 360.0, angleDeg);
  el95->SetFillColor(kOrange-2);
  el95->SetLineWidth(0);
  el95->SetLineColor(kOrange-2);
  el95->Draw();
  //f95->Draw("same");
    
  // check we got the right elipse
  RooDataSet *toys = pdf->generate(RooArgSet(atlas,cms),10000) ;  

  TH1D* h_atlas = new TH1D("h_atlas","h_atlas",100,170.0,175.0);
  TGraph* g_inel = new TGraph();
  int ip=0;
  double total = 0;
  double inelipse = 0;
  const double t=1;
  for(int i=0; i<toys->numEntries(); ++i) {
    auto thisdata = toys->get(i);
    double x = ((RooRealVar*)thisdata->find("ATLAS"))->getVal();
    double y = ((RooRealVar*)thisdata->find("CMS"))->getVal();

    // test if inside elipse
    if( inEllipse(x, y, atlas_mu.getVal(), cms_mu.getVal(), sqrt(t)*sqrt(2.3*eVals[0]), sqrt(t)*sqrt(2.3*eVals[1]), angleDeg) ) {
      inelipse += 1;
      g_inel->SetPoint(ip++, x, y);
    }
    
    h_atlas->Fill(x);
    //if(dist==0) {
      //cout << "Dist 0 " << x << " " << y << endl;
    //  inelipse += 1;
    //}
    total += 1;
	//cout << x << " " << y << " " << dist << endl;
  }
  cout << t << " " << 100.0*inelipse/total << "%" << endl;
  //g_inel->Draw("p same");
  

  
  
  // elipse for 68%
  TEllipse* el68 = new TEllipse(atlas_mu.getVal(), cms_mu.getVal(), sqrt(2.3*eVals[0]), sqrt(2.3*eVals[1]), 0.0, 360.0, angleDeg);
  el68->SetFillColor(kOrange+7);
  el68->SetLineWidth(0);
  el68->SetLineColor(kOrange+7);
  el68->Draw();
  atlas_cms->Draw("pX");

  // Could also have made elipses by fitting the data in RooFit
  /***********************
  RooDataSet* theData = new RooDataSet("theData","theData",RooArgSet(atlas, cms));
  theData->add(RooArgSet(atlas, cms));
  atlas_mu.setConstant(false);
  cms_mu.setConstant(false);
  ProfileLikelihoodCalculator* plc = new ProfileLikelihoodCalculator(*theData, *pdf, RooArgSet(atlas_mu, cms_mu));
  plc->SetConfidenceLevel(0.68);
  LikelihoodInterval *plInt = (LikelihoodInterval *)plc->GetInterval();
  LikelihoodIntervalPlot *plPlot = new LikelihoodIntervalPlot(plInt);
  plPlot->Draw("same");
  ************************/ // end RooFit plotting code

  // Draw x=y
  TF1* f= new TF1("f","x",150,200.0);
  f->SetLineStyle(2);
  f->Draw("same");

  // Draw LHC, error bar is custom
  g_lhc->Draw("p");
  TLine* lhc_errline = new TLine(lhc - lhc_err, lhc - lhc_err, lhc + lhc_err, lhc + lhc_err);
  lhc_errline->SetLineWidth(3);
  lhc_errline->Draw();
  const double tickshift = 0.05;
  TLine* lhc_tick1 = new TLine(lhc - lhc_err - tickshift, lhc - lhc_err + tickshift, lhc - lhc_err + tickshift, lhc - lhc_err - tickshift);
  lhc_tick1->SetLineWidth(3);
  lhc_tick1->Draw();
  TLine* lhc_tick2 = new TLine(lhc + lhc_err - tickshift, lhc + lhc_err + tickshift, lhc + lhc_err + tickshift, lhc + lhc_err - tickshift);
  lhc_tick2->SetLineWidth(3);
  lhc_tick2->Draw();
  
  // redraw some things
  line_atlasI->Draw();
  line_cmsI->Draw();
  
  // Labels
  //TPave* whitebox = new TPave(0.3,0.86,0.7,0.91,0,"NB NDC");
  //whitebox->SetFillColor(kWhite);
  //whitebox->Draw();
  TLatex *LHCLabel1 = new TLatex();
  LHCLabel1->SetNDC(true);
  LHCLabel1->SetTextSize(0.9*text_size);
  LHCLabel1->SetTextColor(1);
  LHCLabel1->SetTextFont(62);
  LHCLabel1->DrawLatex(0.21,0.87, "ATLAS+CMS");
  LHCLabel1->SetTextFont(defaultFont);
  LHCLabel1->SetTextSize(0.6*text_size);
  LHCLabel1->DrawLatex(0.22, 0.64,"LHC#font[52]{#scale[1.2]{top}}WG");
  LHCLabel1->SetTextSize(0.7*text_size);
  LHCLabel1->DrawLatex(0.68,0.87, "#sqrt{s} = 7,8 TeV");

  TLegend* leg = new TLegend(0.21,0.67,0.59,0.86);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.12);
  leg->SetTextSize(0.03);
  leg->SetTextFont(defaultFont);
  leg->AddEntry(atlas_cms,"#splitline{Simultaneous}{combination}","p");
  //leg->AddEntry((TObject*)0,"combination","");
  leg->AddEntry(el68,"68% CL","f");
  leg->AddEntry(el95,"95% CL","f");
  leg->Draw("same");

  TLegend* leg2 = new TLegend(0.65, 0.2, 0.93, 0.4);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  auto g_atlas = new TGraph();
  g_atlas->SetLineColor(kBlue);
  g_atlas->SetFillColor(kBlue);
  g_atlas->SetLineStyle(2);
  g_atlas->SetFillStyle(3005);
  auto g_cms = new TGraph(*g_atlas);
  g_cms->SetLineColor(kRed);
  g_cms->SetFillColor(kRed);
  leg2->SetTextFont(defaultFont);
  leg2->AddEntry(g_atlas, "ATLAS","fl");
  leg2->AddEntry(g_cms, "CMS","fl");
  leg2->AddEntry(g_lhc_err,"ATLAS+CMS","pl");
  leg2->AddEntry(f,"m_{t}^{LHC} = m_{t}^{ATLAS} = m_{t}^{CMS}","l");
  leg2->Draw("same");
  gPad->RedrawAxis();
  
  c2->Print("2dv2.pdf");

  TCanvas* c_atlas = new TCanvas("c_atlas","c_atlas");
  h_atlas->Draw("p");
  h_atlas->Fit("gaus");
}
