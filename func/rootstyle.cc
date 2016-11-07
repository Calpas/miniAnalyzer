void rootstyle(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetTitleFont(102);
  gStyle->SetTitleFont(102,"xyz");
  gStyle->SetStatFont(102);
  gStyle->SetOptFit(111);
  gStyle->SetOptTitle(0);
  gStyle->SetPaintTextFormat(".3f");

  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);  
  gStyle->SetOptStat(0); //don't show statistics box
  gStyle->SetTitleSize(0.05, "xyz"); 
  gStyle->SetLegendBorderSize(1);

  gStyle->SetPadTopMargin(0.1);
  //gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadBottomMargin(0.1);
  //gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);

  gROOT->ForceStyle();
  TH1::SetDefaultSumw2(kTRUE);	
}
