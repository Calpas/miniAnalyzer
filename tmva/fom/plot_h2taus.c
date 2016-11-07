#include "../../func/rootstyle.cc"
#include "../../func/myfunc.cc"
#include <array>
#include <math.h>


bool bdttype[]{0, 0, 1}; // s sv b; s+dy vs wtt, s vs dy;


void plot_h2taus(){

  // pad,axis... style
  rootstyle();

  // input file
  vector<string> vinputfile {
    "GluGluHToTauTauM125.root",
      "VBFHToTauTauM125.root",
      "DYJetsM50.root",
      "WJets.root",
      "TT.root"
  };


  // tree directory
  vector<string> vdir{"mt" ,"et", "em", "tt"};
  vector<string> vdir_{"#mu#tau" ,"e#tau", "e#mu", "#tau#tau"};
  //vector<string>vdir{"et"};
  int nbdir {static_cast<int>(vdir.size())};

  // 1 canvas and 1 pad for each channel, the last is for the legend
  TCanvas *canvas = new TCanvas("fom", "", 1250, 700) ;
  canvas->Divide(4, 3);

  vector<string> vjet{"0j", "1j", "01j"};
  int nbjet {static_cast<int>(vjet.size())};

  TH1D *hfom[nbdir][nbjet];

  vector<string> vcat{"low_pt_t", "high_pt_t", "high_pt_t_boosted", "boosted_tt", "highly_boosted_tt"};
  int nbcat {static_cast<int>(vcat.size())};

  TF1  *fnc[nbdir][nbjet-1][nbcat];

  // legend
  TLegend *legend = new TLegend( 0.15, 0.1, 0.95, 0.9, NULL,"brNDC");
  legend ->SetHeader("FOM: #frac{s}{#sqrt{s + b + (0.1 #times b)^{2}}}");
  legend ->SetTextSize(0.07);
  legend ->SetTextFont(102);
  TLegendEntry *entry;
  TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First(); 
  header ->SetTextSize(.05);
  header ->SetTextFont(102);

  gStyle->SetOptTitle(1); // set hist title
  gStyle->SetTitleSize(0.1,"t");
  gStyle->SetTitleX(0.2); //title X location
  gStyle->SetTitleY(0.9); //title Y location
  gStyle->SetTitleW(0.2); //title width
  gStyle->SetTitleH(0.2); //title height 

  // directories
  int ndir{0};
  for(auto &dir : vdir){ // mt, et, tt, em
    cout<<"\ndir: "<<dir<<endl;

    int njet{0};
    for(auto &jet : vjet){ // 0j, 1j, 01j
      //cout<<"cat: "<<jet<<endl;

      // for each dir. take the 1st hist as model
      string filename{"hist/"+jet+"/"+vinputfile[0]};
      TFile *iifile = new TFile(filename.c_str(), "READ");
      if (!iifile->IsOpen()) { cout<<"can't open "<<filename<<" skipping!!\n"; njet++; continue; }
      string histname{dir+"/bdt"};
      cout<<"histname: "<<histname<<endl;
      string clonename{dir+"_"+jet};
      
      TH1D *test = (TH1D*) iifile->Get(histname.c_str());
      if(!test) {cout<<"dir "<<dir<<" does not exist!!\n"; continue;}

      hfom[ndir][njet] = (TH1D*) iifile->Get(histname.c_str()) -> Clone(clonename.c_str());


      //cout<<"hmodel from bdt name: "<<hfom[ndir][njet]->GetName()<<endl;

      // bdtd steps
      int nbin {hfom[ndir][njet]->GetNbinsX()}; 
      //cout<<"nbin: "<<nbin<<endl;

      // retrieve the fom for each vcat
      vector<double> vfomcut, verrfomcut;


      for(int bin=15; bin>=1; bin--){

	//cout<<"\nbin: "<<bin<<endl;

	double sig{0}, bkg{0}, fom, errfom;

	// input file for BDT fom
	for (auto &inputfile : vinputfile){ 

	  filename = "hist/"+jet+"/"+inputfile;
	  TFile *ifile = new TFile(filename.c_str(), "READ");
	  if(!ifile->IsOpen()) {cout<<"file "<<inputfile<< " not found!!\n"; return;}
	  //cout<<"file: "<<inputfile<<endl;

	  // take the BDT histogram
	  histname = dir+"/bdt";
	  TH1D* hist = (TH1D*)ifile->Get(histname.c_str());
	  if(!hist) {cout<<"hist "<<histname<<" does not exist!!\n"; return;}
	  //cout<<"BDTD hist: "<<histname<<endl;

	  // get nb signal and bkg from max bin to current bin
	  double integral{hist->Integral(bin, 15, "width")};
	  //double integral{hist->Integral(bin, 15)}; 
	  (inputfile.find("125")!=string::npos) ? sig+=integral : bkg+=integral;
	  //cout<<"sig: "<<sig<<endl;
	  //cout<<"bkg: "<<bkg<<endl;

	} //input file


	// compute BDT fom
	FOM(sig, bkg, fom, errfom, 0.1);
	hfom[ndir][njet]->SetBinContent(bin, fom);
	hfom[ndir][njet]->SetBinError( bin, errfom);
	//cout<<"sig, bkg, fom, errfom: "<< sig<<" : "<< bkg<<" : "<< fom<<" : "<< errfom<<endl;


	// input file for CUT fom, only 1 value and for j0 and j1!!
	if(bin==1 && njet<2){
	  //cout<<"\nbin FOR CUT!!: "<<bin<<endl;

	  // take the right hist for each cut
	  for (auto &cat : vcat){ 
            //cout<<"cat: "<<cat<<endl;

	    double sigcut{0}, bkgcut{0};
	    for (auto &inputfile : vinputfile){ 

	      filename = "hist/"+jet+"/"+inputfile;
	      TFile *ifile = new TFile(filename.c_str(), "READ");
	      if(!ifile->IsOpen()) {cout<<"file "<<inputfile<< " not found!!\n"; return;}
	      //cout<<"file: "<<inputfile<<endl;

	      histname = dir+"/"+cat;
	      //cout<<"histname CUT: "<<histname<<endl;
	      TH1D* hist = (TH1D*)ifile->Get(histname.c_str());
	      if(!hist) {cout<<"hist "<<histname<<" does not exist!!\n"; return;}


	      double sum = hist->Integral("width");
	      (inputfile.find("125")!=string::npos) ? sigcut+=sum : bkgcut+=sum;
	      //cout<<"sigcut: "<<sigcut<<endl;
	      //cout<<"bkgcut: "<<bkgcut<<endl;

	      ifile->Close();
	    }// input files

            // fill fom vector for each catgory
	    FOM(sigcut, bkgcut, vfomcut, verrfomcut, 0.1);

	  }// cat

	  //for(auto &i:vfomcut) cout<<"\nfom: "<<i<<endl;
	  //cout<<"vfomcut size(): "<<vfomcut.size()<<endl;

	}// bin=1&&njet<2


      } // bin

      //for(auto &i:vfomcut) {cout<<"\nfom______: "<<i<<endl;}
      //cout<<"vfomcut size()________: "<<vfomcut.size()<<endl;
      //cout<<"verrfomcut size()_______: "<<verrfomcut.size()<<endl;


      // hist settings
      double ymax = hfom[ndir][njet]->GetBinContent(hfom[ndir][njet]->GetMaximumBin());
      double xmin = hfom[ndir][njet]->GetXaxis()->GetBinUpEdge(hfom[ndir][njet]->FindFirstBinAbove(0.));
      double xmax = hfom[ndir][njet]->GetXaxis()->GetBinUpEdge(hfom[ndir][njet]->FindLastBinAbove(0.));
      //hfom[ndir][njet]->GetXaxis()->SetRangeUser(xmin, xmax+0.1 );
      //hfom[ndir][njet]->GetYaxis()->SetRangeUser(0.0001, ymax+0.2*ymax);
      hfom[ndir][njet]->GetXaxis()->SetRangeUser(-0.6, 1 );
      if(dir=="tt")hfom[ndir][njet]->GetYaxis()->SetRangeUser(0.0001, 0.6);
      else if(dir=="em")hfom[ndir][njet]->GetYaxis()->SetRangeUser(0.0001, 0.3);
      else hfom[ndir][njet]->GetYaxis()->SetRangeUser(0.0001, 0.4);
      hfom[ndir][njet]->GetXaxis()->SetTitle("BDT");
      hfom[ndir][njet]->GetYaxis()->SetTitle("FOM");
      hfom[ndir][njet]->GetXaxis()->SetTitleSize(0.06);
      hfom[ndir][njet]->GetYaxis()->SetTitleSize(0.07);
      hfom[ndir][njet]->GetXaxis()->SetTitleOffset(0.8);
      hfom[ndir][njet]->GetYaxis()->SetTitleOffset(1);
      hfom[ndir][njet]->SetLineWidth(2);
      hfom[ndir][njet]->SetLineColor(kBlack);
      hfom[ndir][njet]->SetTitle(vdir_[ndir].c_str());
      hfom[ndir][njet]->SetTitleSize(0.05, "t");


      // draw fom func for each cat
      if(njet<2){
	for(int cat=0; cat<vcat.size(); cat++){

	  string fname {dir+"_"+jet+"_"+vcat[cat]};
	  //cout<<"\nfname: "<<fname<<endl;
	  //if( (dir=="mt" || dir=="et") && vfomcut[cat]>0.1)
	  //  fnc[ndir][njet][cat] = new TF1(fname.c_str(), Form("%g", vfomcut[cat]-0.05), -2, 2);  //!!!!!!
	  //else fnc[ndir][njet][cat] = new TF1(fname.c_str(), Form("%g", vfomcut[cat]), -2, 2);  
	  fnc[ndir][njet][cat] = new TF1(fname.c_str(), Form("%g", vfomcut[cat]), -2, 2);  
	  fnc[ndir][njet][cat] ->SetLineColor(vcolorfom[cat]);
	  if(ndir==0 && njet==0) entry = legend->AddEntry(fnc[ndir][njet][cat], vcat[cat].c_str());
	}
      }
      if(njet==0){
	canvas->cd(ndir+1);
	if(hfom[ndir][njet]->Integral()>0) hfom[ndir][njet]->Draw("hist&same&e");
	for(int cat=0; cat<vcat.size(); cat++) if(fnc[ndir][njet][cat]->Integral(-2, 2)>0) if(!bdttype[1])fnc[ndir][njet][cat] ->Draw("same"); // bdttype[1] is a in between step, no need to draw 
	//cout<<"cat: "<<njet<<endl;
	//cout<<"dir: "<<dir<<endl;
	//for(int cat=0; cat<vcat.size(); cat++){
	//  cout<<"fnc: "<<fnc[ndir][njet][cat]->GetName();
	//  cout<<"fnc int: "<<fnc[ndir][njet][cat]->Integral(-2, 2)<<endl;
	//}

      }
      if(njet==1){
	canvas->cd(ndir+1+4); // draw 1jet on the second row
	if(hfom[ndir][njet]->Integral()>0) hfom[ndir][njet]->Draw("hist&same&e");
	for(int cat=0; cat<vcat.size(); cat++) if(fnc[ndir][njet][cat]->Integral(-2, 2)>0) if(!bdttype[1])fnc[ndir][njet][cat] ->Draw("same"); // bdttype[1] is a in between step, no need to draw
	//cout<<"cat: "<<njet<<endl;
	//cout<<"dir: "<<dir<<endl;
	//for(int cat=0; cat<vcat.size(); cat++){
	//  cout<<"fnc: "<<fnc[ndir][njet][cat]->GetName();
	//  cout<<"fnc int: "<<fnc[ndir][njet][cat]->Integral(-2, 2)<<endl;
	//}

      }
      if(njet==2){
	canvas->cd(ndir+1+8); // draw 1jet on the second row
	if(hfom[ndir][njet]->Integral()>0) hfom[ndir][njet]->Draw("hist&same&e");
      }

      njet++;
    }// njet

    ndir++;
  }// dir

  //canvas->cd(4);
  //legend->Draw();
  string path {"/lstore/cms/calpas/h2taus/CMSSW_8_0_20/src/miniAOD/miniAnalyzer/tmva/fom/plot/"};
  string image {path+"fom.png"};
  canvas->Print(image.c_str()); 

}//void


