//#include "header.h"
#include "../../func/myfunc.cc"
#include "h2taus.h"

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

#include <TKey.h>
#include <iostream>
#include <string>
#include <vector>


using namespace std;

int main ( int argc, char * argv[]) {

  bool bdttype[]{1, 0, 0}; // s sv b; s+dy vs wtt, s vs dy;

  // input file with bdt variables
  string path;
  if     (bdttype[0]) path = "../tmvaclassification/ntuple/";
  else if(bdttype[1]) path = "../tmvaclassification_2steps/ntuple_bdt_sdyVSwtt/";
  else if(bdttype[2]) path = "../tmvaclassification_2steps/ntuple_bdt_sdyVSwtt_sVSdy/";

  vector<string> vinput {
    "GluGluHToTauTauM125.root",
      "VBFHToTauTauM125.root",
      "WJets.root",
      "TT.root",
      "DYJetsM50.root",
  };  

  // fill vector with variables to be plots
  vector<string> vvar;
  ifstream fsvar("variable.txt");
  fillvarname(vvar, fsvar);
  cout<<"main(): nb variable to plot: "<<vvar.size()<<endl;

  // jet category
  vector<string> vcat{"0j", "1j", "01j"};

  // directories
  vector<string>vdir{"mt", "et", "em", "tt"};

  // loop over jet category
  for(auto &cat : vcat){
    // create outdir
    string makepathout{"mkdir -p hist/"+cat};
    const int dir_err = system(makepathout.c_str());
    if (-1 == dir_err){printf("Error creating directory!n"); exit(1);}

    // loop over input files
    for (auto &inp : vinput){ 	

      string input{path+cat+"/"+inp};

      cout << "\nmain(): processing file: " << input << endl;
      TFile* ifile = new TFile(input.c_str(),"READ");
      //if (!ifile->IsOpen()) { cout<<"main(): can't open "<<input<<endl; return 1; }
      if (!ifile->IsOpen()) { cout<<"main(): can't open "<<input<<" skipping!!"<<endl; continue; }

      string outputfile{"hist/"+cat+"/"+inp};
      cout << "main(): output file: " << outputfile << endl;
      TFile* ofile = new TFile(outputfile.c_str(),"RECREATE");
      if (!ofile->IsOpen()) { cout<<"main(): can't create "<<outputfile<<endl; return 1; }

      for(auto &dir : vdir ){

	//write hists in corresponding directory
	TDirectory *ofiledir = ofile->mkdir(dir.c_str());
	ofiledir->cd();  

	string treeName="h2taus/"+dir;
	cout<<"main(): tree name: "<<treeName<<endl;

	TTree *tree = (TTree*) ifile->GetObjectChecked(treeName.c_str(), "TTree");
	if(!tree) {cout<<"tree "<<treeName<<" does not exist!!\n"; continue;}

	// analysis
	string sample{inp};
	sample.replace(sample.find_first_of("."), 5, ""); //remove ".root"
	cout<<"main(): sample: "<<sample<<endl;

	// apply weight to histogram
	TH1I *hnevt = (TH1I*)ifile->Get("h2taus/nevt");
	double nevt{hnevt->Integral("weight")};

	double evtweight{getweight(sample, nevt)}; // cs, evt, lumi
	cout<<"weight: "<<evtweight<<endl;

	// fill histogram
	cout<<"main(): filling hist...\n";
	h2taus an(tree, evtweight, vvar, sample, dir, cat, bdttype);

	an.Loop();
	cout<<"main(): end loop!!\n"<<endl;

	delete tree;
      }// dir
      ifile->Close();
      ofile->Close();

    }// input file

  }// cat

  return 0;
}


