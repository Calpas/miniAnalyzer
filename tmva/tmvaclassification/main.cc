#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include <TMVA/Config.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "TStopwatch.h"

#include "../../func/myfunc.cc"



///////////////
//	     //
// namespace //   
//	     //
///////////////
using namespace TMVA;
using namespace std;



std::map<std::string,int> Use = {{"BDTD", 0}};



/////////////////////
//		   //
// user's function //
//		   //
/////////////////////
void method(TString myMethodList);




//////////////////////
// 		    //
// global variables //
// 		    //
//////////////////////
// root file
string path{"/lstore/cms/calpas/h2taus/CMSSW_8_0_20/src/miniAOD/miniAnalyzer/ntuple/"};
//string path{"/afs/cern.ch/work/c/calpas/ntuple/"};

vector<string> vinput {
  "GluGluHToTauTauM125.root",
    "VBFHToTauTauM125.root",
    "WJets.root",
    "TT.root",
    "DYJetsM50.root",
};  

// jet category
vector<string>vcat{"0j", "1j", "01j"};
//vector<string>vcat{"1j"};

// input file tree directories
vector<string>vdir{"mt", "et", "tt", "em"};
//vector<string>vdir{"mt"};




////////////////////
//                //
// classification //
//                //
////////////////////
int TMVAClassification( TString myMethodList = "" )
{


  cout << "\nTMVAClassification\n";

  method(myMethodList);

  for(auto &cat : vcat){
    cout<<"cat: "<<cat<<endl;

    for(auto &dir : vdir ){
      cout<<"\ndir: "<<dir<<endl;


      // output file 
      string makepathout{"mkdir -p "+cat+"/"+dir+"/"};
      const int dir_err = system(makepathout.c_str());
      if (-1 == dir_err){printf("Error creating directory!n"); exit(1);}


      TString pathout{cat+"/"+dir+"/"};
      TString outfileName{pathout+"TMVA.root"};
      TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

      // factory object.
      TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
	  "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" ); // Transformations=D to make transfo plots

      // weight file directory
      (TMVA::gConfig().GetIONames()).fWeightFileDir = pathout+"weights";


      //////////////////////////////
      //			   //
      // Define factory variables //
      //			   //
      //////////////////////////////


      //// SET1 
      //factory->AddVariable( "variable", "title", "unit", 'F' or 'I' );
      factory->AddVariable( "pt_1",   "pt_1",   "GeV", 'F' ); 
      factory->AddVariable( "pt_2",   "pt_2",   "GeV", 'F' );
      factory->AddVariable( "pt_tt",  "pt_tt",  "GeV", 'F' );
      factory->AddVariable( "eta_1",  "eta_1",  "",    'F' ); 
      factory->AddVariable( "eta_2",  "eta_2",  "",    'F' );
      factory->AddVariable( "m_vis",  "m_vis",  "GeV", 'F' );
      if(cat=="01j")
	factory->AddVariable( "njets",  "njets",  "",    'I' );
      //// SET2 
      factory->AddVariable( "pfmt_1",  "pfmt_1",  "GeV",    'F' ); 
      factory->AddVariable( "pfmt_2",  "pfmt_2",  "GeV",    'F' ); 
      ////SET3  
      if(cat!="0j"){
	factory->AddVariable( "jpt_1",  "jpt_1",  "GeV", 'F' );
	factory->AddVariable( "jeta_1", "jeta_1", "",    'F' );
      }
      //// SET4
      factory->AddVariable( "iso_2", "iso_2", "", 'F' ); 
      //// SET5
      factory->AddVariable( "dphi_tt", "dphi_tt", "", 'F' );
      //factory->AddVariable( "dr_tt",   "dr_tt",   "", 'F' );
      factory->AddVariable( "dr_t1j1", "dr_t1j1", "", 'F' );
      //factory->AddVariable( "dr_t2j1", "dr_t2j1", "", 'F' );
      //// SET6
      //factory->AddVariable( "pzeta := pfpzetamiss-0.85*pzetavis", 'F' );
      //factory->AddVariable( "m_sv",  "m_sv",  "GeV", 'F' );




      // input file
      for (auto &inp : vinput){ 	

	string input{path+inp};

	TFile* ifile = new TFile(input.c_str(),"READ");
	if (!ifile->IsOpen()) { cout<<"TMVAClassification(): can't open "<<ifile->GetName()<<endl; return 1; }

	string treeName{"h2taus/h2taus/"+dir};
	TTree *tree = (TTree*)ifile->Get(treeName.c_str());
	if(!tree) { cout<<"tree "<<tree->GetName()<<" does not exist!!"<<endl; return 1;}

	// sample name
	string sample{input};
	sample.replace(0, sample.find_last_of("/")+1, ""); // remove ".../"
	sample.replace(sample.find_first_of("."), 5, ""); // remove ".root"
	//cout<<"sample: "<<sample<<endl;

	// normalisation
	TH1I *hnevt = (TH1I*)ifile->Get("h2taus/h2taus/nevt");
	double nevt{hnevt->Integral("weight")};

	Double_t weight = getweight(sample, nevt); // (sample, nevt=1, lumi=1)
	//cout<<"weight: "<<weight<<endl;
	
	if(sample.find("125") != string::npos ) factory->AddSignalTree(tree, weight); // Helge: no need weight here nsig will be nbkg
	else factory->AddBackgroundTree(tree, weight);
      } // vinput

      // apply cut
      TCut mycut="";
      if     (cat=="0j")  mycut="nbtag==0 && njets==0";
      else if(cat=="1j")  mycut="nbtag==0 && njets==1";
      else if(cat=="01j") mycut="nbtag==0 && (njets==0 || njets==1)";

      // configure factory
      factory->PrepareTrainingAndTestTree( mycut, "SplitMode=Random:NormMode=NumEvents:!V" ); 

      int ntree{200};
      (dir!="tt") ? ntree = 100 : ntree = 50; 

      if (Use["BDTD"]) // Decorrelation + Adaptive Boost
	factory->BookMethod( TMVA::Types::kBDT, "BDTD",
	    Form("!H:!V:NTrees=%i:MinNodeSize=5%%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate", ntree) );

      else {cout<<"select method!!\n"; return 1;}

      // train/evaluate/compare
      factory->TrainAllMethods();
      factory->TestAllMethods();
      factory->EvaluateAllMethods();

      // Save the output
      outputFile->Close();
      cout << "output: " << outputFile->GetName(); cout << "\ndone!!\n\n";

      delete factory;

    }// vdir

  }// vcat

  return 0;
} // classification




/////////////////
//             //
// application //
//             //
/////////////////
void TMVAClassificationApplication( TString myMethodList = "" ) 
{   

  cout << "\nTMVAClassificationApplication\n"; 

  method(myMethodList);

  for(auto &cat : vcat){
    cout<<"\ncat: "<<cat<<endl;

    // output file 
    string makepathout{"mkdir -p ntuple/"+cat};
    const int dir_err = system(makepathout.c_str());
    if (-1 == dir_err){printf("Error creating directory!n"); exit(1);}

    for (auto &inp : vinput){ 	

      string input{path+inp};

      // input tree 
      TFile* ifile = new TFile(input.c_str(),"READ");
      if (!ifile->IsOpen()) { cout<<"TMVAClassificationApp(): can't open "<<ifile->GetName()<<endl; return; }

      // sample name
      string sample{input};
      sample.replace(0, sample.find_last_of("/")+1, ""); // remove ".../"
      //cout<<"sample "<<sample<<endl;


      TString ofilename{"ntuple/"+cat+"/"+sample};

      TFile* ofile = TFile::Open( ofilename, "RECREATE" );
      if (!ofile->IsOpen()) { cout<<"main(): can't recreate "<<ofile->GetName()<<endl; return; }

      // go to output file directory
      TDirectory *ofiledir = ofile->mkdir("h2taus");
      ofiledir->cd();  

      // copy nevt histogram for evt normalisation
      TH1I *hnevt = (TH1I*)ifile->Get("h2taus/h2taus/nevt");
      hnevt->Write();

      for(auto &dir : vdir ){
	cout<<"dir: "<<dir<<endl;

	// create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    


	/////////////////////////////
	//			  //
	// Add variables to reader //
	//			  //
	/////////////////////////////
	Float_t pt_1;
	Float_t pt_2;
	Float_t pt_tt;
	Float_t eta_1;
	Float_t jpt_1;
	Float_t jeta_1;
	Float_t njets;
	//
	Float_t eta_2;
	Float_t m_vis;
	Float_t mt_tot;
	//
	Float_t iso_2;
	Float_t iso_1;
	//
	Float_t dphi_tt;
	Float_t dr_tt;
	Float_t dr_t1j1;
	Float_t dr_t2j1;
	//
	Float_t pzeta;
	Float_t met;
	Float_t puppimet;
	Float_t pfmt_1;
	Float_t pfmt_2;
	Float_t mt_sv;
	Float_t m_sv;


	//// SET1 
	reader->AddVariable( "pt_1",   &pt_1 ); 
	reader->AddVariable( "pt_2",   &pt_2 );
	reader->AddVariable( "pt_tt",  &pt_tt );
	reader->AddVariable( "eta_1",  &eta_1 ); 
	reader->AddVariable( "eta_2",  &eta_2 );
	reader->AddVariable( "m_vis",  &m_vis );
	if(cat=="01j")
	  reader->AddVariable( "njets",  &njets );
	////SET2  
	reader->AddVariable( "pfmt_1",  &pfmt_1 ); 
	reader->AddVariable( "pfmt_2",  &pfmt_2 ); 
	////SET3  
	if(cat!="0j"){
	  reader->AddVariable( "jpt_1",  &jpt_1 );
	  reader->AddVariable( "jeta_1", &jeta_1 );
	}
	//// SET4 
	reader->AddVariable( "iso_2", &iso_2 ); 
	//// SET5 
        reader->AddVariable( "dphi_tt", &dphi_tt );
	//reader->AddVariable( "dr_tt",   &dr_tt );
	reader->AddVariable( "dr_t1j1", &dr_t1j1 );
	//reader->AddVariable( "dr_t2j1", &dr_t2j1 );
	//// SET6
	//reader->AddVariable( "pfpzetamiss-0.85*pzetavis", &pzeta );
	//reader->AddVariable( "m_sv",  &m_sv );




	// Define the TMVA method and trainning weights files 
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	  if (it->second) {
	    string methodName = it->first+" method";
	    string weightfile = cat+"/"+dir+"/weights/TMVAClassification_"+it->first+".weights.xml";
	    reader->BookMVA( methodName, weightfile ); 
	  }
	}


	// input tree 
	string treename{"h2taus/h2taus/"+dir}; 
	TTree* itree = (TTree*)ifile->Get(treename.c_str());
	if (!itree) {cout<<"tree "<<itree->GetName()<<" not found!!\n"; return;}

	// clone input tree and add the BDT branch
	cout<<"clonning input tree..."<<endl;
	// needeed for h2taus.h in fom the plotting
	itree->SetBranchStatus("*",0);
	itree->SetBranchStatus("d0_1", 1);
	itree->SetBranchStatus("d0_2", 1);
	itree->SetBranchStatus("dZ_1", 1);
	itree->SetBranchStatus("dZ_2", 1);
	itree->SetBranchStatus("pfmt_1", 1);
	itree->SetBranchStatus("pfmt_2", 1);
	itree->SetBranchStatus("met", 1);
	itree->SetBranchStatus("puppimet", 1);
	//
	itree->SetBranchStatus("njets", 1);
	itree->SetBranchStatus("nbtag", 1);
	itree->SetBranchStatus("pt_1",  1);
	itree->SetBranchStatus("pt_2",  1);
	itree->SetBranchStatus("pt_tt", 1);
	itree->SetBranchStatus("eta_1", 1);
	itree->SetBranchStatus("jpt_1", 1);
	itree->SetBranchStatus("jeta_1",1);
	//
	itree->SetBranchStatus("eta_2", 1);
	itree->SetBranchStatus("m_vis", 1);
	itree->SetBranchStatus("mt_tot",1);
	//
	itree->SetBranchStatus("iso_1", 1);
	itree->SetBranchStatus("iso_2", 1);
	//
	itree->SetBranchStatus("dphi_tt",1);
	itree->SetBranchStatus("dr_tt",  1);
	itree->SetBranchStatus("dr_t1j1",1);
	itree->SetBranchStatus("dr_t2j1",1);
	//
	itree->SetBranchStatus("pzetavis",1 );
	itree->SetBranchStatus("pfpzetamiss",1);
	//
	itree->SetBranchStatus("m_sv", 1);


	TTree *otree = itree->CloneTree();
	if(!otree) {cout<<"clone tree failed!!\n"; return;}

	// adding BDT variable to the cloned tree
	Double_t bdt_response;
	TBranch *bdt = otree->Branch("bdt", &bdt_response);

	Int_t njets_; 
	Int_t nbtag_;
	Double_t pt_1_;
	Double_t pt_2_;
	Double_t pt_tt_;
	Double_t eta_1_;
	Double_t jpt_1_;
	Double_t jeta_1_;
	//
	Double_t eta_2_;
	Double_t m_vis_;
	Double_t mt_tot_;
	//
	Double_t iso_1_;
	Double_t iso_2_;
	//
	Double_t dphi_tt_;
	Double_t dr_tt_;
	Double_t dr_t1j1_;
	Double_t dr_t2j1_;
	//
	Double_t pzetavis_;
	Double_t pfpzetamiss_;
	//
	Double_t met_;
	Double_t puppimet_;
	Double_t pfmt_1_;
	Double_t pfmt_2_;
	Double_t mt_sv_;
	Double_t m_sv_;

	itree->SetBranchAddress( "njets",  &njets_  );
	itree->SetBranchAddress( "nbtag",  &nbtag_  );
	itree->SetBranchAddress( "pt_1",   &pt_1_   );
	itree->SetBranchAddress( "pt_2",   &pt_2_   );
	itree->SetBranchAddress( "pt_tt",  &pt_tt_  );
	itree->SetBranchAddress( "eta_1",  &eta_1_  );
	itree->SetBranchAddress( "jpt_1",  &jpt_1_  );
	itree->SetBranchAddress( "jeta_1", &jeta_1_ );
	//
	itree->SetBranchAddress( "eta_2",  &eta_2_  );
	itree->SetBranchAddress( "m_vis",  &m_vis_  );
	itree->SetBranchAddress( "mt_tot", &mt_tot_ );
	//
	itree->SetBranchAddress( "met",    &met_ );
	itree->SetBranchAddress( "puppimet", &puppimet_ );
	itree->SetBranchAddress( "pfmt_1", &pfmt_1_ );
	itree->SetBranchAddress( "pfmt_2", &pfmt_2_ );
	itree->SetBranchAddress( "mt_sv",  &mt_sv_ );
	itree->SetBranchAddress( "m_sv", &m_sv_ );
	//
	itree->SetBranchAddress( "iso_1",  &iso_1_  );
	itree->SetBranchAddress( "iso_2",  &iso_2_  );
	//
	itree->SetBranchAddress( "dphi_tt", &dphi_tt_ );
	itree->SetBranchAddress( "dr_tt",   &dr_tt_   );
	itree->SetBranchAddress( "dr_t1j1", &dr_t1j1_ );
	itree->SetBranchAddress( "dr_t2j1", &dr_t2j1_ );
	//
	itree->SetBranchAddress( "pzetavis", &pzetavis_ );
	itree->SetBranchAddress( "pfpzetamiss", &pfpzetamiss_ );
	//

	TStopwatch sw;
	sw.Start();
	for (Long64_t ievt=0; ievt<itree->GetEntries();ievt++) {

	  if (ievt%100000 == 0) cout << "processing event: " << ievt << endl;

	  itree->GetEntry(ievt);

	  // set the reader variables
	  pt_1    = float(pt_1_);
	  pt_2    = float(pt_2_);
	  pt_tt   = float(pt_tt_);
	  eta_1   = float(eta_1_);
	  jpt_1   = float(jpt_1_);
	  jeta_1  = float(jeta_1_);
	  eta_2   = float(eta_2_);
	  njets   = float(njets_);
	  m_vis   = float(m_vis_);
	  mt_tot  = float(mt_tot_);
	  iso_1   = float(iso_1_);
	  iso_2   = float(iso_2_);
	  dphi_tt = float(dphi_tt_);
	  dr_tt   = float(dr_tt_);
	  dr_t1j1 = float(dr_t1j1_);
	  dr_t2j1 = float(dr_t2j1_);
	  pzeta   = float(pfpzetamiss_-0.85*pzetavis_ );
	  met     = float(met_);
	  puppimet = float(puppimet_);
	  pfmt_1  = float(pfmt_1_);
	  pfmt_2  = float(pfmt_2_);
	  mt_sv   = float(mt_sv);
	  m_sv    = float(m_sv_);

	  bdt_response = -10;
	  if( (cat=="0j"  &&  njets_==0 && nbtag_==0) ||
	      (cat=="1j"  &&  njets_==1 && nbtag_==0) ||
	      (cat=="01j" && (njets_==0 || njets_==1) && nbtag_==0)
	    )
	    bdt_response = reader->EvaluateMVA("BDTD method");

	  bdt->Fill();
	}
	sw.Stop();

	delete reader;

	cout << "output: " << ofile->GetName(); cout <<"\ndone!!\n\n";

      }// vdir
      ofile->Write();
      ofile->Close();
    }// file
  }// cat

}// application






int main( int argc, char** argv )
{


  // loads the library
  TMVA::Tools::Instance();

  // Select methods 
  TString methodList; 
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(","); 
    methodList += regMethod;
  }

  TMVAClassification(methodList); 
  TMVAClassificationApplication(methodList); 

  return 0;
}




/////////////////////
//		   //
// user's function //
//		   //
/////////////////////

void method(TString myMethodList){

  // default MVA methods to be trained + tested
  std::cout << "checking option..." << std::endl;

  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	std::cout << std::endl;
	return;
      }
      Use[regMethod] = 1;
    }
  }

}




// MVA methods
// http://tmva.sourceforge.net/optionRef.html#MVA::BDT
// H: Print method-specific help message
// V: Verbose output
// NTrees: Number of trees in the forest
// MinNodeSize: Minimum percentage of training events required in a leaf node
// MaxDepth: Max depth of the decision tree allowed
// BoostType: Boosting type for the trees in the forest (AdaBoost, RealAdaBoost, Bagging, AdaBoostR2, Grad)
// SeparationType: Separation criterion for node splitting (CrossEntropy, GiniIndex, GiniIndexWithLaplace, MisClassificationError, SDivSqrtSPlusB, RegressionVariance)
// nCuts: Number of grid points in variable range used in finding optimal cut in node splitting
// VarTransform: List of variable transformations performed before training, e.g., D_Background,P_Signal,G,N_AllClasses for: Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is assumed)     

//if (Use["BDTD"]) // Decorrelation + Adaptive Boost
//  factory->BookMethod( TMVA::Types::kBDT, "BDTD",
//	   "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );


