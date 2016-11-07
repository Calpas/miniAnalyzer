#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <locale> 

#include "TH1.h"
#include "TMath.h"

#include <boost/algorithm/string.hpp>


using namespace std;
using namespace boost;


////////////////////////
//		      // 	
// create directories //
//		      // 	
////////////////////////
void makedir(string path){
  string makepath{"mkdir -p "+path};
  const int dir_err = system(makepath.c_str());
  if (-1 == dir_err){printf("Error creating directory!n"); exit(1);}
}


/////////////////////
//		   //
// figure of merit //
//		   //
/////////////////////
void FOM(double sig, double bkg, double &fom, double &errfom, double f=0.1){

  double den;

  den = sig+bkg+pow(f*bkg, 2);

  den > 0  ? fom = sig/sqrt(den) : fom = 0; 

  den = pow(sig+bkg+pow(f*bkg, 2), 3);

  double dFdS; 
  den != 0 ? dFdS = pow((sig/2)+bkg+pow(f*bkg, 2), 2)/den  : dFdS = 0;

  double dFdB;
  den != 0 ? dFdB = pow(sig,2)*pow(1+(2*f*f*bkg), 2)/den : dFdB = 0;

  errfom = sqrt((dFdS*sig) + (dFdB*bkg));
}


void FOM(double sig, double bkg, vector<double> &vfom, vector<double> &verrfom, double f=0.1){

  double den;

  den = sig+bkg+pow(f*bkg, 2);

  den > 0  ? vfom.push_back(sig/sqrt(den)) : vfom.push_back(0); 

  den = pow(sig+bkg+pow(f*bkg, 2), 3);

  double dFdS; 
  den != 0 ? dFdS = pow((sig/2)+bkg+pow(f*bkg, 2), 2)/den  : dFdS = 0;

  double dFdB;
  den != 0 ? dFdB = pow(sig, 2)*pow(1+(2*f*f*bkg), 2)/den : dFdB = 0;

  verrfom.push_back(sqrt((dFdS*sig) + (dFdB*bkg)));
}



/////////////////////
// 		   //
// vector of color //
// 		   //
/////////////////////
vector<int> vcolorfom{
  kRed,
    kBlue, 	
    kGreen, 	
    kViolet, 	
    kOrange, 
    kMagenta, 	
    kCyan, 	
    kGray+3, 	
    kBlack, 
};

vector<int> vcolor{
  kRed,
    kBlue, 	
    kOrange, 	
    kGreen, 	
    kGray+3, 	
    kAzure, 	
    kMagenta, 	
    kGray+1, 	
    kPink, 
    kOrange, 
};


///////////////////////////////////////////////////////
//						     //	
// fill a vector with variable taken from a txt file //
//						     //	
///////////////////////////////////////////////////////
void fillvarname(vector<string> &vvar, ifstream &fsvar){
  string varname;
  if ( fsvar.is_open() ) {
    while ( getline (fsvar, varname)){ 
      //remove white space
      trim(varname);
      if(!varname.empty()) vvar.push_back(varname); 
    }
    fsvar.close();
  }   
}


//////////////////////////////////////
//				    //
// return sample weight(cs/ngenevt) //
//				    //
//////////////////////////////////////
double getweight(const string &sample, const double &nevt=1, const double &lumi=40000){ // lumi=40000 pb-1
//double getweight(const string &sample, const double &nevt=1, const double &lumi=24600){ 
  double cs{0};
  if     (sample=="GluGluHToTauTauM125") cs = 44.14*0.0627;
  else if(sample=="VBFHToTauTauM125")    cs = 3.782*0.0627;
  else if(sample=="Signal_VBF_GG")       cs = (3.782+44.14)*0.0627;
  else if(sample=="DYJetsM50")           cs = 6025.2;
  else if(sample=="WJets")               cs = 61526.7;
  else if(sample=="TT")                  cs = 831.76;
  return lumi*cs/nevt;
}


//////////////////////////////////////////////////////
//						    //
// set the xaxis hist to the max btw signal and bkg //
//						    //
//////////////////////////////////////////////////////

using My3DTH1DArray1 = std::array<std::array<std::array<std::array<TH1D*, 5>, 20>, 10>, 200>; // cat, file, dir, var
void setaxis(My3DTH1DArray1& h, 
	     const std::vector<std::string> &vcat,
	     const std::vector<std::string> &vfile,
	     const std::vector<std::string> &vdir,
	     const std::vector<std::string> &vvar,
	     const int &rebin
    )
{ 

  for(unsigned int cat=0; cat<vcat.size(); cat++){ 

    for(unsigned int dir=0; dir<vdir.size(); dir++){ 

      for(unsigned int var=0; var<vvar.size(); var++) { 

	int xtempmax;
	int xtempmin;
	double ytemp;
	int xlastbin{0};
	int xfirstbin{10000};
	double xmax{0};
	double xmin{0};
	double ymax{0};
	for(unsigned int file=0; file<vfile.size(); file++){  

	  h[cat][file][dir][var]->Rebin(rebin); 

	  // keep the highest xaxis 
	  double limit{0.001*h[cat][file][dir][var]->Integral("width")}; // 1 per mille
	  xtempmax = h[cat][file][dir][var] -> FindLastBinAbove(limit); //get last bin on x axis with data
	  xtempmin = h[cat][file][dir][var] -> FindFirstBinAbove(limit); //get last bin on x axis with data
	  if(xtempmax > xlastbin )  xlastbin  = xtempmax;
	  if(xtempmin < xfirstbin ) xfirstbin = xtempmin;

	  // keep the max yaxis 
	  ytemp = h[cat][file][dir][var] -> GetBinContent(h[cat][file][dir][var] -> GetMaximumBin());
	  if(ytemp > ymax ) ymax = ytemp;
	}

	// set best axis in hist
	for(unsigned int file=0; file<vfile.size(); file++) { 
	  xmax = h[cat][file][dir][var]-> GetXaxis()->GetBinUpEdge(xlastbin);
	  xmin = h[cat][file][dir][var]-> GetXaxis()->GetBinLowEdge(xfirstbin);
	  h[cat][file][dir][var] -> GetXaxis()->SetRangeUser(xmin, xmax);
	  h[cat][file][dir][var] -> GetYaxis()->SetRangeUser(0., ymax+(0.1*ymax));
	  h[cat][file][dir][var] -> SetLineWidth(2);
	  h[cat][file][dir][var] -> GetXaxis()->SetTitleOffset(0.9);
	  h[cat][file][dir][var] -> GetYaxis()->SetTitleOffset(0.9);
	  h[cat][file][dir][var] -> GetYaxis()->SetTitle("1/dm[1/GeV]");
	  h[cat][file][dir][var] -> GetXaxis()->SetTitleSize(0.055);
	  h[cat][file][dir][var] -> GetYaxis()->SetTitleSize(0.07);
	  h[cat][file][dir][var] -> SetLineColor(vcolor[file]); 
	}

      } // var
    } // dir
  } // cat
} // setaxis



using My3DTH1DArray2 = std::array<std::array<std::array<TH1D*, 20>, 10>, 200>;
void setaxis(My3DTH1DArray2& h, 
    const std::vector<std::string> &vfile,
    const std::vector<std::string> &vdir,
    const std::vector<std::string> &vvar,
    const int &rebin
    )
{ 

  for(unsigned int dir=0; dir<vdir.size(); dir++){ 

    for(unsigned int var=0; var<vvar.size(); var++) { 

      int xtempmax;
      int xtempmin;
      double ytemp;
      int xlastbin{0};
      int xfirstbin{10000};
      double xmax{0};
      double xmin{0};
      double ymax{0};
      for(unsigned int file=0; file<vfile.size(); file++){  

	h[file][dir][var]->Rebin(rebin); 

	// keep the highest xaxis 
	double limit{0.001*h[file][dir][var]->Integral("width")}; // 1 per mille
	xtempmax = h[file][dir][var] -> FindLastBinAbove(limit); //get last bin on x axis with data
	xtempmin = h[file][dir][var] -> FindFirstBinAbove(limit); //get last bin on x axis with data
	if(xtempmax > xlastbin )  xlastbin  = xtempmax;
	if(xtempmin < xfirstbin ) xfirstbin = xtempmin;

	// keep the max yaxis 
	ytemp = h[file][dir][var] -> GetBinContent(h[file][dir][var] -> GetMaximumBin());
	if(ytemp > ymax ) ymax = ytemp;
      }

      // set best axis in hist
      for(unsigned int file=0; file<vfile.size(); file++) { 
	xmax = h[file][dir][var]-> GetXaxis()->GetBinUpEdge(xlastbin);
	xmin = h[file][dir][var]-> GetXaxis()->GetBinLowEdge(xfirstbin);
	h[file][dir][var] -> GetXaxis()->SetRangeUser(xmin, xmax);
	h[file][dir][var] -> GetYaxis()->SetRangeUser(0., ymax+(0.1*ymax));
	h[file][dir][var] -> SetLineWidth(2);
	h[file][dir][var] -> GetXaxis()->SetTitleOffset(0.8);
	h[file][dir][var] -> GetYaxis()->SetTitleOffset(1.1);
	h[file][dir][var] -> GetYaxis()->SetTitle("1/dm[1/GeV]");
	h[file][dir][var] -> GetXaxis()->SetTitleSize(0.06);
	h[file][dir][var] -> GetYaxis()->SetTitleSize(0.06);
	h[file][dir][var] -> SetLineColor(vcolor[file]); 
	h[file][dir][var] -> GetYaxis()->SetTitleOffset(1.4);; 
      }

    } // var
  } // dir
} // setaxis


//////////////////////////////////////////////////////
//						    //
// set the xaxis hist to the max btw signal and bkg //
//						    //
//////////////////////////////////////////////////////
void setaxis(TH1D *hsig, TH1D *hbkg){

  int xsiglastbin{hsig -> FindLastBinAbove(0.0005)}; //get last bin on x axis with data
  double xsigmax = hsig -> GetXaxis()->GetBinUpEdge(xsiglastbin);

  int xbkglastbin{hbkg -> FindLastBinAbove(0.0005)}; //get last bin on x axis with data
  double xbkgmax = hbkg -> GetXaxis()->GetBinUpEdge(xbkglastbin);

  double xmax{max(xsigmax, xbkgmax)};
  hsig -> GetXaxis()->SetRangeUser(0., xmax);
  hbkg -> GetXaxis()->SetRangeUser(0., xmax);

  double ysigmax{hsig -> GetBinContent(hsig -> GetMaximumBin())};
  double ybkgmax{hbkg -> GetBinContent(hbkg -> GetMaximumBin())};
  double ymax{max(ysigmax, ybkgmax)+0.002};
  hsig -> GetYaxis()->SetRangeUser(0., ymax);
  hbkg -> GetYaxis()->SetRangeUser(0., ymax);

  hsig -> SetLineWidth(2);
  hsig -> SetLineColor(kRed);
  hsig -> GetXaxis()->SetTitleOffset(0.8);
  hsig -> GetYaxis()->SetTitleOffset(1.1);
  hsig -> GetYaxis()->SetTitle("1/dm[1/GeV]");

  hsig -> GetXaxis()->SetTitleSize(0.06);
  hsig -> GetYaxis()->SetTitleSize(0.06);

  hbkg -> SetLineWidth(2);

}


/*
//return median of a histogram
double Median(const TH1D * h1) {
int n = h1->GetXaxis()->GetNbins(); 
std::vector<double>  x(n);
h1->GetXaxis()->GetCenter( &x[0] );
const double * y = h1->GetArray();
// exclude underflow/overflows from bin content array y
return TMath::Median(n, &x[0], &y[1]); 
}


//round a value
double rounded(double val)
{
double decimal;
double resultat;

val*=1000;

decimal=val-floor(val);
if (decimal< 0.5) resultat=floor(val);
else resultat=ceil(val);

resultat/=1000;

return resultat;
}
*/
