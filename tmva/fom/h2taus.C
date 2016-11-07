#define h2taus_cxx
#include "h2taus.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void h2taus::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  cout<<"h2taus::Loop(): nentries "<<nentries<<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if(jentry%100000==0)cout<<"h2taus::Loop(): processing entry "<<jentry<<endl;

    if( nbtag!=0) continue;


    //////////////////
    // 		    //
    // BDT cut hist //
    // 		    //
    //////////////////
    
    if( !bdttype[2] || (bdttype[2] && bdt_sVSdy>0.1)) {
    //if( !bdttype[2] || (bdttype[2] && bdt_sVSdy>0.2)) {
    //if( !bdttype[2] || (bdttype[2] && bdt_sVSdy>0.3)) {
    //if( !bdttype[2] || (bdttype[2] && bdt_sVSdy>0.4)) {

      if( (cat=="0j"  &&  njets==0) ||    
	  (cat=="1j"  &&  njets==1) ||
	  (cat=="01j" && (njets==0 || njets==1))
	)
	
	// fill hist with the bdt you want
	if     (bdttype[0]) hist_bdt->Fill(bdt);
	else if(bdttype[1]) hist_bdt->Fill(bdt_sdyVSwtt);
	else if(bdttype[2]) hist_bdt->Fill(bdt_sVSdy);
    }

    //////////////
    //	        //
    // cut hist //
    //	        //
    //////////////
    bool baseline = d0_1<0.045 && dZ_1<0.2 && d0_2<0.045 && dZ_2<0.2;

    bool ltcut = pt_1>20 && fabs(eta_1)<2.1 && iso_1<0.1 &&
      pt_2>30 && fabs(eta_2)<2.4 && iso_2<1.5;

    bool emcut = pt_1>10 && fabs(eta_1)<2.3 && iso_1<0.1 &&
      pt_2>20 && fabs(eta_2)<2.1 && iso_2<0.1;

    bool ttcut = pt_1>45 && fabs(eta_1)<2.1 && iso_1<1 &&
      pt_2>45 && fabs(eta_2)<2.1 && iso_2<1;

    if( baseline ){ 

      if ( dir=="mt" && ltcut && pfmt_1<30){
	if(cat=="0j" && njets==0){
	  if  ( pt_2<=45) hist[0]->Fill(pt_1); 
	  else            hist[1]->Fill(pt_1); 
	}
	if(cat=="1j" && njets==1){
	  if( pt_2<=45)         	    hist[0]->Fill(pt_1); 
	  if( pt_2>45 && fabs(pt_tt<=100) ) hist[1]->Fill(pt_1); 
	  if( pt_2>45 && fabs(pt_tt> 100) ) hist[2]->Fill(pt_1); 
	}
      }// mt

      if ( dir=="et" && ltcut && pfmt_1<30){
	if(cat=="0j" && njets==0){
	  if  ( pt_2<=45) hist[0]->Fill(pt_1); 
	  else            hist[1]->Fill(pt_1); 
	}
	if(cat=="1j" && njets==1 && met>30){
	  if( pt_2<=45 )       	            hist[0]->Fill(pt_1); 
	  if( pt_2>45 && fabs(pt_tt>100) )  hist[2]->Fill(pt_1);
	}
      }// et

      if ( dir=="em" && emcut ){ // missing BDT cut!!
	if(cat=="0j" && njets==0){
	  if  ( pt_2<=35) hist[0]->Fill(pt_1); 
	  else            hist[1]->Fill(pt_1); 
	}
	if(cat=="1j" && njets==1){
	  if( pt_2<=35 )  hist[0]->Fill(pt_1); 
	  else            hist[1]->Fill(pt_1); 
	}
      }// em

      if ( (cat=="0j" || cat=="1j") && dir=="tt" && ttcut ){
	if( pt_tt>100 && pt_tt<=170 ) hist[3]->Fill(pt_1); 
	if( pt_tt>170 )               hist[4]->Fill(pt_1); 
      }// tt

    }// baseline

  }// nentries

  // scale and write hist
  hist_bdt->Scale(evtweight); hist_bdt->Write();  

  for(int i=0; i<5; i++) {hist[i]->Scale(evtweight); hist[i]->Write();}
    
}
