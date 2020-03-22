#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include <math.h> 
#include "TMath.h" 
#include <limits>
#include "TSystem.h"
#include <vector>
#include <string>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "SubmitSVFit/ROOT/interface/scenario.h"
//If doES       0 does not apply any ES shifts
//              1 applies ES shifts to TT channel, no effect on other channels
//
//If metType    1 use mvamet
//        -1 use pf met

ClassicSVfit svfitAlgorithm;

int copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser,  char TreeToUse[], int doES, std::string unc = "") ;
int CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);
double tesUncertainties(unsigned int year, float decaymode, float pT); 
double tesUncertainties_MC(unsigned int year, float decaymode);
double tesUncertainties_MChighpT(unsigned int year, float decaymode);
double pt_shifted(float pt, double tesUnc, bool isDM, int updown);
double metcorr_shifted(double metcorr, 
		       float pt1, float phi1, bool isDM1, double tesUnc1, 
		       float pt2, float phi2, bool isDM2, double tesUnc2, 
		       int xory, int updown);
void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
              double measuredMETx, double measuredMETy,
              TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta,
              float &svFitPhi, float &svFitMET, float &svFitTransverseMass,
              TLorentzVector& tau1, TLorentzVector& tau2);

int main (int argc, char* argv[]) 
{
  optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
  parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
  parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
  parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
  parser.addOption("doRecoil",optutl::CommandLineParser::kDouble,"doRecoil",0.0);
  parser.addOption("doMET",optutl::CommandLineParser::kDouble,"doMET",0.0);
  parser.addOption("print",optutl::CommandLineParser::kDouble,"print",0.0);
  parser.parseArguments (argc, argv);
  
  std::cout << "EXTRA COMMANDS:"
	    << "\n --- print: " << parser.doubleValue("print") 
	    << "\n --- doRecoil: " << parser.doubleValue("doRecoil") 
	    << "\n --- doMET: " << parser.doubleValue("doMET") 
	    << "\n --- doES: " << parser.doubleValue("doES") << std::endl;
  
  // No DiTauMass constraint
  svfitAlgorithm.setDiTauMassConstraint(-1.0);
  
  char TreeToUse[80]="first";
  TFile *fProduce;//= new TFile(parser.stringValue("newFile").c_str(),"UPDATE");
  
  TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"READ");
  std::cout<<"Creating new outputfile"<<std::endl;
  std::string newFileName = parser.stringValue("newFile");
  
  fProduce = new TFile(newFileName.c_str(),"RECREATE");
  if ( copyFiles(parser, f, fProduce) == 0 ) return -1;

  fProduce = new TFile(newFileName.c_str(),"UPDATE");
  std::cout<<"listing the directories================="<<std::endl;
  fProduce->ls();
  std::string unc[2] = {
    "nominal","UESUp"
  };
  
  // Add systematics here! 
  std::string uncJES[24] = { // All MC
    "JERUp", "JERDown",
    "JetAbsoluteUp","JetAbsoluteDown",
    "JetAbsoluteyearUp","JetAbsoluteyearDown",
    "JetBBEC1Up","JetBBEC1Down",
    "JetBBEC1yearUp","JetBBEC1yearDown",
    "JetEC2Up","JetEC2Down",
    "JetEC2yearUp","JetEC2yearDown",
    "JetFlavorQCDUp","JetFlavorQCDDown",
    "JetHFUp","JetHFDown",
    "JetHFyearUp","JetHFyearDown",
    "JetRelativeSampleUp","JetRelativeSampleDown",
    "JetRelativeBalUp","JetRelativeBalDown",
  };
  std::string uncTES[8] = { // All MC
    "DM0_Up","DM0_Down",
    "DM1_Up","DM1_Down",
    "DM10_Up","DM10_Down",
    "DM11_Up","DM11_Down"
  };
  std::string uncRecoil[4] = { // qqH, ggH, DY
    "resolutionUp","resolutionDown",
    "responseUp","responseDown"
  };
  std::string uncMET[2]={ // All MC but for qqH, ggH, DY
    "UESUp","UESDown"
  };
  

  // Nominal
  readdir(fProduce,parser,TreeToUse,0,"nominal");
  // Test Uncertainties
  /*
  for (unsigned int i=0; i != sizeof(unc)/sizeof(std::string); ++i) {
    //if ( i != 0 ) break; // Un comment for deugging - won't run over uncertainties
    std::cout << "\n\n\033[1;31mRUN "<< unc[i] << "\033[0m" <<std::endl;
    readdir(fProduce,parser,TreeToUse,0,unc[i]);
    } */     

  // Uncertainties  
  if (parser.doubleValue("doRecoil")) {
    for (unsigned int i=0; i != sizeof(uncRecoil)/sizeof(std::string); ++i) {
      //if ( i != 0 ) break; // Un comment for deugging - won't run over uncertainties
      std::cout << "\n\n\033[1;31mRUN "<< uncRecoil[i] << "\033[0m" <<std::endl;
      readdir(fProduce,parser,TreeToUse,0,uncRecoil[i]);
    }      
  }  
  if (parser.doubleValue("doMET")) {
    for (unsigned int i=0; i != sizeof(uncMET)/sizeof(std::string); ++i) {
      //if ( i != 0 ) break; // Un comment for deugging - won't run over uncertainties
      std::cout << "\n\n\033[1;31mRUN "<< uncMET[i] << "\033[0m" <<std::endl;
      readdir(fProduce,parser,TreeToUse,0,uncMET[i]);
    }      
  }
  if (parser.doubleValue("doES")) {
    // JES
    for (unsigned int i=0; i != sizeof(uncJES)/sizeof(std::string); ++i) {
      std::cout << "\n\n\033[1;31mRUN "<< uncJES[i] << "\033[0m" <<std::endl;
      readdir(fProduce,parser,TreeToUse,0,uncJES[i]);
    }

    // TES
    for (unsigned int i=0; i != sizeof(uncTES)/sizeof(std::string); ++i) {
      std::cout << "\n\n\033[1;31mRUN "<< uncTES[i] << "\033[0m" <<std::endl;
      readdir(fProduce,parser,TreeToUse,1,uncTES[i]);
    }
  }


  fProduce->Close();
  f->Close();
} 


void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES, std::string unc) 
{
  
  TLorentzVector tau1, tau2;
  
  classic_svFit::MeasuredTauLepton::kDecayType decayType1 = classic_svFit::MeasuredTauLepton::kUndefinedDecayType;
  classic_svFit::MeasuredTauLepton::kDecayType decayType2 = classic_svFit::MeasuredTauLepton::kUndefinedDecayType; 
  
  // Both masses should depend on decay mode and particle?
  float mass1;
  float mass2;
  std::string channel = "x";
  
  TDirectory *dirsav = gDirectory;
  TKey *key;
  dir->cd();      
  
  std::vector<TString> processedNames;
  
  TIter next(dir->GetListOfKeys());
  while ((key = (TKey*)next())) {
    printf("Found key=%s \n",key->GetName());
    
    TObject *obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      std::cout << "This is a directory, diving in!" << std::endl;
      // zero the processedNames vector, to allow processing trees with duplicate names in separate directories
      processedNames.clear();
      
      dir->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      sprintf(TreeToUse,"%s",key->GetName());
      readdir(subdir, parser, TreeToUse, doES, unc);
      dirsav->cd();
    }
    else if(obj->IsA()->InheritsFrom(TTree::Class())) {
      // check  if this tree was already processed
      std::vector<TString>::const_iterator it = find(processedNames.begin(), processedNames.end(), key->GetName());
      if ( it != processedNames.end() ) {
        std::cout << "This tree was already processed, skipping..." <<  std::endl;
        continue;
      }
      std::cout << "This is the tree! Start processing" << std::endl;
      processedNames.push_back(key->GetName());
      
      // Identify the process
      if ( std::string(key->GetName()).find("tt") != std::string::npos )  {
        decayType1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        mass1 = 0.13957;
        mass2 = 0.13957;
        channel = "tt";
        std::cout << "Identified channel tt and using kappa = 5" << std::endl;
        svfitAlgorithm.addLogM_fixed(true, 5);
      } 
      else {
        std::cout<<"Tree "<< key->GetName() <<" does not match ... Skipping!!"<<std::endl;
        return;
      }

    
      TTree *t = (TTree*)obj;
      scenario s(t, unc);
      std::string postfix = unc!="nominal"?"_"+unc:"";
      std::string m_sv_unc = "m_sv"+postfix;
      std::string pt_sv_unc = "pt_sv"+postfix;
      std::string eta_sv_unc = "eta_sv"+postfix;
      std::string phi_sv_unc = "phi_sv"+postfix;
      std::string met_sv_unc = "met_sv"+postfix;
      std::string mt_sv_unc = "mt_sv"+postfix;
      std::string metcorr_ex_sv_unc = "metcorr_ex_sv"+postfix;
      std::string metcorr_ey_sv_unc = "metcorr_ey_sv"+postfix;
      std::string metcor_sv_unc = "metcor_sv"+postfix;
      std::string metcorphi_sv_unc = "metcorphi_sv"+postfix;

      TBranch *newBranch1 = t->Branch(m_sv_unc.c_str(), &s.svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch(pt_sv_unc.c_str(), &s.svFitPt, "pt_sv/F");
      TBranch *newBranch3 = t->Branch(eta_sv_unc.c_str(), &s.svFitEta, "eta_sv/F");
      TBranch *newBranch4 = t->Branch(phi_sv_unc.c_str(), &s.svFitPhi, "phi_sv/F");
      TBranch *newBranch5 = t->Branch(met_sv_unc.c_str(), &s.svFitMET, "met_sv/F");
      TBranch *newBranch6 = t->Branch(mt_sv_unc.c_str(), &s.svFitTransverseMass, "mt_sv/F");      
      TBranch *newBranch7 = t->Branch(metcorr_ex_sv_unc.c_str(), &s.metcorr_ex, "metcorr_ex_sv/F");
      TBranch *newBranch8 = t->Branch(metcorr_ey_sv_unc.c_str(), &s.metcorr_ey, "metcorr_ey_sv/F");
      TBranch *newBranch9 = t->Branch(metcor_sv_unc.c_str(), &s.metcor, "metcor_sv/F");
      TBranch *newBranch10 = t->Branch(metcorphi_sv_unc.c_str(), &s.metcorphi, "metcorphi_sv/F");
      
      // adding tau-related branches
      std::vector<TBranch*> tau4VectorBranches;
      std::string tau1_pt_unc = "tau1_pt"+postfix;
      std::string tau1_eta_unc = "tau1_eta"+postfix;
      std::string tau1_phi_unc = "tau1_phi"+postfix;
      std::string tau1_m_unc = "tau1_m"+postfix;
      std::string tau2_pt_unc = "tau2_pt"+postfix;
      std::string tau2_eta_unc = "tau2_eta"+postfix;
      std::string tau2_phi_unc = "tau2_phi"+postfix;
      std::string tau2_m_unc = "tau2_m"+postfix;
      tau4VectorBranches.push_back(t->Branch(tau1_pt_unc.c_str(),  &s.tau1_pt,  "tau1_pt/F"));
      tau4VectorBranches.push_back(t->Branch(tau1_eta_unc.c_str(), &s.tau1_eta, "tau1_eta/F"));
      tau4VectorBranches.push_back(t->Branch(tau1_phi_unc.c_str(), &s.tau1_phi, "tau1_phi/F"));
      tau4VectorBranches.push_back(t->Branch(tau1_m_unc.c_str(),   &s.tau1_m,   "tau1_m/F"));
      tau4VectorBranches.push_back(t->Branch(tau2_pt_unc.c_str(),  &s.tau2_pt,  "tau2_pt/F"));
      tau4VectorBranches.push_back(t->Branch(tau2_eta_unc.c_str(), &s.tau2_eta, "tau2_eta/F"));
      tau4VectorBranches.push_back(t->Branch(tau2_phi_unc.c_str(), &s.tau2_phi, "tau2_phi/F"));
      tau4VectorBranches.push_back(t->Branch(tau2_m_unc.c_str(),   &s.tau2_m,   "tau2_m/F"));
      std::cout << "That's a lot of tau 4-vector branches! N = " << tau4VectorBranches.size() << std::endl;
    
      //unsigned long long evt;
      unsigned int evt, run, lumi, NtupleVer;
      float eta1, phi1, m1;
      int gen_match_1;
      float eta2, phi2, m2;
      int gen_match_2;
      float decayMode=-999.;
      float decayMode2;
      float pfCovMatrix00, pfCovMatrix01, pfCovMatrix10, pfCovMatrix11;

      // define MET
      double measuredMETx = 0.;
      double measuredMETy = 0.;
      TLorentzVector TMet(0,0,0,0);
      // define MET covariance
      TMatrixD covMET(2, 2);
      // Branches from input - no variation
      t->SetBranchAddress("NtupleVer",&NtupleVer);
      t->SetBranchAddress("evt",&evt);
      t->SetBranchAddress("run",&run);
      t->SetBranchAddress("lumi",&lumi);      
      t->SetBranchAddress("gen_match_1",&gen_match_1);
      t->SetBranchAddress("gen_match_2",&gen_match_2);
      t->SetBranchAddress("t1_decayMode", &decayMode);
      t->SetBranchAddress("t2_decayMode", &decayMode2);
      t->SetBranchAddress("eta_1",&eta1);
      t->SetBranchAddress("phi_1",&phi1);
      t->SetBranchAddress("m_1",&m1);
      t->SetBranchAddress("eta_2",&eta2);
      t->SetBranchAddress("phi_2",&phi2);
      t->SetBranchAddress("m_2",&m2);
      t->SetBranchAddress("metcov00",&pfCovMatrix00);
      t->SetBranchAddress("metcov01",&pfCovMatrix01);
      t->SetBranchAddress("metcov10",&pfCovMatrix10);
      t->SetBranchAddress("metcov11",&pfCovMatrix11);

      printf("Found tree -> weighting\n");
    
      int nevents = t->GetEntries();
      //if ( parser.integerValue("numEvents") != -1 ) nevents = parser.integerValue("numEvents");
      for(Int_t i=0;i<nevents;++i){
        t->GetEntry(i);

	covMET[0][0] =  pfCovMatrix00;
	covMET[1][0] =  pfCovMatrix10;
	covMET[0][1] =  pfCovMatrix01;
	covMET[1][1] =  pfCovMatrix11;

	// diff from unc by unc
	float pt1 = s.get_pt_1();
	float pt2 = s.get_pt_2();
	float pfmet = s.get_met();
	float pfmetphi = s.get_metphi();       
	TMet.SetPtEtaPhiM(pfmet,0,pfmetphi,0);
	measuredMETx = pfmet*TMath::Cos(pfmetphi);
	measuredMETy = pfmet*TMath::Sin(pfmetphi);     
        s.metcorr_ex = measuredMETx;
        s.metcorr_ey = measuredMETy;
        s.metcor = TMath::Sqrt( s.metcorr_ex*s.metcorr_ex + s.metcorr_ey*s.metcorr_ey);
        s.metcorphi = TMath::ATan2( s.metcorr_ey, s.metcorr_ex );

	// TES
	if (doES) {
	  bool isDM0_1 = false, isDM1_1 = false, isDM10_1 = false, isDM11_1 = false, isDM0_2 = false, isDM1_2 = false, isDM10_2 = false, isDM11_2 = false;
	  if (gen_match_1==5) {
	    if (decayMode==0) isDM0_1 = true;
	    else if (decayMode==1) isDM1_1 = true;
	    else if (decayMode==10) isDM10_1 = true;
	    else if (decayMode==11) isDM11_1 = true;
	  }
	  if (gen_match_2==5) {
	    if (decayMode2==0) isDM0_2 = true;
	    else if (decayMode2==1) isDM1_2 = true;
	    else if (decayMode2==10) isDM10_2 = true;
	    else if (decayMode2==11) isDM11_2 = true;
	  }
	  std::cout.precision(11);
	  /*
	  std::cout << std::endl << "Before : "<< std::endl <<
	    pt1 << "\t" << decayMode << "\t" << gen_match_1 <<  std::endl <<
	    pt2 << "\t" << decayMode2 << "\t" << gen_match_2 <<std::endl <<
	    s.metcorr_ex << "\t" << s.metcorr_ey << std::endl;
	  */
	  int updown = (unc.find("Up")!= std::string::npos)?+1:-1;
	  if (unc.find("DM0_")!= std::string::npos) {
	    pt1 = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode, pt1), isDM0_1, updown);
	    pt2 = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2, pt2), isDM0_2, updown);
	    s.metcorr_ex = metcorr_shifted(s.metcorr_ex, 
					   pt1, phi1, isDM0_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM0_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   +1, updown);
	    s.metcorr_ey = metcorr_shifted(s.metcorr_ey, 
					   pt1, phi1, isDM0_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM0_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   -1, updown);
	  }
	  /*
	  std::cout << std::endl << "After : "<< std::endl <<
	    pt1 << "\t" << decayMode <<"\t" << gen_match_1<< std::endl <<
	    pt2 << "\t" << decayMode2 <<"\t" << gen_match_2 <<std::endl <<
	    s.metcorr_ex << "\t" << s.metcorr_ey << std::endl;
	  */
	  if (unc.find("DM1_")!= std::string::npos) {
	    pt1 = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode, pt1), isDM1_1, updown);
	    pt2 = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2, pt2), isDM1_2, updown);
	    s.metcorr_ex = metcorr_shifted(s.metcorr_ex, 
					   pt1, phi1, isDM1_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM1_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   +1, updown);
	    s.metcorr_ey = metcorr_shifted(s.metcorr_ey, 
					   pt1, phi1, isDM1_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM1_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   -1, updown);
	  }	  
	  if (unc.find("DM10_")!= std::string::npos) {
	    pt1 = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode, pt1), isDM10_1, updown);
	    pt2 = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2, pt2), isDM10_2, updown);
	    s.metcorr_ex = metcorr_shifted(s.metcorr_ex, 
					   pt1, phi1, isDM10_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM10_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   +1, updown);
	    s.metcorr_ey = metcorr_shifted(s.metcorr_ey, 
					   pt1, phi1, isDM10_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM10_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   -1, updown);
	  }	  
	  if (unc.find("DM11_")!= std::string::npos) {
	    pt1 = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode, pt1), isDM11_1, updown);
	    pt2 = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2, pt2), isDM11_2, updown);
	    s.metcorr_ex = metcorr_shifted(s.metcorr_ex, 
					   pt1, phi1, isDM11_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM11_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   +1, updown);
	    s.metcorr_ey = metcorr_shifted(s.metcorr_ey, 
					   pt1, phi1, isDM11_1, tesUncertainties(NtupleVer, decayMode, pt1),
					   pt2, phi2, isDM11_2, tesUncertainties(NtupleVer, decayMode2, pt2),
					   -1, updown);
	  }	  
	}

        if (parser.doubleValue("print")) std::cout << " - metcor "<< s.metcor << " metcorphi " << s.metcorphi << std::endl;        
	// define lepton four vectors
       std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
       // Add Tau of higest Pt first
       if (pt1 > pt2) {
         measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1, decayMode));
         measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2));
       }
       else {
         measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2));
         measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1, decayMode));
       }
       
       if (parser.doubleValue("print")) std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        

       // run svfit
       runSVFit(measuredTauLeptons, s.metcorr_ex, s.metcorr_ey, covMET, 0, 
		s.svFitMass, s.svFitPt, s.svFitEta, s.svFitPhi, s.svFitMET, s.svFitTransverseMass, tau1, tau2);
       if (parser.doubleValue("print")) std::cout<<"finished running SVFit for "<< unc << std::endl<<std::endl<<std::endl;       
       // fill the tau 4-vector parameters for branch-filling
       s.tau1_pt  = tau1.Pt();
       s.tau1_eta = tau1.Eta();
       s.tau1_phi = tau1.Phi();
       s.tau1_m   = tau1.M();
       s.tau2_pt  = tau2.Pt();
       s.tau2_eta = tau2.Eta();
       s.tau2_phi = tau2.Phi();
       s.tau2_m   = tau2.M();       
       newBranch1->Fill();
       newBranch2->Fill();
       newBranch3->Fill();
       newBranch4->Fill();
       newBranch5->Fill();
       newBranch6->Fill();
       newBranch7->Fill();
       newBranch8->Fill();
       newBranch9->Fill();
       newBranch10->Fill();
       
       
       for(unsigned int i = 0; i != tau4VectorBranches.size(); ++i)
	 (tau4VectorBranches[i])->Fill();       
      }
      dir->cd();
      t->Write("",TObject::kOverwrite);
      delete  t;
    } // if the iterator of the key is a TTree
  }
}

void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons, 
	      double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, 
	      float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass, 
	      TLorentzVector& tau1, TLorentzVector& tau2){
  
  svfitAlgorithm.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  if ( svfitAlgorithm.isValidSolution()) {
    
    svFitMass = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getMass();
    svFitPt = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getPt();
    svFitEta = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getEta();
    svFitPhi = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getPhi();
    svFitTransverseMass = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getTransverseMass();

    classic_svFit::HistogramAdapterTau* h_tau1 = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->tau1();
    classic_svFit::HistogramAdapterTau* h_tau2 = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->tau2();
    const classic_svFit::LorentzVector tau1_p4 = h_tau1->getP4();
    const classic_svFit::LorentzVector tau2_p4 = h_tau2->getP4();
    tau1.SetPtEtaPhiM(tau1_p4.Pt(), tau1_p4.Eta(), tau1_p4.Phi(), tau1_p4.M());
    tau2.SetPtEtaPhiM(tau2_p4.Pt(), tau2_p4.Eta(), tau2_p4.Phi(), tau2_p4.M());

    const classic_svFit::LorentzVector ditau = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getP4();
    svFitMET = (ditau - (measuredTauLeptons[0].p4() + measuredTauLeptons[1].p4())).Pt();

    //TLorentzVector testTau1, testTau2;
    //testTau1.SetPtEtaPhiM(measuredTauLeptons[0].p4().Pt(), measuredTauLeptons[0].p4().Eta(), measuredTauLeptons[0].p4().Phi(), measuredTauLeptons[0].p4().M());
    //testTau2.SetPtEtaPhiM(measuredTauLeptons[1].p4().Pt(), measuredTauLeptons[1].p4().Eta(), measuredTauLeptons[1].p4().Phi(), measuredTauLeptons[1].p4().M());

    //TLorentzVector ditau;
    //ditau.SetPtEtaPhiM(svFitPt, svFitEta, svFitPhi, svFitMass);
  }

}

//Thank you Renee Brun :)
void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir; 
  if(source->GetName()!=parser.stringValue("inputFile")){
    adir = savdir->mkdir(source->GetName());
    std::cout<<"Source name is not outputfile name"<<std::endl;
    adir->cd();    
  }
  else{
    //adir = savdir->mkdir("input");
    adir->cd();    
  }

  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir,parser);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}
int CopyFile(const char *fname, optutl::CommandLineParser parser) {
  //Copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return 0;
  }
  target->cd();
  CopyDir(f,parser);
  delete f;
  target->cd();
  return 1;
}
int copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) 
{
  //prepare files to be copied
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  if ( CopyFile(parser.stringValue("inputFile").c_str(),parser) == 0 ) return 0;
  fNew->ls();
  fNew->Close();
  return 1;

}
/*
double tesUncertainties_MC(unsigned int year, float decaymode) {
  // https://gitlab.cern.ch/doyeong/myelog/blob/master/HiggsToTauTau/TES_deeptau_Feb28.txt
  // https://github.com/cms-tau-pog/TauIDSFs/blob/master/python/TauIDSFTool.py#L162-L187
  double tesSize = -1000;
  if (year==2016) {
    if (decaymode == 0) tesSize = 0.008;
    else if (decaymode == 1) tesSize = 0.006;
    else if (decaymode == 10) tesSize = 0.008;
    else if (decaymode == 11) tesSize = 0.011;
  }
  if (year == 2017) {
    if (decaymode == 0) tesSize = 0.01; 
    else if (decaymode == 1) tesSize = 0.006;
    else if (decaymode == 10) tesSize = 0.007;
    else if (decaymode == 11) tesSize = 0.014;
  }
  if (year == 2018) {
    if (decaymode == 0) tesSize = 0.009; 
    else if (decaymode == 1) tesSize = 0.006; 
    else if (decaymode == 10) tesSize = 0.007;
    else if (decaymode == 11) tesSize = 0.012;
  }
  
  return tesSize;
}

double tesUncertainties_MChighpT(unsigned int year, float decaymode) {
  // https://gitlab.cern.ch/doyeong/myelog/blob/master/HiggsToTauTau/TES_deeptau_Feb28.txt
  // https://github.com/cms-tau-pog/TauIDSFs/blob/master/python/TauIDSFTool.py#L162-L187
  double tesSize = -1000;
  if (year==2016) {
    if (decaymode == 0) tesSize = 0.030;
    else if (decaymode == 1) tesSize = 0.020;
    else if (decaymode == 10) tesSize = 0.012;
    else if (decaymode == 11) tesSize = 0.027;
  }
  if (year == 2017) {
    if (decaymode == 0) tesSize = 0.030; 
    else if (decaymode == 1) tesSize = 0.027;
    else if (decaymode == 10) tesSize = 0.017;
    else if (decaymode == 11) tesSize = 0.040;
  }
  if (year == 2018) {
    if (decaymode == 0) tesSize = 0.030; 
    else if (decaymode == 1) tesSize = 0.020; 
    else if (decaymode == 10) tesSize = 0.011;
    else if (decaymode == 11) tesSize = 0.039;
  }
  
  return tesSize;
}

double tesUncertainties(unsigned int year, float decaymode, float pT) {
  float lowBound=34.0, highBound=170.0;
  double tesSize = -1000;
  float err_high = tesUncertainties_MChighpT(year, decaymode);
  float err_low = tesUncertainties_MC(year, decaymode);

  if (pT>=highBound) 
    tesSize=err_high;
  else if (pT>lowBound) 
    tesSize = err_low + (err_high-err_low)/(highBound-lowBound)*(lowBound);
  else
    tesSize=err_low;
  
  return tesSize;
}
*/

double tesUncertainties(unsigned int year, float decaymode, float pT) {
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingLegacyRun2#Tau_energy_scale_uncertainty
  double tesSize = -1000;
  if (year == 2016) {
    if (decaymode == 0) tesSize = 0.010; 
    else if (decaymode == 1) tesSize =0.009; 
    else if (decaymode == 10) tesSize =0.011;
    else if (decaymode == 11) tesSize =0.011;
  }
  if (year == 2017) {
    if (decaymode == 0) tesSize = 0.008; 
    else if (decaymode == 1) tesSize = 0.008; 
    else if (decaymode == 10) tesSize = 0.009;
    else if (decaymode == 11) tesSize = 0.009;
  }
  if (year == 2018) {
    if (decaymode == 0) tesSize = 0.011; 
    else if (decaymode == 1) tesSize = 0.008; 
    else if (decaymode == 10) tesSize = 0.009;
    else if (decaymode == 11) tesSize = 0.009;
  }
  
  return tesSize;
}


double pt_shifted(float pt, double tesUnc, bool isDM, int updown) {
  double shifted_pT = pt;
  if (isDM) {
    if (updown>0) shifted_pT = pt*(1.0+tesUnc);
    else if (updown<0) shifted_pT = pt*(1.0-tesUnc);
  }
  return shifted_pT;
}

double metcorr_shifted(double metcorr, float pt1, float phi1, bool isDM1, double tesUnc1, float pt2, float phi2, bool isDM2, double tesUnc2, int xory, int updown) {
  double dx1 = 0.0, dx2 = 0.0;
  double shifted_metcorr = metcorr;
  if (isDM1 || isDM2) {
    double tesScale1 = 0.0, tesScale2 = 0.0;
    if (updown>0) {
      tesScale1 = 1.0+tesUnc1;
      tesScale2 = 1.0+tesUnc2;
    }
    else if (updown<0) {
      tesScale1 = 1.0-tesUnc1;  
      tesScale2 = 1.0-tesUnc2;  
    }
    if (xory>0) {
      dx1 = pt1 * TMath::Cos( phi1 ) * (( 1. / tesScale1 ) - 1.);
      dx2 = pt2 * TMath::Cos( phi2 ) * (( 1. / tesScale2 ) - 1.);
    }
    else if (xory<0) {
      dx1 = pt1 * TMath::Sin( phi1 ) * (( 1. / tesScale1 ) - 1.);
      dx2 = pt2 * TMath::Sin( phi2 ) * (( 1. / tesScale2 ) - 1.);
    }
    shifted_metcorr = metcorr + dx1 + dx2;
  }
  return shifted_metcorr;
}
