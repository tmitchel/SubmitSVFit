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

//If doES       0 does not apply any ES shifts
//              1 applies ES shifts to TT channel, no effect on other channels
//
//If isWJets    0 no shift in number of jets used for recoil corrections
//              1 shifts njets + 1 for recoil corrections
//
//If metType    1 use mvamet
//        -1 use pf met

ClassicSVfit svfitAlgorithm;

int copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser,  char TreeToUse[], int doES, int isWJets, int metType, double tesSize) ;
int CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);
double tesUncertainties(unsigned int year, float decaymode); 
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
  parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__svFit__");
  parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
  parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
  parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
  parser.addOption("isWJets",optutl::CommandLineParser::kDouble,"isWJets",0.0);
  parser.addOption("metType",optutl::CommandLineParser::kDouble,"metType",-1.0); // 1 = mvamet, -1 = pf met
  parser.addOption("tesSize",optutl::CommandLineParser::kDouble,"tesSize",0.012); // Default TES = 1.2%
  parser.addOption("numEvents",optutl::CommandLineParser::kInteger,"numEvents",-1);
  parser.parseArguments (argc, argv);
  
  std::cout << "EXTRA COMMANDS:"
        << "\n --- numEvents: " << parser.integerValue("numEvents")
        << "\n --- doES: " << parser.doubleValue("doES")
        << "\n --- isWJets: " << parser.doubleValue("isWJets")
        << "\n --- metType: " << parser.doubleValue("metType")
        << "\n --- tesSize: " << parser.doubleValue("tesSize") << std::endl;
  
  // Make sure a proper Met Type is chosen
  assert (parser.doubleValue("metType") == 1.0 || parser.doubleValue("metType") == -1.0);
  
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
  readdir(fProduce,parser,TreeToUse,parser.doubleValue("doES"),parser.doubleValue("isWJets"),
          parser.doubleValue("metType"),parser.doubleValue("tesSize"));
  
  fProduce->Close();
  f->Close();
} 


void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES, int isWJets, int metType, double tesSize) 
{
  
  TLorentzVector tau1, tau2;
  // up systematics
  TLorentzVector tau1_Up, tau1_DM0_Up, tau1_DM1_Up, tau1_DM10_Up, tau1_UncMet_Up, tau1_ClusteredMet_Up;
  TLorentzVector tau1_JetEC2_Up, tau1_JetEta0to3_Up, tau1_JetEta0to5_Up, tau1_JetEta3to5_Up, tau1_JetRelativeBal_Up, tau1_JetRelativeSample_Up;
  TLorentzVector tau2_Up, tau2_DM0_Up, tau2_DM1_Up, tau2_DM10_Up, tau2_UncMet_Up, tau2_ClusteredMet_Up;
  TLorentzVector tau2_JetEC2_Up, tau2_JetEta0to3_Up, tau2_JetEta0to5_Up, tau2_JetEta3to5_Up, tau2_JetRelativeBal_Up, tau2_JetRelativeSample_Up;
  // down systematics
  TLorentzVector tau1_Down, tau1_DM0_Down, tau1_DM1_Down, tau1_DM10_Down, tau1_UncMet_Down, tau1_ClusteredMet_Down;
  TLorentzVector tau1_JetEC2_Down, tau1_JetEta0to3_Down, tau1_JetEta0to5_Down, tau1_JetEta3to5_Down, tau1_JetRelativeBal_Down, tau1_JetRelativeSample_Down;
  TLorentzVector tau2_Down, tau2_DM0_Down, tau2_DM1_Down, tau2_DM10_Down, tau2_UncMet_Down, tau2_ClusteredMet_Down;
  TLorentzVector tau2_JetEC2_Down, tau2_JetEta0to3_Down, tau2_JetEta0to5_Down, tau2_JetEta3to5_Down, tau2_JetRelativeBal_Down, tau2_JetRelativeSample_Down;
  
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
      readdir(subdir,parser,TreeToUse,parser.doubleValue("doES"),parser.doubleValue("isWJets"),
          parser.doubleValue("metType"),parser.doubleValue("tesSize"));
      
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
      if ( std::string(key->GetName()).find("tautau") != std::string::npos )  {
        decayType1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        mass1 = 0.13957;
        mass2 = 0.13957;
        channel = "tt";
        std::cout << "Identified channel tt and using kappa = 5" << std::endl;
        svfitAlgorithm.addLogM_fixed(true, 5);
      } 
      else if ( std::string(key->GetName()).find("em") != std::string::npos )  {
        std::cout<< "EMu sample" <<std::endl;
        decayType1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
        decayType2 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
        mass1 = 0.00051100;
        mass2 = 0.105658;
        channel = "em";
        std::cout << "Identified channel em and using kappa = 3" << std::endl;
        svfitAlgorithm.addLogM_fixed(true, 3);
      }
      else if ( std::string(key->GetName()).find("et") != std::string::npos ) {
        std::cout<<"eleTauTree"<<std::endl;
        decayType1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
        decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        mass1 = 0.00051100;
        mass2 = 0;
        channel = "et";
        std::cout << "Identified channel et and using kappa = 4" << std::endl;
        svfitAlgorithm.addLogM_fixed(true, 4);
      } 
      else if ( std::string(key->GetName()).find("mt") != std::string::npos || std::string(key->GetName()).find("mutau") != std::string::npos ) {

        std::cout << "muTauEvent" << std::endl;
        decayType1 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
        decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        mass1 = 0.105658;
        mass2 = 0;
        channel = "mt";
        std::cout << "Identified channel mt and using kappa = 4" << std::endl;
        svfitAlgorithm.addLogM_fixed(true, 4);
      } else {
        std::cout<<"Tree "<< key->GetName() <<" does not match ... Skipping!!"<<std::endl;
        return;
      }
      
      TTree *t = (TTree*)obj;
      float svFitMass = -10;
      float svFitPt = -10;
      float svFitEta = -10;
      float svFitPhi = -10;
      float svFitMET = -10;
      float svFitTransverseMass = -10;
      
      float metcorr_ex = -10; // corrected met px (float)
      float metcorr_ey = -10;  // corrected met py (float)

      // For saving
      float metcor = -10; // corrected metcor
      float metcorphi = -10; // corrected metcorphi

      TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
      TBranch *newBranch3 = t->Branch("eta_sv", &svFitEta, "eta_sv/F");
      TBranch *newBranch4 = t->Branch("phi_sv", &svFitPhi, "phi_sv/F");
      TBranch *newBranch5 = t->Branch("met_sv", &svFitMET, "met_sv/F");
      TBranch *newBranch6 = t->Branch("mt_sv", &svFitTransverseMass, "mt_sv/F");
      
      TBranch *newBranch7 = t->Branch("metcorr_ex_sv", &metcorr_ex, "metcorr_ex_sv/F");
      TBranch *newBranch8 = t->Branch("metcorr_ey_sv", &metcorr_ey, "metcorr_ey_sv/F");
      TBranch *newBranch9 = t->Branch("metcor_sv", &metcor, "metcor_sv/F");
      TBranch *newBranch10 = t->Branch("metcorphi_sv", &metcorphi, "metcorphi_sv/F");
      
      // tau leptons                                                                                                  
      // nominal ========================
      float tau1_pt  = -10;
      float tau1_eta = -10;
      float tau1_phi = -10;
      float tau1_m   = -10;
      float tau2_pt  = -10;
      float tau2_eta = -10;
      float tau2_phi = -10;
      float tau2_m   = -10;
    
      // adding tau-related branches
      std::vector<TBranch*> tau4VectorBranches;
      tau4VectorBranches.push_back(t->Branch("tau1_pt",  &tau1_pt,  "tau1_pt/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta", &tau1_eta, "tau1_eta/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi", &tau1_phi, "tau1_phi/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m",   &tau1_m,   "tau1_m/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt",  &tau2_pt,  "tau2_pt/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta", &tau2_eta, "tau2_eta/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi", &tau2_phi, "tau2_phi/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m",   &tau2_m,   "tau2_m/F"));
      std::cout << "That's a lot of tau 4-vector branches! N = " << tau4VectorBranches.size() << std::endl;
    
      unsigned long long evt;
      unsigned int run, lumi, NtupleVer;
      float pt1;
      float eta1;
      float phi1;
      int gen_match_1;
      float pt2;
      float eta2;
      float phi2;
      float m1;
      float m2;
      int gen_match_2;
      float decayMode=-999.;
      float decayMode2;

      float pfCovMatrix00;
      float pfCovMatrix10;
      float pfCovMatrix01;
      float pfCovMatrix11;
      //float mvamet_ex, // uncorrected mva met px (float)
      //  mvamet_ey, // uncorrected mva met py (float)

      int njets =-999.   ;  // number of jets (hadronic jet multiplicity) (int)

      // define MET
      double measuredMETx = 0.;
      double measuredMETy = 0.;
      float pfmet;
      float pfmetphi;
      TLorentzVector TMet(0,0,0,0);
      // define MET covariance
      TMatrixD covMET(2, 2);

      //ele/mu variables
      TBranch *pt1branch;
       
      t->SetBranchAddress("NtupleVer",&NtupleVer);
      t->SetBranchAddress("evt",&evt);
      t->SetBranchAddress("run",&run);
      t->SetBranchAddress("lumi",&lumi);
      t->SetBranchAddress("gen_match_1",&gen_match_1);
      t->SetBranchAddress("gen_match_2",&gen_match_2);
      t->SetBranchAddress("t1_decayMode", &decayMode);
      t->SetBranchAddress("t2_decayMode", &decayMode2);
      t->SetBranchAddress("t1_pt",&pt1,&pt1branch);
      t->SetBranchAddress("t1_eta",&eta1);
      t->SetBranchAddress("t1_phi",&phi1);
      t->SetBranchAddress("t2_pt",&pt2);
      t->SetBranchAddress("t2_eta",&eta2);
      t->SetBranchAddress("t2_phi",&phi2);
      t->SetBranchAddress("t1_mass",&m1);
      t->SetBranchAddress("t2_mass",&m2);

      t->SetBranchAddress("njets", &njets);
      t->SetBranchAddress("met",&pfmet);
      t->SetBranchAddress("metphi",&pfmetphi);
      // FOR PF MET ANALYSIS
      t->SetBranchAddress("metcov00",&pfCovMatrix00);
      t->SetBranchAddress("metcov01",&pfCovMatrix01);
      t->SetBranchAddress("metcov10",&pfCovMatrix10);
      t->SetBranchAddress("metcov11",&pfCovMatrix11);

      printf("Found tree -> weighting\n");
    
      //double tesUp = 1.0 + tesSize;
      //double tesDown = 1.0 - tesSize;

      int nevents = t->GetEntries();
      if ( parser.integerValue("numEvents") != -1 ) nevents = parser.integerValue("numEvents");
      for(Int_t i=0;i<nevents;++i){
        t->GetEntry(i);

        // Using PF Met or Mva Met?
        if (metType == 1) { // 1 = Mva Met
          std::cerr << "Only PF Met is currently supported. You're in for a world full of problems if you continue." << std::endl;
        } // mva met
        if (metType == -1) { // -1 = PF Met
          TMet.SetPtEtaPhiM(pfmet,0,pfmetphi,0);
          measuredMETx = pfmet*TMath::Cos(pfmetphi);
          measuredMETy = pfmet*TMath::Sin(pfmetphi);

          covMET[0][0] =  pfCovMatrix00;
          covMET[1][0] =  pfCovMatrix10;
          covMET[0][1] =  pfCovMatrix01;
          covMET[1][1] =  pfCovMatrix11;
        } // pf met
     
        metcorr_ex = measuredMETx;
        metcorr_ey = measuredMETy;

        metcor = TMath::Sqrt( metcorr_ex*metcorr_ex + metcorr_ey*metcorr_ey);
        metcorphi = TMath::ATan2( metcorr_ey, metcorr_ex );
        std::cout << " - metcor "<<metcor<<" metcorphi "<<metcorphi<<std::endl;
        
       
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
       
       std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
       runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, 
		svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
       std::cout<<"finished running non-TES SVFit in TT"<<std::endl;
       

       // fill the tau 4-vector parameters for branch-filling
       tau1_pt  = tau1.Pt();
       tau1_eta = tau1.Eta();
       tau1_phi = tau1.Phi();
       tau1_m   = tau1.M();
       tau2_pt  = tau2.Pt();
       tau2_eta = tau2.Eta();
       tau2_phi = tau2.Phi();
       tau2_m   = tau2.M();
       
       std::cout << "\n\n" << std::endl;
       //std::cout << "\n\nex: " << metcorr_ex << "   ey: " << metcorr_ey <<  " phi: " << metcorphi<<"\n"<<std::endl; 
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

double tesUncertainties(unsigned int year, float decaymode) {
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingLegacyRun2#Tau_energy_scale_uncertainty
  double tesSize = -1000;
  if (year == 2016) {
    if (decaymode == 0) tesSize = 0.010; 
    else if (decaymode == 1) tesSize =0.009; 
    else if (decaymode == 10) tesSize =0.011;
  }
  if (year == 2017) {
    if (decaymode == 0) tesSize = 0.008; 
    else if (decaymode == 1) tesSize = 0.008; 
    else if (decaymode == 10) tesSize = 0.009;
  }
  if (year == 2018) {
    if (decaymode == 0) tesSize = 0.011; 
    else if (decaymode == 1) tesSize = 0.008; 
    else if (decaymode == 10) tesSize = 0.009;
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
