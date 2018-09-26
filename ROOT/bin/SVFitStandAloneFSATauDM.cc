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

void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser,  char TreeToUse[], int doES, int isWJets, int metType, double tesSize) ;
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);

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
  copyFiles(parser, f, fProduce);//new TFile(parser.stringValue("inputFile").c_str()+"SVFit","UPDATE");
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
  TLorentzVector tau2_Up, tau2_DM0_Up, tau2_DM1_Up, tau2_DM10_Up, tau2_UncMet_Up, tau2_ClusteredMet_Up;
  // down systematics
  TLorentzVector tau1_Down, tau1_DM0_Down, tau1_DM1_Down, tau1_DM10_Down, tau1_UncMet_Down, tau1_ClusteredMet_Down;
  TLorentzVector tau2_Down, tau2_DM0_Down, tau2_DM1_Down, tau2_DM10_Down, tau2_UncMet_Down, tau2_ClusteredMet_Down;
  
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
      if ( std::string(key->GetName()).find("tt") != std::string::npos )  {
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
      float metcorrUncUp_ex = -10; // corrUncUpected met px (float)
      float metcorrUncUp_ey = -10;  // corrUncUpected met py (float)
      float metcorrUncDown_ex = -10; // corrUncDownected met px (float)
      float metcorrUncDown_ey = -10;  // corrUncDownected met py (float)
      float metcorrClusteredUp_ex = -10; // corrClusteredUpected met px (float)
      float metcorrClusteredUp_ey = -10;  // corrClusteredUpected met py (float)
      float metcorrClusteredDown_ex = -10; // corrClusteredDownected met px (float)
      float metcorrClusteredDown_ey = -10;  // corrClusteredDownected met py (float)
      
      // For saving
      float metcor = -10; // corrected metcor
      float metcorphi = -10; // corrected metcorphi
      float metcorClusteredDown = -10;   
      float metcorphiClusteredDown = -10;
      float metcorClusteredUp = -10;     
      float metcorphiClusteredUp = -10;  
      float metcorUncDown = -10;         
      float metcorphiUncDown = -10;      
      float metcorUncUp = -10;           
      float metcorphiUncUp = -10;        
      
      TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
      TBranch *newBranch3 = t->Branch("eta_sv", &svFitEta, "eta_sv/F");
      TBranch *newBranch4 = t->Branch("phi_sv", &svFitPhi, "phi_sv/F");
      TBranch *newBranch5 = t->Branch("met_sv", &svFitMET, "met_sv/F");
      TBranch *newBranch6 = t->Branch("mt_sv", &svFitTransverseMass, "mt_sv/F");
      
      TBranch *newBranch7 = t->Branch("metcorr_ex", &metcorr_ex, "metcorr_ex/F");
      TBranch *newBranch8 = t->Branch("metcorr_ey", &metcorr_ey, "metcorr_ey/F");
      TBranch *newBranch9 = t->Branch("metcor", &metcor, "metcor/F");
      TBranch *newBranch10 = t->Branch("metcorphi", &metcorphi, "metcorphi/F");
      
      // If doing ES shifts, we need extra ouput branches
      //
      float svFitMass_Up = -10;
      float svFitPt_Up = -10;
      float svFitEta_Up = -10;
      float svFitPhi_Up = -10;
      float svFitMET_Up = -10;
      float svFitTransverseMass_Up = -10;
      float svFitMass_Down = -10;
      float svFitPt_Down = -10;
      float svFitEta_Down = -10;
      float svFitPhi_Down = -10;
      float svFitMET_Down = -10;
      float svFitTransverseMass_Down = -10;
      
      float svFitMass_DM0_Up = -10;
      float svFitPt_DM0_Up = -10;
      float svFitEta_DM0_Up = -10;
      float svFitPhi_DM0_Up = -10;
      float svFitMET_DM0_Up = -10;
      float svFitTransverseMass_DM0_Up = -10;
      float svFitMass_DM0_Down = -10;
      float svFitPt_DM0_Down = -10;
      float svFitEta_DM0_Down = -10;
      float svFitPhi_DM0_Down = -10;
      float svFitMET_DM0_Down = -10;
      float svFitTransverseMass_DM0_Down = -10;
      
      float svFitMass_DM1_Up = -10;
      float svFitPt_DM1_Up = -10;
      float svFitEta_DM1_Up = -10;
      float svFitPhi_DM1_Up = -10;
      float svFitMET_DM1_Up = -10;
      float svFitTransverseMass_DM1_Up = -10;
      float svFitMass_DM1_Down = -10;
      float svFitPt_DM1_Down = -10;
      float svFitEta_DM1_Down = -10;
      float svFitPhi_DM1_Down = -10;
      float svFitMET_DM1_Down = -10;
      float svFitTransverseMass_DM1_Down = -10;

      float svFitMass_DM10_Up = -10;
      float svFitPt_DM10_Up = -10;
      float svFitEta_DM10_Up = -10;
      float svFitPhi_DM10_Up = -10;
      float svFitMET_DM10_Up = -10;
      float svFitTransverseMass_DM10_Up = -10;
      float svFitMass_DM10_Down = -10;
      float svFitPt_DM10_Down = -10;
      float svFitEta_DM10_Down = -10;
      float svFitPhi_DM10_Down = -10;
      float svFitMET_DM10_Down = -10;
      float svFitTransverseMass_DM10_Down = -10;

      float svFitMass_UncMet_Up = -10;
      float svFitPt_UncMet_Up = -10;
      float svFitEta_UncMet_Up = -10;
      float svFitPhi_UncMet_Up = -10;
      float svFitMET_UncMet_Up = -10;
      float svFitTransverseMass_UncMet_Up = -10;
      float svFitMass_UncMet_Down = -10;
      float svFitPt_UncMet_Down = -10;
      float svFitEta_UncMet_Down = -10;
      float svFitPhi_UncMet_Down = -10;
      float svFitMET_UncMet_Down = -10;
      float svFitTransverseMass_UncMet_Down = -10;

      float svFitMass_ClusteredMet_Up = -10;
      float svFitPt_ClusteredMet_Up = -10;
      float svFitEta_ClusteredMet_Up = -10;
      float svFitPhi_ClusteredMet_Up = -10;
      float svFitMET_ClusteredMet_Up = -10;
      float svFitTransverseMass_ClusteredMet_Up = -10;
      float svFitMass_ClusteredMet_Down = -10;
      float svFitPt_ClusteredMet_Down = -10;
      float svFitEta_ClusteredMet_Down = -10;
      float svFitPhi_ClusteredMet_Down = -10;
      float svFitMET_ClusteredMet_Down = -10;
      float svFitTransverseMass_ClusteredMet_Down = -10;
      
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
      // up (whatever it is) ============
      float tau1_pt_Up  = -10;
      float tau1_eta_Up = -10;
      float tau1_phi_Up = -10;
      float tau1_m_Up   = -10;
      float tau2_pt_Up  = -10;
      float tau2_eta_Up = -10;
      float tau2_phi_Up = -10;
      float tau2_m_Up   = -10;
      // down 
      float tau1_pt_Down  = -10;
      float tau1_eta_Down = -10;
      float tau1_phi_Down = -10;
      float tau1_m_Down   = -10;
      float tau2_pt_Down  = -10;
      float tau2_eta_Down = -10;
      float tau2_phi_Down = -10;
      float tau2_m_Down   = -10;
      // up DM0 =========================
      float tau1_pt_DM0_Up  = -10;
      float tau1_eta_DM0_Up = -10;
      float tau1_phi_DM0_Up = -10;
      float tau1_m_DM0_Up   = -10;
      float tau2_pt_DM0_Up  = -10;
      float tau2_eta_DM0_Up = -10;
      float tau2_phi_DM0_Up = -10;
      float tau2_m_DM0_Up   = -10;
      // down 
      float tau1_pt_DM0_Down  = -10;
      float tau1_eta_DM0_Down = -10;
      float tau1_phi_DM0_Down = -10;
      float tau1_m_DM0_Down   = -10;
      float tau2_pt_DM0_Down  = -10;
      float tau2_eta_DM0_Down = -10;
      float tau2_phi_DM0_Down = -10;
      float tau2_m_DM0_Down   = -10;
      // up DM1 =========================
      float tau1_pt_DM1_Up  = -10;
      float tau1_eta_DM1_Up = -10;
      float tau1_phi_DM1_Up = -10;
      float tau1_m_DM1_Up   = -10;
      float tau2_pt_DM1_Up  = -10;
      float tau2_eta_DM1_Up = -10;
      float tau2_phi_DM1_Up = -10;
      float tau2_m_DM1_Up   = -10;
      // down 
      float tau1_pt_DM1_Down  = -10;
      float tau1_eta_DM1_Down = -10;
      float tau1_phi_DM1_Down = -10;
      float tau1_m_DM1_Down   = -10;
      float tau2_pt_DM1_Down  = -10;
      float tau2_eta_DM1_Down = -10;
      float tau2_phi_DM1_Down = -10;
      float tau2_m_DM1_Down   = -10;
      // up DM10 =========================
      float tau1_pt_DM10_Up  = -10;
      float tau1_eta_DM10_Up = -10;
      float tau1_phi_DM10_Up = -10;
      float tau1_m_DM10_Up   = -10;
      float tau2_pt_DM10_Up  = -10;
      float tau2_eta_DM10_Up = -10;
      float tau2_phi_DM10_Up = -10;
      float tau2_m_DM10_Up   = -10;
      // down 
      float tau1_pt_DM10_Down  = -10;
      float tau1_eta_DM10_Down = -10;
      float tau1_phi_DM10_Down = -10;
      float tau1_m_DM10_Down   = -10;
      float tau2_pt_DM10_Down  = -10;
      float tau2_eta_DM10_Down = -10;
      float tau2_phi_DM10_Down = -10;
      float tau2_m_DM10_Down   = -10;
      // up UncMet ========================= 
      float tau1_pt_UncMet_Up  = -10;
      float tau1_eta_UncMet_Up = -10;
      float tau1_phi_UncMet_Up = -10;
      float tau1_m_UncMet_Up   = -10;
      float tau2_pt_UncMet_Up  = -10;
      float tau2_eta_UncMet_Up = -10;
      float tau2_phi_UncMet_Up = -10;
      float tau2_m_UncMet_Up   = -10;
      // down 
      float tau1_pt_UncMet_Down  = -10;
      float tau1_eta_UncMet_Down = -10;
      float tau1_phi_UncMet_Down = -10;
      float tau1_m_UncMet_Down   = -10;
      float tau2_pt_UncMet_Down  = -10;
      float tau2_eta_UncMet_Down = -10;
      float tau2_phi_UncMet_Down = -10;
      float tau2_m_UncMet_Down   = -10;
      // up ClusteredMet ========================= 
      float tau1_pt_ClusteredMet_Up  = -10;
      float tau1_eta_ClusteredMet_Up = -10;
      float tau1_phi_ClusteredMet_Up = -10;
      float tau1_m_ClusteredMet_Up   = -10;
      float tau2_pt_ClusteredMet_Up  = -10;
      float tau2_eta_ClusteredMet_Up = -10;
      float tau2_phi_ClusteredMet_Up = -10;
      float tau2_m_ClusteredMet_Up   = -10;
      // down 
      float tau1_pt_ClusteredMet_Down  = -10;
      float tau1_eta_ClusteredMet_Down = -10;
      float tau1_phi_ClusteredMet_Down = -10;
      float tau1_m_ClusteredMet_Down   = -10;
      float tau2_pt_ClusteredMet_Down  = -10;
      float tau2_eta_ClusteredMet_Down = -10;
      float tau2_phi_ClusteredMet_Down = -10;
      float tau2_m_ClusteredMet_Down   = -10;
                                                                                  
      TBranch *newBranch11 = t->Branch("m_sv_Up", &svFitMass_Up, "m_sv_Up/F");
      TBranch *newBranch12 = t->Branch("pt_sv_Up", &svFitPt_Up, "pt_sv_Up/F");
      TBranch *newBranch13 = t->Branch("eta_sv_Up", &svFitEta_Up, "eta_sv_Up/F");
      TBranch *newBranch14 = t->Branch("phi_sv_Up", &svFitPhi_Up, "phi_sv_Up/F");
      TBranch *newBranch15 = t->Branch("met_sv_Up", &svFitMET_Up, "met_sv_Up/F");
      TBranch *newBranch16 = t->Branch("mt_sv_Up", &svFitTransverseMass_Up, "mt_sv_Up/F");

      TBranch *newBranch17 = t->Branch("m_sv_Down", &svFitMass_Down, "m_sv_Down/F");
      TBranch *newBranch18 = t->Branch("pt_sv_Down", &svFitPt_Down, "pt_sv_Down/F");
      TBranch *newBranch19 = t->Branch("eta_sv_Down", &svFitEta_Down, "eta_sv_Down/F");
      TBranch *newBranch20 = t->Branch("phi_sv_Down", &svFitPhi_Down, "phi_sv_Down/F");
      TBranch *newBranch21 = t->Branch("met_sv_Down", &svFitMET_Down, "met_sv_Down/F");
      TBranch *newBranch22 = t->Branch("mt_sv_Down", &svFitTransverseMass_Down, "mt_sv_Down/F");

      TBranch *newBranch23 = t->Branch("m_sv_DM0_Up", &svFitMass_DM0_Up, "m_sv_DM0_Up/F");
      TBranch *newBranch24 = t->Branch("pt_sv_DM0_Up", &svFitPt_DM0_Up, "pt_sv_DM0_Up/F");
      TBranch *newBranch25 = t->Branch("eta_sv_DM0_Up", &svFitEta_DM0_Up, "eta_sv_DM0_Up/F");
      TBranch *newBranch26 = t->Branch("phi_sv_DM0_Up", &svFitPhi_DM0_Up, "phi_sv_DM0_Up/F");
      TBranch *newBranch27 = t->Branch("met_sv_DM0_Up", &svFitMET_DM0_Up, "met_sv_DM0_Up/F");
      TBranch *newBranch28 = t->Branch("mt_sv_DM0_Up", &svFitTransverseMass_DM0_Up, "mt_sv_DM0_Up/F");

      TBranch *newBranch29 = t->Branch("m_sv_DM0_Down", &svFitMass_DM0_Down, "m_sv_DM0_Down/F");
      TBranch *newBranch30 = t->Branch("pt_sv_DM0_Down", &svFitPt_DM0_Down, "pt_sv_DM0_Down/F");
      TBranch *newBranch31 = t->Branch("eta_sv_DM0_Down", &svFitEta_DM0_Down, "eta_sv_DM0_Down/F");
      TBranch *newBranch32 = t->Branch("phi_sv_DM0_Down", &svFitPhi_DM0_Down, "phi_sv_DM0_Down/F");
      TBranch *newBranch33 = t->Branch("met_sv_DM0_Down", &svFitMET_DM0_Down, "met_sv_DM0_Down/F");
      TBranch *newBranch34 = t->Branch("mt_sv_DM0_Down", &svFitTransverseMass_DM0_Down, "mt_sv_DM0_Down/F");

      TBranch *newBranch35 = t->Branch("m_sv_DM1_Up", &svFitMass_DM1_Up, "m_sv_DM1_Up/F");
      TBranch *newBranch36 = t->Branch("pt_sv_DM1_Up", &svFitPt_DM1_Up, "pt_sv_DM1_Up/F");
      TBranch *newBranch37 = t->Branch("eta_sv_DM1_Up", &svFitEta_DM1_Up, "eta_sv_DM1_Up/F");
      TBranch *newBranch38 = t->Branch("phi_sv_DM1_Up", &svFitPhi_DM1_Up, "phi_sv_DM1_Up/F");
      TBranch *newBranch39 = t->Branch("met_sv_DM1_Up", &svFitMET_DM1_Up, "met_sv_DM1_Up/F");
      TBranch *newBranch40 = t->Branch("mt_sv_DM1_Up", &svFitTransverseMass_DM1_Up, "mt_sv_DM1_Up/F");

      TBranch *newBranch41 = t->Branch("m_sv_DM1_Down", &svFitMass_DM1_Down, "m_sv_DM1_Down/F");
      TBranch *newBranch42 = t->Branch("pt_sv_DM1_Down", &svFitPt_DM1_Down, "pt_sv_DM1_Down/F");
      TBranch *newBranch43 = t->Branch("eta_sv_DM1_Down", &svFitEta_DM1_Down, "eta_sv_DM1_Down/F");
      TBranch *newBranch44 = t->Branch("phi_sv_DM1_Down", &svFitPhi_DM1_Down, "phi_sv_DM1_Down/F");
      TBranch *newBranch45 = t->Branch("met_sv_DM1_Down", &svFitMET_DM1_Down, "met_sv_DM1_Down/F");
      TBranch *newBranch46 = t->Branch("mt_sv_DM1_Down", &svFitTransverseMass_DM1_Down, "mt_sv_DM1_Down/F");

      TBranch *newBranch47 = t->Branch("m_sv_DM10_Up", &svFitMass_DM10_Up, "m_sv_DM10_Up/F");
      TBranch *newBranch48 = t->Branch("pt_sv_DM10_Up", &svFitPt_DM10_Up, "pt_sv_DM10_Up/F");
      TBranch *newBranch49 = t->Branch("eta_sv_DM10_Up", &svFitEta_DM10_Up, "eta_sv_DM10_Up/F");
      TBranch *newBranch50 = t->Branch("phi_sv_DM10_Up", &svFitPhi_DM10_Up, "phi_sv_DM10_Up/F");
      TBranch *newBranch51 = t->Branch("met_sv_DM10_Up", &svFitMET_DM10_Up, "met_sv_DM10_Up/F");
      TBranch *newBranch52 = t->Branch("mt_sv_DM10_Up", &svFitTransverseMass_DM10_Up, "mt_sv_DM10_Up/F");

      TBranch *newBranch53 = t->Branch("m_sv_DM10_Down", &svFitMass_DM10_Down, "m_sv_DM10_Down/F");
      TBranch *newBranch54 = t->Branch("pt_sv_DM10_Down", &svFitPt_DM10_Down, "pt_sv_DM10_Down/F");
      TBranch *newBranch55 = t->Branch("eta_sv_DM10_Down", &svFitEta_DM10_Down, "eta_sv_DM10_Down/F");
      TBranch *newBranch56 = t->Branch("phi_sv_DM10_Down", &svFitPhi_DM10_Down, "phi_sv_DM10_Down/F");
      TBranch *newBranch57 = t->Branch("met_sv_DM10_Down", &svFitMET_DM10_Down, "met_sv_DM10_Down/F");
      TBranch *newBranch58 = t->Branch("mt_sv_DM10_Down", &svFitTransverseMass_DM10_Down, "mt_sv_DM10_Down/F");

      TBranch *newBranch59 = t->Branch("m_sv_UncMet_Up", &svFitMass_UncMet_Up, "m_sv_UncMet_Up/F");
      TBranch *newBranch60 = t->Branch("pt_sv_UncMet_Up", &svFitPt_UncMet_Up, "pt_sv_UncMet_Up/F");
      TBranch *newBranch61 = t->Branch("eta_sv_UncMet_Up", &svFitEta_UncMet_Up, "eta_sv_UncMet_Up/F");
      TBranch *newBranch62 = t->Branch("phi_sv_UncMet_Up", &svFitPhi_UncMet_Up, "phi_sv_UncMet_Up/F");
      TBranch *newBranch63 = t->Branch("met_sv_UncMet_Up", &svFitMET_UncMet_Up, "met_sv_UncMet_Up/F");
      TBranch *newBranch64 = t->Branch("mt_sv_UncMet_Up", &svFitTransverseMass_UncMet_Up, "mt_sv_UncMet_Up/F");

      TBranch *newBranch65 = t->Branch("m_sv_UncMet_Down", &svFitMass_UncMet_Down, "m_sv_UncMet_Down/F");
      TBranch *newBranch66 = t->Branch("pt_sv_UncMet_Down", &svFitPt_UncMet_Down, "pt_sv_UncMet_Down/F");
      TBranch *newBranch67 = t->Branch("eta_sv_UncMet_Down", &svFitEta_UncMet_Down, "eta_sv_UncMet_Down/F");
      TBranch *newBranch68 = t->Branch("phi_sv_UncMet_Down", &svFitPhi_UncMet_Down, "phi_sv_UncMet_Down/F");
      TBranch *newBranch69 = t->Branch("met_sv_UncMet_Down", &svFitMET_UncMet_Down, "met_sv_UncMet_Down/F");
      TBranch *newBranch70 = t->Branch("mt_sv_UncMet_Down", &svFitTransverseMass_UncMet_Down, "mt_sv_UncMet_Down/F");

      TBranch *newBranch71 = t->Branch("m_sv_ClusteredMet_Up", &svFitMass_ClusteredMet_Up, "m_sv_ClusteredMet_Up/F");
      TBranch *newBranch72 = t->Branch("pt_sv_ClusteredMet_Up", &svFitPt_ClusteredMet_Up, "pt_sv_ClusteredMet_Up/F");
      TBranch *newBranch73 = t->Branch("eta_sv_ClusteredMet_Up", &svFitEta_ClusteredMet_Up, "eta_sv_ClusteredMet_Up/F");
      TBranch *newBranch74 = t->Branch("phi_sv_ClusteredMet_Up", &svFitPhi_ClusteredMet_Up, "phi_sv_ClusteredMet_Up/F");
      TBranch *newBranch75 = t->Branch("met_sv_ClusteredMet_Up", &svFitMET_ClusteredMet_Up, "met_sv_ClusteredMet_Up/F");
      TBranch *newBranch76 = t->Branch("mt_sv_ClusteredMet_Up", &svFitTransverseMass_ClusteredMet_Up, "mt_sv_ClusteredMet_Up/F");

      TBranch *newBranch77 = t->Branch("m_sv_ClusteredMet_Down", &svFitMass_ClusteredMet_Down, "m_sv_ClusteredMet_Down/F");
      TBranch *newBranch78 = t->Branch("pt_sv_ClusteredMet_Down", &svFitPt_ClusteredMet_Down, "pt_sv_ClusteredMet_Down/F");
      TBranch *newBranch79 = t->Branch("eta_sv_ClusteredMet_Down", &svFitEta_ClusteredMet_Down, "eta_sv_ClusteredMet_Down/F");
      TBranch *newBranch80 = t->Branch("phi_sv_ClusteredMet_Down", &svFitPhi_ClusteredMet_Down, "phi_sv_ClusteredMet_Down/F");
      TBranch *newBranch81 = t->Branch("met_sv_ClusteredMet_Down", &svFitMET_ClusteredMet_Down, "met_sv_ClusteredMet_Down/F");
      TBranch *newBranch82 = t->Branch("mt_sv_ClusteredMet_Down", &svFitTransverseMass_ClusteredMet_Down, "mt_sv_ClusteredMet_Down/F");
    
      TBranch *newBranch83 = t->Branch("metcorClusteredDown",    &metcorClusteredDown,   "metcorClusteredDown/F");
      TBranch *newBranch84 = t->Branch("metcorphiClusteredDown", &metcorphiClusteredDown,"metcorphiClusteredDown/F");
      TBranch *newBranch85 = t->Branch("metcorClusteredUp",      &metcorClusteredUp,     "metcorClusteredUp/F");
      TBranch *newBranch86 = t->Branch("metcorphiClusteredUp",   &metcorphiClusteredUp,  "metcorphiClusteredUp/F");
      TBranch *newBranch87 = t->Branch("metcorUncDown",          &metcorUncDown,         "metcorUncDown/F");
      TBranch *newBranch89 = t->Branch("metcorphiUncDown",       &metcorphiUncDown,      "metcorphiUncDown/F");
      TBranch *newBranch90 = t->Branch("metcorUncUp",            &metcorUncUp,           "metcorUncUp/F");
      TBranch *newBranch91 = t->Branch("metcorphiUncUp",         &metcorphiUncUp,        "metcorphiUncUp/F");
    
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
      // up
      tau4VectorBranches.push_back(t->Branch("tau1_pt_Up",  &tau1_pt_Up,  "tau1_pt_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_Up", &tau1_eta_Up, "tau1_eta_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_Up", &tau1_phi_Up, "tau1_phi_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_Up",   &tau1_m_Up,   "tau1_m_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_Up",  &tau2_pt_Up,  "tau2_pt_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_Up", &tau2_eta_Up, "tau2_eta_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_Up", &tau2_phi_Up, "tau2_phi_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_Up",   &tau2_m_Up,   "tau2_m_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_Down",  &tau1_pt_Down,  "tau1_pt_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_Down", &tau1_eta_Down, "tau1_eta_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_Down", &tau1_phi_Down, "tau1_phi_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_Down",   &tau1_m_Down,   "tau1_m_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_Down",  &tau2_pt_Down,  "tau2_pt_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_Down", &tau2_eta_Down, "tau2_eta_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_Down", &tau2_phi_Down, "tau2_phi_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_Down",   &tau2_m_Down,   "tau2_m_Down/F"));
      // up DM0
      tau4VectorBranches.push_back(t->Branch("tau1_pt_DM0_Up",  &tau1_pt_DM0_Up,  "tau1_pt_DM0_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_DM0_Up", &tau1_eta_DM0_Up, "tau1_eta_DM0_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_DM0_Up", &tau1_phi_DM0_Up, "tau1_phi_DM0_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_DM0_Up",   &tau1_m_DM0_Up,   "tau1_m_DM0_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_DM0_Up",  &tau2_pt_DM0_Up,  "tau2_pt_DM0_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_DM0_Up", &tau2_eta_DM0_Up, "tau2_eta_DM0_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_DM0_Up", &tau2_phi_DM0_Up, "tau2_phi_DM0_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_DM0_Up",   &tau2_m_DM0_Up,   "tau2_m_DM0_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_DM0_Down",  &tau1_pt_DM0_Down,  "tau1_pt_DM0_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_DM0_Down", &tau1_eta_DM0_Down, "tau1_eta_DM0_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_DM0_Down", &tau1_phi_DM0_Down, "tau1_phi_DM0_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_DM0_Down",   &tau1_m_DM0_Down,   "tau1_m_DM0_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_DM0_Down",  &tau2_pt_DM0_Down,  "tau2_pt_DM0_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_DM0_Down", &tau2_eta_DM0_Down, "tau2_eta_DM0_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_DM0_Down", &tau2_phi_DM0_Down, "tau2_phi_DM0_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_DM0_Down",   &tau2_m_DM0_Down,   "tau2_m_DM0_Down/F"));
      // up DM1
      tau4VectorBranches.push_back(t->Branch("tau1_pt_DM1_Up",  &tau1_pt_DM1_Up,  "tau1_pt_DM1_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_DM1_Up", &tau1_eta_DM1_Up, "tau1_eta_DM1_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_DM1_Up", &tau1_phi_DM1_Up, "tau1_phi_DM1_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_DM1_Up",   &tau1_m_DM1_Up,   "tau1_m_DM1_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_DM1_Up",  &tau2_pt_DM1_Up,  "tau2_pt_DM1_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_DM1_Up", &tau2_eta_DM1_Up, "tau2_eta_DM1_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_DM1_Up", &tau2_phi_DM1_Up, "tau2_phi_DM1_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_DM1_Up",   &tau2_m_DM1_Up,   "tau2_m_DM1_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_DM1_Down",  &tau1_pt_DM1_Down,  "tau1_pt_DM1_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_DM1_Down", &tau1_eta_DM1_Down, "tau1_eta_DM1_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_DM1_Down", &tau1_phi_DM1_Down, "tau1_phi_DM1_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_DM1_Down",   &tau1_m_DM1_Down,   "tau1_m_DM1_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_DM1_Down",  &tau2_pt_DM1_Down,  "tau2_pt_DM1_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_DM1_Down", &tau2_eta_DM1_Down, "tau2_eta_DM1_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_DM1_Down", &tau2_phi_DM1_Down, "tau2_phi_DM1_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_DM1_Down",   &tau2_m_DM1_Down,   "tau2_m_DM1_Down/F"));
      // up DM10
    
      tau4VectorBranches.push_back(t->Branch("tau1_pt_DM10_Up",  &tau1_pt_DM10_Up,  "tau1_pt_DM10_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_DM10_Up", &tau1_eta_DM10_Up, "tau1_eta_DM10_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_DM10_Up", &tau1_phi_DM10_Up, "tau1_phi_DM10_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_DM10_Up",   &tau1_m_DM10_Up,   "tau1_m_DM10_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_DM10_Up",  &tau2_pt_DM10_Up,  "tau2_pt_DM10_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_DM10_Up", &tau2_eta_DM10_Up, "tau2_eta_DM10_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_DM10_Up", &tau2_phi_DM10_Up, "tau2_phi_DM10_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_DM10_Up",   &tau2_m_DM10_Up,   "tau2_m_DM10_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_DM10_Down",  &tau1_pt_DM10_Down,  "tau1_pt_DM10_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_DM10_Down", &tau1_eta_DM10_Down, "tau1_eta_DM10_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_DM10_Down", &tau1_phi_DM10_Down, "tau1_phi_DM10_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_DM10_Down",   &tau1_m_DM10_Down,   "tau1_m_DM10_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_DM10_Down",  &tau2_pt_DM10_Down,  "tau2_pt_DM10_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_DM10_Down", &tau2_eta_DM10_Down, "tau2_eta_DM10_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_DM10_Down", &tau2_phi_DM10_Down, "tau2_phi_DM10_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_DM10_Down",   &tau2_m_DM10_Down,   "tau2_m_DM10_Down/F"));
      // up UncMet
      tau4VectorBranches.push_back(t->Branch("tau1_pt_UncMet_Up",  &tau1_pt_UncMet_Up,  "tau1_pt_UncMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_UncMet_Up", &tau1_eta_UncMet_Up, "tau1_eta_UncMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_UncMet_Up", &tau1_phi_UncMet_Up, "tau1_phi_UncMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_UncMet_Up",   &tau1_m_UncMet_Up,   "tau1_m_UncMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_UncMet_Up",  &tau2_pt_UncMet_Up,  "tau2_pt_UncMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_UncMet_Up", &tau2_eta_UncMet_Up, "tau2_eta_UncMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_UncMet_Up", &tau2_phi_UncMet_Up, "tau2_phi_UncMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_UncMet_Up",   &tau2_m_UncMet_Up,   "tau2_m_UncMet_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_UncMet_Down",  &tau1_pt_UncMet_Down,  "tau1_pt_UncMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_UncMet_Down", &tau1_eta_UncMet_Down, "tau1_eta_UncMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_UncMet_Down", &tau1_phi_UncMet_Down, "tau1_phi_UncMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_UncMet_Down",   &tau1_m_UncMet_Down,   "tau1_m_UncMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_UncMet_Down",  &tau2_pt_UncMet_Down,  "tau2_pt_UncMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_UncMet_Down", &tau2_eta_UncMet_Down, "tau2_eta_UncMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_UncMet_Down", &tau2_phi_UncMet_Down, "tau2_phi_UncMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_UncMet_Down",   &tau2_m_UncMet_Down,   "tau2_m_UncMet_Down/F"));
      // up ClusteredMet
      tau4VectorBranches.push_back(t->Branch("tau1_pt_ClusteredMet_Up",  &tau1_pt_ClusteredMet_Up,  "tau1_pt_ClusteredMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_ClusteredMet_Up", &tau1_eta_ClusteredMet_Up, "tau1_eta_ClusteredMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_ClusteredMet_Up", &tau1_phi_ClusteredMet_Up, "tau1_phi_ClusteredMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_ClusteredMet_Up",   &tau1_m_ClusteredMet_Up,   "tau1_m_ClusteredMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_ClusteredMet_Up",  &tau2_pt_ClusteredMet_Up,  "tau2_pt_ClusteredMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_ClusteredMet_Up", &tau2_eta_ClusteredMet_Up, "tau2_eta_ClusteredMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_ClusteredMet_Up", &tau2_phi_ClusteredMet_Up, "tau2_phi_ClusteredMet_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_ClusteredMet_Up",   &tau2_m_ClusteredMet_Up,   "tau2_m_ClusteredMet_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_ClusteredMet_Down",  &tau1_pt_ClusteredMet_Down,  "tau1_pt_ClusteredMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_ClusteredMet_Down", &tau1_eta_ClusteredMet_Down, "tau1_eta_ClusteredMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_ClusteredMet_Down", &tau1_phi_ClusteredMet_Down, "tau1_phi_ClusteredMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_ClusteredMet_Down",   &tau1_m_ClusteredMet_Down,   "tau1_m_ClusteredMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_ClusteredMet_Down",  &tau2_pt_ClusteredMet_Down,  "tau2_pt_ClusteredMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_ClusteredMet_Down", &tau2_eta_ClusteredMet_Down, "tau2_eta_ClusteredMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_ClusteredMet_Down", &tau2_phi_ClusteredMet_Down, "tau2_phi_ClusteredMet_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_ClusteredMet_Down",   &tau2_m_ClusteredMet_Down,   "tau2_m_ClusteredMet_Down/F"));
      std::cout << "That's a lot of tau 4-vector branches! N = " << tau4VectorBranches.size() << std::endl;
    
      unsigned long long evt;
      int run, lumi;
      float pt1;
      float eta1;
      float phi1;
      int gen_match_1;
      float pt2;
      float eta2;
      float phi2;
      float m2;
      int gen_match_2;
      float decayMode=-999.;
      float decayMode2;
      float mvaCovMatrix00;
      float mvaCovMatrix10;
      float mvaCovMatrix01;
      float mvaCovMatrix11;
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
      float mvamet;
      float mvametphi;
      float pfmet;
      float pfmetphi;
      TLorentzVector TMet(0,0,0,0);
      // define MET covariance
      TMatrixD covMET(2, 2);

      // MET Uncertainties
      float uncMetPtUp;
      float uncMetPtDown;
      float uncMetPhiUp;
      float uncMetPhiDown;
      float clusteredMetPtUp;
      float clusteredMetPtDown;
      float clusteredMetPhiUp;
      float clusteredMetPhiDown;
      double uncMetUpMETx = 0.;
      double uncMetUpMETy = 0.;
      double uncMetDownMETx = 0.;
      double uncMetDownMETy = 0.;
      double clusteredMetUpMETx = 0.;
      double clusteredMetUpMETy = 0.;
      double clusteredMetDownMETx = 0.;
      double clusteredMetDownMETy = 0.;
      
      //ele/mu variables
      TBranch *pt1branch;
      
      t->SetBranchAddress("evt",&evt);
      t->SetBranchAddress("run",&run);
      t->SetBranchAddress("lumi",&lumi);
      t->SetBranchAddress("gen_match_1",&gen_match_1);
      t->SetBranchAddress("gen_match_2",&gen_match_2);
      if ( channel == "tt" ) t->SetBranchAddress("t1_decayMode", &decayMode);
      if ( channel == "tt" ) t->SetBranchAddress("t2_decayMode", &decayMode2);
      t->SetBranchAddress("pt_1",&pt1,&pt1branch);
      t->SetBranchAddress("eta_1",&eta1);
      t->SetBranchAddress("phi_1",&phi1);
      t->SetBranchAddress("pt_2",&pt2);
      t->SetBranchAddress("eta_2",&eta2);
      t->SetBranchAddress("phi_2",&phi2);
      t->SetBranchAddress("m_2",&m2);
      //t->SetBranchAddress("l1_decayMode",&decayMode);
      if ( channel != "tt" ) t->SetBranchAddress("l2_decayMode",&decayMode2);
      t->SetBranchAddress("mvacov00",&mvaCovMatrix00);  // branch not stored in et/mt trees
      t->SetBranchAddress("mvacov01",&mvaCovMatrix01);  // branch not stored in et/mt trees
      t->SetBranchAddress("mvacov10",&mvaCovMatrix10);  // branch not stored in et/mt trees
      t->SetBranchAddress("mvacov11",&mvaCovMatrix11);  // branch not stored in et/mt trees
      t->SetBranchAddress("mvamet",&mvamet);            // branch not stored in et/mt trees
      t->SetBranchAddress("mvametphi",&mvametphi);      // branch not stored in et/mt trees
      t->SetBranchAddress("njets", &njets);
      t->SetBranchAddress("met",&pfmet);
      t->SetBranchAddress("metphi",&pfmetphi);
      // FOR PF MET ANALYSIS
      t->SetBranchAddress("metcov00",&pfCovMatrix00);
      t->SetBranchAddress("metcov01",&pfCovMatrix01);
      t->SetBranchAddress("metcov10",&pfCovMatrix10);
      t->SetBranchAddress("metcov11",&pfCovMatrix11);
      // Met Unc
      if ( channel ==  "tt" ) {
        t->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnUp",&uncMetPtUp);
        t->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnDown",&uncMetPtDown);
        t->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnUp",&uncMetPhiUp);
        t->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnDown",&uncMetPhiDown);
        t->SetBranchAddress("type1_pfMet_shiftedPt_JetEnUp",&clusteredMetPtUp);
        t->SetBranchAddress("type1_pfMet_shiftedPt_JetEnDown",&clusteredMetPtDown);
        t->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnUp",&clusteredMetPhiUp);
        t->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnDown",&clusteredMetPhiDown);
      } else  {
        t->SetBranchAddress("met_UESUp", &uncMetPtUp);
        t->SetBranchAddress("met_UESDown", &uncMetPtDown);
        t->SetBranchAddress("metphi_UESUp", &uncMetPhiUp);
        t->SetBranchAddress("metphi_UESDown", &uncMetPhiDown);

        t->SetBranchAddress("met_JESUp", &clusteredMetPtUp);
        t->SetBranchAddress("met_JESDown", &clusteredMetPtDown);
        t->SetBranchAddress("metphi_JESUp", &clusteredMetPhiUp);
        t->SetBranchAddress("metphi_JESDown", &clusteredMetPhiDown);
      }
      
      printf("Found tree -> weighting\n");
    
      double tesUp = 1.0 + tesSize;
      double tesDown = 1.0 - tesSize;

      int nevents = t->GetEntries();
      if ( parser.integerValue("numEvents") != -1 ) nevents = parser.integerValue("numEvents");
      for(Int_t i=0;i<nevents;++i){
        t->GetEntry(i);

        // Using PF Met or Mva Met?
        if (metType == 1) { // 1 = Mva Met
          std::cerr << "Only PF Met is currently supported. You're in for a world full of problems if you continue." << std::endl;
          TMet.SetPtEtaPhiM(mvamet,0,mvametphi,0);
          measuredMETx = mvamet*TMath::Cos(mvametphi);
          measuredMETy = mvamet*TMath::Sin(mvametphi);
          
          covMET[0][0] =  mvaCovMatrix00;
          covMET[1][0] =  mvaCovMatrix10;
          covMET[0][1] =  mvaCovMatrix01;
          covMET[1][1] =  mvaCovMatrix11;
        } // mva met
        if (metType == -1) { // -1 = PF Met
          TMet.SetPtEtaPhiM(pfmet,0,pfmetphi,0);
          measuredMETx = pfmet*TMath::Cos(pfmetphi);
          measuredMETy = pfmet*TMath::Sin(pfmetphi);
          // Shifted METs
          uncMetUpMETx         = uncMetPtUp*TMath::Cos(uncMetPhiUp);
          uncMetUpMETy         = uncMetPtUp*TMath::Sin(uncMetPhiUp);
          uncMetDownMETx       = uncMetPtDown*TMath::Cos(uncMetPhiDown);
          uncMetDownMETy       = uncMetPtDown*TMath::Sin(uncMetPhiDown);
          clusteredMetUpMETx   = clusteredMetPtUp*TMath::Cos(clusteredMetPhiUp);
          clusteredMetUpMETy   = clusteredMetPtUp*TMath::Sin(clusteredMetPhiUp);
          clusteredMetDownMETx = clusteredMetPtDown*TMath::Cos(clusteredMetPhiDown);
          clusteredMetDownMETy = clusteredMetPtDown*TMath::Sin(clusteredMetPhiDown);
          
          covMET[0][0] =  pfCovMatrix00;
          covMET[1][0] =  pfCovMatrix10;
          covMET[0][1] =  pfCovMatrix01;
          covMET[1][1] =  pfCovMatrix11;
        } // pf met
     
        metcorr_ex = measuredMETx;
        metcorr_ey = measuredMETy;
        metcorrUncUp_ex = uncMetUpMETx;
        metcorrUncUp_ey = uncMetUpMETy;
        metcorrUncDown_ex = uncMetDownMETx;
        metcorrUncDown_ey = uncMetDownMETy;
        metcorrClusteredUp_ex = clusteredMetUpMETx;
        metcorrClusteredUp_ey = clusteredMetUpMETy;
        metcorrClusteredDown_ex = clusteredMetDownMETx;
        metcorrClusteredDown_ey = clusteredMetDownMETy;
        
        metcor = TMath::Sqrt( metcorr_ex*metcorr_ex + metcorr_ey*metcorr_ey);
        metcorphi = TMath::ATan2( metcorr_ey, metcorr_ex );
        std::cout << " - metcor "<<metcor<<" metcorphi "<<metcorphi<<std::endl;
        
        // Corrected MET values for saving
        // Will be re-corrected with TEC if running tautau channel
        metcorClusteredDown = TMath::Sqrt( metcorrClusteredDown_ex*metcorrClusteredDown_ex + metcorrClusteredDown_ey*metcorrClusteredDown_ey);
        metcorphiClusteredDown = TMath::ATan2( metcorrClusteredDown_ey, metcorrClusteredDown_ex );
        
        metcorClusteredUp = TMath::Sqrt( metcorrClusteredUp_ex*metcorrClusteredUp_ex + metcorrClusteredUp_ey*metcorrClusteredUp_ey);
        metcorphiClusteredUp = TMath::ATan2( metcorrClusteredUp_ey, metcorrClusteredUp_ex );
        
        metcorUncDown = TMath::Sqrt( metcorrUncDown_ex*metcorrUncDown_ex + metcorrUncDown_ey*metcorrUncDown_ey);
        metcorphiUncDown = TMath::ATan2( metcorrUncDown_ey, metcorrUncDown_ex );
        
        metcorUncUp = TMath::Sqrt( metcorrUncUp_ex*metcorrUncUp_ex + metcorrUncUp_ey*metcorrUncUp_ey);
        metcorphiUncUp = TMath::ATan2( metcorrUncUp_ey, metcorrUncUp_ex );
     
        if( channel == "mt" || channel == "et" ) {
          
          mass2 = m2;
          std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons{
              classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1),
              classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2)
          };
 
          std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
          //modified
          runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, 
               svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
          std::cout<<"finished runningSVFit"<<std::endl;

          // MET systematics
          std::cout << "MET Unclustered Energy Up   ---  ";
          runSVFit(measuredTauLeptons, metcorrUncUp_ex, metcorrUncUp_ey, covMET, 0, svFitMass_UncMet_Up, svFitPt_UncMet_Up, svFitEta_UncMet_Up, 
              svFitPhi_UncMet_Up, svFitMET_UncMet_Up, svFitTransverseMass_UncMet_Up, tau1_UncMet_Up, tau2_UncMet_Up);

          std::cout << "MET Unclustered Energy Down ---  ";
          runSVFit(measuredTauLeptons, metcorrUncDown_ex, metcorrUncDown_ey, covMET, 0, svFitMass_UncMet_Down, svFitPt_UncMet_Down, svFitEta_UncMet_Down, 
              svFitPhi_UncMet_Down, svFitMET_UncMet_Down, svFitTransverseMass_UncMet_Down, tau1_UncMet_Down, tau2_UncMet_Down);

          std::cout << "MET Clustered Energy Up     ---  ";
          runSVFit(measuredTauLeptons, metcorrClusteredUp_ex, metcorrClusteredUp_ey, covMET, 0, svFitMass_ClusteredMet_Up, svFitPt_ClusteredMet_Up, 
              svFitEta_ClusteredMet_Up, svFitPhi_ClusteredMet_Up, svFitMET_ClusteredMet_Up, svFitTransverseMass_ClusteredMet_Up, tau1_ClusteredMet_Up, tau2_ClusteredMet_Up);

          std::cout << "MET Clustered Energy Down   ---  ";
          runSVFit(measuredTauLeptons, metcorrClusteredDown_ex, metcorrClusteredDown_ey, covMET, 0, svFitMass_ClusteredMet_Down, svFitPt_ClusteredMet_Down, 
              svFitEta_ClusteredMet_Down, svFitPhi_ClusteredMet_Down, svFitMET_ClusteredMet_Down, svFitTransverseMass_ClusteredMet_Down, tau1_ClusteredMet_Down, tau2_ClusteredMet_Down);

          std::cout<< "Shifted MET Summary:\n\tmetcorr_ex " << metcorr_ex << " --- metcorrUncUp_ex " << metcorrUncUp_ex << " metcorrUncDown_ex " << metcorrUncDown_ex
              << " metcorrClusteredUp_ex " << metcorrClusteredUp_ex << " metcorrClusteredDown_ex " << metcorrClusteredDown_ex << std::endl;
          std::cout<< "\tmetcorr_ey " << metcorr_ey << " --- metcorrUncUp_ey " << metcorrUncUp_ey << " metcorrUncDown_ey " << metcorrUncDown_ey
              << " metcorrClusteredUp_ey " << metcorrClusteredUp_ey << " metcorrClusteredDown_ey " << metcorrClusteredDown_ey << std::endl;

          // new met to store 
          metcor = TMath::Sqrt( metcorr_ex*metcorr_ex + metcorr_ey*metcorr_ey);
          metcorphi = TMath::ATan2( metcorr_ey, metcorr_ex );
          metcorClusteredDown = TMath::Sqrt( metcorrClusteredDown_ex*metcorrClusteredDown_ex + metcorrClusteredDown_ey*metcorrClusteredDown_ey);
          metcorphiClusteredDown = TMath::ATan2( metcorrClusteredDown_ey, metcorrClusteredDown_ex );
          metcorClusteredUp = TMath::Sqrt( metcorrClusteredUp_ex*metcorrClusteredUp_ex + metcorrClusteredUp_ey*metcorrClusteredUp_ey);
          metcorphiClusteredUp = TMath::ATan2( metcorrClusteredUp_ey, metcorrClusteredUp_ex );
          metcorUncDown = TMath::Sqrt( metcorrUncDown_ex*metcorrUncDown_ex + metcorrUncDown_ey*metcorrUncDown_ey);
          metcorphiUncDown = TMath::ATan2( metcorrUncDown_ey, metcorrUncDown_ex );
          metcorUncUp = TMath::Sqrt( metcorrUncUp_ex*metcorrUncUp_ex + metcorrUncUp_ey*metcorrUncUp_ey);
          metcorphiUncUp = TMath::ATan2( metcorrUncUp_ey, metcorrUncUp_ex );
          
          if(doES) {
            
            // corrections only need to be done once
            float ES_Up( 1. ), ES_Down( 1. ); // shift TES
            if (gen_match_2 == 5) { // 0.6% uncertainty on hadronic tau
              ES_Up = 1.012;
              ES_Down = 0.988;
            } else if (gen_match_2 < 5) { // 1.0% uncertainty on gen matched ele/mu
              ES_Up = 1.01;
              ES_Down = 0.99;
            }
            
            double pt_Up  ( pt2 * ES_Up ), pt_Down( pt2 * ES_Down ); // shift tau pT by energy scale
            double dx_Up  ( pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Up ) - 1.) ),
                   dy_Up  ( pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Up ) - 1.) ),
                   dx_Down( pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Down ) - 1.) ),
                   dy_Down( pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Down ) - 1.) ); 
            double metcorr_ex_Up  ( metcorr_ex + dx_Up ), 
                   metcorr_ey_Up  ( metcorr_ey + dy_Up ),
                   metcorr_ex_Down( metcorr_ex + dx_Down ),
                   metcorr_ey_Down( metcorr_ey + dy_Down );

            // leptons shifted up
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp{
                classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1),
                classic_svFit::MeasuredTauLepton(decayType2,  pt_Up, eta2, phi2,  mass2, decayMode2)
            };

            // leptons shifted down
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown{
                classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1),
                classic_svFit::MeasuredTauLepton(decayType2,  pt_Down, eta2, phi2,  mass2, decayMode2)
            };
          
            /////////////////////////////
            // All upward shifts below //
            /////////////////////////////

            // all tau shift up
            if (gen_match_2 < 6) {
              std::cout << "Tau shift up" << std::endl;
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_Up, svFitPt_Up, 
                    svFitEta_Up, svFitPhi_Up, svFitMET_Up, svFitTransverseMass_Up, tau1_Up, tau2_Up);
            } else {
              svFitMass_Up=svFitMass;
              svFitPt_Up=svFitPt;
              svFitEta_Up=svFitEta;
              svFitPhi_Up=svFitPhi;
              svFitMET_Up=svFitMET;
              svFitTransverseMass_Up=svFitTransverseMass;
              tau1_Up = tau1;
              tau2_Up = tau2;
            }

            // tau DM0 shifted up
            if (gen_match_2 == 5 && decayMode2 == 0) {
              std::cout << "DM0 shift up" << std::endl;
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM0_Up, svFitPt_DM0_Up, 
                    svFitEta_DM0_Up, svFitPhi_DM0_Up, svFitMET_DM0_Up, svFitTransverseMass_DM0_Up, tau1_DM0_Up, tau2_DM0_Up);

            } else {
              svFitMass_DM0_Up=svFitMass;
              svFitPt_DM0_Up=svFitPt;
              svFitEta_DM0_Up=svFitEta;
              svFitPhi_DM0_Up=svFitPhi;
              svFitMET_DM0_Up=svFitMET;
              svFitTransverseMass_DM0_Up=svFitTransverseMass;
              tau1_DM0_Up = tau1;
              tau2_DM0_Up = tau2;
            }

            // tau DM1 shifted up
            if (gen_match_2 == 5 && decayMode2 == 1) {
              std::cout << "DM1 shift up" << std::endl;
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM1_Up, svFitPt_DM1_Up, 
                    svFitEta_DM1_Up, svFitPhi_DM1_Up, svFitMET_DM1_Up, svFitTransverseMass_DM1_Up, tau1_DM1_Up, tau2_DM1_Up);
            } else {
              svFitMass_DM1_Up=svFitMass;
              svFitPt_DM1_Up=svFitPt;
              svFitEta_DM1_Up=svFitEta;
              svFitPhi_DM1_Up=svFitPhi;
              svFitMET_DM1_Up=svFitMET;
              svFitTransverseMass_DM1_Up=svFitTransverseMass;
              tau1_DM1_Up = tau1;
              tau2_DM1_Up = tau2;
            }

            // tau DM10 shifted up
            if (gen_match_2 == 5 && decayMode2 == 10) {
              std::cout << "DM10 shift up" << std::endl;
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM10_Up, svFitPt_DM10_Up, 
                    svFitEta_DM10_Up, svFitPhi_DM10_Up, svFitMET_DM10_Up, svFitTransverseMass_DM10_Up, tau1_DM10_Up, tau2_DM10_Up);
            } else {
              svFitMass_DM10_Up=svFitMass;
              svFitPt_DM10_Up=svFitPt;
              svFitEta_DM10_Up=svFitEta;
              svFitPhi_DM10_Up=svFitPhi;
              svFitMET_DM10_Up=svFitMET;
              svFitTransverseMass_DM10_Up=svFitTransverseMass;
              tau1_DM10_Up = tau1;
              tau2_DM10_Up = tau2;
            }

            ///////////////////////////////
            // All downward shifts below //
            ///////////////////////////////

            // all tau shift down
            if (gen_match_2 < 6) {
              std::cout << "Tau shift down" << std::endl;
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_Down, svFitPt_Down, 
                    svFitEta_Down, svFitPhi_Down, svFitMET_Down, svFitTransverseMass_Down, tau1_Down, tau2_Down);
            } else {
              svFitMass_Down=svFitMass;
              svFitPt_Down=svFitPt;
              svFitEta_Down=svFitEta;
              svFitPhi_Down=svFitPhi;
              svFitMET_Down=svFitMET;
              svFitTransverseMass_Down=svFitTransverseMass;
              tau1_Down = tau1;
              tau2_Down = tau2;
            }

            // tau DM0 shifted down
            if (gen_match_2 == 5 && decayMode2 == 0) {
              std::cout << "DM0 shift Down" << std::endl;
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM0_Down, svFitPt_DM0_Down, 
                    svFitEta_DM0_Down, svFitPhi_DM0_Down, svFitMET_DM0_Down, svFitTransverseMass_DM0_Down, tau1_DM0_Down, tau2_DM0_Down);
            } else {
              svFitMass_DM0_Down=svFitMass;
              svFitPt_DM0_Down=svFitPt;
              svFitEta_DM0_Down=svFitEta;
              svFitPhi_DM0_Down=svFitPhi;
              svFitMET_DM0_Down=svFitMET;
              svFitTransverseMass_DM0_Down=svFitTransverseMass;
              tau1_DM0_Down = tau1;
              tau2_DM0_Down = tau2;
            }

            // tau DM1 shifted down
            if (gen_match_2 == 5 && decayMode2 == 1) {
              std::cout << "DM1 shift down" << std::endl;
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM1_Down, svFitPt_DM1_Down, 
                    svFitEta_DM1_Down, svFitPhi_DM1_Down, svFitMET_DM1_Down, svFitTransverseMass_DM1_Down, tau1_DM1_Down, tau2_DM1_Down);
            } else {
              svFitMass_DM1_Down=svFitMass;
              svFitPt_DM1_Down=svFitPt;
              svFitEta_DM1_Down=svFitEta;
              svFitPhi_DM1_Down=svFitPhi;
              svFitMET_DM1_Down=svFitMET;
              svFitTransverseMass_DM1_Down=svFitTransverseMass;
              tau1_DM1_Down = tau1;
              tau2_DM1_Down = tau2;
            }

            // tau DM10 shifted down
            if (gen_match_2 == 5 && decayMode2 == 10) {
              std::cout << "DM10 shift down" << std::endl;
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM10_Down, svFitPt_DM10_Down, 
                    svFitEta_DM10_Down, svFitPhi_DM10_Down, svFitMET_DM10_Down, svFitTransverseMass_DM10_Down, tau1_DM10_Down, tau2_DM10_Down);
            } else {
              svFitMass_DM10_Down=svFitMass;
              svFitPt_DM10_Down=svFitPt;
              svFitEta_DM10_Down=svFitEta;
              svFitPhi_DM10_Down=svFitPhi;
              svFitMET_DM10_Down=svFitMET;
              svFitTransverseMass_DM10_Down=svFitTransverseMass;
              tau1_DM10_Down = tau1;
              tau2_DM10_Down = tau2;
            }
          } // end doES
        } // eTau / muTau
     
     else if(channel=="em"){
       // define lepton four vectors
       std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
       measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
       measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2, mass2));
       
       std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
       runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
       std::cout<<"finished running non-EES SVFit in EMu"<<std::endl;
       if(doES) {
         
         // Only shift Electrons
         // 1% shift in barrel and 2.5% shift in endcap
         // applied to electrons in emu channel
         float etaBarrelEndcap  = 1.479;
         float ES_Up_scale;
         if (abs(eta1) < etaBarrelEndcap) ES_Up_scale = 1.01;
         else ES_Up_scale = 1.025;
         double pt1_Up = pt1 * ES_Up_scale;
         std::cout << "E eta: " << eta1 << " ees SF: " << ES_Up_scale << std::endl;
         double metcorr_ex_Up, metcorr_ey_Up;
         double dx1_Up, dy1_Up;
         dx1_Up = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Up_scale ) - 1.);
         dy1_Up = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Up_scale ) - 1.);
         metcorr_ex_Up = metcorr_ex + dx1_Up;
         metcorr_ey_Up = metcorr_ey + dy1_Up;
         std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_Up <<std::endl;
         std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_Up <<std::endl;
         std::cout << "pt1 " << pt1 << "  pt1_Up " << pt1_Up <<std::endl;
         std::cout << "metcor_ex " << metcorr_ex << " ees: " << metcorr_ex_Up << std::endl;
         std::cout << "metcor_ey " << metcorr_ey << " ees: " << metcorr_ey_Up << std::endl;
         
         std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1));
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2, mass2));
         
         std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1_Up << " mass1 " << mass1 << " pt2: "<< pt2 << " mass2: "<< mass2 <<std::endl;        
         runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_Up, svFitPt_Up, svFitEta_Up, svFitPhi_Up, svFitMET_Up, svFitTransverseMass_Up, tau1_Up, tau2_Up);
         std::cout<<"finished runningSVFit in EMu EES Up"<<std::endl;
         
         
         // EES Down
         float ES_Down_scale;
         if (abs(eta1) < etaBarrelEndcap) ES_Down_scale = 0.99;
         else ES_Down_scale = 0.975;
         std::cout << "E eta: " << eta1 << " ees SF: " << ES_Down_scale << std::endl;
         double pt1_Down;
         pt1_Down = pt1 * ES_Down_scale;
         double metcorr_ex_Down, metcorr_ey_Down;
         double dx1_Down, dy1_Down;
         dx1_Down = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Down_scale ) - 1.);
         dy1_Down = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Down_scale ) - 1.);
         metcorr_ex_Down = metcorr_ex + dx1_Down;
         metcorr_ey_Down = metcorr_ey + dy1_Down;
         std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_Down <<std::endl;
         std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_Down <<std::endl;
         std::cout << "metcor_ex " << metcorr_ex << " ees: " << metcorr_ex_Down << std::endl;
         std::cout << "metcor_ey " << metcorr_ey << " ees: " << metcorr_ey_Down << std::endl;
         
         std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2));
         
         runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_Down, svFitPt_Down, svFitEta_Down, svFitPhi_Down, svFitMET_Down, svFitTransverseMass_Down, tau1_Down, tau2_Down);
         
         
       } // end doES
     } // eMu
     
     
     
     else if(channel=="tt"){
       
       float shiftTau1=1.0;
       float shiftTau2=1.0;
       if (gen_match_1==5 && decayMode==0) shiftTau1=0.982;
       if (gen_match_1==5 && decayMode==1) shiftTau1=1.010;
       if (gen_match_1==5 && decayMode==10) shiftTau1=1.004;
       if (gen_match_2==5 && decayMode2==0) shiftTau2=0.982;
       if (gen_match_2==5 && decayMode2==1) shiftTau2=1.010;
       if (gen_match_2==5 && decayMode2==10) shiftTau2=1.004;
       pt1 = pt1 * shiftTau1;
       pt2 = pt2 * shiftTau2;
       double dx1, dy1;
       double dx2, dy2;
       dx1 = pt1 * TMath::Cos( phi1 ) * (( 1. / shiftTau1 ) - 1.);
       dy1 = pt1 * TMath::Sin( phi1 ) * (( 1. / shiftTau1 ) - 1.);
       dx2 = pt2 * TMath::Cos( phi2 ) * (( 1. / shiftTau2 ) - 1.);
       dy2 = pt2 * TMath::Sin( phi2 ) * (( 1. / shiftTau2 ) - 1.);
       metcorr_ex = metcorr_ex + dx1;
       metcorr_ey = metcorr_ey + dy1;
       metcorr_ex = metcorr_ex + dx2;
       metcorr_ey = metcorr_ey + dy2;
       
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
       runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
       std::cout<<"finished running non-TES SVFit in TT"<<std::endl;
       
       
       //*****************************************************
       // MET SYSTEMATICS
       // Taus Pt have been corrected (TEC)
       // Still need TEC propagated to shifted METs
       //*****************************************************
       
       std::cout << "MET Unclustered Energy Up   ---  ";
       metcorrUncUp_ex = metcorrUncUp_ex + dx1 + dx2;
       metcorrUncUp_ey = metcorrUncUp_ey + dy1 + dy2;
       runSVFit(measuredTauLeptons, metcorrUncUp_ex, metcorrUncUp_ey, covMET, 0, svFitMass_UncMet_Up, svFitPt_UncMet_Up, svFitEta_UncMet_Up, svFitPhi_UncMet_Up, svFitMET_UncMet_Up, svFitTransverseMass_UncMet_Up, tau1_UncMet_Up, tau2_UncMet_Up);
       std::cout << "MET Unclustered Energy Down ---  ";
       metcorrUncDown_ex = metcorrUncDown_ex + dx1 + dx2;
       metcorrUncDown_ey = metcorrUncDown_ey + dy1 + dy2;
       runSVFit(measuredTauLeptons, metcorrUncDown_ex, metcorrUncDown_ey, covMET, 0, svFitMass_UncMet_Down, svFitPt_UncMet_Down, svFitEta_UncMet_Down, svFitPhi_UncMet_Down, svFitMET_UncMet_Down, svFitTransverseMass_UncMet_Down, tau1_UncMet_Down, tau2_UncMet_Down);
       std::cout << "MET Clustered Energy Up     ---  ";
       metcorrClusteredUp_ex = metcorrClusteredUp_ex + dx1 + dx2;
       metcorrClusteredUp_ey = metcorrClusteredUp_ey + dy1 + dy2;
       runSVFit(measuredTauLeptons, metcorrClusteredUp_ex, metcorrClusteredUp_ey, covMET, 0, svFitMass_ClusteredMet_Up, svFitPt_ClusteredMet_Up, svFitEta_ClusteredMet_Up, svFitPhi_ClusteredMet_Up, svFitMET_ClusteredMet_Up, svFitTransverseMass_ClusteredMet_Up, tau1_ClusteredMet_Up, tau2_ClusteredMet_Up);
       std::cout << "MET Clustered Energy Down   ---  ";
       metcorrClusteredDown_ex = metcorrClusteredDown_ex + dx1 + dx2;
       metcorrClusteredDown_ey = metcorrClusteredDown_ey + dy1 + dy2;
       runSVFit(measuredTauLeptons, metcorrClusteredDown_ex, metcorrClusteredDown_ey, covMET, 0, svFitMass_ClusteredMet_Down, svFitPt_ClusteredMet_Down, svFitEta_ClusteredMet_Down, svFitPhi_ClusteredMet_Down, svFitMET_ClusteredMet_Down, svFitTransverseMass_ClusteredMet_Down, tau1_ClusteredMet_Down, tau2_ClusteredMet_Down);
       std::cout<< "Shifted MET Summary:\nmetcorr_ex " << metcorr_ex << "\n --- metcorrUncUp_ex " << metcorrUncUp_ex << " metcorrUncDown_ex " << metcorrUncDown_ex
            << " metcorrClusteredUp_ex " << metcorrClusteredUp_ex << " metcorrClusteredDown_ex " << metcorrClusteredDown_ex << std::endl;
       std::cout<< "metcorr_ey " << metcorr_ey << "\n --- metcorrUncUp_ey " << metcorrUncUp_ey << " metcorrUncDown_ey " << metcorrUncDown_ey
            << " metcorrClusteredUp_ey " << metcorrClusteredUp_ey << " metcorrClusteredDown_ey " << metcorrClusteredDown_ey << std::endl;
       
       
       // Corrected MET values for saving
       metcor = TMath::Sqrt( metcorr_ex*metcorr_ex + metcorr_ey*metcorr_ey);
       metcorphi = TMath::ATan2( metcorr_ey, metcorr_ex );
       metcorClusteredDown = TMath::Sqrt( metcorrClusteredDown_ex*metcorrClusteredDown_ex + metcorrClusteredDown_ey*metcorrClusteredDown_ey);
       metcorphiClusteredDown = TMath::ATan2( metcorrClusteredDown_ey, metcorrClusteredDown_ex );
       metcorClusteredUp = TMath::Sqrt( metcorrClusteredUp_ex*metcorrClusteredUp_ex + metcorrClusteredUp_ey*metcorrClusteredUp_ey);
       metcorphiClusteredUp = TMath::ATan2( metcorrClusteredUp_ey, metcorrClusteredUp_ex );
       metcorUncDown = TMath::Sqrt( metcorrUncDown_ex*metcorrUncDown_ex + metcorrUncDown_ey*metcorrUncDown_ey);
       metcorphiUncDown = TMath::ATan2( metcorrUncDown_ey, metcorrUncDown_ex );
       metcorUncUp = TMath::Sqrt( metcorrUncUp_ex*metcorrUncUp_ex + metcorrUncUp_ey*metcorrUncUp_ey);
       metcorphiUncUp = TMath::ATan2( metcorrUncUp_ey, metcorrUncUp_ex );
       
       if(doES) {
         
         //***************************************************************************
         //********************* Two taus shifted up *********************************
         //***************************************************************************
         
         if (gen_match_2==5 or gen_match_1==5){
           std::cout << "Two Up    ---  ";
           float ES_Up_scale1 = 1.0;
           float ES_Up_scale2 = 1.0;
           if(gen_match_1==5) ES_Up_scale1 = tesUp;
           if(gen_match_2==5) ES_Up_scale2 = tesUp;
           std::cout << "TES values: gen1: " << gen_match_1 << "   dm_1: " << decayMode;
           std::cout << "   tes1: " << ES_Up_scale1;
           std::cout << "   gen2: " << gen_match_2 << "   dm_2: " << decayMode2;
           std::cout << "   tes2: " << ES_Up_scale2 << std::endl;
           double pt1_Up, pt2_Up;
           pt1_Up = pt1 * ES_Up_scale1;
           pt2_Up = pt2 * ES_Up_scale2;
           double metcorr_ex_Up, metcorr_ey_Up;
           double dx1_Up, dy1_Up, dx2_Up, dy2_Up;
           dx1_Up = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dy1_Up = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dx2_Up = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           dy2_Up = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           metcorr_ex_Up = metcorr_ex + dx1_Up + dx2_Up;
           metcorr_ey_Up = metcorr_ey + dy1_Up + dy2_Up;
           
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
           
           // Add Tau of higest Pt first
           if (pt1 > pt2) {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
         
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
         
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_Up, svFitPt_Up, svFitEta_Up, svFitPhi_Up, svFitMET_Up, svFitTransverseMass_Up, tau1_Up, tau2_Up);
         }
         else {
           svFitMass_Up=svFitMass;
           svFitPt_Up=svFitPt;
           svFitEta_Up=svFitEta;
           svFitPhi_Up=svFitPhi;
           svFitMET_Up=svFitMET;
           svFitTransverseMass_Up=svFitTransverseMass;
           tau1_Up = tau1;
           tau2_Up = tau2;
         }
         
         //***************************************************************************
         //********************** Tau DM0 shifted up *********************************
         //***************************************************************************
         
         if ((gen_match_2==5 && decayMode2==0) or (gen_match_1==5 && decayMode==0)){
           std::cout << "DM0 Up    ---  ";
           float ES_Up_scale1 = 1.0;
           float ES_Up_scale2 = 1.0;
           if(gen_match_1==5 && decayMode==0) ES_Up_scale1 = tesUp;
           if(gen_match_2==5 && decayMode2==0) ES_Up_scale2 = tesUp;
           double pt1_Up, pt2_Up;
           pt1_Up = pt1 * ES_Up_scale1;
           pt2_Up = pt2 * ES_Up_scale2;
           double metcorr_ex_Up, metcorr_ey_Up;
           double dx1_Up, dy1_Up, dx2_Up, dy2_Up;
           dx1_Up = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dy1_Up = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dx2_Up = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           dy2_Up = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           metcorr_ex_Up = metcorr_ex + dx1_Up + dx2_Up;
           metcorr_ey_Up = metcorr_ey + dy1_Up + dy2_Up;
           
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
           // Add Tau of higest Pt first
           if (pt1 > pt2) {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM0_Up, svFitPt_DM0_Up, svFitEta_DM0_Up, svFitPhi_DM0_Up, svFitMET_DM0_Up, svFitTransverseMass_DM0_Up, tau1_DM0_Up, tau2_DM0_Down);
         }
         else {
           svFitMass_DM0_Up=svFitMass;
           svFitPt_DM0_Up=svFitPt;
           svFitEta_DM0_Up=svFitEta;
           svFitPhi_DM0_Up=svFitPhi;
           svFitMET_DM0_Up=svFitMET;
           svFitTransverseMass_DM0_Up=svFitTransverseMass;
           tau1_DM0_Up = tau1;
           tau2_DM0_Down = tau2;
         }
         
         //***************************************************************************
         //********************** Tau DM1 shifted up *********************************
         //***************************************************************************
         
         if ((decayMode==1 && gen_match_1==5) or (decayMode2==1 && gen_match_2==5)){
           std::cout << "DM1 Up    ---  ";
           float ES_Up_scale1 = 1.0;
           float ES_Up_scale2 = 1.0;
           if (decayMode==1 && gen_match_1==5) ES_Up_scale1 = tesUp;
           if (decayMode2==1 && gen_match_2==5) ES_Up_scale2 = tesUp;
           double pt1_Up, pt2_Up;
           pt1_Up = pt1 * ES_Up_scale1;
           pt2_Up = pt2 * ES_Up_scale2;
           double metcorr_ex_Up, metcorr_ey_Up;
           double dx1_Up, dy1_Up, dx2_Up, dy2_Up;
           dx1_Up = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dy1_Up = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dx2_Up = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           dy2_Up = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           metcorr_ex_Up = metcorr_ex + dx1_Up + dx2_Up;
           metcorr_ey_Up = metcorr_ey + dy1_Up + dy2_Up;
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
           // Add Tau of higest Pt first
           if (pt1 > pt2) {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM1_Up, svFitPt_DM1_Up, svFitEta_DM1_Up, svFitPhi_DM1_Up, svFitMET_DM1_Up, svFitTransverseMass_DM1_Up, tau1_DM1_Up, tau2_DM1_Up);
         }
         else {
           svFitMass_DM1_Up=svFitMass;
           svFitPt_DM1_Up=svFitPt;
           svFitEta_DM1_Up=svFitEta;
           svFitPhi_DM1_Up=svFitPhi;
           svFitMET_DM1_Up=svFitMET;
           svFitTransverseMass_DM1_Up=svFitTransverseMass;
           tau1_DM1_Up = tau1;
           tau2_DM1_Up = tau2;
         }
         
         //***************************************************************************
         //********************* Tau DM10 shifted up *********************************
         //***************************************************************************
         
         if ((decayMode2==10 && gen_match_2==5) or (decayMode==10 && gen_match_1==5)){
           std::cout << "DM10 Up    ---  ";
           float ES_Up_scale1 = 1.0;
           float ES_Up_scale2 = 1.0;
           if(decayMode==10 && gen_match_1==5) ES_Up_scale1 = tesUp;
           if(decayMode2==10 && gen_match_2==5) ES_Up_scale2 = tesUp;
           double pt1_Up, pt2_Up;
           pt1_Up = pt1 * ES_Up_scale1;
           pt2_Up = pt2 * ES_Up_scale2;
           double metcorr_ex_Up, metcorr_ey_Up;
           double dx1_Up, dy1_Up, dx2_Up, dy2_Up;
           dx1_Up = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dy1_Up = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Up_scale1 ) - 1.);
           dx2_Up = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           dy2_Up = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Up_scale2 ) - 1.);
           metcorr_ex_Up = metcorr_ex + dx1_Up + dx2_Up;
           metcorr_ey_Up = metcorr_ey + dy1_Up + dy2_Up;
           
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
           
           // Add Tau of higest Pt first
           if (pt1 > pt2) {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Up, eta2, phi2,  mass2, decayMode2));
         measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM10_Up, svFitPt_DM10_Up, svFitEta_DM10_Up, svFitPhi_DM10_Up, svFitMET_DM10_Up, svFitTransverseMass_DM10_Up, tau1_DM10_Up, tau2_DM10_Up);
         }
         else {
           svFitMass_DM10_Up=svFitMass;
           svFitPt_DM10_Up=svFitPt;
           svFitEta_DM10_Up=svFitEta;
           svFitPhi_DM10_Up=svFitPhi;
           svFitMET_DM10_Up=svFitMET;
           svFitTransverseMass_DM10_Up=svFitTransverseMass;
           tau1_DM10_Up = tau1;
           tau2_DM10_Up = tau2;
         }
         
         //*****************************************************
         //************* Two taus shifted down *****************
         //*****************************************************
         
         if (gen_match_1==5 or gen_match_2==5){
           std::cout << "Two Down  ---  ";
           float ES_Down_scale1 = 1.0;
           float ES_Down_scale2 = 1.0;
           if (gen_match_1==5) ES_Down_scale1 = tesDown;
           if (gen_match_2==5) ES_Down_scale2 = tesDown;
           double pt1_Down, pt2_Down;
           pt1_Down = pt1 * ES_Down_scale1;
           pt2_Down = pt2 * ES_Down_scale2;
           double metcorr_ex_Down, metcorr_ey_Down;
           double dx1_Down, dy1_Down, dx2_Down, dy2_Down;
           dx1_Down = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dy1_Down = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dx2_Down = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           dy2_Down = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           metcorr_ex_Down = metcorr_ex + dx1_Down + dx2_Down;
           metcorr_ey_Down = metcorr_ey + dy1_Down + dy2_Down;
           
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
           
           if (pt1 > pt2) {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Down, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Down, eta2, phi2,  mass2, decayMode2));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_Down, svFitPt_Down, svFitEta_Down, svFitPhi_Down, svFitMET_Down, svFitTransverseMass_Down, tau1_Down, tau2_Down);
         }
         else {
           svFitMass_Down=svFitMass;
           svFitPt_Down=svFitPt;
           svFitEta_Down=svFitEta;
           svFitPhi_Down=svFitPhi;
           svFitMET_Down=svFitMET;
           svFitTransverseMass_Down=svFitTransverseMass;
           tau1_Down = tau1;
           tau2_Down = tau2;
         }
         
         //*****************************************************
         //************* Tau DM0 shifted down  *****************
         //*****************************************************
         
         if ((decayMode==0 && gen_match_1==5) or (decayMode2==0 && gen_match_2==5)){
           std::cout << "DM0 Down  ---  ";
           float ES_Down_scale1 = 1.0;
           float ES_Down_scale2 = 1.0;
           if (decayMode==0 && gen_match_1==5) ES_Down_scale1 = tesDown;
           if (decayMode2==0 && gen_match_2==5) ES_Down_scale2 = tesDown;
           double pt1_Down, pt2_Down;
           pt1_Down = pt1 * ES_Down_scale1;
           pt2_Down = pt2 * ES_Down_scale2;
           double metcorr_ex_Down, metcorr_ey_Down;
           double dx1_Down, dy1_Down, dx2_Down, dy2_Down;
           dx1_Down = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dy1_Down = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dx2_Down = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           dy2_Down = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           metcorr_ex_Down = metcorr_ex + dx1_Down + dx2_Down;
           metcorr_ey_Down = metcorr_ey + dy1_Down + dy2_Down;
           
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
           if (pt1 > pt2) {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Down, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Down, eta2, phi2,  mass2, decayMode2));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM0_Down, svFitPt_DM0_Down, svFitEta_DM0_Down, svFitPhi_DM0_Down, svFitMET_DM0_Down, svFitTransverseMass_DM0_Down, tau1_DM0_Down, tau2_DM0_Down);
         }
         else {
           svFitMass_DM0_Down=svFitMass;
           svFitPt_DM0_Down=svFitPt;
           svFitEta_DM0_Down=svFitEta;
           svFitPhi_DM0_Down=svFitPhi;
           svFitMET_DM0_Down=svFitMET;
           svFitTransverseMass_DM0_Down=svFitTransverseMass;
           tau1_DM0_Down = tau1;
           tau2_DM0_Down = tau2;
         }
         
         //*****************************************************
         //************** Tau DM1 shifted down *****************
         //*****************************************************
         
         if ((decayMode==1 && gen_match_1==5) or (decayMode2==1 && gen_match_2==5)){
           std::cout << "DM1 Down  ---  ";
           float ES_Down_scale1 = 1.0;
           float ES_Down_scale2 = 1.0;
           if (decayMode==1 && gen_match_1==5) ES_Down_scale1 = tesDown;
           if (decayMode2==1 && gen_match_2==5) ES_Down_scale2 = tesDown;
           double pt1_Down, pt2_Down;
           pt1_Down = pt1 * ES_Down_scale1;
           pt2_Down = pt2 * ES_Down_scale2;
           double metcorr_ex_Down, metcorr_ey_Down;
           double dx1_Down, dy1_Down, dx2_Down, dy2_Down;
           dx1_Down = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dy1_Down = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dx2_Down = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           dy2_Down = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           metcorr_ex_Down = metcorr_ex + dx1_Down + dx2_Down;
           metcorr_ey_Down = metcorr_ey + dy1_Down + dy2_Down;
           
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
           
           if (pt1 > pt2) {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Down, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Down, eta2, phi2,  mass2, decayMode2));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM1_Down, svFitPt_DM1_Down, svFitEta_DM1_Down, svFitPhi_DM1_Down, svFitMET_DM1_Down, svFitTransverseMass_DM1_Down, tau1_DM1_Down, tau2_DM1_Down);
         }
         else {
           svFitMass_DM1_Down=svFitMass;
           svFitPt_DM1_Down=svFitPt;
           svFitEta_DM1_Down=svFitEta;
           svFitPhi_DM1_Down=svFitPhi;
           svFitMET_DM1_Down=svFitMET;
           svFitTransverseMass_DM1_Down=svFitTransverseMass;
           tau1_DM1_Down =  tau1;
           tau2_DM1_Down =  tau2;
         }
         
         //*****************************************************
         //************* Tau DM10 shifted down *****************
         //*****************************************************
         
         if ((decayMode==10 && gen_match_1==5) or (decayMode2==10 && gen_match_2==5)){
           std::cout << "DM10 Down  ---  ";
           float ES_Down_scale1 = 1.0;
           float ES_Down_scale2 = 1.0;
           if (decayMode==10 && gen_match_1==5) ES_Down_scale1 = tesDown;
           if (decayMode2==10 && gen_match_2==5) ES_Down_scale2 = tesDown;
           double pt1_Down, pt2_Down;
           pt1_Down = pt1 * ES_Down_scale1;
           pt2_Down = pt2 * ES_Down_scale2;
           double metcorr_ex_Down, metcorr_ey_Down;
           double dx1_Down, dy1_Down, dx2_Down, dy2_Down;
           dx1_Down = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dy1_Down = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_Down_scale1 ) - 1.);
           dx2_Down = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           dy2_Down = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_Down_scale2 ) - 1.);
           metcorr_ex_Down = metcorr_ex + dx1_Down + dx2_Down;
           metcorr_ey_Down = metcorr_ey + dy1_Down + dy2_Down;
           
           std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
           if (pt1 > pt2) {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2,  pt2_Down, eta2, phi2,  mass2, decayMode2));
           }
           else {
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2, pt2_Down, eta2, phi2, mass2, decayMode2));
         measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1, phi1, mass1, decayMode));
           }
           
           runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM10_Down, svFitPt_DM10_Down, svFitEta_DM10_Down, svFitPhi_DM10_Down, svFitMET_DM10_Down, svFitTransverseMass_DM10_Down, tau1_DM10_Down, tau2_DM10_Down);
         }
         else {
           svFitMass_DM10_Down=svFitMass;
           svFitPt_DM10_Down=svFitPt;
           svFitEta_DM10_Down=svFitEta;
           svFitPhi_DM10_Down=svFitPhi;
           svFitMET_DM10_Down=svFitMET;
           svFitTransverseMass_DM10_Down=svFitTransverseMass;
           tau1_DM10_Down = tau1;
           tau2_DM10_Down = tau2;
         }
       }// Do ES (TT)
       
     } // Double Hadronic (TT)
     
     else {
       svFitMass = -100;
       svFitPt = -100;
       svFitEta = -100;
       svFitPhi = -100;
       svFitMET = -100;
       svFitTransverseMass = -100;
       tau1.SetPtEtaPhiM(0,0,0,0);
       tau2.SetPtEtaPhiM(0,0,0,0); 
     }
     // fill the tau 4-vector parameters for branch-filling
     tau1_pt  = tau1.Pt();
     tau1_eta = tau1.Eta();
     tau1_phi = tau1.Phi();
     tau1_m   = tau1.M();
     tau2_pt  = tau2.Pt();
     tau2_eta = tau2.Eta();
     tau2_phi = tau2.Phi();
     tau2_m   = tau2.M();
     // up
     tau1_pt_Up  = tau1_Up.Pt();
     tau1_eta_Up = tau1_Up.Eta();
     tau1_phi_Up = tau1_Up.Phi();
     tau1_m_Up   = tau1_Up.M();
     tau2_pt_Up  = tau2_Up.Pt();
     tau2_eta_Up = tau2_Up.Eta();
     tau2_phi_Up = tau2_Up.Phi();
     tau2_m_Up   = tau2_Up.M();
     // down
     tau1_pt_Down  = tau1_Down.Pt();
     tau1_eta_Down = tau1_Down.Eta();
     tau1_phi_Down = tau1_Down.Phi();
     tau1_m_Down   = tau1_Down.M();
     tau2_pt_Down  = tau2_Down.Pt();
     tau2_eta_Down = tau2_Down.Eta();
     tau2_phi_Down = tau2_Down.Phi();
     tau2_m_Down   = tau2_Down.M();
     // DM0 up
     tau1_pt_DM0_Up  = tau1_DM0_Up.Pt();
     tau1_eta_DM0_Up = tau1_DM0_Up.Eta();
     tau1_phi_DM0_Up = tau1_DM0_Up.Phi();
     tau1_m_DM0_Up   = tau1_DM0_Up.M();
     tau2_pt_DM0_Up  = tau2_DM0_Up.Pt();
     tau2_eta_DM0_Up = tau2_DM0_Up.Eta();
     tau2_phi_DM0_Up = tau2_DM0_Up.Phi();
     tau2_m_DM0_Up   = tau2_DM0_Up.M();
     // down
     tau1_pt_DM0_Down  = tau1_DM0_Down.Pt();
     tau1_eta_DM0_Down = tau1_DM0_Down.Eta();
     tau1_phi_DM0_Down = tau1_DM0_Down.Phi();
     tau1_m_DM0_Down   = tau1_DM0_Down.M();
     tau2_pt_DM0_Down  = tau2_DM0_Down.Pt();
     tau2_eta_DM0_Down = tau2_DM0_Down.Eta();
     tau2_phi_DM0_Down = tau2_DM0_Down.Phi();
     tau2_m_DM0_Down   = tau2_DM0_Down.M();
     // DM1 up
     tau1_pt_DM1_Up  = tau1_DM1_Up.Pt();
     tau1_eta_DM1_Up = tau1_DM1_Up.Eta();
     tau1_phi_DM1_Up = tau1_DM1_Up.Phi();
     tau1_m_DM1_Up   = tau1_DM1_Up.M();
     tau2_pt_DM1_Up  = tau2_DM1_Up.Pt();
     tau2_eta_DM1_Up = tau2_DM1_Up.Eta();
     tau2_phi_DM1_Up = tau2_DM1_Up.Phi();
     tau2_m_DM1_Up   = tau2_DM1_Up.M();
     // down
     tau1_pt_DM1_Down  = tau1_DM1_Down.Pt();
     tau1_eta_DM1_Down = tau1_DM1_Down.Eta();
     tau1_phi_DM1_Down = tau1_DM1_Down.Phi();
     tau1_m_DM1_Down   = tau1_DM1_Down.M();
     tau2_pt_DM1_Down  = tau2_DM1_Down.Pt();
     tau2_eta_DM1_Down = tau2_DM1_Down.Eta();
     tau2_phi_DM1_Down = tau2_DM1_Down.Phi();
     tau2_m_DM1_Down   = tau2_DM1_Down.M();
     // DM10 up
     tau1_pt_DM10_Up  = tau1_DM10_Up.Pt();
     tau1_eta_DM10_Up = tau1_DM10_Up.Eta();
     tau1_phi_DM10_Up = tau1_DM10_Up.Phi();
     tau1_m_DM10_Up   = tau1_DM10_Up.M();
     tau2_pt_DM10_Up  = tau2_DM10_Up.Pt();
     tau2_eta_DM10_Up = tau2_DM10_Up.Eta();
     tau2_phi_DM10_Up = tau2_DM10_Up.Phi();
     tau2_m_DM10_Up   = tau2_DM10_Up.M();
     // down
     tau1_pt_DM10_Down  = tau1_DM10_Down.Pt();
     tau1_eta_DM10_Down = tau1_DM10_Down.Eta();
     tau1_phi_DM10_Down = tau1_DM10_Down.Phi();
     tau1_m_DM10_Down   = tau1_DM10_Down.M();
     tau2_pt_DM10_Down  = tau2_DM10_Down.Pt();
     tau2_eta_DM10_Down = tau2_DM10_Down.Eta();
     tau2_phi_DM10_Down = tau2_DM10_Down.Phi();
     tau2_m_DM10_Down   = tau2_DM10_Down.M();
     // UncMet up
     tau1_pt_UncMet_Up  = tau1_UncMet_Up.Pt();
     tau1_eta_UncMet_Up = tau1_UncMet_Up.Eta();
     tau1_phi_UncMet_Up = tau1_UncMet_Up.Phi();
     tau1_m_UncMet_Up   = tau1_UncMet_Up.M();
     tau2_pt_UncMet_Up  = tau2_UncMet_Up.Pt();
     tau2_eta_UncMet_Up = tau2_UncMet_Up.Eta();
     tau2_phi_UncMet_Up = tau2_UncMet_Up.Phi();
     tau2_m_UncMet_Up   = tau2_UncMet_Up.M();
     // down
     tau1_pt_UncMet_Down  = tau1_UncMet_Down.Pt();
     tau1_eta_UncMet_Down = tau1_UncMet_Down.Eta();
     tau1_phi_UncMet_Down = tau1_UncMet_Down.Phi();
     tau1_m_UncMet_Down   = tau1_UncMet_Down.M();
     tau2_pt_UncMet_Down  = tau2_UncMet_Down.Pt();
     tau2_eta_UncMet_Down = tau2_UncMet_Down.Eta();
     tau2_phi_UncMet_Down = tau2_UncMet_Down.Phi();
     tau2_m_UncMet_Down   = tau2_UncMet_Down.M();
     // UnclusteredMet up
     tau1_pt_ClusteredMet_Up  = tau1_ClusteredMet_Up.Pt();
     tau1_eta_ClusteredMet_Up = tau1_ClusteredMet_Up.Eta();
     tau1_phi_ClusteredMet_Up = tau1_ClusteredMet_Up.Phi();
     tau1_m_ClusteredMet_Up   = tau1_ClusteredMet_Up.M();
     tau2_pt_ClusteredMet_Up  = tau2_ClusteredMet_Up.Pt();
     tau2_eta_ClusteredMet_Up = tau2_ClusteredMet_Up.Eta();
     tau2_phi_ClusteredMet_Up = tau2_ClusteredMet_Up.Phi();
     tau2_m_ClusteredMet_Up   = tau2_ClusteredMet_Up.M();
     // down
     tau1_pt_ClusteredMet_Down  = tau1_ClusteredMet_Down.Pt();
     tau1_eta_ClusteredMet_Down = tau1_ClusteredMet_Down.Eta();
     tau1_phi_ClusteredMet_Down = tau1_ClusteredMet_Down.Phi();
     tau1_m_ClusteredMet_Down   = tau1_ClusteredMet_Down.M();
     tau2_pt_ClusteredMet_Down  = tau2_ClusteredMet_Down.Pt();
     tau2_eta_ClusteredMet_Down = tau2_ClusteredMet_Down.Eta();
     tau2_phi_ClusteredMet_Down = tau2_ClusteredMet_Down.Phi();
     tau2_m_ClusteredMet_Down   = tau2_ClusteredMet_Down.M();
     
     
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
     
     newBranch11->Fill();
     newBranch12->Fill();
     newBranch13->Fill();
     newBranch14->Fill();
     newBranch15->Fill();
     newBranch16->Fill();
     newBranch17->Fill();
     newBranch18->Fill();
     newBranch19->Fill();
     newBranch20->Fill();
     newBranch21->Fill();
     newBranch22->Fill();
     
     newBranch23->Fill();
     newBranch24->Fill();
     newBranch25->Fill();
     newBranch26->Fill();
     newBranch27->Fill();
     newBranch28->Fill();
     newBranch29->Fill();
     newBranch30->Fill();
     newBranch31->Fill();
     newBranch32->Fill();
     newBranch33->Fill();
     newBranch34->Fill();
     
     newBranch35->Fill();
     newBranch36->Fill();
     newBranch37->Fill();
     newBranch38->Fill();
     newBranch39->Fill();
     newBranch40->Fill();
     newBranch41->Fill();
     newBranch42->Fill();
     newBranch43->Fill();
     newBranch44->Fill();
     newBranch45->Fill();
     newBranch46->Fill();
     
     newBranch47->Fill();
     newBranch48->Fill();
     newBranch49->Fill();
     newBranch50->Fill();
     newBranch51->Fill();
     newBranch52->Fill();
     newBranch53->Fill();
     newBranch54->Fill();
     newBranch55->Fill();
     newBranch56->Fill();
     newBranch57->Fill();
     newBranch58->Fill();
     
     newBranch59->Fill();
     newBranch60->Fill();
     newBranch61->Fill();
     newBranch62->Fill();
     newBranch63->Fill();
     newBranch64->Fill();
     newBranch65->Fill();
     newBranch66->Fill();
     newBranch67->Fill();
     newBranch68->Fill();
     newBranch69->Fill();
     newBranch70->Fill();
     
     newBranch71->Fill();
     newBranch72->Fill();
     newBranch73->Fill();
     newBranch74->Fill();
     newBranch75->Fill();
     newBranch76->Fill();
     newBranch77->Fill();
     newBranch78->Fill();
     newBranch79->Fill();
     newBranch80->Fill();
     newBranch81->Fill();
     newBranch82->Fill();
     
     newBranch83->Fill();
     newBranch84->Fill();
     newBranch85->Fill();
     newBranch86->Fill();
     newBranch87->Fill();
     newBranch89->Fill();
     newBranch90->Fill();
     newBranch91->Fill();
     
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
void CopyFile(const char *fname, optutl::CommandLineParser parser) {
  //Copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return;
  }
  target->cd();
  CopyDir(f,parser);
  delete f;
  target->cd();
}
void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) 
{
  //prepare files to be copied
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("inputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();

}

