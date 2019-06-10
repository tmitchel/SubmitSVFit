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

      float metcorrJetEC2Up_ex = -10;
      float metcorrJetEC2Up_ey = -10;
      float metcorrJetEC2Down_ex = -10;
      float metcorrJetEC2Down_ey = -10;
      float metcorrJetEta0to3Up_ex = -10;
      float metcorrJetEta0to3Up_ey = -10;
      float metcorrJetEta0to3Down_ex = -10;
      float metcorrJetEta0to3Down_ey = -10;
      float metcorrJetEta0to5Up_ex = -10;
      float metcorrJetEta0to5Up_ey = -10;
      float metcorrJetEta0to5Down_ex = -10;
      float metcorrJetEta0to5Down_ey = -10;
      float metcorrJetEta3to5Up_ex = -10;
      float metcorrJetEta3to5Up_ey = -10;
      float metcorrJetEta3to5Down_ex = -10;
      float metcorrJetEta3to5Down_ey = -10;
      float metcorrJetRelativeBalUp_ex = -10;
      float metcorrJetRelativeBalUp_ey = -10;
      float metcorrJetRelativeBalDown_ex = -10;
      float metcorrJetRelativeBalDown_ey = -10;
      float metcorrJetRelativeSampleUp_ex = -10;
      float metcorrJetRelativeSampleUp_ey = -10;
      float metcorrJetRelativeSampleDown_ex = -10;
      float metcorrJetRelativeSampleDown_ey = -10;


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
      float metcorJetEC2Down = -10;
      float metcorphiJetEC2Down = -10;
      float metcorJetEC2Up = -10;
      float metcorphiJetEC2Up = -10;
      float metcorJetEta0to3Down = -10;
      float metcorphiJetEta0to3Down = -10;
      float metcorJetEta0to3Up = -10;
      float metcorphiJetEta0to3Up = -10;
      float metcorJetEta0to5Down = -10;
      float metcorphiJetEta0to5Down = -10;
      float metcorJetEta0to5Up = -10;
      float metcorphiJetEta0to5Up = -10;
      float metcorJetEta3to5Down = -10;
      float metcorphiJetEta3to5Down = -10;
      float metcorJetEta3to5Up = -10;
      float metcorphiJetEta3to5Up = -10;
      float metcorJetRelativeBalDown = -10;
      float metcorphiJetRelativeBalDown = -10;
      float metcorJetRelativeBalUp = -10;
      float metcorphiJetRelativeBalUp = -10;
      float metcorJetRelativeSampleDown = -10;
      float metcorphiJetRelativeSampleDown = -10;
      float metcorJetRelativeSampleUp = -10;
      float metcorphiJetRelativeSampleUp = -10;

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
      
      float svFitMass_JetEC2_Up = -10;
      float svFitPt_JetEC2_Up = -10;
      float svFitEta_JetEC2_Up = -10;
      float svFitPhi_JetEC2_Up = -10;
      float svFitMET_JetEC2_Up = -10;
      float svFitTransverseMass_JetEC2_Up = -10;

      float svFitMass_JetEC2_Down = -10;
      float svFitPt_JetEC2_Down = -10;
      float svFitEta_JetEC2_Down = -10;
      float svFitPhi_JetEC2_Down = -10;
      float svFitMET_JetEC2_Down = -10;
      float svFitTransverseMass_JetEC2_Down = -10;

      float svFitMass_JetEta0to3_Up = -10;
      float svFitPt_JetEta0to3_Up = -10;
      float svFitEta_JetEta0to3_Up = -10;
      float svFitPhi_JetEta0to3_Up = -10;
      float svFitMET_JetEta0to3_Up = -10;
      float svFitTransverseMass_JetEta0to3_Up = -10;

      float svFitMass_JetEta0to3_Down = -10;
      float svFitPt_JetEta0to3_Down = -10;
      float svFitEta_JetEta0to3_Down = -10;
      float svFitPhi_JetEta0to3_Down = -10;
      float svFitMET_JetEta0to3_Down = -10;
      float svFitTransverseMass_JetEta0to3_Down = -10;

      float svFitMass_JetEta0to5_Up = -10;
      float svFitPt_JetEta0to5_Up = -10;
      float svFitEta_JetEta0to5_Up = -10;
      float svFitPhi_JetEta0to5_Up = -10;
      float svFitMET_JetEta0to5_Up = -10;
      float svFitTransverseMass_JetEta0to5_Up = -10;

      float svFitMass_JetEta0to5_Down = -10;
      float svFitPt_JetEta0to5_Down = -10;
      float svFitEta_JetEta0to5_Down = -10;
      float svFitPhi_JetEta0to5_Down = -10;
      float svFitMET_JetEta0to5_Down = -10;
      float svFitTransverseMass_JetEta0to5_Down = -10;

      float svFitMass_JetEta3to5_Up = -10;
      float svFitPt_JetEta3to5_Up = -10;
      float svFitEta_JetEta3to5_Up = -10;
      float svFitPhi_JetEta3to5_Up = -10;
      float svFitMET_JetEta3to5_Up = -10;
      float svFitTransverseMass_JetEta3to5_Up = -10;

      float svFitMass_JetEta3to5_Down = -10;
      float svFitPt_JetEta3to5_Down = -10;
      float svFitEta_JetEta3to5_Down = -10;
      float svFitPhi_JetEta3to5_Down = -10;
      float svFitMET_JetEta3to5_Down = -10;
      float svFitTransverseMass_JetEta3to5_Down = -10;

      float svFitMass_JetRelativeBal_Up = -10;
      float svFitPt_JetRelativeBal_Up = -10;
      float svFitEta_JetRelativeBal_Up = -10;
      float svFitPhi_JetRelativeBal_Up = -10;
      float svFitMET_JetRelativeBal_Up = -10;
      float svFitTransverseMass_JetRelativeBal_Up = -10;

      float svFitMass_JetRelativeBal_Down = -10;
      float svFitPt_JetRelativeBal_Down = -10;
      float svFitEta_JetRelativeBal_Down = -10;
      float svFitPhi_JetRelativeBal_Down = -10;
      float svFitMET_JetRelativeBal_Down = -10;
      float svFitTransverseMass_JetRelativeBal_Down = -10;

      float svFitMass_JetRelativeSample_Up = -10;
      float svFitPt_JetRelativeSample_Up = -10;
      float svFitEta_JetRelativeSample_Up = -10;
      float svFitPhi_JetRelativeSample_Up = -10;
      float svFitMET_JetRelativeSample_Up = -10;
      float svFitTransverseMass_JetRelativeSample_Up = -10;

      float svFitMass_JetRelativeSample_Down = -10;
      float svFitPt_JetRelativeSample_Down = -10;
      float svFitEta_JetRelativeSample_Down = -10;
      float svFitPhi_JetRelativeSample_Down = -10;
      float svFitMET_JetRelativeSample_Down = -10;
      float svFitTransverseMass_JetRelativeSample_Down = -10;


	

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
      // up JES JetEC2 ========================= 
      float tau1_pt_JetEC2_Up = -10;
      float tau1_eta_JetEC2_Up = -10;
      float tau1_phi_JetEC2_Up = -10;
      float tau1_m_JetEC2_Up = -10;
      float tau2_pt_JetEC2_Up = -10;
      float tau2_eta_JetEC2_Up = -10;
      float tau2_phi_JetEC2_Up = -10;
      float tau2_m_JetEC2_Up = -10;
      // down
      float tau1_pt_JetEC2_Down = -10;
      float tau1_eta_JetEC2_Down = -10;
      float tau1_phi_JetEC2_Down = -10;
      float tau1_m_JetEC2_Down = -10;
      float tau2_pt_JetEC2_Down = -10;
      float tau2_eta_JetEC2_Down = -10;
      float tau2_phi_JetEC2_Down = -10;
      float tau2_m_JetEC2_Down = -10;
      // up JES JetEta0to3 ========================= 
      float tau1_pt_JetEta0to3_Up = -10;
      float tau1_eta_JetEta0to3_Up = -10;
      float tau1_phi_JetEta0to3_Up = -10;
      float tau1_m_JetEta0to3_Up = -10;
      float tau2_pt_JetEta0to3_Up = -10;
      float tau2_eta_JetEta0to3_Up = -10;
      float tau2_phi_JetEta0to3_Up = -10;
      float tau2_m_JetEta0to3_Up = -10;
      // down
      float tau1_pt_JetEta0to3_Down = -10;
      float tau1_eta_JetEta0to3_Down = -10;
      float tau1_phi_JetEta0to3_Down = -10;
      float tau1_m_JetEta0to3_Down = -10;
      float tau2_pt_JetEta0to3_Down = -10;
      float tau2_eta_JetEta0to3_Down = -10;
      float tau2_phi_JetEta0to3_Down = -10;
      float tau2_m_JetEta0to3_Down = -10;
      // up JES JetEta0to5 ========================= 
      float tau1_pt_JetEta0to5_Up = -10;
      float tau1_eta_JetEta0to5_Up = -10;
      float tau1_phi_JetEta0to5_Up = -10;
      float tau1_m_JetEta0to5_Up = -10;
      float tau2_pt_JetEta0to5_Up = -10;
      float tau2_eta_JetEta0to5_Up = -10;
      float tau2_phi_JetEta0to5_Up = -10;
      float tau2_m_JetEta0to5_Up = -10;
      // down
      float tau1_pt_JetEta0to5_Down = -10;
      float tau1_eta_JetEta0to5_Down = -10;
      float tau1_phi_JetEta0to5_Down = -10;
      float tau1_m_JetEta0to5_Down = -10;
      float tau2_pt_JetEta0to5_Down = -10;
      float tau2_eta_JetEta0to5_Down = -10;
      float tau2_phi_JetEta0to5_Down = -10;
      float tau2_m_JetEta0to5_Down = -10;
      // up JES JetEta3to5 ========================= 
      float tau1_pt_JetEta3to5_Up = -10;
      float tau1_eta_JetEta3to5_Up = -10;
      float tau1_phi_JetEta3to5_Up = -10;
      float tau1_m_JetEta3to5_Up = -10;
      float tau2_pt_JetEta3to5_Up = -10;
      float tau2_eta_JetEta3to5_Up = -10;
      float tau2_phi_JetEta3to5_Up = -10;
      float tau2_m_JetEta3to5_Up = -10;
      // down
      float tau1_pt_JetEta3to5_Down = -10;
      float tau1_eta_JetEta3to5_Down = -10;
      float tau1_phi_JetEta3to5_Down = -10;
      float tau1_m_JetEta3to5_Down = -10;
      float tau2_pt_JetEta3to5_Down = -10;
      float tau2_eta_JetEta3to5_Down = -10;
      float tau2_phi_JetEta3to5_Down = -10;
      float tau2_m_JetEta3to5_Down = -10;
      // up JES JetRelativeBal ========================= 
      float tau1_pt_JetRelativeBal_Up = -10;
      float tau1_eta_JetRelativeBal_Up = -10;
      float tau1_phi_JetRelativeBal_Up = -10;
      float tau1_m_JetRelativeBal_Up = -10;
      float tau2_pt_JetRelativeBal_Up = -10;
      float tau2_eta_JetRelativeBal_Up = -10;
      float tau2_phi_JetRelativeBal_Up = -10;
      float tau2_m_JetRelativeBal_Up = -10;
      // down
      float tau1_pt_JetRelativeBal_Down = -10;
      float tau1_eta_JetRelativeBal_Down = -10;
      float tau1_phi_JetRelativeBal_Down = -10;
      float tau1_m_JetRelativeBal_Down = -10;
      float tau2_pt_JetRelativeBal_Down = -10;
      float tau2_eta_JetRelativeBal_Down = -10;
      float tau2_phi_JetRelativeBal_Down = -10;
      float tau2_m_JetRelativeBal_Down = -10;
      // up JES JetRelativeSample ========================= 
      float tau1_pt_JetRelativeSample_Up = -10;
      float tau1_eta_JetRelativeSample_Up = -10;
      float tau1_phi_JetRelativeSample_Up = -10;
      float tau1_m_JetRelativeSample_Up = -10;
      float tau2_pt_JetRelativeSample_Up = -10;
      float tau2_eta_JetRelativeSample_Up = -10;
      float tau2_phi_JetRelativeSample_Up = -10;
      float tau2_m_JetRelativeSample_Up = -10;
      // down
      float tau1_pt_JetRelativeSample_Down = -10;
      float tau1_eta_JetRelativeSample_Down = -10;
      float tau1_phi_JetRelativeSample_Down = -10;
      float tau1_m_JetRelativeSample_Down = -10;
      float tau2_pt_JetRelativeSample_Down = -10;
      float tau2_eta_JetRelativeSample_Down = -10;
      float tau2_phi_JetRelativeSample_Down = -10;
      float tau2_m_JetRelativeSample_Down = -10;

                               
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

      TBranch *newBranch83 = t->Branch("m_sv_JetEC2_Up", &svFitMass_JetEC2_Up, "m_sv_JetEC2_Up/F");
      TBranch *newBranch84 = t->Branch("pt_sv_JetEC2_Up", &svFitPt_JetEC2_Up, "pt_sv_JetEC2_Up/F");
      TBranch *newBranch85 = t->Branch("eta_sv_JetEC2_Up", &svFitEta_JetEC2_Up, "eta_sv_JetEC2_Up/F");
      TBranch *newBranch86 = t->Branch("phi_sv_JetEC2_Up", &svFitPhi_JetEC2_Up, "phi_sv_JetEC2_Up/F");
      TBranch *newBranch87 = t->Branch("met_sv_JetEC2_Up", &svFitMET_JetEC2_Up, "met_sv_JetEC2_Up/F");
      TBranch *newBranch88 = t->Branch("mt_sv_JetEC2_Up", &svFitTransverseMass_JetEC2_Up, "mt_sv_JetEC2_Up/F");

      TBranch *newBranch89 = t->Branch("m_sv_JetEC2_Down", &svFitMass_JetEC2_Down, "m_sv_JetEC2_Down/F");
      TBranch *newBranch90 = t->Branch("pt_sv_JetEC2_Down", &svFitPt_JetEC2_Down, "pt_sv_JetEC2_Down/F");
      TBranch *newBranch91 = t->Branch("eta_sv_JetEC2_Down", &svFitEta_JetEC2_Down, "eta_sv_JetEC2_Down/F");
      TBranch *newBranch92 = t->Branch("phi_sv_JetEC2_Down", &svFitPhi_JetEC2_Down, "phi_sv_JetEC2_Down/F");
      TBranch *newBranch93 = t->Branch("met_sv_JetEC2_Down", &svFitMET_JetEC2_Down, "met_sv_JetEC2_Down/F");
      TBranch *newBranch94 = t->Branch("mt_sv_JetEC2_Down", &svFitTransverseMass_JetEC2_Down, "mt_sv_JetEC2_Down/F");

      TBranch *newBranch95 = t->Branch("m_sv_JetEta0to3_Up", &svFitMass_JetEta0to3_Up, "m_sv_JetEta0to3_Up/F");
      TBranch *newBranch96 = t->Branch("pt_sv_JetEta0to3_Up", &svFitPt_JetEta0to3_Up, "pt_sv_JetEta0to3_Up/F");
      TBranch *newBranch97 = t->Branch("eta_sv_JetEta0to3_Up", &svFitEta_JetEta0to3_Up, "eta_sv_JetEta0to3_Up/F");
      TBranch *newBranch98 = t->Branch("phi_sv_JetEta0to3_Up", &svFitPhi_JetEta0to3_Up, "phi_sv_JetEta0to3_Up/F");
      TBranch *newBranch99 = t->Branch("met_sv_JetEta0to3_Up", &svFitMET_JetEta0to3_Up, "met_sv_JetEta0to3_Up/F");
      TBranch *newBranch100 = t->Branch("mt_sv_JetEta0to3_Up", &svFitTransverseMass_JetEta0to3_Up, "mt_sv_JetEta0to3_Up/F");

      TBranch *newBranch101 = t->Branch("m_sv_JetEta0to3_Down", &svFitMass_JetEta0to3_Down, "m_sv_JetEta0to3_Down/F");
      TBranch *newBranch102 = t->Branch("pt_sv_JetEta0to3_Down", &svFitPt_JetEta0to3_Down, "pt_sv_JetEta0to3_Down/F");
      TBranch *newBranch103 = t->Branch("eta_sv_JetEta0to3_Down", &svFitEta_JetEta0to3_Down, "eta_sv_JetEta0to3_Down/F");
      TBranch *newBranch104 = t->Branch("phi_sv_JetEta0to3_Down", &svFitPhi_JetEta0to3_Down, "phi_sv_JetEta0to3_Down/F");
      TBranch *newBranch105 = t->Branch("met_sv_JetEta0to3_Down", &svFitMET_JetEta0to3_Down, "met_sv_JetEta0to3_Down/F");
      TBranch *newBranch106 = t->Branch("mt_sv_JetEta0to3_Down", &svFitTransverseMass_JetEta0to3_Down, "mt_sv_JetEta0to3_Down/F");

      TBranch *newBranch107 = t->Branch("m_sv_JetEta0to5_Up", &svFitMass_JetEta0to5_Up, "m_sv_JetEta0to5_Up/F");
      TBranch *newBranch108 = t->Branch("pt_sv_JetEta0to5_Up", &svFitPt_JetEta0to5_Up, "pt_sv_JetEta0to5_Up/F");
      TBranch *newBranch109 = t->Branch("eta_sv_JetEta0to5_Up", &svFitEta_JetEta0to5_Up, "eta_sv_JetEta0to5_Up/F");
      TBranch *newBranch110 = t->Branch("phi_sv_JetEta0to5_Up", &svFitPhi_JetEta0to5_Up, "phi_sv_JetEta0to5_Up/F");
      TBranch *newBranch111 = t->Branch("met_sv_JetEta0to5_Up", &svFitMET_JetEta0to5_Up, "met_sv_JetEta0to5_Up/F");
      TBranch *newBranch112 = t->Branch("mt_sv_JetEta0to5_Up", &svFitTransverseMass_JetEta0to5_Up, "mt_sv_JetEta0to5_Up/F");

      TBranch *newBranch113 = t->Branch("m_sv_JetEta0to5_Down", &svFitMass_JetEta0to5_Down, "m_sv_JetEta0to5_Down/F");
      TBranch *newBranch114 = t->Branch("pt_sv_JetEta0to5_Down", &svFitPt_JetEta0to5_Down, "pt_sv_JetEta0to5_Down/F");
      TBranch *newBranch115 = t->Branch("eta_sv_JetEta0to5_Down", &svFitEta_JetEta0to5_Down, "eta_sv_JetEta0to5_Down/F");
      TBranch *newBranch116 = t->Branch("phi_sv_JetEta0to5_Down", &svFitPhi_JetEta0to5_Down, "phi_sv_JetEta0to5_Down/F");
      TBranch *newBranch117 = t->Branch("met_sv_JetEta0to5_Down", &svFitMET_JetEta0to5_Down, "met_sv_JetEta0to5_Down/F");
      TBranch *newBranch118 = t->Branch("mt_sv_JetEta0to5_Down", &svFitTransverseMass_JetEta0to5_Down, "mt_sv_JetEta0to5_Down/F");

      TBranch *newBranch119 = t->Branch("m_sv_JetEta3to5_Up", &svFitMass_JetEta3to5_Up, "m_sv_JetEta3to5_Up/F");
      TBranch *newBranch120 = t->Branch("pt_sv_JetEta3to5_Up", &svFitPt_JetEta3to5_Up, "pt_sv_JetEta3to5_Up/F");
      TBranch *newBranch121 = t->Branch("eta_sv_JetEta3to5_Up", &svFitEta_JetEta3to5_Up, "eta_sv_JetEta3to5_Up/F");
      TBranch *newBranch122 = t->Branch("phi_sv_JetEta3to5_Up", &svFitPhi_JetEta3to5_Up, "phi_sv_JetEta3to5_Up/F");
      TBranch *newBranch123 = t->Branch("met_sv_JetEta3to5_Up", &svFitMET_JetEta3to5_Up, "met_sv_JetEta3to5_Up/F");
      TBranch *newBranch124 = t->Branch("mt_sv_JetEta3to5_Up", &svFitTransverseMass_JetEta3to5_Up, "mt_sv_JetEta3to5_Up/F");

      TBranch *newBranch125 = t->Branch("m_sv_JetEta3to5_Down", &svFitMass_JetEta3to5_Down, "m_sv_JetEta3to5_Down/F");
      TBranch *newBranch126 = t->Branch("pt_sv_JetEta3to5_Down", &svFitPt_JetEta3to5_Down, "pt_sv_JetEta3to5_Down/F");
      TBranch *newBranch127 = t->Branch("eta_sv_JetEta3to5_Down", &svFitEta_JetEta3to5_Down, "eta_sv_JetEta3to5_Down/F");
      TBranch *newBranch128 = t->Branch("phi_sv_JetEta3to5_Down", &svFitPhi_JetEta3to5_Down, "phi_sv_JetEta3to5_Down/F");
      TBranch *newBranch129 = t->Branch("met_sv_JetEta3to5_Down", &svFitMET_JetEta3to5_Down, "met_sv_JetEta3to5_Down/F");
      TBranch *newBranch130 = t->Branch("mt_sv_JetEta3to5_Down", &svFitTransverseMass_JetEta3to5_Down, "mt_sv_JetEta3to5_Down/F");

      TBranch *newBranch131 = t->Branch("m_sv_JetRelativeBal_Up", &svFitMass_JetRelativeBal_Up, "m_sv_JetRelativeBal_Up/F");
      TBranch *newBranch132 = t->Branch("pt_sv_JetRelativeBal_Up", &svFitPt_JetRelativeBal_Up, "pt_sv_JetRelativeBal_Up/F");
      TBranch *newBranch133 = t->Branch("eta_sv_JetRelativeBal_Up", &svFitEta_JetRelativeBal_Up, "eta_sv_JetRelativeBal_Up/F");
      TBranch *newBranch134 = t->Branch("phi_sv_JetRelativeBal_Up", &svFitPhi_JetRelativeBal_Up, "phi_sv_JetRelativeBal_Up/F");
      TBranch *newBranch135 = t->Branch("met_sv_JetRelativeBal_Up", &svFitMET_JetRelativeBal_Up, "met_sv_JetRelativeBal_Up/F");
      TBranch *newBranch136 = t->Branch("mt_sv_JetRelativeBal_Up", &svFitTransverseMass_JetRelativeBal_Up, "mt_sv_JetRelativeBal_Up/F");

      TBranch *newBranch137 = t->Branch("m_sv_JetRelativeBal_Down", &svFitMass_JetRelativeBal_Down, "m_sv_JetRelativeBal_Down/F");
      TBranch *newBranch138 = t->Branch("pt_sv_JetRelativeBal_Down", &svFitPt_JetRelativeBal_Down, "pt_sv_JetRelativeBal_Down/F");
      TBranch *newBranch139 = t->Branch("eta_sv_JetRelativeBal_Down", &svFitEta_JetRelativeBal_Down, "eta_sv_JetRelativeBal_Down/F");
      TBranch *newBranch140 = t->Branch("phi_sv_JetRelativeBal_Down", &svFitPhi_JetRelativeBal_Down, "phi_sv_JetRelativeBal_Down/F");
      TBranch *newBranch141 = t->Branch("met_sv_JetRelativeBal_Down", &svFitMET_JetRelativeBal_Down, "met_sv_JetRelativeBal_Down/F");
      TBranch *newBranch142 = t->Branch("mt_sv_JetRelativeBal_Down", &svFitTransverseMass_JetRelativeBal_Down, "mt_sv_JetRelativeBal_Down/F");

      TBranch *newBranch143 = t->Branch("m_sv_JetRelativeSample_Up", &svFitMass_JetRelativeSample_Up, "m_sv_JetRelativeSample_Up/F");
      TBranch *newBranch144 = t->Branch("pt_sv_JetRelativeSample_Up", &svFitPt_JetRelativeSample_Up, "pt_sv_JetRelativeSample_Up/F");
      TBranch *newBranch145 = t->Branch("eta_sv_JetRelativeSample_Up", &svFitEta_JetRelativeSample_Up, "eta_sv_JetRelativeSample_Up/F");
      TBranch *newBranch146 = t->Branch("phi_sv_JetRelativeSample_Up", &svFitPhi_JetRelativeSample_Up, "phi_sv_JetRelativeSample_Up/F");
      TBranch *newBranch147 = t->Branch("met_sv_JetRelativeSample_Up", &svFitMET_JetRelativeSample_Up, "met_sv_JetRelativeSample_Up/F");
      TBranch *newBranch148 = t->Branch("mt_sv_JetRelativeSample_Up", &svFitTransverseMass_JetRelativeSample_Up, "mt_sv_JetRelativeSample_Up/F");

      TBranch *newBranch149 = t->Branch("m_sv_JetRelativeSample_Down", &svFitMass_JetRelativeSample_Down, "m_sv_JetRelativeSample_Down/F");
      TBranch *newBranch150 = t->Branch("pt_sv_JetRelativeSample_Down", &svFitPt_JetRelativeSample_Down, "pt_sv_JetRelativeSample_Down/F");
      TBranch *newBranch151 = t->Branch("eta_sv_JetRelativeSample_Down", &svFitEta_JetRelativeSample_Down, "eta_sv_JetRelativeSample_Down/F");
      TBranch *newBranch152 = t->Branch("phi_sv_JetRelativeSample_Down", &svFitPhi_JetRelativeSample_Down, "phi_sv_JetRelativeSample_Down/F");
      TBranch *newBranch153 = t->Branch("met_sv_JetRelativeSample_Down", &svFitMET_JetRelativeSample_Down, "met_sv_JetRelativeSample_Down/F");
      TBranch *newBranch154 = t->Branch("mt_sv_JetRelativeSample_Down", &svFitTransverseMass_JetRelativeSample_Down, "mt_sv_JetRelativeSample_Down/F");


      TBranch *newBranch155 = t->Branch("metcorClusteredDown",    &metcorClusteredDown,   "metcorClusteredDown/F");
      TBranch *newBranch156 = t->Branch("metcorphiClusteredDown", &metcorphiClusteredDown,"metcorphiClusteredDown/F");
      TBranch *newBranch157 = t->Branch("metcorClusteredUp",      &metcorClusteredUp,     "metcorClusteredUp/F");
      TBranch *newBranch158 = t->Branch("metcorphiClusteredUp",   &metcorphiClusteredUp,  "metcorphiClusteredUp/F");
      TBranch *newBranch159 = t->Branch("metcorUncDown",          &metcorUncDown,         "metcorUncDown/F");
      TBranch *newBranch160 = t->Branch("metcorphiUncDown",       &metcorphiUncDown,      "metcorphiUncDown/F");
      TBranch *newBranch161 = t->Branch("metcorUncUp",            &metcorUncUp,           "metcorUncUp/F");
      TBranch *newBranch162 = t->Branch("metcorphiUncUp",         &metcorphiUncUp,        "metcorphiUncUp/F");

      TBranch *newBranch163 = t->Branch("metcorJetEC2Down",          &metcorJetEC2Down,         "metcorJetEC2Down/F");
      TBranch *newBranch164 = t->Branch("metcorphiJetEC2Down",       &metcorphiJetEC2Down,      "metcorphiJetEC2Down/F");
      TBranch *newBranch165 = t->Branch("metcorJetEC2Up",            &metcorJetEC2Up,           "metcorJetEC2Up/F");
      TBranch *newBranch166 = t->Branch("metcorphiJetEC2Up",         &metcorphiJetEC2Up,        "metcorphiJetEC2Up/F");
      TBranch *newBranch167 = t->Branch("metcorJetEta0to3Down",          &metcorJetEta0to3Down,         "metcorJetEta0to3Down/F");
      TBranch *newBranch168 = t->Branch("metcorphiJetEta0to3Down",       &metcorphiJetEta0to3Down,      "metcorphiJetEta0to3Down/F");
      TBranch *newBranch169 = t->Branch("metcorJetEta0to3Up",            &metcorJetEta0to3Up,           "metcorJetEta0to3Up/F");
      TBranch *newBranch170 = t->Branch("metcorphiJetEta0to3Up",         &metcorphiJetEta0to3Up,        "metcorphiJetEta0to3Up/F");
      TBranch *newBranch171 = t->Branch("metcorJetEta0to5Down",          &metcorJetEta0to5Down,         "metcorJetEta0to5Down/F");
      TBranch *newBranch172 = t->Branch("metcorphiJetEta0to5Down",       &metcorphiJetEta0to5Down,      "metcorphiJetEta0to5Down/F");
      TBranch *newBranch173 = t->Branch("metcorJetEta0to5Up",            &metcorJetEta0to5Up,           "metcorJetEta0to5Up/F");
      TBranch *newBranch174 = t->Branch("metcorphiJetEta0to5Up",         &metcorphiJetEta0to5Up,        "metcorphiJetEta0to5Up/F");
      TBranch *newBranch175 = t->Branch("metcorJetEta3to5Down",          &metcorJetEta3to5Down,         "metcorJetEta3to5Down/F");
      TBranch *newBranch176 = t->Branch("metcorphiJetEta3to5Down",       &metcorphiJetEta3to5Down,      "metcorphiJetEta3to5Down/F");
      TBranch *newBranch177 = t->Branch("metcorJetEta3to5Up",            &metcorJetEta3to5Up,           "metcorJetEta3to5Up/F");
      TBranch *newBranch178 = t->Branch("metcorphiJetEta3to5Up",         &metcorphiJetEta3to5Up,        "metcorphiJetEta3to5Up/F");
      TBranch *newBranch179 = t->Branch("metcorJetRelativeBalDown",          &metcorJetRelativeBalDown,         "metcorJetRelativeBalDown/F");
      TBranch *newBranch180 = t->Branch("metcorphiJetRelativeBalDown",       &metcorphiJetRelativeBalDown,      "metcorphiJetRelativeBalDown/F");
      TBranch *newBranch181 = t->Branch("metcorJetRelativeBalUp",            &metcorJetRelativeBalUp,           "metcorJetRelativeBalUp/F");
      TBranch *newBranch182 = t->Branch("metcorphiJetRelativeBalUp",         &metcorphiJetRelativeBalUp,        "metcorphiJetRelativeBalUp/F");
      TBranch *newBranch183 = t->Branch("metcorJetRelativeSampleDown",          &metcorJetRelativeSampleDown,         "metcorJetRelativeSampleDown/F");
      TBranch *newBranch184 = t->Branch("metcorphiJetRelativeSampleDown",       &metcorphiJetRelativeSampleDown,      "metcorphiJetRelativeSampleDown/F");
      TBranch *newBranch185 = t->Branch("metcorJetRelativeSampleUp",            &metcorJetRelativeSampleUp,           "metcorJetRelativeSampleUp/F");
      TBranch *newBranch186 = t->Branch("metcorphiJetRelativeSampleUp",         &metcorphiJetRelativeSampleUp,        "metcorphiJetRelativeSampleUp/F");

    
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

      // up JetEC2
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEC2_Up",  &tau1_pt_JetEC2_Up,  "tau1_pt_JetEC2_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEC2_Up", &tau1_eta_JetEC2_Up, "tau1_eta_JetEC2_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEC2_Up", &tau1_phi_JetEC2_Up, "tau1_phi_JetEC2_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEC2_Up",   &tau1_m_JetEC2_Up,   "tau1_m_JetEC2_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEC2_Up",  &tau2_pt_JetEC2_Up,  "tau2_pt_JetEC2_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEC2_Up", &tau2_eta_JetEC2_Up, "tau2_eta_JetEC2_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEC2_Up", &tau2_phi_JetEC2_Up, "tau2_phi_JetEC2_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEC2_Up",   &tau2_m_JetEC2_Up,   "tau2_m_JetEC2_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEC2_Down",  &tau1_pt_JetEC2_Down,  "tau1_pt_JetEC2_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEC2_Down", &tau1_eta_JetEC2_Down, "tau1_eta_JetEC2_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEC2_Down", &tau1_phi_JetEC2_Down, "tau1_phi_JetEC2_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEC2_Down",   &tau1_m_JetEC2_Down,   "tau1_m_JetEC2_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEC2_Down",  &tau2_pt_JetEC2_Down,  "tau2_pt_JetEC2_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEC2_Down", &tau2_eta_JetEC2_Down, "tau2_eta_JetEC2_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEC2_Down", &tau2_phi_JetEC2_Down, "tau2_phi_JetEC2_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEC2_Down",   &tau2_m_JetEC2_Down,   "tau2_m_JetEC2_Down/F"));
      // up JetEta0to3
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEta0to3_Up",  &tau1_pt_JetEta0to3_Up,  "tau1_pt_JetEta0to3_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEta0to3_Up", &tau1_eta_JetEta0to3_Up, "tau1_eta_JetEta0to3_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEta0to3_Up", &tau1_phi_JetEta0to3_Up, "tau1_phi_JetEta0to3_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEta0to3_Up",   &tau1_m_JetEta0to3_Up,   "tau1_m_JetEta0to3_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEta0to3_Up",  &tau2_pt_JetEta0to3_Up,  "tau2_pt_JetEta0to3_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEta0to3_Up", &tau2_eta_JetEta0to3_Up, "tau2_eta_JetEta0to3_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEta0to3_Up", &tau2_phi_JetEta0to3_Up, "tau2_phi_JetEta0to3_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEta0to3_Up",   &tau2_m_JetEta0to3_Up,   "tau2_m_JetEta0to3_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEta0to3_Down",  &tau1_pt_JetEta0to3_Down,  "tau1_pt_JetEta0to3_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEta0to3_Down", &tau1_eta_JetEta0to3_Down, "tau1_eta_JetEta0to3_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEta0to3_Down", &tau1_phi_JetEta0to3_Down, "tau1_phi_JetEta0to3_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEta0to3_Down",   &tau1_m_JetEta0to3_Down,   "tau1_m_JetEta0to3_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEta0to3_Down",  &tau2_pt_JetEta0to3_Down,  "tau2_pt_JetEta0to3_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEta0to3_Down", &tau2_eta_JetEta0to3_Down, "tau2_eta_JetEta0to3_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEta0to3_Down", &tau2_phi_JetEta0to3_Down, "tau2_phi_JetEta0to3_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEta0to3_Down",   &tau2_m_JetEta0to3_Down,   "tau2_m_JetEta0to3_Down/F"));
      // up JetEta0to5
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEta0to5_Up",  &tau1_pt_JetEta0to5_Up,  "tau1_pt_JetEta0to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEta0to5_Up", &tau1_eta_JetEta0to5_Up, "tau1_eta_JetEta0to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEta0to5_Up", &tau1_phi_JetEta0to5_Up, "tau1_phi_JetEta0to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEta0to5_Up",   &tau1_m_JetEta0to5_Up,   "tau1_m_JetEta0to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEta0to5_Up",  &tau2_pt_JetEta0to5_Up,  "tau2_pt_JetEta0to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEta0to5_Up", &tau2_eta_JetEta0to5_Up, "tau2_eta_JetEta0to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEta0to5_Up", &tau2_phi_JetEta0to5_Up, "tau2_phi_JetEta0to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEta0to5_Up",   &tau2_m_JetEta0to5_Up,   "tau2_m_JetEta0to5_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEta0to5_Down",  &tau1_pt_JetEta0to5_Down,  "tau1_pt_JetEta0to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEta0to5_Down", &tau1_eta_JetEta0to5_Down, "tau1_eta_JetEta0to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEta0to5_Down", &tau1_phi_JetEta0to5_Down, "tau1_phi_JetEta0to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEta0to5_Down",   &tau1_m_JetEta0to5_Down,   "tau1_m_JetEta0to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEta0to5_Down",  &tau2_pt_JetEta0to5_Down,  "tau2_pt_JetEta0to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEta0to5_Down", &tau2_eta_JetEta0to5_Down, "tau2_eta_JetEta0to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEta0to5_Down", &tau2_phi_JetEta0to5_Down, "tau2_phi_JetEta0to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEta0to5_Down",   &tau2_m_JetEta0to5_Down,   "tau2_m_JetEta0to5_Down/F"));
      // up JetEta3to5
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEta3to5_Up",  &tau1_pt_JetEta3to5_Up,  "tau1_pt_JetEta3to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEta3to5_Up", &tau1_eta_JetEta3to5_Up, "tau1_eta_JetEta3to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEta3to5_Up", &tau1_phi_JetEta3to5_Up, "tau1_phi_JetEta3to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEta3to5_Up",   &tau1_m_JetEta3to5_Up,   "tau1_m_JetEta3to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEta3to5_Up",  &tau2_pt_JetEta3to5_Up,  "tau2_pt_JetEta3to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEta3to5_Up", &tau2_eta_JetEta3to5_Up, "tau2_eta_JetEta3to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEta3to5_Up", &tau2_phi_JetEta3to5_Up, "tau2_phi_JetEta3to5_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEta3to5_Up",   &tau2_m_JetEta3to5_Up,   "tau2_m_JetEta3to5_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetEta3to5_Down",  &tau1_pt_JetEta3to5_Down,  "tau1_pt_JetEta3to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetEta3to5_Down", &tau1_eta_JetEta3to5_Down, "tau1_eta_JetEta3to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetEta3to5_Down", &tau1_phi_JetEta3to5_Down, "tau1_phi_JetEta3to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetEta3to5_Down",   &tau1_m_JetEta3to5_Down,   "tau1_m_JetEta3to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetEta3to5_Down",  &tau2_pt_JetEta3to5_Down,  "tau2_pt_JetEta3to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetEta3to5_Down", &tau2_eta_JetEta3to5_Down, "tau2_eta_JetEta3to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetEta3to5_Down", &tau2_phi_JetEta3to5_Down, "tau2_phi_JetEta3to5_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetEta3to5_Down",   &tau2_m_JetEta3to5_Down,   "tau2_m_JetEta3to5_Down/F"));
      // up JetRelativeBal
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetRelativeBal_Up",  &tau1_pt_JetRelativeBal_Up,  "tau1_pt_JetRelativeBal_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetRelativeBal_Up", &tau1_eta_JetRelativeBal_Up, "tau1_eta_JetRelativeBal_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetRelativeBal_Up", &tau1_phi_JetRelativeBal_Up, "tau1_phi_JetRelativeBal_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetRelativeBal_Up",   &tau1_m_JetRelativeBal_Up,   "tau1_m_JetRelativeBal_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetRelativeBal_Up",  &tau2_pt_JetRelativeBal_Up,  "tau2_pt_JetRelativeBal_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetRelativeBal_Up", &tau2_eta_JetRelativeBal_Up, "tau2_eta_JetRelativeBal_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetRelativeBal_Up", &tau2_phi_JetRelativeBal_Up, "tau2_phi_JetRelativeBal_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetRelativeBal_Up",   &tau2_m_JetRelativeBal_Up,   "tau2_m_JetRelativeBal_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetRelativeBal_Down",  &tau1_pt_JetRelativeBal_Down,  "tau1_pt_JetRelativeBal_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetRelativeBal_Down", &tau1_eta_JetRelativeBal_Down, "tau1_eta_JetRelativeBal_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetRelativeBal_Down", &tau1_phi_JetRelativeBal_Down, "tau1_phi_JetRelativeBal_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetRelativeBal_Down",   &tau1_m_JetRelativeBal_Down,   "tau1_m_JetRelativeBal_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetRelativeBal_Down",  &tau2_pt_JetRelativeBal_Down,  "tau2_pt_JetRelativeBal_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetRelativeBal_Down", &tau2_eta_JetRelativeBal_Down, "tau2_eta_JetRelativeBal_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetRelativeBal_Down", &tau2_phi_JetRelativeBal_Down, "tau2_phi_JetRelativeBal_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetRelativeBal_Down",   &tau2_m_JetRelativeBal_Down,   "tau2_m_JetRelativeBal_Down/F"));
      // up  JetRelativeSample
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetRelativeSample_Up",  &tau1_pt_JetRelativeSample_Up,  "tau1_pt_JetRelativeSample_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetRelativeSample_Up", &tau1_eta_JetRelativeSample_Up, "tau1_eta_JetRelativeSample_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetRelativeSample_Up", &tau1_phi_JetRelativeSample_Up, "tau1_phi_JetRelativeSample_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetRelativeSample_Up",   &tau1_m_JetRelativeSample_Up,   "tau1_m_JetRelativeSample_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetRelativeSample_Up",  &tau2_pt_JetRelativeSample_Up,  "tau2_pt_JetRelativeSample_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetRelativeSample_Up", &tau2_eta_JetRelativeSample_Up, "tau2_eta_JetRelativeSample_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetRelativeSample_Up", &tau2_phi_JetRelativeSample_Up, "tau2_phi_JetRelativeSample_Up/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetRelativeSample_Up",   &tau2_m_JetRelativeSample_Up,   "tau2_m_JetRelativeSample_Up/F"));
      // down
      tau4VectorBranches.push_back(t->Branch("tau1_pt_JetRelativeSample_Down",  &tau1_pt_JetRelativeSample_Down,  "tau1_pt_JetRelativeSample_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_eta_JetRelativeSample_Down", &tau1_eta_JetRelativeSample_Down, "tau1_eta_JetRelativeSample_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_phi_JetRelativeSample_Down", &tau1_phi_JetRelativeSample_Down, "tau1_phi_JetRelativeSample_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau1_m_JetRelativeSample_Down",   &tau1_m_JetRelativeSample_Down,   "tau1_m_JetRelativeSample_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_pt_JetRelativeSample_Down",  &tau2_pt_JetRelativeSample_Down,  "tau2_pt_JetRelativeSample_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_eta_JetRelativeSample_Down", &tau2_eta_JetRelativeSample_Down, "tau2_eta_JetRelativeSample_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_phi_JetRelativeSample_Down", &tau2_phi_JetRelativeSample_Down, "tau2_phi_JetRelativeSample_Down/F"));
      tau4VectorBranches.push_back(t->Branch("tau2_m_JetRelativeSample_Down",   &tau2_m_JetRelativeSample_Down,   "tau2_m_JetRelativeSample_Down/F"));
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

      // JES Uncertainties
      float met_JetEC2Up, met_JetEC2Down, metphi_JetEC2Up, metphi_JetEC2Down;
      float JetEC2UpMETx = 0., JetEC2UpMETy = 0., JetEC2DownMETx = 0., JetEC2DownMETy = 0.;
      float met_JetEta0to3Up, met_JetEta0to3Down, metphi_JetEta0to3Up, metphi_JetEta0to3Down;
      float JetEta0to3UpMETx = 0., JetEta0to3UpMETy = 0., JetEta0to3DownMETx = 0., JetEta0to3DownMETy = 0.;
      float met_JetEta0to5Up, met_JetEta0to5Down, metphi_JetEta0to5Up, metphi_JetEta0to5Down;
      float JetEta0to5UpMETx = 0., JetEta0to5UpMETy = 0., JetEta0to5DownMETx = 0., JetEta0to5DownMETy = 0.;
      float met_JetEta3to5Up, met_JetEta3to5Down, metphi_JetEta3to5Up, metphi_JetEta3to5Down;
      float JetEta3to5UpMETx = 0., JetEta3to5UpMETy = 0., JetEta3to5DownMETx = 0., JetEta3to5DownMETy = 0.;
      float met_JetRelativeBalUp, met_JetRelativeBalDown, metphi_JetRelativeBalUp, metphi_JetRelativeBalDown;
      float JetRelativeBalUpMETx = 0., JetRelativeBalUpMETy = 0., JetRelativeBalDownMETx = 0., JetRelativeBalDownMETy = 0.;
      float met_JetRelativeSampleUp, met_JetRelativeSampleDown, metphi_JetRelativeSampleUp, metphi_JetRelativeSampleDown;
      float JetRelativeSampleUpMETx = 0., JetRelativeSampleUpMETy = 0., JetRelativeSampleDownMETx = 0., JetRelativeSampleDownMETy = 0.;

      //ele/mu variables
      TBranch *pt1branch;
       
      t->SetBranchAddress("NtupleVer",&NtupleVer);
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
      t->SetBranchAddress("met_UESUp", &uncMetPtUp);
      t->SetBranchAddress("met_UESDown", &uncMetPtDown);
      t->SetBranchAddress("metphi_UESUp", &uncMetPhiUp);
      t->SetBranchAddress("metphi_UESDown", &uncMetPhiDown);
      
      t->SetBranchAddress("met_JESUp", &clusteredMetPtUp);
      t->SetBranchAddress("met_JESDown", &clusteredMetPtDown);
      t->SetBranchAddress("metphi_JESUp", &clusteredMetPhiUp);
      t->SetBranchAddress("metphi_JESDown", &clusteredMetPhiDown);

      // JES Unc - JetEC2
      t->SetBranchAddress("met_JetEC2Up", &met_JetEC2Up);
      t->SetBranchAddress("met_JetEC2Down", &met_JetEC2Down);
      t->SetBranchAddress("metphi_JetEC2Up", &metphi_JetEC2Up);
      t->SetBranchAddress("metphi_JetEC2Down", &metphi_JetEC2Down);

      // JES Unc - JetEta0to3
      t->SetBranchAddress("met_JetEta0to3Up", &met_JetEta0to3Up);
      t->SetBranchAddress("met_JetEta0to3Down", &met_JetEta0to3Down);
      t->SetBranchAddress("metphi_JetEta0to3Up", &metphi_JetEta0to3Up);
      t->SetBranchAddress("metphi_JetEta0to3Down", &metphi_JetEta0to3Down);

      // JES Unc - JetEta0to5
      t->SetBranchAddress("met_JetEta0to5Up", &met_JetEta0to5Up);
      t->SetBranchAddress("met_JetEta0to5Down", &met_JetEta0to5Down);
      t->SetBranchAddress("metphi_JetEta0to5Up", &metphi_JetEta0to5Up);
      t->SetBranchAddress("metphi_JetEta0to5Down", &metphi_JetEta0to5Down);

      // JES Unc - JetEta3to5
      t->SetBranchAddress("met_JetEta3to5Up", &met_JetEta3to5Up);
      t->SetBranchAddress("met_JetEta3to5Down", &met_JetEta3to5Down);
      t->SetBranchAddress("metphi_JetEta3to5Up", &metphi_JetEta3to5Up);
      t->SetBranchAddress("metphi_JetEta3to5Down", &metphi_JetEta3to5Down);

      // JES Unc - JetRelativeBal
      t->SetBranchAddress("met_JetRelativeBalUp", &met_JetRelativeBalUp);
      t->SetBranchAddress("met_JetRelativeBalDown", &met_JetRelativeBalDown);
      t->SetBranchAddress("metphi_JetRelativeBalUp", &metphi_JetRelativeBalUp);
      t->SetBranchAddress("metphi_JetRelativeBalDown", &metphi_JetRelativeBalDown);

      // JES Unc - JetRelativeSample
      t->SetBranchAddress("met_JetRelativeSampleUp", &met_JetRelativeSampleUp);
      t->SetBranchAddress("met_JetRelativeSampleDown", &met_JetRelativeSampleDown);
      t->SetBranchAddress("metphi_JetRelativeSampleUp", &metphi_JetRelativeSampleUp);
      t->SetBranchAddress("metphi_JetRelativeSampleDown", &metphi_JetRelativeSampleDown);

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
	  JetEC2UpMETx         = met_JetEC2Up*TMath::Cos(metphi_JetEC2Up);
	  JetEC2UpMETy         = met_JetEC2Up*TMath::Sin(metphi_JetEC2Up);
	  JetEC2DownMETx       = met_JetEC2Down*TMath::Cos(metphi_JetEC2Down);
	  JetEC2DownMETy       = met_JetEC2Down*TMath::Sin(metphi_JetEC2Down);
	  JetEta0to3UpMETx         = met_JetEta0to3Up*TMath::Cos(metphi_JetEta0to3Up);
	  JetEta0to3UpMETy         = met_JetEta0to3Up*TMath::Sin(metphi_JetEta0to3Up);
	  JetEta0to3DownMETx       = met_JetEta0to3Down*TMath::Cos(metphi_JetEta0to3Down);
	  JetEta0to3DownMETy       = met_JetEta0to3Down*TMath::Sin(metphi_JetEta0to3Down);
	  JetEta0to5UpMETx         = met_JetEta0to5Up*TMath::Cos(metphi_JetEta0to5Up);
	  JetEta0to5UpMETy         = met_JetEta0to5Up*TMath::Sin(metphi_JetEta0to5Up);
	  JetEta0to5DownMETx       = met_JetEta0to5Down*TMath::Cos(metphi_JetEta0to5Down);
	  JetEta0to5DownMETy       = met_JetEta0to5Down*TMath::Sin(metphi_JetEta0to5Down);
	  JetEta3to5UpMETx         = met_JetEta3to5Up*TMath::Cos(metphi_JetEta3to5Up);
	  JetEta3to5UpMETy         = met_JetEta3to5Up*TMath::Sin(metphi_JetEta3to5Up);
	  JetEta3to5DownMETx       = met_JetEta3to5Down*TMath::Cos(metphi_JetEta3to5Down);
	  JetEta3to5DownMETy       = met_JetEta3to5Down*TMath::Sin(metphi_JetEta3to5Down);
	  JetRelativeBalUpMETx         = met_JetRelativeBalUp*TMath::Cos(metphi_JetRelativeBalUp);
	  JetRelativeBalUpMETy         = met_JetRelativeBalUp*TMath::Sin(metphi_JetRelativeBalUp);
	  JetRelativeBalDownMETx       = met_JetRelativeBalDown*TMath::Cos(metphi_JetRelativeBalDown);
	  JetRelativeBalDownMETy       = met_JetRelativeBalDown*TMath::Sin(metphi_JetRelativeBalDown);
	  JetRelativeSampleUpMETx         = met_JetRelativeSampleUp*TMath::Cos(metphi_JetRelativeSampleUp);
	  JetRelativeSampleUpMETy         = met_JetRelativeSampleUp*TMath::Sin(metphi_JetRelativeSampleUp);
	  JetRelativeSampleDownMETx       = met_JetRelativeSampleDown*TMath::Cos(metphi_JetRelativeSampleDown);
	  JetRelativeSampleDownMETy       = met_JetRelativeSampleDown*TMath::Sin(metphi_JetRelativeSampleDown);

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
	// JES
	metcorrJetEC2Up_ex = JetEC2UpMETx;
	metcorrJetEC2Up_ey = JetEC2UpMETy;
	metcorrJetEC2Down_ex = JetEC2DownMETx;
	metcorrJetEC2Down_ey = JetEC2DownMETy;
	metcorrJetEta0to3Up_ex = JetEta0to3UpMETx;
	metcorrJetEta0to3Up_ey = JetEta0to3UpMETy;
	metcorrJetEta0to3Down_ex = JetEta0to3DownMETx;
	metcorrJetEta0to3Down_ey = JetEta0to3DownMETy;
	metcorrJetEta0to5Up_ex = JetEta0to5UpMETx;
	metcorrJetEta0to5Up_ey = JetEta0to5UpMETy;
	metcorrJetEta0to5Down_ex = JetEta0to5DownMETx;
	metcorrJetEta0to5Down_ey = JetEta0to5DownMETy;
	metcorrJetEta3to5Up_ex = JetEta3to5UpMETx;
	metcorrJetEta3to5Up_ey = JetEta3to5UpMETy;
	metcorrJetEta3to5Down_ex = JetEta3to5DownMETx;
	metcorrJetEta3to5Down_ey = JetEta3to5DownMETy;
	metcorrJetRelativeBalUp_ex = JetRelativeBalUpMETx;
	metcorrJetRelativeBalUp_ey = JetRelativeBalUpMETy;
	metcorrJetRelativeBalDown_ex = JetRelativeBalDownMETx;
	metcorrJetRelativeBalDown_ey = JetRelativeBalDownMETy;
	metcorrJetRelativeSampleUp_ex = JetRelativeSampleUpMETx;
	metcorrJetRelativeSampleUp_ey = JetRelativeSampleUpMETy;
	metcorrJetRelativeSampleDown_ex = JetRelativeSampleDownMETx;
	metcorrJetRelativeSampleDown_ey = JetRelativeSampleDownMETy;
	
        
        metcor = TMath::Sqrt( metcorr_ex*metcorr_ex + metcorr_ey*metcorr_ey);
        metcorphi = TMath::ATan2( metcorr_ey, metcorr_ex );
        std::cout << " - metcor "<<metcor<<" metcorphi "<<metcorphi<<std::endl;
        
        // Corrected MET values for saving
        metcorClusteredDown = TMath::Sqrt( metcorrClusteredDown_ex*metcorrClusteredDown_ex + metcorrClusteredDown_ey*metcorrClusteredDown_ey);
        metcorphiClusteredDown = TMath::ATan2( metcorrClusteredDown_ey, metcorrClusteredDown_ex );
        
        metcorClusteredUp = TMath::Sqrt( metcorrClusteredUp_ex*metcorrClusteredUp_ex + metcorrClusteredUp_ey*metcorrClusteredUp_ey);
        metcorphiClusteredUp = TMath::ATan2( metcorrClusteredUp_ey, metcorrClusteredUp_ex );
        
        metcorUncDown = TMath::Sqrt( metcorrUncDown_ex*metcorrUncDown_ex + metcorrUncDown_ey*metcorrUncDown_ey);
        metcorphiUncDown = TMath::ATan2( metcorrUncDown_ey, metcorrUncDown_ex );
        
        metcorUncUp = TMath::Sqrt( metcorrUncUp_ex*metcorrUncUp_ex + metcorrUncUp_ey*metcorrUncUp_ey);
        metcorphiUncUp = TMath::ATan2( metcorrUncUp_ey, metcorrUncUp_ex );
     
	metcorJetEC2Down = TMath::Sqrt( metcorrJetEC2Down_ex*metcorrJetEC2Down_ex + metcorrJetEC2Down_ey*metcorrJetEC2Down_ey);
	metcorphiJetEC2Down = TMath::ATan2( metcorrJetEC2Down_ey, metcorrJetEC2Down_ex );

	metcorJetEC2Up = TMath::Sqrt( metcorrJetEC2Up_ex*metcorrJetEC2Up_ex + metcorrJetEC2Up_ey*metcorrJetEC2Up_ey);
	metcorphiJetEC2Up = TMath::ATan2( metcorrJetEC2Up_ey, metcorrJetEC2Up_ex );

	metcorJetEta0to3Down = TMath::Sqrt( metcorrJetEta0to3Down_ex*metcorrJetEta0to3Down_ex + metcorrJetEta0to3Down_ey*metcorrJetEta0to3Down_ey);
	metcorphiJetEta0to3Down = TMath::ATan2( metcorrJetEta0to3Down_ey, metcorrJetEta0to3Down_ex );

	metcorJetEta0to3Up = TMath::Sqrt( metcorrJetEta0to3Up_ex*metcorrJetEta0to3Up_ex + metcorrJetEta0to3Up_ey*metcorrJetEta0to3Up_ey);
	metcorphiJetEta0to3Up = TMath::ATan2( metcorrJetEta0to3Up_ey, metcorrJetEta0to3Up_ex );

	metcorJetEta0to5Down = TMath::Sqrt( metcorrJetEta0to5Down_ex*metcorrJetEta0to5Down_ex + metcorrJetEta0to5Down_ey*metcorrJetEta0to5Down_ey);
	metcorphiJetEta0to5Down = TMath::ATan2( metcorrJetEta0to5Down_ey, metcorrJetEta0to5Down_ex );

	metcorJetEta0to5Up = TMath::Sqrt( metcorrJetEta0to5Up_ex*metcorrJetEta0to5Up_ex + metcorrJetEta0to5Up_ey*metcorrJetEta0to5Up_ey);
	metcorphiJetEta0to5Up = TMath::ATan2( metcorrJetEta0to5Up_ey, metcorrJetEta0to5Up_ex );

	metcorJetEta3to5Down = TMath::Sqrt( metcorrJetEta3to5Down_ex*metcorrJetEta3to5Down_ex + metcorrJetEta3to5Down_ey*metcorrJetEta3to5Down_ey);
	metcorphiJetEta3to5Down = TMath::ATan2( metcorrJetEta3to5Down_ey, metcorrJetEta3to5Down_ex );

	metcorJetEta3to5Up = TMath::Sqrt( metcorrJetEta3to5Up_ex*metcorrJetEta3to5Up_ex + metcorrJetEta3to5Up_ey*metcorrJetEta3to5Up_ey);
	metcorphiJetEta3to5Up = TMath::ATan2( metcorrJetEta3to5Up_ey, metcorrJetEta3to5Up_ex );

	metcorJetRelativeSampleDown = TMath::Sqrt( metcorrJetRelativeSampleDown_ex*metcorrJetRelativeSampleDown_ex + metcorrJetRelativeSampleDown_ey*metcorrJetRelativeSampleDown_ey);
	metcorphiJetRelativeSampleDown = TMath::ATan2( metcorrJetRelativeSampleDown_ey, metcorrJetRelativeSampleDown_ex );

	metcorJetRelativeSampleUp = TMath::Sqrt( metcorrJetRelativeSampleUp_ex*metcorrJetRelativeSampleUp_ex + metcorrJetRelativeSampleUp_ey*metcorrJetRelativeSampleUp_ey);
	metcorphiJetRelativeSampleUp = TMath::ATan2( metcorrJetRelativeSampleUp_ey, metcorrJetRelativeSampleUp_ex );

	metcorJetRelativeSampleDown = TMath::Sqrt( metcorrJetRelativeSampleDown_ex*metcorrJetRelativeSampleDown_ex + metcorrJetRelativeSampleDown_ey*metcorrJetRelativeSampleDown_ey);
	metcorphiJetRelativeSampleDown = TMath::ATan2( metcorrJetRelativeSampleDown_ey, metcorrJetRelativeSampleDown_ex );

	metcorJetRelativeSampleUp = TMath::Sqrt( metcorrJetRelativeSampleUp_ex*metcorrJetRelativeSampleUp_ex + metcorrJetRelativeSampleUp_ey*metcorrJetRelativeSampleUp_ey);
	metcorphiJetRelativeSampleUp = TMath::ATan2( metcorrJetRelativeSampleUp_ey, metcorrJetRelativeSampleUp_ex );


	
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
  
// Added by Abdollah
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
       
       
       //*****************************************************
       // MET SYSTEMATICS
       // Taus Pt have been corrected (TEC)
       // Still need TEC propagated to shifted METs
       //*****************************************************
       
       std::cout << "MET Unclustered Energy Up   ---  ";
       runSVFit(measuredTauLeptons, metcorrUncUp_ex, metcorrUncUp_ey, covMET, 0, 
		svFitMass_UncMet_Up, svFitPt_UncMet_Up, svFitEta_UncMet_Up, svFitPhi_UncMet_Up, svFitMET_UncMet_Up, svFitTransverseMass_UncMet_Up, tau1_UncMet_Up, tau2_UncMet_Up);
       std::cout << "MET Unclustered Energy Down ---  ";
       runSVFit(measuredTauLeptons, metcorrUncDown_ex, metcorrUncDown_ey, covMET, 0, 
		svFitMass_UncMet_Down, svFitPt_UncMet_Down, svFitEta_UncMet_Down, svFitPhi_UncMet_Down, svFitMET_UncMet_Down, svFitTransverseMass_UncMet_Down, tau1_UncMet_Down, tau2_UncMet_Down);
       std::cout << "JES JetEC2 Energy Up   ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEC2Up_ex, metcorrJetEC2Up_ey, covMET, 0, 
		svFitMass_JetEC2_Up, svFitPt_JetEC2_Up, svFitEta_JetEC2_Up, svFitPhi_JetEC2_Up, svFitMET_JetEC2_Up, svFitTransverseMass_JetEC2_Up, tau1_JetEC2_Up, tau2_JetEC2_Up);
       std::cout << "JES JetEC2 Energy Down ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEC2Down_ex, metcorrJetEC2Down_ey, covMET, 0, 
		svFitMass_JetEC2_Down, svFitPt_JetEC2_Down, svFitEta_JetEC2_Down, svFitPhi_JetEC2_Down, svFitMET_JetEC2_Down, svFitTransverseMass_JetEC2_Down, tau1_JetEC2_Down, tau2_JetEC2_Down);
       std::cout << "JES JetEta0to3 Energy Up   ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEta0to3Up_ex, metcorrJetEta0to3Up_ey, covMET, 0, 
		svFitMass_JetEta0to3_Up, svFitPt_JetEta0to3_Up, svFitEta_JetEta0to3_Up, svFitPhi_JetEta0to3_Up, svFitMET_JetEta0to3_Up, svFitTransverseMass_JetEta0to3_Up, tau1_JetEta0to3_Up, tau2_JetEta0to3_Up);
       std::cout << "JES JetEta0to3 Energy Down ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEta0to3Down_ex, metcorrJetEta0to3Down_ey, covMET, 0, 
		svFitMass_JetEta0to3_Down, svFitPt_JetEta0to3_Down, svFitEta_JetEta0to3_Down, svFitPhi_JetEta0to3_Down, svFitMET_JetEta0to3_Down, svFitTransverseMass_JetEta0to3_Down, tau1_JetEta0to3_Down, tau2_JetEta0to3_Down);
       std::cout << "JES JetEta0to5 Energy Up   ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEta0to5Up_ex, metcorrJetEta0to5Up_ey, covMET, 0, 
		svFitMass_JetEta0to5_Up, svFitPt_JetEta0to5_Up, svFitEta_JetEta0to5_Up, svFitPhi_JetEta0to5_Up, svFitMET_JetEta0to5_Up, svFitTransverseMass_JetEta0to5_Up, tau1_JetEta0to5_Up, tau2_JetEta0to5_Up);
       std::cout << "JES JetEta0to5 Energy Down ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEta0to5Down_ex, metcorrJetEta0to5Down_ey, covMET, 0, 
		svFitMass_JetEta0to5_Down, svFitPt_JetEta0to5_Down, svFitEta_JetEta0to5_Down, svFitPhi_JetEta0to5_Down, svFitMET_JetEta0to5_Down, svFitTransverseMass_JetEta0to5_Down, tau1_JetEta0to5_Down, tau2_JetEta0to5_Down);
       std::cout << "JES JetEta3to5 Energy Up   ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEta3to5Up_ex, metcorrJetEta3to5Up_ey, covMET, 0, 
		svFitMass_JetEta3to5_Up, svFitPt_JetEta3to5_Up, svFitEta_JetEta3to5_Up, svFitPhi_JetEta3to5_Up, svFitMET_JetEta3to5_Up, svFitTransverseMass_JetEta3to5_Up, tau1_JetEta3to5_Up, tau2_JetEta3to5_Up);
       std::cout << "JES JetEta3to5 Energy Down ---  ";
       runSVFit(measuredTauLeptons, metcorrJetEta3to5Down_ex, metcorrJetEta3to5Down_ey, covMET, 0, 
		svFitMass_JetEta3to5_Down, svFitPt_JetEta3to5_Down, svFitEta_JetEta3to5_Down, svFitPhi_JetEta3to5_Down, svFitMET_JetEta3to5_Down, svFitTransverseMass_JetEta3to5_Down, tau1_JetEta3to5_Down, tau2_JetEta3to5_Down);
       std::cout << "JES JetRelativeBal Energy Up   ---  ";
       runSVFit(measuredTauLeptons, metcorrJetRelativeBalUp_ex, metcorrJetRelativeBalUp_ey, covMET, 0, 
		svFitMass_JetRelativeBal_Up, svFitPt_JetRelativeBal_Up, svFitEta_JetRelativeBal_Up, svFitPhi_JetRelativeBal_Up, svFitMET_JetRelativeBal_Up, svFitTransverseMass_JetRelativeBal_Up, tau1_JetRelativeBal_Up, tau2_JetRelativeBal_Up);
       std::cout << "JES JetRelativeBal Energy Down ---  ";
       runSVFit(measuredTauLeptons, metcorrJetRelativeBalDown_ex, metcorrJetRelativeBalDown_ey, covMET, 0, 
		svFitMass_JetRelativeBal_Down, svFitPt_JetRelativeBal_Down, svFitEta_JetRelativeBal_Down, svFitPhi_JetRelativeBal_Down, svFitMET_JetRelativeBal_Down, svFitTransverseMass_JetRelativeBal_Down, tau1_JetRelativeBal_Down, tau2_JetRelativeBal_Down);
       std::cout << "JES JetRelativeSample Energy Up   ---  ";
       runSVFit(measuredTauLeptons, metcorrJetRelativeSampleUp_ex, metcorrJetRelativeSampleUp_ey, covMET, 0, 
		svFitMass_JetRelativeSample_Up, svFitPt_JetRelativeSample_Up, svFitEta_JetRelativeSample_Up, svFitPhi_JetRelativeSample_Up, svFitMET_JetRelativeSample_Up, svFitTransverseMass_JetRelativeSample_Up, tau1_JetRelativeSample_Up, tau2_JetRelativeSample_Up);
       std::cout << "JES JetRelativeSample Energy Down ---  ";
       runSVFit(measuredTauLeptons, metcorrJetRelativeSampleDown_ex, metcorrJetRelativeSampleDown_ey, covMET, 0, 
		svFitMass_JetRelativeSample_Down, svFitPt_JetRelativeSample_Down, svFitEta_JetRelativeSample_Down, svFitPhi_JetRelativeSample_Down, svFitMET_JetRelativeSample_Down, svFitTransverseMass_JetRelativeSample_Down, tau1_JetRelativeSample_Down, tau2_JetRelativeSample_Down);


       /*
       std::cout << "MET Clustered Energy Up     ---  ";
       runSVFit(measuredTauLeptons, metcorrClusteredUp_ex, metcorrClusteredUp_ey, covMET, 0, 
		svFitMass_ClusteredMet_Up, svFitPt_ClusteredMet_Up, svFitEta_ClusteredMet_Up, svFitPhi_ClusteredMet_Up, svFitMET_ClusteredMet_Up, svFitTransverseMass_ClusteredMet_Up, tau1_ClusteredMet_Up, tau2_ClusteredMet_Up);
       std::cout << "MET Clustered Energy Down   ---  ";
       runSVFit(measuredTauLeptons, metcorrClusteredDown_ex, metcorrClusteredDown_ey, covMET, 0, 
		svFitMass_ClusteredMet_Down, svFitPt_ClusteredMet_Down, svFitEta_ClusteredMet_Down, svFitPhi_ClusteredMet_Down, svFitMET_ClusteredMet_Down, svFitTransverseMass_ClusteredMet_Down, tau1_ClusteredMet_Down, tau2_ClusteredMet_Down);
       std::cout<< "Shifted MET Summary:\nmetcorr_ex " << metcorr_ex << "\n --- metcorrUncUp_ex " << metcorrUncUp_ex << " metcorrUncDown_ex " << metcorrUncDown_ex
            << " metcorrClusteredUp_ex " << metcorrClusteredUp_ex << " metcorrClusteredDown_ex " << metcorrClusteredDown_ex << std::endl;
       std::cout<< "metcorr_ey " << metcorr_ey << "\n --- metcorrUncUp_ey " << metcorrUncUp_ey << " metcorrUncDown_ey " << metcorrUncDown_ey
            << " metcorrClusteredUp_ey " << metcorrClusteredUp_ey << " metcorrClusteredDown_ey " << metcorrClusteredDown_ey << std::endl;
       */
       
       // Corrected MET values for saving
       /*
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
       */

       if(doES) {
	 bool isDM0_1 = false, isDM1_1 = false, isDM10_1 = false, isDM0_2 = false, isDM1_2 = false, isDM10_2 = false;
	 if (gen_match_1==5) {
	   if (decayMode==0) isDM0_1 =true;
	   else if (decayMode==1) isDM1_1 =true;
	   else if (decayMode==10) isDM10_1 =true;
	 }
	 if (gen_match_2==5) {
	   if (decayMode2==0) isDM0_2 =true;
	   else if (decayMode2==1) isDM1_2 =true;
	   else if (decayMode2==10) isDM10_2 =true;
	 }

	 //***********************************//
	 //****  Tau DM0 shifted up/down  ****//
	 //***********************************//
	 if (isDM0_1 || isDM0_2) {
	   std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
	   std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
	   double pt1_Up, pt2_Up, pt1_Down, pt2_Down;
	   pt1_Up = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode), isDM0_1, +1);
	   pt2_Up = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2), isDM0_2, +1);
	   pt1_Down = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode), isDM0_1,-1);
	   pt2_Down = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2), isDM0_2, -1);
	   
	   measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
	   measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2, pt2_Up, eta2, phi2,  mass2, decayMode2));	   
	   measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
	   measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2, pt2_Down, eta2, phi2,  mass2, decayMode2));
	   
	   double metcorr_ex_Up, metcorr_ey_Up, metcorr_ex_Down, metcorr_ey_Down;
	   metcorr_ex_Up = metcorr_shifted(metcorr_ex, 
					   pt1, phi1, isDM0_1, tesUncertainties(NtupleVer, decayMode),
					   pt2, phi2, isDM0_2, tesUncertainties(NtupleVer, decayMode2),
					   +1, +1);
	   metcorr_ey_Up = metcorr_shifted(metcorr_ey, 
					   pt1, phi1, isDM0_1, tesUncertainties(NtupleVer, decayMode),
					   pt2, phi2, isDM0_2, tesUncertainties(NtupleVer, decayMode2),
					   -1, +1);
	   metcorr_ex_Down = metcorr_shifted(metcorr_ex, 
					     pt1, phi1, isDM0_1, tesUncertainties(NtupleVer, decayMode),
					     pt2, phi2, isDM0_2, tesUncertainties(NtupleVer, decayMode2),
					     +1, -1);
	   metcorr_ey_Down = metcorr_shifted(metcorr_ey, 
					     pt1, phi1, isDM0_1, tesUncertainties(NtupleVer, decayMode),
					     pt2, phi2, isDM0_2, tesUncertainties(NtupleVer, decayMode2),
					     -1, -1);
	   
	   runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, 
		    svFitMass_DM0_Up, svFitPt_DM0_Up, svFitEta_DM0_Up, svFitPhi_DM0_Up, svFitMET_DM0_Up, svFitTransverseMass_DM0_Up, tau1_DM0_Up, tau2_DM0_Up);
	   
	   runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, 
		    svFitMass_DM0_Down, svFitPt_DM0_Down, svFitEta_DM0_Down, svFitPhi_DM0_Down, svFitMET_DM0_Down, svFitTransverseMass_DM0_Down, tau1_DM0_Down, tau2_DM0_Down);
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

           svFitMass_DM0_Down=svFitMass;
           svFitPt_DM0_Down=svFitPt;
           svFitEta_DM0_Down=svFitEta;
           svFitPhi_DM0_Down=svFitPhi;
           svFitMET_DM0_Down=svFitMET;
           svFitTransverseMass_DM0_Down=svFitTransverseMass;
           tau1_DM0_Down = tau1;
           tau2_DM0_Down = tau2;
	 }



	 //***********************************//
	 //****  Tau DM1 shifted up/down  ****//
	 //***********************************//
	 if (isDM1_1 || isDM1_2) {
	   std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
	   std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
	   double pt1_Up, pt2_Up, pt1_Down, pt2_Down;
	   pt1_Up = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode), isDM1_1, +1);
	   pt2_Up = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2), isDM1_2, +1);
	   pt1_Down = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode), isDM1_1, -1);
	   pt2_Down = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2), isDM1_2, -1);
	   
	   measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
	   measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2, pt2_Up, eta2, phi2,  mass2, decayMode2));	   
	   measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
	   measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2, pt2_Down, eta2, phi2,  mass2, decayMode2));
	   
	   double metcorr_ex_Up, metcorr_ey_Up, metcorr_ex_Down, metcorr_ey_Down;
	   metcorr_ex_Up = metcorr_shifted(metcorr_ex, 
					   pt1, phi1, isDM1_1, tesUncertainties(NtupleVer, decayMode),
					   pt2, phi2, isDM1_2, tesUncertainties(NtupleVer, decayMode2),
					   +1, +1);
	   metcorr_ey_Up = metcorr_shifted(metcorr_ey, 
					   pt1, phi1, isDM1_1, tesUncertainties(NtupleVer, decayMode),
					   pt2, phi2, isDM1_2, tesUncertainties(NtupleVer, decayMode2),
					   -1, +1);
	   metcorr_ex_Down = metcorr_shifted(metcorr_ex, 
					     pt1, phi1, isDM1_1, tesUncertainties(NtupleVer, decayMode),
					     pt2, phi2, isDM1_2, tesUncertainties(NtupleVer, decayMode2),
					     +1, -1);
	   metcorr_ey_Down = metcorr_shifted(metcorr_ey, 
					     pt1, phi1, isDM1_1, tesUncertainties(NtupleVer, decayMode),
					     pt2, phi2, isDM1_2, tesUncertainties(NtupleVer, decayMode2),
					     -1, -1);
	   
	   runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, 
		    svFitMass_DM1_Up, svFitPt_DM1_Up, svFitEta_DM1_Up, svFitPhi_DM1_Up, svFitMET_DM1_Up, svFitTransverseMass_DM1_Up, tau1_DM1_Up, tau2_DM1_Up);
	   
	   runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, 
		    svFitMass_DM1_Down, svFitPt_DM1_Down, svFitEta_DM1_Down, svFitPhi_DM1_Down, svFitMET_DM1_Down, svFitTransverseMass_DM1_Down, tau1_DM1_Down, tau2_DM1_Down);
	 }
	 else {
	   svFitMass_DM1_Up=svFitMass;
           svFitPt_DM1_Up=svFitPt;
           svFitEta_DM1_Up=svFitEta;
           svFitPhi_DM1_Up=svFitPhi;
           svFitMET_DM1_Up=svFitMET;
           svFitTransverseMass_DM1_Up=svFitTransverseMass;
           tau1_DM1_Up = tau1;
           tau2_DM1_Down = tau2;

           svFitMass_DM1_Down=svFitMass;
           svFitPt_DM1_Down=svFitPt;
           svFitEta_DM1_Down=svFitEta;
           svFitPhi_DM1_Down=svFitPhi;
           svFitMET_DM1_Down=svFitMET;
           svFitTransverseMass_DM1_Down=svFitTransverseMass;
           tau1_DM1_Down = tau1;
           tau2_DM1_Down = tau2;
	 }




	 //***********************************//
	 //****  Tau DM10 shifted up/down  ****//
	 //***********************************//
	 if (isDM10_1 || isDM10_2) {
	   std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp;
	   std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown;
	   double pt1_Up, pt2_Up, pt1_Down, pt2_Down;
	   pt1_Up = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode), isDM10_1, +1);
	   pt2_Up = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2), isDM10_2, +1);
	   pt1_Down = pt_shifted(pt1, tesUncertainties(NtupleVer, decayMode), isDM10_1, -1);
	   pt2_Down = pt_shifted(pt2, tesUncertainties(NtupleVer, decayMode2), isDM10_2, -1);
	   
	   measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Up, eta1,  phi1, mass1, decayMode));
	   measuredTauLeptonsUp.push_back(classic_svFit::MeasuredTauLepton(decayType2, pt2_Up, eta2, phi2,  mass2, decayMode2));	   
	   measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType1, pt1_Down, eta1,  phi1, mass1, decayMode));
	   measuredTauLeptonsDown.push_back(classic_svFit::MeasuredTauLepton(decayType2, pt2_Down, eta2, phi2,  mass2, decayMode2));
	   
	   double metcorr_ex_Up, metcorr_ey_Up, metcorr_ex_Down, metcorr_ey_Down;
	   metcorr_ex_Up = metcorr_shifted(metcorr_ex, 
					   pt1, phi1, isDM10_1, tesUncertainties(NtupleVer, decayMode),
					   pt2, phi2, isDM10_2, tesUncertainties(NtupleVer, decayMode2),
					   +1, +1);
	   metcorr_ey_Up = metcorr_shifted(metcorr_ey, 
					   pt1, phi1, isDM10_1, tesUncertainties(NtupleVer, decayMode),
					   pt2, phi2, isDM10_2, tesUncertainties(NtupleVer, decayMode2),
					   -1, +1);
	   metcorr_ex_Down = metcorr_shifted(metcorr_ex, 
					     pt1, phi1, isDM10_1, tesUncertainties(NtupleVer, decayMode),
					     pt2, phi2, isDM10_2, tesUncertainties(NtupleVer, decayMode2),
					     +1, -1);
	   metcorr_ey_Down = metcorr_shifted(metcorr_ey, 
					     pt1, phi1, isDM10_1, tesUncertainties(NtupleVer, decayMode),
					     pt2, phi2, isDM10_2, tesUncertainties(NtupleVer, decayMode2),
					     -1, -1);
	   
	   runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, 
		    svFitMass_DM10_Up, svFitPt_DM10_Up, svFitEta_DM10_Up, svFitPhi_DM10_Up, svFitMET_DM10_Up, svFitTransverseMass_DM10_Up, tau1_DM10_Up, tau2_DM10_Up);
	   
	   runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, 
		    svFitMass_DM10_Down, svFitPt_DM10_Down, svFitEta_DM10_Down, svFitPhi_DM10_Down, svFitMET_DM10_Down, svFitTransverseMass_DM10_Down, tau1_DM10_Down, tau2_DM10_Down);
	 }
	 else {
	   svFitMass_DM10_Up=svFitMass;
           svFitPt_DM10_Up=svFitPt;
           svFitEta_DM10_Up=svFitEta;
           svFitPhi_DM10_Up=svFitPhi;
           svFitMET_DM10_Up=svFitMET;
           svFitTransverseMass_DM10_Up=svFitTransverseMass;
           tau1_DM10_Up = tau1;
           tau2_DM10_Down = tau2;

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

     // _JetEC2_ up
     tau1_pt_JetEC2_Up  = tau1_JetEC2_Up.Pt();
     tau1_eta_JetEC2_Up = tau1_JetEC2_Up.Eta();
     tau1_phi_JetEC2_Up = tau1_JetEC2_Up.Phi();
     tau1_m_JetEC2_Up   = tau1_JetEC2_Up.M();
     tau2_pt_JetEC2_Up  = tau2_JetEC2_Up.Pt();
     tau2_eta_JetEC2_Up = tau2_JetEC2_Up.Eta();
     tau2_phi_JetEC2_Up = tau2_JetEC2_Up.Phi();
     tau2_m_JetEC2_Up   = tau2_JetEC2_Up.M();
     // down
     tau1_pt_JetEC2_Down  = tau1_JetEC2_Down.Pt();
     tau1_eta_JetEC2_Down = tau1_JetEC2_Down.Eta();
     tau1_phi_JetEC2_Down = tau1_JetEC2_Down.Phi();
     tau1_m_JetEC2_Down   = tau1_JetEC2_Down.M();
     tau2_pt_JetEC2_Down  = tau2_JetEC2_Down.Pt();
     tau2_eta_JetEC2_Down = tau2_JetEC2_Down.Eta();
     tau2_phi_JetEC2_Down = tau2_JetEC2_Down.Phi();
     tau2_m_JetEC2_Down   = tau2_JetEC2_Down.M();

     // _JetEta0to3_ up
     tau1_pt_JetEta0to3_Up  = tau1_JetEta0to3_Up.Pt();
     tau1_eta_JetEta0to3_Up = tau1_JetEta0to3_Up.Eta();
     tau1_phi_JetEta0to3_Up = tau1_JetEta0to3_Up.Phi();
     tau1_m_JetEta0to3_Up   = tau1_JetEta0to3_Up.M();
     tau2_pt_JetEta0to3_Up  = tau2_JetEta0to3_Up.Pt();
     tau2_eta_JetEta0to3_Up = tau2_JetEta0to3_Up.Eta();
     tau2_phi_JetEta0to3_Up = tau2_JetEta0to3_Up.Phi();
     tau2_m_JetEta0to3_Up   = tau2_JetEta0to3_Up.M();
     // down
     tau1_pt_JetEta0to3_Down  = tau1_JetEta0to3_Down.Pt();
     tau1_eta_JetEta0to3_Down = tau1_JetEta0to3_Down.Eta();
     tau1_phi_JetEta0to3_Down = tau1_JetEta0to3_Down.Phi();
     tau1_m_JetEta0to3_Down   = tau1_JetEta0to3_Down.M();
     tau2_pt_JetEta0to3_Down  = tau2_JetEta0to3_Down.Pt();
     tau2_eta_JetEta0to3_Down = tau2_JetEta0to3_Down.Eta();
     tau2_phi_JetEta0to3_Down = tau2_JetEta0to3_Down.Phi();
     tau2_m_JetEta0to3_Down   = tau2_JetEta0to3_Down.M();

     // _JetEta0to5_ up
     tau1_pt_JetEta0to5_Up  = tau1_JetEta0to5_Up.Pt();
     tau1_eta_JetEta0to5_Up = tau1_JetEta0to5_Up.Eta();
     tau1_phi_JetEta0to5_Up = tau1_JetEta0to5_Up.Phi();
     tau1_m_JetEta0to5_Up   = tau1_JetEta0to5_Up.M();
     tau2_pt_JetEta0to5_Up  = tau2_JetEta0to5_Up.Pt();
     tau2_eta_JetEta0to5_Up = tau2_JetEta0to5_Up.Eta();
     tau2_phi_JetEta0to5_Up = tau2_JetEta0to5_Up.Phi();
     tau2_m_JetEta0to5_Up   = tau2_JetEta0to5_Up.M();
     // down
     tau1_pt_JetEta0to5_Down  = tau1_JetEta0to5_Down.Pt();
     tau1_eta_JetEta0to5_Down = tau1_JetEta0to5_Down.Eta();
     tau1_phi_JetEta0to5_Down = tau1_JetEta0to5_Down.Phi();
     tau1_m_JetEta0to5_Down   = tau1_JetEta0to5_Down.M();
     tau2_pt_JetEta0to5_Down  = tau2_JetEta0to5_Down.Pt();
     tau2_eta_JetEta0to5_Down = tau2_JetEta0to5_Down.Eta();
     tau2_phi_JetEta0to5_Down = tau2_JetEta0to5_Down.Phi();
     tau2_m_JetEta0to5_Down   = tau2_JetEta0to5_Down.M();

     // _JetEta3to5_ up
     tau1_pt_JetEta3to5_Up  = tau1_JetEta3to5_Up.Pt();
     tau1_eta_JetEta3to5_Up = tau1_JetEta3to5_Up.Eta();
     tau1_phi_JetEta3to5_Up = tau1_JetEta3to5_Up.Phi();
     tau1_m_JetEta3to5_Up   = tau1_JetEta3to5_Up.M();
     tau2_pt_JetEta3to5_Up  = tau2_JetEta3to5_Up.Pt();
     tau2_eta_JetEta3to5_Up = tau2_JetEta3to5_Up.Eta();
     tau2_phi_JetEta3to5_Up = tau2_JetEta3to5_Up.Phi();
     tau2_m_JetEta3to5_Up   = tau2_JetEta3to5_Up.M();
     // down
     tau1_pt_JetEta3to5_Down  = tau1_JetEta3to5_Down.Pt();
     tau1_eta_JetEta3to5_Down = tau1_JetEta3to5_Down.Eta();
     tau1_phi_JetEta3to5_Down = tau1_JetEta3to5_Down.Phi();
     tau1_m_JetEta3to5_Down   = tau1_JetEta3to5_Down.M();
     tau2_pt_JetEta3to5_Down  = tau2_JetEta3to5_Down.Pt();
     tau2_eta_JetEta3to5_Down = tau2_JetEta3to5_Down.Eta();
     tau2_phi_JetEta3to5_Down = tau2_JetEta3to5_Down.Phi();
     tau2_m_JetEta3to5_Down   = tau2_JetEta3to5_Down.M();

     // _JetRelativeBal_ up
     tau1_pt_JetRelativeBal_Up  = tau1_JetRelativeBal_Up.Pt();
     tau1_eta_JetRelativeBal_Up = tau1_JetRelativeBal_Up.Eta();
     tau1_phi_JetRelativeBal_Up = tau1_JetRelativeBal_Up.Phi();
     tau1_m_JetRelativeBal_Up   = tau1_JetRelativeBal_Up.M();
     tau2_pt_JetRelativeBal_Up  = tau2_JetRelativeBal_Up.Pt();
     tau2_eta_JetRelativeBal_Up = tau2_JetRelativeBal_Up.Eta();
     tau2_phi_JetRelativeBal_Up = tau2_JetRelativeBal_Up.Phi();
     tau2_m_JetRelativeBal_Up   = tau2_JetRelativeBal_Up.M();
     // down
     tau1_pt_JetRelativeBal_Down  = tau1_JetRelativeBal_Down.Pt();
     tau1_eta_JetRelativeBal_Down = tau1_JetRelativeBal_Down.Eta();
     tau1_phi_JetRelativeBal_Down = tau1_JetRelativeBal_Down.Phi();
     tau1_m_JetRelativeBal_Down   = tau1_JetRelativeBal_Down.M();
     tau2_pt_JetRelativeBal_Down  = tau2_JetRelativeBal_Down.Pt();
     tau2_eta_JetRelativeBal_Down = tau2_JetRelativeBal_Down.Eta();
     tau2_phi_JetRelativeBal_Down = tau2_JetRelativeBal_Down.Phi();
     tau2_m_JetRelativeBal_Down   = tau2_JetRelativeBal_Down.M();

     // _JetRelativeSample_ up
     tau1_pt_JetRelativeSample_Up  = tau1_JetRelativeSample_Up.Pt();
     tau1_eta_JetRelativeSample_Up = tau1_JetRelativeSample_Up.Eta();
     tau1_phi_JetRelativeSample_Up = tau1_JetRelativeSample_Up.Phi();
     tau1_m_JetRelativeSample_Up   = tau1_JetRelativeSample_Up.M();
     tau2_pt_JetRelativeSample_Up  = tau2_JetRelativeSample_Up.Pt();
     tau2_eta_JetRelativeSample_Up = tau2_JetRelativeSample_Up.Eta();
     tau2_phi_JetRelativeSample_Up = tau2_JetRelativeSample_Up.Phi();
     tau2_m_JetRelativeSample_Up   = tau2_JetRelativeSample_Up.M();
     // down
     tau1_pt_JetRelativeSample_Down  = tau1_JetRelativeSample_Down.Pt();
     tau1_eta_JetRelativeSample_Down = tau1_JetRelativeSample_Down.Eta();
     tau1_phi_JetRelativeSample_Down = tau1_JetRelativeSample_Down.Phi();
     tau1_m_JetRelativeSample_Down   = tau1_JetRelativeSample_Down.M();
     tau2_pt_JetRelativeSample_Down  = tau2_JetRelativeSample_Down.Pt();
     tau2_eta_JetRelativeSample_Down = tau2_JetRelativeSample_Down.Eta();
     tau2_phi_JetRelativeSample_Down = tau2_JetRelativeSample_Down.Phi();
     tau2_m_JetRelativeSample_Down   = tau2_JetRelativeSample_Down.M();
     
     
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
     newBranch88->Fill();
     newBranch89->Fill();
     newBranch90->Fill();
     newBranch91->Fill();
     newBranch92->Fill();
     newBranch93->Fill();
     newBranch94->Fill();
     newBranch95->Fill();
     newBranch96->Fill();
     newBranch97->Fill();
     newBranch98->Fill();
     newBranch99->Fill();
     newBranch100->Fill();
     newBranch101->Fill();
     newBranch102->Fill();
     newBranch103->Fill();
     newBranch104->Fill();
     newBranch105->Fill();
     newBranch106->Fill();
     newBranch107->Fill();
     newBranch108->Fill();
     newBranch109->Fill();
     newBranch110->Fill();
     newBranch111->Fill();
     newBranch112->Fill();
     newBranch113->Fill();
     newBranch114->Fill();
     newBranch115->Fill();
     newBranch116->Fill();
     newBranch117->Fill();
     newBranch118->Fill();
     newBranch119->Fill();
     newBranch120->Fill();
     newBranch121->Fill();
     newBranch122->Fill();
     newBranch123->Fill();
     newBranch124->Fill();
     newBranch125->Fill();
     newBranch126->Fill();
     newBranch127->Fill();
     newBranch128->Fill();
     newBranch129->Fill();
     newBranch130->Fill();
     newBranch131->Fill();
     newBranch132->Fill();
     newBranch133->Fill();
     newBranch134->Fill();
     newBranch135->Fill();
     newBranch136->Fill();
     newBranch137->Fill();
     newBranch138->Fill();
     newBranch139->Fill();
     newBranch140->Fill();
     newBranch141->Fill();
     newBranch142->Fill();
     newBranch143->Fill();
     newBranch144->Fill();
     newBranch145->Fill();
     newBranch146->Fill();
     newBranch147->Fill();
     newBranch148->Fill();
     newBranch149->Fill();
     newBranch150->Fill();
     newBranch151->Fill();
     newBranch152->Fill();
     newBranch153->Fill();
     newBranch154->Fill();
     newBranch155->Fill();
     newBranch156->Fill();
     newBranch157->Fill();
     newBranch158->Fill();
     newBranch159->Fill();
     newBranch160->Fill();
     newBranch161->Fill();
     newBranch162->Fill();
     newBranch163->Fill();
     newBranch164->Fill();
     newBranch165->Fill();
     newBranch166->Fill();
     newBranch167->Fill();
     newBranch168->Fill();
     newBranch169->Fill();
     newBranch170->Fill();
     newBranch171->Fill();
     newBranch172->Fill();
     newBranch173->Fill();
     newBranch174->Fill();
     newBranch175->Fill();
     newBranch176->Fill();
     newBranch177->Fill();
     newBranch178->Fill();
     newBranch179->Fill();
     newBranch180->Fill();
     newBranch181->Fill();
     newBranch182->Fill();
     newBranch183->Fill();
     newBranch184->Fill();
     newBranch185->Fill();
     newBranch186->Fill();


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
