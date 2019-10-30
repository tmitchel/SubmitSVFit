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

      float metcorrRecoilResoUp_ex = -10;
      float metcorrRecoilResoUp_ey = -10;
      float metcorrRecoilResoDown_ex = -10;
      float metcorrRecoilResoDown_ey = -10;
      float metcorrRecoilRespUp_ex = -10;
      float metcorrRecoilRespUp_ey = -10;
      float metcorrRecoilRespDown_ex = -10;
      float metcorrRecoilRespDown_ey = -10;

      float metcorrJetEC2UpMET_ex = -10;
      float metcorrJetEC2UpMET_ey = -10;
      float metcorrJetEC2DownMET_ex = -10;
      float metcorrJetEC2DownMET_ey = -10;
      float metcorrJetEta0to3UpMET_ex = -10;
      float metcorrJetEta0to3UpMET_ey = -10;
      float metcorrJetEta0to3DownMET_ex = -10;
      float metcorrJetEta0to3DownMET_ey = -10;
      float metcorrJetEta0to5UpMET_ex = -10;
      float metcorrJetEta0to5UpMET_ey = -10;
      float metcorrJetEta0to5DownMET_ex = -10;
      float metcorrJetEta0to5DownMET_ey = -10;
      float metcorrJetEta3to5UpMET_ex = -10;
      float metcorrJetEta3to5UpMET_ey = -10;
      float metcorrJetEta3to5DownMET_ex = -10;
      float metcorrJetEta3to5DownMET_ey = -10;
      float metcorrJetRelativeBalUpMET_ex = -10;
      float metcorrJetRelativeBalUpMET_ey = -10;
      float metcorrJetRelativeBalDownMET_ex = -10;
      float metcorrJetRelativeBalDownMET_ey = -10;
      float metcorrJetRelativeSampleUpMET_ex = -10;
      float metcorrJetRelativeSampleUpMET_ey = -10;
      float metcorrJetRelativeSampleDownMET_ex = -10;
      float metcorrJetRelativeSampleDownMET_ey = -10;

      TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
      
      // If doing ES shifts, we need extra ouput branches
      //
      float svFitMass_Up = -10;
      float svFitPt_Up = -10;
      float svFitMass_Down = -10;
      float svFitPt_Down = -10;
      
      float svFitMass_DM0_Up = -10;
      float svFitPt_DM0_Up = -10;
      float svFitMass_DM0_Down = -10;
      float svFitPt_DM0_Down = -10;
      
      float svFitMass_DM1_Up = -10;
      float svFitPt_DM1_Up = -10;
      float svFitMass_DM1_Down = -10;
      float svFitPt_DM1_Down = -10;

      float svFitMass_DM10_Up = -10;
      float svFitPt_DM10_Up = -10;
      float svFitMass_DM10_Down = -10;
      float svFitPt_DM10_Down = -10;

      float svFitMass_LES_DM0_Up = -10;
      float svFitPt_LES_DM0_Up = -10;
      float svFitMass_LES_DM0_Down = -10;
      float svFitPt_LES_DM0_Down = -10;

      float svFitMass_LES_DM1_Up = -10;
      float svFitPt_LES_DM1_Up = -10;
      float svFitMass_LES_DM1_Down = -10;
      float svFitPt_LES_DM1_Down = -10;

      float svFitMass_EEScale_Up = -10;
      float svFitPt_EEScale_Up = -10;
      float svFitMass_EEScale_Down = -10;
      float svFitPt_EEScale_Down = -10;

      float svFitMass_EESigma_Up = -10;
      float svFitPt_EESigma_Up = -10;
      float svFitMass_EESigma_Down = -10;
      float svFitPt_EESigma_Down = -10;

      float svFitMass_MESbin1_Up = -10;
      float svFitPt_MESbin1_Up = -10;
      float svFitMass_MESbin1_Down = -10;
      float svFitPt_MESbin1_Down = -10;

      float svFitMass_MESbin2_Up = -10;
      float svFitPt_MESbin2_Up = -10;
      float svFitMass_MESbin2_Down = -10;
      float svFitPt_MESbin2_Down = -10;

      float svFitMass_MESbin3_Up = -10;
      float svFitPt_MESbin3_Up = -10;
      float svFitMass_MESbin3_Down = -10;
      float svFitPt_MESbin3_Down = -10;

      float svFitMass_MESbin4_Up = -10;
      float svFitPt_MESbin4_Up = -10;
      float svFitMass_MESbin4_Down = -10;
      float svFitPt_MESbin4_Down = -10;

      float svFitMass_MESbin5_Up = -10;
      float svFitPt_MESbin5_Up = -10;
      float svFitMass_MESbin5_Down = -10;
      float svFitPt_MESbin5_Down = -10;

      float svFitMass_UncMet_Up = -10;
      float svFitPt_UncMet_Up = -10;
      float svFitMass_UncMet_Down = -10;
      float svFitPt_UncMet_Down = -10;

      float svFitMass_ClusteredMet_Up = -10;
      float svFitPt_ClusteredMet_Up = -10;
      float svFitMass_ClusteredMet_Down = -10;
      float svFitPt_ClusteredMet_Down = -10;
      
      float svFitMass_JetEC2_Up = -10;
      float svFitPt_JetEC2_Up = -10;
      float svFitMass_JetEC2_Down = -10;
      float svFitPt_JetEC2_Down = -10;

      float svFitMass_JetEta0to3_Up = -10;
      float svFitPt_JetEta0to3_Up = -10;
      float svFitMass_JetEta0to3_Down = -10;
      float svFitPt_JetEta0to3_Down = -10;

      float svFitMass_JetEta0to5_Up = -10;
      float svFitPt_JetEta0to5_Up = -10;
      float svFitMass_JetEta0to5_Down = -10;
      float svFitPt_JetEta0to5_Down = -10;

      float svFitMass_JetEta3to5_Up = -10;
      float svFitPt_JetEta3to5_Up = -10;
      float svFitMass_JetEta3to5_Down = -10;
      float svFitPt_JetEta3to5_Down = -10;

      float svFitMass_JetRelativeBal_Up = -10;
      float svFitPt_JetRelativeBal_Up = -10;
      float svFitMass_JetRelativeBal_Down = -10;
      float svFitPt_JetRelativeBal_Down = -10;

      float svFitMass_JetRelativeSample_Up = -10;
      float svFitPt_JetRelativeSample_Up = -10;
      float svFitMass_JetRelativeSample_Down = -10;
      float svFitPt_JetRelativeSample_Down = -10;

      float svFitMass_RecoilReso_Up = -10;
      float svFitPt_RecoilReso_Up = -10;
      float svFitMass_RecoilReso_Down = -10;
      float svFitPt_RecoilReso_Down = -10;

      float svFitMass_RecoilResp_Up = -10;
      float svFitPt_RecoilResp_Up = -10;
      float svFitMass_RecoilResp_Down = -10;
      float svFitPt_RecoilResp_Down = -10;
	

      // tau leptons                                                                                                  

      TBranch *newDMBranch1 = t->Branch("m_sv_Up", &svFitMass_Up, "m_sv_Up/F");
      TBranch *newDMBranch2 = t->Branch("pt_sv_Up", &svFitPt_Up, "pt_sv_Up/F");
      TBranch *newDMBranch3 = t->Branch("m_sv_Down", &svFitMass_Down, "m_sv_Down/F");
      TBranch *newDMBranch4 = t->Branch("pt_sv_Down", &svFitPt_Down, "pt_sv_Down/F");

      TBranch *newDM0Branch1 = t->Branch("m_sv_DM0_Up", &svFitMass_DM0_Up, "m_sv_DM0_Up/F");
      TBranch *newDM0Branch2 = t->Branch("pt_sv_DM0_Up", &svFitPt_DM0_Up, "pt_sv_DM0_Up/F");
      TBranch *newDM0Branch3 = t->Branch("m_sv_DM0_Down", &svFitMass_DM0_Down, "m_sv_DM0_Down/F");
      TBranch *newDM0Branch4 = t->Branch("pt_sv_DM0_Down", &svFitPt_DM0_Down, "pt_sv_DM0_Down/F");

      TBranch *newEEScaleBranch1 = t->Branch("m_sv_EEScale_Up", &svFitMass_EEScale_Up, "m_sv_EEScale_Up/F");
      TBranch *newEEScaleBranch2 = t->Branch("pt_sv_EEScale_Up", &svFitPt_EEScale_Up, "pt_sv_EEScale_Up/F");
      TBranch *newEEScaleBranch3 = t->Branch("m_sv_EEScale_Down", &svFitMass_EEScale_Down, "m_sv_EEScale_Down/F");
      TBranch *newEEScaleBranch4 = t->Branch("pt_sv_EEScale_Down", &svFitPt_EEScale_Down, "pt_sv_EEScale_Down/F");

      TBranch *newEESigmaBranch1 = t->Branch("m_sv_EESigma_Up", &svFitMass_EESigma_Up, "m_sv_EESigma_Up/F");
      TBranch *newEESigmaBranch2 = t->Branch("pt_sv_EESigma_Up", &svFitPt_EESigma_Up, "pt_sv_EESigma_Up/F");
      TBranch *newEESigmaBranch3 = t->Branch("m_sv_EESigma_Down", &svFitMass_EESigma_Down, "m_sv_EESigma_Down/F");
      TBranch *newEESigmaBranch4 = t->Branch("pt_sv_EESigma_Down", &svFitPt_EESigma_Down, "pt_sv_EESigma_Down/F");

      TBranch *newMESbin1Branch1 = t->Branch("m_sv_MESbin1_Up", &svFitMass_MESbin1_Up, "m_sv_MESbin1_Up/F");
      TBranch *newMESbin1Branch2 = t->Branch("pt_sv_MESbin1_Up", &svFitPt_MESbin1_Up, "pt_sv_MESbin1_Up/F");
      TBranch *newMESbin1Branch3 = t->Branch("m_sv_MESbin1_Down", &svFitMass_MESbin1_Down, "m_sv_MESbin1_Down/F");
      TBranch *newMESbin1Branch4 = t->Branch("pt_sv_MESbin1_Down", &svFitPt_MESbin1_Down, "pt_sv_MESbin1_Down/F");

      TBranch *newMESbin2Branch1 = t->Branch("m_sv_MESbin2_Up", &svFitMass_MESbin2_Up, "m_sv_MESbin2_Up/F");
      TBranch *newMESbin2Branch2 = t->Branch("pt_sv_MESbin2_Up", &svFitPt_MESbin2_Up, "pt_sv_MESbin2_Up/F");
      TBranch *newMESbin2Branch3 = t->Branch("m_sv_MESbin2_Down", &svFitMass_MESbin2_Down, "m_sv_MESbin2_Down/F");
      TBranch *newMESbin2Branch4 = t->Branch("pt_sv_MESbin2_Down", &svFitPt_MESbin2_Down, "pt_sv_MESbin2_Down/F");

      TBranch *newMESbin3Branch1 = t->Branch("m_sv_MESbin3_Up", &svFitMass_MESbin3_Up, "m_sv_MESbin3_Up/F");
      TBranch *newMESbin3Branch2 = t->Branch("pt_sv_MESbin3_Up", &svFitPt_MESbin3_Up, "pt_sv_MESbin3_Up/F");
      TBranch *newMESbin3Branch3 = t->Branch("m_sv_MESbin3_Down", &svFitMass_MESbin3_Down, "m_sv_MESbin3_Down/F");
      TBranch *newMESbin3Branch4 = t->Branch("pt_sv_MESbin3_Down", &svFitPt_MESbin3_Down, "pt_sv_MESbin3_Down/F");

      TBranch *newMESbin4Branch1 = t->Branch("m_sv_MESbin4_Up", &svFitMass_MESbin4_Up, "m_sv_MESbin4_Up/F");
      TBranch *newMESbin4Branch2 = t->Branch("pt_sv_MESbin4_Up", &svFitPt_MESbin4_Up, "pt_sv_MESbin4_Up/F");
      TBranch *newMESbin4Branch3 = t->Branch("m_sv_MESbin4_Down", &svFitMass_MESbin4_Down, "m_sv_MESbin4_Down/F");
      TBranch *newMESbin4Branch4 = t->Branch("pt_sv_MESbin4_Down", &svFitPt_MESbin4_Down, "pt_sv_MESbin4_Down/F");

      TBranch *newMESbin5Branch1 = t->Branch("m_sv_MESbin5_Up", &svFitMass_MESbin5_Up, "m_sv_MESbin5_Up/F");
      TBranch *newMESbin5Branch2 = t->Branch("pt_sv_MESbin5_Up", &svFitPt_MESbin5_Up, "pt_sv_MESbin5_Up/F");
      TBranch *newMESbin5Branch3 = t->Branch("m_sv_MESbin5_Down", &svFitMass_MESbin5_Down, "m_sv_MESbin5_Down/F");
      TBranch *newMESbin5Branch4 = t->Branch("pt_sv_MESbin5_Down", &svFitPt_MESbin5_Down, "pt_sv_MESbin5_Down/F");

      TBranch *newDM1Branch1 = t->Branch("m_sv_DM1_Up", &svFitMass_DM1_Up, "m_sv_DM1_Up/F");
      TBranch *newDM1Branch2 = t->Branch("pt_sv_DM1_Up", &svFitPt_DM1_Up, "pt_sv_DM1_Up/F");
      TBranch *newDM1Branch3 = t->Branch("m_sv_DM1_Down", &svFitMass_DM1_Down, "m_sv_DM1_Down/F");
      TBranch *newDM1Branch4 = t->Branch("pt_sv_DM1_Down", &svFitPt_DM1_Down, "pt_sv_DM1_Down/F");

      TBranch *newDM10Branch1 = t->Branch("m_sv_DM10_Up", &svFitMass_DM10_Up, "m_sv_DM10_Up/F");
      TBranch *newDM10Branch2 = t->Branch("pt_sv_DM10_Up", &svFitPt_DM10_Up, "pt_sv_DM10_Up/F");
      TBranch *newDM10Branch3 = t->Branch("m_sv_DM10_Down", &svFitMass_DM10_Down, "m_sv_DM10_Down/F");
      TBranch *newDM10Branch4 = t->Branch("pt_sv_DM10_Down", &svFitPt_DM10_Down, "pt_sv_DM10_Down/F");

      TBranch *newLESDM0Branch1 = t->Branch("m_sv_LES_DM0_Up", &svFitMass_LES_DM0_Up, "m_sv_LES_DM0_Up/F");
      TBranch *newLESDM0Branch2 = t->Branch("pt_sv_LES_DM0_Up", &svFitPt_LES_DM0_Up, "pt_sv_LES_DM0_Up/F");
      TBranch *newLESDM0Branch3 = t->Branch("m_sv_LES_DM0_Down", &svFitMass_LES_DM0_Down, "m_sv_LES_DM0_Down/F");
      TBranch *newLESDM0Branch4 = t->Branch("pt_sv_LES_DM0_Down", &svFitPt_LES_DM0_Down, "pt_sv_LES_DM0_Down/F");

      TBranch *newLESDM1Branch1 = t->Branch("m_sv_LES_DM1_Up", &svFitMass_LES_DM1_Up, "m_sv_LES_DM1_Up/F");
      TBranch *newLESDM1Branch2 = t->Branch("pt_sv_LES_DM1_Up", &svFitPt_LES_DM1_Up, "pt_sv_LES_DM1_Up/F");
      TBranch *newLESDM1Branch3 = t->Branch("m_sv_LES_DM1_Down", &svFitMass_LES_DM1_Down, "m_sv_LES_DM1_Down/F");
      TBranch *newLESDM1Branch4 = t->Branch("pt_sv_LES_DM1_Down", &svFitPt_LES_DM1_Down, "pt_sv_LES_DM1_Down/F");

      TBranch *newUncMetBranch1 = t->Branch("m_sv_UncMet_Up", &svFitMass_UncMet_Up, "m_sv_UncMet_Up/F");
      TBranch *newUncMetBranch2 = t->Branch("pt_sv_UncMet_Up", &svFitPt_UncMet_Up, "pt_sv_UncMet_Up/F");
      TBranch *newUncMetBranch3 = t->Branch("m_sv_UncMet_Down", &svFitMass_UncMet_Down, "m_sv_UncMet_Down/F");
      TBranch *newUncMetBranch4 = t->Branch("pt_sv_UncMet_Down", &svFitPt_UncMet_Down, "pt_sv_UncMet_Down/F");

      TBranch *newClusteredMetBranch1 = t->Branch("m_sv_ClusteredMet_Up", &svFitMass_ClusteredMet_Up, "m_sv_ClusteredMet_Up/F");
      TBranch *newClusteredMetBranch2 = t->Branch("pt_sv_ClusteredMet_Up", &svFitPt_ClusteredMet_Up, "pt_sv_ClusteredMet_Up/F");
      TBranch *newClusteredMetBranch3 = t->Branch("m_sv_ClusteredMet_Down", &svFitMass_ClusteredMet_Down, "m_sv_ClusteredMet_Down/F");
      TBranch *newClusteredMetBranch4 = t->Branch("pt_sv_ClusteredMet_Down", &svFitPt_ClusteredMet_Down, "pt_sv_ClusteredMet_Down/F");

      TBranch *newJetEC2Branch1 = t->Branch("m_sv_JetEC2_Up", &svFitMass_JetEC2_Up, "m_sv_JetEC2_Up/F");
      TBranch *newJetEC2Branch2 = t->Branch("pt_sv_JetEC2_Up", &svFitPt_JetEC2_Up, "pt_sv_JetEC2_Up/F");
      TBranch *newJetEC2Branch3 = t->Branch("m_sv_JetEC2_Down", &svFitMass_JetEC2_Down, "m_sv_JetEC2_Down/F");
      TBranch *newJetEC2Branch4 = t->Branch("pt_sv_JetEC2_Down", &svFitPt_JetEC2_Down, "pt_sv_JetEC2_Down/F");

      TBranch *newJetEta0to3Branch1 = t->Branch("m_sv_JetEta0to3_Up", &svFitMass_JetEta0to3_Up, "m_sv_JetEta0to3_Up/F");
      TBranch *newJetEta0to3Branch2 = t->Branch("pt_sv_JetEta0to3_Up", &svFitPt_JetEta0to3_Up, "pt_sv_JetEta0to3_Up/F");
      TBranch *newJetEta0to3Branch3 = t->Branch("m_sv_JetEta0to3_Down", &svFitMass_JetEta0to3_Down, "m_sv_JetEta0to3_Down/F");
      TBranch *newJetEta0to3Branch4 = t->Branch("pt_sv_JetEta0to3_Down", &svFitPt_JetEta0to3_Down, "pt_sv_JetEta0to3_Down/F");

      TBranch *newJetEta0to5Branch1 = t->Branch("m_sv_JetEta0to5_Up", &svFitMass_JetEta0to5_Up, "m_sv_JetEta0to5_Up/F");
      TBranch *newJetEta0to5Branch2 = t->Branch("pt_sv_JetEta0to5_Up", &svFitPt_JetEta0to5_Up, "pt_sv_JetEta0to5_Up/F");
      TBranch *newJetEta0to5Branch3 = t->Branch("m_sv_JetEta0to5_Down", &svFitMass_JetEta0to5_Down, "m_sv_JetEta0to5_Down/F");
      TBranch *newJetEta0to5Branch4 = t->Branch("pt_sv_JetEta0to5_Down", &svFitPt_JetEta0to5_Down, "pt_sv_JetEta0to5_Down/F");

      TBranch *newJetEta3to5Branch1 = t->Branch("m_sv_JetEta3to5_Up", &svFitMass_JetEta3to5_Up, "m_sv_JetEta3to5_Up/F");
      TBranch *newJetEta3to5Branch2 = t->Branch("pt_sv_JetEta3to5_Up", &svFitPt_JetEta3to5_Up, "pt_sv_JetEta3to5_Up/F");
      TBranch *newJetEta3to5Branch3 = t->Branch("m_sv_JetEta3to5_Down", &svFitMass_JetEta3to5_Down, "m_sv_JetEta3to5_Down/F");
      TBranch *newJetEta3to5Branch4 = t->Branch("pt_sv_JetEta3to5_Down", &svFitPt_JetEta3to5_Down, "pt_sv_JetEta3to5_Down/F");

      TBranch *newJetRelativeBalBranch1 = t->Branch("m_sv_JetRelativeBal_Up", &svFitMass_JetRelativeBal_Up, "m_sv_JetRelativeBal_Up/F");
      TBranch *newJetRelativeBalBranch2 = t->Branch("pt_sv_JetRelativeBal_Up", &svFitPt_JetRelativeBal_Up, "pt_sv_JetRelativeBal_Up/F");
      TBranch *newJetRelativeBalBranch3 = t->Branch("m_sv_JetRelativeBal_Down", &svFitMass_JetRelativeBal_Down, "m_sv_JetRelativeBal_Down/F");
      TBranch *newJetRelativeBalBranch4 = t->Branch("pt_sv_JetRelativeBal_Down", &svFitPt_JetRelativeBal_Down, "pt_sv_JetRelativeBal_Down/F");

      TBranch *newJetRelativeSampleBranch1 = t->Branch("m_sv_JetRelativeSample_Up", &svFitMass_JetRelativeSample_Up, "m_sv_JetRelativeSample_Up/F");
      TBranch *newJetRelativeSampleBranch2 = t->Branch("pt_sv_JetRelativeSample_Up", &svFitPt_JetRelativeSample_Up, "pt_sv_JetRelativeSample_Up/F");
      TBranch *newJetRelativeSampleBranch3 = t->Branch("m_sv_JetRelativeSample_Down", &svFitMass_JetRelativeSample_Down, "m_sv_JetRelativeSample_Down/F");
      TBranch *newJetRelativeSampleBranch4 = t->Branch("pt_sv_JetRelativeSample_Down", &svFitPt_JetRelativeSample_Down, "pt_sv_JetRelativeSample_Down/F");

      TBranch *newRecoilResoBranch1 = t->Branch("m_sv_RecoilReso_Up", &svFitMass_RecoilReso_Up, "m_sv_RecoilReso_Up/F");
      TBranch *newRecoilResoBranch2 = t->Branch("pt_sv_RecoilReso_Up", &svFitPt_RecoilReso_Up, "pt_sv_RecoilReso_Up/F");
      TBranch *newRecoilResoBranch3 = t->Branch("m_sv_RecoilReso_Down", &svFitMass_RecoilReso_Down, "m_sv_RecoilReso_Down/F");
      TBranch *newRecoilResoBranch4 = t->Branch("pt_sv_RecoilReso_Down", &svFitPt_RecoilReso_Down, "pt_sv_RecoilReso_Down/F");

      TBranch *newRecoilRespBranch1 = t->Branch("m_sv_RecoilResp_Up", &svFitMass_RecoilResp_Up, "m_sv_RecoilResp_Up/F");
      TBranch *newRecoilRespBranch2 = t->Branch("pt_sv_RecoilResp_Up", &svFitPt_RecoilResp_Up, "pt_sv_RecoilResp_Up/F");
      TBranch *newRecoilRespBranch3 = t->Branch("m_sv_RecoilResp_Down", &svFitMass_RecoilResp_Down, "m_sv_RecoilResp_Down/F");
      TBranch *newRecoilRespBranch4 = t->Branch("pt_sv_RecoilResp_Down", &svFitPt_RecoilResp_Down, "pt_sv_RecoilResp_Down/F");
    
      // adding tau-related branches
      std::vector<TBranch*> tau4VectorBranches;
    
      Int_t era;
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

      // JES Uncertainties
      float met_JetEC2Up, met_JetEC2Down, metphi_JetEC2Up, metphi_JetEC2Down;
      float met_JetEta0to3Up, met_JetEta0to3Down, metphi_JetEta0to3Up, metphi_JetEta0to3Down;
      float met_JetEta0to5Up, met_JetEta0to5Down, metphi_JetEta0to5Up, metphi_JetEta0to5Down;
      float met_JetEta3to5Up, met_JetEta3to5Down, metphi_JetEta3to5Up, metphi_JetEta3to5Down;
      float met_JetRelativeBalUp, met_JetRelativeBalDown, metphi_JetRelativeBalUp, metphi_JetRelativeBalDown;
      float met_JetRelativeSampleUp, met_JetRelativeSampleDown, metphi_JetRelativeSampleUp, metphi_JetRelativeSampleDown;
      float met_reso_Up, met_reso_Down, met_resp_Up, met_resp_Down;
      float metphi_reso_Up, metphi_reso_Down, metphi_resp_Up, metphi_resp_Down;

      float eCorrectedEt = 0., eEnergyScaleUp = 0., eEnergyScaleDown = 0., eEnergySigmaUp = 0., eEnergySigmaDown = 0.;

      //ele/mu variables
      TBranch *pt1branch;
       
      t->SetBranchAddress("era",&era);
      t->SetBranchAddress("NtupleVer",&NtupleVer);
      t->SetBranchAddress("evt",&evt);
      t->SetBranchAddress("run",&run);
      t->SetBranchAddress("lumi",&lumi);
      t->SetBranchAddress("gen_match_1",&gen_match_1);
      t->SetBranchAddress("gen_match_2",&gen_match_2);
      if ( channel == "tt" ) t->SetBranchAddress("t1_decayMode", &decayMode);
      if ( channel == "tt" ) t->SetBranchAddress("t2_decayMode", &decayMode2);
      if ( channel == "et") {
        t->SetBranchAddress("eCorrectedEt", &eCorrectedEt);
        t->SetBranchAddress("eEnergyScaleUp", &eEnergyScaleUp);
        t->SetBranchAddress("eEnergyScaleDown", &eEnergyScaleDown);
        t->SetBranchAddress("eEnergySigmaUp", &eEnergySigmaUp);
        t->SetBranchAddress("eEnergySigmaDown", &eEnergySigmaDown);
      }
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

      t->SetBranchAddress("met_reso_Up", &met_reso_Up);
      t->SetBranchAddress("met_reso_Down", &met_reso_Down);
      t->SetBranchAddress("metphi_reso_Up", &metphi_reso_Up);
      t->SetBranchAddress("metphi_reso_Down", &metphi_reso_Down);

      t->SetBranchAddress("met_resp_Up", &met_resp_Up);
      t->SetBranchAddress("met_resp_Down", &met_resp_Down);
      t->SetBranchAddress("metphi_resp_Up", &metphi_resp_Up);
      t->SetBranchAddress("metphi_resp_Down", &metphi_resp_Down);


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
          
          covMET[0][0] =  mvaCovMatrix00;
          covMET[1][0] =  mvaCovMatrix10;
          covMET[0][1] =  mvaCovMatrix01;
          covMET[1][1] =  mvaCovMatrix11;
        } // mva met
        if (metType == -1) { // -1 = PF Met
          TMet.SetPtEtaPhiM(pfmet,0,pfmetphi,0);
          metcorr_ex = pfmet*TMath::Cos(pfmetphi);
          metcorr_ey = pfmet*TMath::Sin(pfmetphi);
          // Shifted METs
          
          metcorrUncUp_ex         = uncMetPtUp*TMath::Cos(uncMetPhiUp);
          metcorrUncUp_ey         = uncMetPtUp*TMath::Sin(uncMetPhiUp);
          metcorrUncDown_ex       = uncMetPtDown*TMath::Cos(uncMetPhiDown);
          metcorrUncDown_ey       = uncMetPtDown*TMath::Sin(uncMetPhiDown);
          metcorrClusteredUp_ex   = clusteredMetPtUp*TMath::Cos(clusteredMetPhiUp);
          metcorrClusteredUp_ey   = clusteredMetPtUp*TMath::Sin(clusteredMetPhiUp);
          metcorrClusteredDown_ex = clusteredMetPtDown*TMath::Cos(clusteredMetPhiDown);
          metcorrClusteredDown_ey = clusteredMetPtDown*TMath::Sin(clusteredMetPhiDown);
          metcorrJetEC2UpMET_ex         = met_JetEC2Up*TMath::Cos(metphi_JetEC2Up);
          metcorrJetEC2UpMET_ey         = met_JetEC2Up*TMath::Sin(metphi_JetEC2Up);
          metcorrJetEC2DownMET_ex       = met_JetEC2Down*TMath::Cos(metphi_JetEC2Down);
          metcorrJetEC2DownMET_ey       = met_JetEC2Down*TMath::Sin(metphi_JetEC2Down);
          metcorrJetEta0to3UpMET_ex         = met_JetEta0to3Up*TMath::Cos(metphi_JetEta0to3Up);
          metcorrJetEta0to3UpMET_ey         = met_JetEta0to3Up*TMath::Sin(metphi_JetEta0to3Up);
          metcorrJetEta0to3DownMET_ex       = met_JetEta0to3Down*TMath::Cos(metphi_JetEta0to3Down);
          metcorrJetEta0to3DownMET_ey       = met_JetEta0to3Down*TMath::Sin(metphi_JetEta0to3Down);
          metcorrJetEta0to5UpMET_ex         = met_JetEta0to5Up*TMath::Cos(metphi_JetEta0to5Up);
          metcorrJetEta0to5UpMET_ey         = met_JetEta0to5Up*TMath::Sin(metphi_JetEta0to5Up);
          metcorrJetEta0to5DownMET_ex       = met_JetEta0to5Down*TMath::Cos(metphi_JetEta0to5Down);
          metcorrJetEta0to5DownMET_ey       = met_JetEta0to5Down*TMath::Sin(metphi_JetEta0to5Down);
          metcorrJetEta3to5UpMET_ex         = met_JetEta3to5Up*TMath::Cos(metphi_JetEta3to5Up);
          metcorrJetEta3to5UpMET_ey         = met_JetEta3to5Up*TMath::Sin(metphi_JetEta3to5Up);
          metcorrJetEta3to5DownMET_ex       = met_JetEta3to5Down*TMath::Cos(metphi_JetEta3to5Down);
          metcorrJetEta3to5DownMET_ey       = met_JetEta3to5Down*TMath::Sin(metphi_JetEta3to5Down);
          metcorrJetRelativeBalUpMET_ex         = met_JetRelativeBalUp*TMath::Cos(metphi_JetRelativeBalUp);
          metcorrJetRelativeBalUpMET_ey         = met_JetRelativeBalUp*TMath::Sin(metphi_JetRelativeBalUp);
          metcorrJetRelativeBalDownMET_ex       = met_JetRelativeBalDown*TMath::Cos(metphi_JetRelativeBalDown);
          metcorrJetRelativeBalDownMET_ey       = met_JetRelativeBalDown*TMath::Sin(metphi_JetRelativeBalDown);
          metcorrJetRelativeSampleUpMET_ex         = met_JetRelativeSampleUp*TMath::Cos(metphi_JetRelativeSampleUp);
          metcorrJetRelativeSampleUpMET_ey         = met_JetRelativeSampleUp*TMath::Sin(metphi_JetRelativeSampleUp);
          metcorrJetRelativeSampleDownMET_ex       = met_JetRelativeSampleDown*TMath::Cos(metphi_JetRelativeSampleDown);
          metcorrJetRelativeSampleDownMET_ey       = met_JetRelativeSampleDown*TMath::Sin(metphi_JetRelativeSampleDown);
          metcorrRecoilResoUp_ex = met_reso_Up*TMath::Cos(metphi_reso_Up);
          metcorrRecoilResoUp_ey = met_reso_Down*TMath::Cos(metphi_reso_Down);
          metcorrRecoilResoDown_ex = met_reso_Up*TMath::Cos(metphi_reso_Up);
          metcorrRecoilResoDown_ey = met_reso_Down*TMath::Cos(metphi_reso_Down);
          metcorrRecoilRespUp_ex = met_resp_Up*TMath::Cos(metphi_resp_Up);
          metcorrRecoilRespUp_ey = met_resp_Down*TMath::Cos(metphi_resp_Down);
          metcorrRecoilRespDown_ex = met_resp_Up*TMath::Cos(metphi_resp_Up);
          metcorrRecoilRespDown_ey = met_resp_Down*TMath::Cos(metphi_resp_Down);

          covMET[0][0] =  pfCovMatrix00;
          covMET[1][0] =  pfCovMatrix10;
          covMET[0][1] =  pfCovMatrix01;
          covMET[1][1] =  pfCovMatrix11;
        } // pf met


        if (channel == "mt" || channel == "et") {
          mass2 = m2;
          std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons{
              classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1),
              classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2)
          };

          std::cout << "era: " << era << "evt: " << evt << " run: " << run << " lumi: " << lumi << " pt1 " << pt1
                    << " mass1 " << mass1 << " pt2: " << pt2 << " mass2: " << mass2 << std::endl;

          runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0,
               svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          // MET systematics
          runSVFit(measuredTauLeptons, metcorrUncUp_ex, metcorrUncUp_ey, covMET, 0, svFitMass_UncMet_Up, svFitPt_UncMet_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrUncDown_ex, metcorrUncDown_ey, covMET, 0, svFitMass_UncMet_Down, svFitPt_UncMet_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrClusteredUp_ex, metcorrClusteredUp_ey, covMET, 0, svFitMass_ClusteredMet_Up,
              svFitPt_ClusteredMet_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrClusteredDown_ex, metcorrClusteredDown_ey, covMET, 0, svFitMass_ClusteredMet_Down,
              svFitPt_ClusteredMet_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEC2UpMET_ex, metcorrJetEC2UpMET_ey, covMET, 0, svFitMass_JetEC2_Up, svFitPt_JetEC2_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEC2DownMET_ex, metcorrJetEC2DownMET_ey, covMET, 0, svFitMass_JetEC2_Down, svFitPt_JetEC2_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEta0to3UpMET_ex, metcorrJetEta0to3UpMET_ey, covMET, 0, svFitMass_JetEta0to3_Up, svFitPt_JetEta0to3_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEta0to3DownMET_ex, metcorrJetEta0to3DownMET_ey, covMET, 0, svFitMass_JetEta0to3_Down,
              svFitPt_JetEta0to3_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEta0to5UpMET_ex, metcorrJetEta0to5UpMET_ey, covMET, 0, svFitMass_JetEta0to5_Up,
              svFitPt_JetEta0to5_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEta0to5DownMET_ex, metcorrJetEta0to5DownMET_ey, covMET, 0, svFitMass_JetEta0to5_Down,
              svFitPt_JetEta0to5_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEta3to5UpMET_ex, metcorrJetEta3to5UpMET_ey, covMET, 0, svFitMass_JetEta3to5_Up, svFitPt_JetEta3to5_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetEta3to5DownMET_ex, metcorrJetEta3to5DownMET_ey, covMET, 0, svFitMass_JetEta3to5_Down,
              svFitPt_JetEta3to5_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetRelativeBalUpMET_ex, metcorrJetRelativeBalUpMET_ey, covMET, 0, svFitMass_JetRelativeBal_Up,
              svFitPt_JetRelativeBal_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetRelativeBalDownMET_ex, metcorrJetRelativeBalDownMET_ey, covMET, 0, svFitMass_JetRelativeBal_Down,
              svFitPt_JetRelativeBal_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetRelativeSampleUpMET_ex, metcorrJetRelativeSampleUpMET_ey, covMET, 0, svFitMass_JetRelativeSample_Up,
              svFitPt_JetRelativeSample_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJetRelativeSampleDownMET_ex, metcorrJetRelativeSampleDownMET_ey, covMET, 0, svFitMass_JetRelativeSample_Down,
              svFitPt_JetRelativeSample_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRecoilResoUp_ex, metcorrRecoilResoUp_ey, covMET, 0, svFitMass_RecoilReso_Up,
              svFitPt_RecoilReso_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRecoilResoDown_ex, metcorrRecoilResoDown_ey, covMET, 0, svFitMass_RecoilReso_Down,
              svFitPt_RecoilReso_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRecoilRespUp_ex, metcorrRecoilRespUp_ey, covMET, 0, svFitMass_RecoilResp_Up,
              svFitPt_RecoilResp_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRecoilRespDown_ex, metcorrRecoilRespDown_ey, covMET, 0, svFitMass_RecoilResp_Down,
              svFitPt_RecoilResp_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          if (doES) {
            // corrections only need to be done once
            float ES_Up(1.), ES_Down(1.);  // shift TES
            if (gen_match_2 == 5) {  // 0.6% uncertainty on hadronic tau
              ES_Up = 1 + tesUncertainties(era, decayMode2);
              ES_Down = 1 - tesUncertainties(era, decayMode2);
            } else if (gen_match_2 < 5) {  // flat 2% on el/mu -> tau energy scale systematics
              ES_Up = 1.02;
              ES_Down = 0.98;
            }

            if (channel == "et") {
              TLorentzVector orig_el;
              orig_el.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
              TLorentzVector scale_up = orig_el * (eEnergyScaleUp / eCorrectedEt);
              TLorentzVector scale_dn = orig_el * (eEnergyScaleDown / eCorrectedEt);
              TLorentzVector sigma_up = orig_el * (eEnergySigmaUp / eCorrectedEt);
              TLorentzVector sigma_dn = orig_el * (eEnergySigmaDown / eCorrectedEt);
              std::vector<classic_svFit::MeasuredTauLepton> measuredTauScaleUp{
                classic_svFit::MeasuredTauLepton(decayType1, scale_up.Pt(), scale_up.Eta(), scale_up.Phi(), scale_up.M()),
                classic_svFit::MeasuredTauLepton(decayType2, pt2, eta2, phi2, mass2, decayMode2)
              };
              runSVFit(measuredTauScaleUp, metcorr_ex, metcorr_ey, covMET, 0,
               svFitMass_EEScale_Up, svFitPt_EEScale_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

              std::vector<classic_svFit::MeasuredTauLepton> measuredTauScaleDn{
                classic_svFit::MeasuredTauLepton(decayType1, scale_dn.Pt(), scale_dn.Eta(), scale_dn.Phi(), scale_dn.M()),
                classic_svFit::MeasuredTauLepton(decayType2, pt2, eta2, phi2, mass2, decayMode2)
              };
              runSVFit(measuredTauScaleDn, metcorr_ex, metcorr_ey, covMET, 0,
               svFitMass_EEScale_Down, svFitPt_EEScale_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

              std::vector<classic_svFit::MeasuredTauLepton> measuredTauSigmaUp{
                classic_svFit::MeasuredTauLepton(decayType1, sigma_up.Pt(), sigma_up.Eta(), sigma_up.Phi(), sigma_up.M()),
                classic_svFit::MeasuredTauLepton(decayType2, pt2, eta2, phi2, mass2, decayMode2)
              };
              runSVFit(measuredTauSigmaUp, metcorr_ex, metcorr_ey, covMET, 0,
               svFitMass_EESigma_Up, svFitPt_EESigma_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

              std::vector<classic_svFit::MeasuredTauLepton> measuredTauSigmaDn{
                classic_svFit::MeasuredTauLepton(decayType1, sigma_dn.Pt(), sigma_dn.Eta(), sigma_dn.Phi(), sigma_dn.M()),
                classic_svFit::MeasuredTauLepton(decayType2, pt2, eta2, phi2, mass2, decayMode2)
              };
              runSVFit(measuredTauSigmaDn, metcorr_ex, metcorr_ey, covMET, 0,
               svFitMass_EESigma_Down, svFitPt_EESigma_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else if (channel == "mt") {
              TLorentzVector orig_mu;
              orig_mu.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
              TLorentzVector scaled_mu_up = orig_mu;
              TLorentzVector scaled_mu_dn = orig_mu;
              if (orig_mu.Eta() > -2.4 && orig_mu.Eta() < -2.1) {
                scaled_mu_up *= 1 + 0.027;
                scaled_mu_dn *= 1 - 0.027;
              } else if (orig_mu.Eta() > -2.1 && orig_mu.Eta() < -1.2)  {
                scaled_mu_up *= 1 + 0.009;
                scaled_mu_dn *= 1 - 0.009;
              } else if (orig_mu.Eta() > -1.2 && orig_mu.Eta() < 1.2)  {
                scaled_mu_up *= 1 + 0.004;
                scaled_mu_dn *= 1 - 0.004;
              } else if (orig_mu.Eta() > 1.2 && orig_mu.Eta() < 2.1)  {
                scaled_mu_up *= 1 + 0.009;
                scaled_mu_dn *= 1 - 0.009;
              } else if (orig_mu.Eta() > 2.1 && orig_mu.Eta() < 2.4)  {
                scaled_mu_up *= 1 + 0.017;
                scaled_mu_dn *= 1 - 0.017;
              }

              std::vector<classic_svFit::MeasuredTauLepton> measuredTauUp{
                classic_svFit::MeasuredTauLepton(decayType1, scaled_mu_up.Pt(), scaled_mu_up.Eta(), scaled_mu_up.Phi(), scaled_mu_up.M()),
                classic_svFit::MeasuredTauLepton(decayType2, pt2, eta2, phi2, mass2, decayMode2)
              };
              std::vector<classic_svFit::MeasuredTauLepton> measuredTauDn{
                classic_svFit::MeasuredTauLepton(decayType1, scaled_mu_dn.Pt(), scaled_mu_dn.Eta(), scaled_mu_dn.Phi(), scaled_mu_dn.M()),
                classic_svFit::MeasuredTauLepton(decayType2, pt2, eta2, phi2, mass2, decayMode2)
              };

              if (orig_mu.Eta() > -2.4 && orig_mu.Eta() < -2.1) {
                runSVFit(measuredTauUp, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin1_Up, svFitPt_MESbin1_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
                runSVFit(measuredTauDn, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin1_Down, svFitPt_MESbin1_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
              } else {
                svFitMass_MESbin1_Up = svFitMass;
                svFitPt_MESbin1_Up = svFitPt;
                svFitMass_MESbin1_Down = svFitMass;
                svFitPt_MESbin1_Down = svFitPt;
              }
              if (orig_mu.Eta() > -2.1 && orig_mu.Eta() < -1.2) {
                runSVFit(measuredTauUp, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin2_Up, svFitPt_MESbin2_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
                runSVFit(measuredTauDn, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin2_Down, svFitPt_MESbin2_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
              } else {
                svFitMass_MESbin2_Up = svFitMass;
                svFitPt_MESbin2_Up = svFitPt;
                svFitMass_MESbin2_Down = svFitMass;
                svFitPt_MESbin2_Down = svFitPt;
              }
              if (orig_mu.Eta() > -1.2 && orig_mu.Eta() < 1.2) {
                runSVFit(measuredTauUp, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin3_Up, svFitPt_MESbin3_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
                runSVFit(measuredTauDn, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin3_Down, svFitPt_MESbin3_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
              } else {
                svFitMass_MESbin3_Up = svFitMass;
                svFitPt_MESbin3_Up = svFitPt;
                svFitMass_MESbin3_Down = svFitMass;
                svFitPt_MESbin3_Down = svFitPt;
              }
              if (orig_mu.Eta() > 1.2 && orig_mu.Eta() < 2.1) {
                runSVFit(measuredTauUp, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin4_Up, svFitPt_MESbin4_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
                runSVFit(measuredTauDn, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin4_Down, svFitPt_MESbin4_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
              } else {
                svFitMass_MESbin4_Up = svFitMass;
                svFitPt_MESbin4_Up = svFitPt;
                svFitMass_MESbin4_Down = svFitMass;
                svFitPt_MESbin4_Down = svFitPt;
              }
              if (orig_mu.Eta() > 2.1 && orig_mu.Eta() < 2.4) {
                runSVFit(measuredTauUp, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin5_Up, svFitPt_MESbin5_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
                runSVFit(measuredTauDn, metcorr_ex, metcorr_ey, covMET, 0,
                  svFitMass_MESbin5_Down, svFitPt_MESbin5_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
              } else {
                svFitMass_MESbin5_Up = svFitMass;
                svFitPt_MESbin5_Up = svFitPt;
                svFitMass_MESbin5_Down = svFitMass;
                svFitPt_MESbin5_Down = svFitPt;
              }
            }

            double pt_Up(pt2 * ES_Up), pt_Down(pt2 * ES_Down);  // shift tau pT by energy scale
            double dx_Up(pt2 * TMath::Cos(phi2) * ((1. / ES_Up) - 1.)),
                   dy_Up(pt2 * TMath::Sin(phi2) * ((1. / ES_Up) - 1.)),
                   dx_Down(pt2 * TMath::Cos(phi2) * ((1. / ES_Down) - 1.)),
                   dy_Down(pt2 * TMath::Sin(phi2) * ((1. / ES_Down) - 1.));
            double metcorr_ex_Up(metcorr_ex + dx_Up),
                   metcorr_ey_Up(metcorr_ey + dy_Up),
                   metcorr_ex_Down(metcorr_ex + dx_Down),
                   metcorr_ey_Down(metcorr_ey + dy_Down);

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
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_Up, svFitPt_Up,
                  svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_Up = svFitMass;
              svFitPt_Up = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // tau DM0 shifted up
            if (gen_match_2 == 5 && decayMode2 == 0) {
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM0_Up, svFitPt_DM0_Up,
                  svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

            } else {
              svFitMass_DM0_Up = svFitMass;
              svFitPt_DM0_Up = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // tau DM1 shifted up
            if (gen_match_2 == 5 && decayMode2 == 1) {
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM1_Up, svFitPt_DM1_Up,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_DM1_Up = svFitMass;
              svFitPt_DM1_Up = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // tau DM10 shifted up
            if (gen_match_2 == 5 && decayMode2 == 10) {
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM10_Up, svFitPt_DM10_Up,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_DM10_Up = svFitMass;
              svFitPt_DM10_Up = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // lep->tau DM0 shifted up
            if (gen_match_2 < 5 && decayMode2 == 0) {
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_LES_DM0_Up, svFitPt_LES_DM0_Up,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

            } else {
              svFitMass_LES_DM0_Up = svFitMass;
              svFitPt_LES_DM0_Up = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // lep->tau DM1 shifted up
            if (gen_match_2 < 5 && decayMode2 == 1) {
              runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_LES_DM1_Up, svFitPt_LES_DM1_Up,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_LES_DM1_Up = svFitMass;
              svFitPt_LES_DM1_Up = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            ///////////////////////////////
            // All downward shifts below //
            ///////////////////////////////

            // all tau shift down
            if (gen_match_2 < 6) {
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_Down, svFitPt_Down,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_Down = svFitMass;
              svFitPt_Down = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // tau DM0 shifted down
            if (gen_match_2 == 5 && decayMode2 == 0) {
              std::cout << "DM0 shift Down" << std::endl;
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM0_Down, svFitPt_DM0_Down,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_DM0_Down = svFitMass;
              svFitPt_DM0_Down = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // tau DM1 shifted down
            if (gen_match_2 == 5 && decayMode2 == 1) {
              std::cout << "DM1 shift down" << std::endl;
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM1_Down, svFitPt_DM1_Down,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_DM1_Down = svFitMass;
              svFitPt_DM1_Down = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // tau DM10 shifted down
            if (gen_match_2 == 5 && decayMode2 == 10) {
              std::cout << "DM10 shift down" << std::endl;
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM10_Down, svFitPt_DM10_Down,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_DM10_Down = svFitMass;
              svFitPt_DM10_Down = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // lep->tau DM0 shifted down
            if (gen_match_2 < 5 && decayMode2 == 0) {
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_LES_DM0_Down, svFitPt_LES_DM0_Down,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_LES_DM0_Down = svFitMass;
              svFitPt_LES_DM0_Down = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }

            // lep->tau DM1 shifted down
            if (gen_match_2 < 5 && decayMode2 == 1) {
              runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_LES_DM1_Down, svFitPt_LES_DM1_Down,
                    svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
            } else {
              svFitMass_LES_DM1_Down = svFitMass;
              svFitPt_LES_DM1_Down = svFitPt;
              tau1 = tau1;
              tau2 = tau2;
            }
          }  // end doES
        }  // eTau / muTau     
     else {
       svFitMass = -100;
       svFitPt = -100;
     }
     
     std::cout << "\n\n" << std::endl;
     //std::cout << "\n\nex: " << metcorr_ex << "   ey: " << metcorr_ey <<  " phi: " << metcorphi<<"\n"<<std::endl; 
    newBranch1->Fill();
    newBranch2->Fill();
    newDMBranch1->Fill();
    newDMBranch2->Fill();
    newDMBranch3->Fill();
    newDMBranch4->Fill();
    newDM0Branch1->Fill();
    newDM0Branch2->Fill();
    newDM0Branch3->Fill();
    newDM0Branch4->Fill();
    newEEScaleBranch1->Fill();
    newEEScaleBranch2->Fill();
    newEEScaleBranch3->Fill();
    newEEScaleBranch4->Fill();
    newEESigmaBranch1->Fill();
    newEESigmaBranch2->Fill();
    newEESigmaBranch3->Fill();
    newEESigmaBranch4->Fill();
    newMESbin1Branch1->Fill();
    newMESbin1Branch2->Fill();
    newMESbin1Branch3->Fill();
    newMESbin1Branch4->Fill();
    newMESbin2Branch1->Fill();
    newMESbin2Branch2->Fill();
    newMESbin2Branch3->Fill();
    newMESbin2Branch4->Fill();
    newMESbin3Branch1->Fill();
    newMESbin3Branch2->Fill();
    newMESbin3Branch3->Fill();
    newMESbin3Branch4->Fill();
    newMESbin4Branch1->Fill();
    newMESbin4Branch2->Fill();
    newMESbin4Branch3->Fill();
    newMESbin4Branch4->Fill();
    newMESbin5Branch1->Fill();
    newMESbin5Branch2->Fill();
    newMESbin5Branch3->Fill();
    newMESbin5Branch4->Fill();
    newDM1Branch1->Fill();
    newDM1Branch2->Fill();
    newDM1Branch3->Fill();
    newDM1Branch4->Fill();
    newDM10Branch1->Fill();
    newDM10Branch2->Fill();
    newDM10Branch3->Fill();
    newDM10Branch4->Fill();
    newLESDM0Branch1->Fill();
    newLESDM0Branch2->Fill();
    newLESDM0Branch3->Fill();
    newLESDM0Branch4->Fill();
    newLESDM1Branch1->Fill();
    newLESDM1Branch2->Fill();
    newLESDM1Branch3->Fill();
    newLESDM1Branch4->Fill();
    newUncMetBranch1->Fill();
    newUncMetBranch2->Fill();
    newUncMetBranch3->Fill();
    newUncMetBranch4->Fill();
    newClusteredMetBranch1->Fill();
    newClusteredMetBranch2->Fill();
    newClusteredMetBranch3->Fill();
    newClusteredMetBranch4->Fill();
    newJetEC2Branch1->Fill();
    newJetEC2Branch2->Fill();
    newJetEC2Branch3->Fill();
    newJetEC2Branch4->Fill();
    newJetEta0to3Branch1->Fill();
    newJetEta0to3Branch2->Fill();
    newJetEta0to3Branch3->Fill();
    newJetEta0to3Branch4->Fill();
    newJetEta0to5Branch1->Fill();
    newJetEta0to5Branch2->Fill();
    newJetEta0to5Branch3->Fill();
    newJetEta0to5Branch4->Fill();
    newJetEta3to5Branch1->Fill();
    newJetEta3to5Branch2->Fill();
    newJetEta3to5Branch3->Fill();
    newJetEta3to5Branch4->Fill();
    newJetRelativeBalBranch1->Fill();
    newJetRelativeBalBranch2->Fill();
    newJetRelativeBalBranch3->Fill();
    newJetRelativeBalBranch4->Fill();
    newJetRelativeSampleBranch1->Fill();
    newJetRelativeSampleBranch2->Fill();
    newJetRelativeSampleBranch3->Fill();
    newJetRelativeSampleBranch4->Fill();
    newRecoilResoBranch1->Fill();
    newRecoilResoBranch2->Fill();
    newRecoilResoBranch3->Fill();
    newRecoilResoBranch4->Fill();
    newRecoilRespBranch1->Fill();
    newRecoilRespBranch2->Fill();
    newRecoilRespBranch3->Fill();
    newRecoilRespBranch4->Fill();
     
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
