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
#include "SubmitSVFit/ROOT/interface/TauFESTool.h"

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
  parser.addOption("idyear", optutl::CommandLineParser::kString, "year", "2018");
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

      // get TES/FES weighter
      TauFESTool* tfes = new TauFESTool(parser.stringValue("year"), "DeepTau2017v2p1VSe", "SubmitSVFit/ROOT/data");
      
      TTree *t = (TTree*)obj;
      float svFitMass = -10;
      float svFitPt = -10;
      float svFitEta = -10;
      float svFitPhi = -10;
      float svFitMET = -10;
      float svFitTransverseMass = -10;
      
      float metcorr_ex = -10; // corrected met px (float)
      float metcorr_ey = -10;  // corrected met py (float)

      // recoil
      float metcorrRecoilResoUp_ex = -10;
      float metcorrRecoilResoUp_ey = -10;
      float metcorrRecoilResoDown_ex = -10;
      float metcorrRecoilResoDown_ey = -10;
      float metcorrRecoilRespUp_ex = -10;
      float metcorrRecoilRespUp_ey = -10;
      float metcorrRecoilRespDown_ex = -10;
      float metcorrRecoilRespDown_ey = -10;

      // JEC/JER
      float metcorrJERDown_ex = -10;
      float metcorrJERUp_ex = -10;
      float metcorrAbsoluteDown_ex = -10;
      float metcorrAbsoluteUp_ex = -10;
      float metcorrAbsoluteyearDown_ex = -10;
      float metcorrAbsoluteyearUp_ex = -10;
      float metcorrBBEC1Down_ex = -10;
      float metcorrBBEC1Up_ex = -10;
      float metcorrBBEC1yearDown_ex = -10;
      float metcorrBBEC1yearUp_ex = -10;
      float metcorrEC2Down_ex = -10;
      float metcorrEC2Up_ex = -10;
      float metcorrEC2yearDown_ex = -10;
      float metcorrEC2yearUp_ex = -10;
      float metcorrFlavorQCDDown_ex = -10;
      float metcorrFlavorQCDUp_ex = -10;
      float metcorrHFDown_ex = -10;
      float metcorrHFUp_ex = -10;
      float metcorrHFyearDown_ex = -10;
      float metcorrHFyearUp_ex = -10;
      float metcorrRelBalDown_ex = -10;
      float metcorrRelBalUp_ex = -10;
      float metcorrRelSamDown_ex = -10;
      float metcorrRelSamUp_ex = -10;
      float metcorrResDown_ex = -10;
      float metcorrResUp_ex = -10;
      float metcorrUESDown_ex = -10;
      float metcorrUESUp_ex = -10;

      float metcorrJERDown_ey = -10;
      float metcorrJERUp_ey = -10;
      float metcorrAbsoluteDown_ey = -10;
      float metcorrAbsoluteUp_ey = -10;
      float metcorrAbsoluteyearDown_ey = -10;
      float metcorrAbsoluteyearUp_ey = -10;
      float metcorrBBEC1Down_ey = -10;
      float metcorrBBEC1Up_ey = -10;
      float metcorrBBEC1yearDown_ey = -10;
      float metcorrBBEC1yearUp_ey = -10;
      float metcorrEC2Down_ey = -10;
      float metcorrEC2Up_ey = -10;
      float metcorrEC2yearDown_ey = -10;
      float metcorrEC2yearUp_ey = -10;
      float metcorrFlavorQCDDown_ey = -10;
      float metcorrFlavorQCDUp_ey = -10;
      float metcorrHFDown_ey = -10;
      float metcorrHFUp_ey = -10;
      float metcorrHFyearDown_ey = -10;
      float metcorrHFyearUp_ey = -10;
      float metcorrRelBalDown_ey = -10;
      float metcorrRelBalUp_ey = -10;
      float metcorrRelSamDown_ey = -10;
      float metcorrRelSamUp_ey = -10;
      float metcorrResDown_ey = -10;
      float metcorrResUp_ey = -10;
      float metcorrUESDown_ey = -10;
      float metcorrUESUp_ey = -10;


      TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
      
      // If doing ES shifts, we need extra ouput branches
      
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

      float svFitMass_MES_Up = -10;
      float svFitPt_MES_Up = -10;
      float svFitMass_MES_Down = -10;
      float svFitPt_MES_Down = -10;

      float svFitMass_JER_Down = -10;
      float svFitPt_JER_Down = -10;
      float svFitMass_JER_Up = -10;
      float svFitPt_JER_Up = -10;

      float svFitMass_Absolute_Down = -10;
      float svFitPt_Absolute_Down = -10;
      float svFitMass_Absolute_Up = -10;
      float svFitPt_Absolute_Up = -10;

      float svFitMass_Absoluteyear_Down = -10;
      float svFitPt_Absoluteyear_Down = -10;
      float svFitMass_Absoluteyear_Up = -10;
      float svFitPt_Absoluteyear_Up = -10;

      float svFitMass_BBEC1_Down = -10;
      float svFitPt_BBEC1_Down = -10;
      float svFitMass_BBEC1_Up = -10;
      float svFitPt_BBEC1_Up = -10;

      float svFitMass_BBEC1year_Down = -10;
      float svFitPt_BBEC1year_Down = -10;
      float svFitMass_BBEC1year_Up = -10;
      float svFitPt_BBEC1year_Up = -10;

      float svFitMass_EC2_Down = -10;
      float svFitPt_EC2_Down = -10;
      float svFitMass_EC2_Up = -10;
      float svFitPt_EC2_Up = -10;

      float svFitMass_EC2year_Down = -10;
      float svFitPt_EC2year_Down = -10;
      float svFitMass_EC2year_Up = -10;
      float svFitPt_EC2year_Up = -10;

      float svFitMass_FlavorQCD_Down = -10;
      float svFitPt_FlavorQCD_Down = -10;
      float svFitMass_FlavorQCD_Up = -10;
      float svFitPt_FlavorQCD_Up = -10;

      float svFitMass_HF_Down = -10;
      float svFitPt_HF_Down = -10;
      float svFitMass_HF_Up = -10;
      float svFitPt_HF_Up = -10;

      float svFitMass_HFyear_Down = -10;
      float svFitPt_HFyear_Down = -10;
      float svFitMass_HFyear_Up = -10;
      float svFitPt_HFyear_Up = -10;

      float svFitMass_RelBal_Down = -10;
      float svFitPt_RelBal_Down = -10;
      float svFitMass_RelBal_Up = -10;
      float svFitPt_RelBal_Up = -10;

      float svFitMass_RelSam_Down = -10;
      float svFitPt_RelSam_Down = -10;
      float svFitMass_RelSam_Up = -10;
      float svFitPt_RelSam_Up = -10;

      float svFitMass_Res_Down = -10;
      float svFitPt_Res_Down = -10;
      float svFitMass_Res_Up = -10;
      float svFitPt_Res_Up = -10;

      float svFitMass_UES_Down = -10;
      float svFitPt_UES_Down = -10;
      float svFitMass_UES_Up = -10;
      float svFitPt_UES_Up = -10;

      float svFitMass_RecoilReso_Up = -10;
      float svFitPt_RecoilReso_Up = -10;
      float svFitMass_RecoilReso_Down = -10;
      float svFitPt_RecoilReso_Down = -10;

      float svFitMass_RecoilResp_Up = -10;
      float svFitPt_RecoilResp_Up = -10;
      float svFitMass_RecoilResp_Down = -10;
      float svFitPt_RecoilResp_Down = -10;

      // tau leptons           

      float tau1_pt  = -10;
      float tau1_eta = -10;
      float tau1_phi = -10;
      float tau1_m   = -10;
      float tau2_pt  = -10;
      float tau2_eta = -10;
      float tau2_phi = -10;
      float tau2_m   = -10;

      TBranch *tauBranch1 = t->Branch("tau1_pt", &tau1_pt, "tau1_pt/F");
      TBranch *tauBranch2 = t->Branch("tau1_eta", &tau1_eta, "tau1_eta/F");
      TBranch *tauBranch3 = t->Branch("tau1_phi", &tau1_phi, "tau1_phi/F");
      TBranch *tauBranch4 = t->Branch("tau1_m", &tau1_m, "tau1_m/F");
      TBranch *tauBranch5 = t->Branch("tau2_pt", &tau2_pt, "tau2_pt/F");
      TBranch *tauBranch6 = t->Branch("tau2_eta", &tau2_eta, "tau2_eta/F");
      TBranch *tauBranch7 = t->Branch("tau2_phi", &tau2_phi, "tau2_phi/F");
      TBranch *tauBranch8 = t->Branch("tau2_m", &tau2_m, "tau2_m/F");

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

      TBranch *newMESBranch1 = t->Branch("m_sv_MES_Up", &svFitMass_MES_Up, "m_sv_MES_Up/F");
      TBranch *newMESBranch2 = t->Branch("pt_sv_MES_Up", &svFitPt_MES_Up, "pt_sv_MES_Up/F");
      TBranch *newMESBranch3 = t->Branch("m_sv_MES_Down", &svFitMass_MES_Down, "m_sv_MES_Down/F");
      TBranch *newMESBranch4 = t->Branch("pt_sv_MES_Down", &svFitPt_MES_Down, "pt_sv_MES_Down/F");

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

      TBranch *newBranchsvFitMass_JER_Down = t->Branch("m_sv_JetJER_Down", &svFitMass_JER_Down, "m_sv_JetJER_Down/F");
      TBranch *newBranchsvFitPt_JER_Down = t->Branch("pt_sv_JetJER_Down", &svFitPt_JER_Down, "pt_sv_JetJER_Down/F");
      TBranch *newBranchsvFitMass_JER_Up = t->Branch("m_sv_JetJER_Up", &svFitMass_JER_Up, "m_sv_JetJER_Up/F");
      TBranch *newBranchsvFitPt_JER_Up = t->Branch("pt_sv_JetJER_Up", &svFitPt_JER_Up, "pt_sv_JetJER_Up/F");

      TBranch *newBranchsvFitMass_Absolute_Down = t->Branch("m_sv_JetAbsolute_Down", &svFitMass_Absolute_Down, "m_sv_JetAbsolute_Down/F");
      TBranch *newBranchsvFitPt_Absolute_Down = t->Branch("pt_sv_JetAbsolute_Down", &svFitPt_Absolute_Down, "pt_sv_JetAbsolute_Down/F");
      TBranch *newBranchsvFitMass_Absolute_Up = t->Branch("m_sv_JetAbsolute_Up", &svFitMass_Absolute_Up, "m_sv_JetAbsolute_Up/F");
      TBranch *newBranchsvFitPt_Absolute_Up = t->Branch("pt_sv_JetAbsolute_Up", &svFitPt_Absolute_Up, "pt_sv_JetAbsolute_Up/F");

      TBranch *newBranchsvFitMass_Absoluteyear_Down = t->Branch("m_sv_JetAbsoluteyear_Down", &svFitMass_Absoluteyear_Down, "m_sv_JetAbsoluteyear_Down/F");
      TBranch *newBranchsvFitPt_Absoluteyear_Down = t->Branch("pt_sv_JetAbsoluteyear_Down", &svFitPt_Absoluteyear_Down, "pt_sv_JetAbsoluteyear_Down/F");
      TBranch *newBranchsvFitMass_Absoluteyear_Up = t->Branch("m_sv_JetAbsoluteyear_Up", &svFitMass_Absoluteyear_Up, "m_sv_JetAbsoluteyear_Up/F");
      TBranch *newBranchsvFitPt_Absoluteyear_Up = t->Branch("pt_sv_JetAbsoluteyear_Up", &svFitPt_Absoluteyear_Up, "pt_sv_JetAbsoluteyear_Up/F");

      TBranch *newBranchsvFitMass_BBEC1_Down = t->Branch("m_sv_JetBBEC1_Down", &svFitMass_BBEC1_Down, "m_sv_JetBBEC1_Down/F");
      TBranch *newBranchsvFitPt_BBEC1_Down = t->Branch("pt_sv_JetBBEC1_Down", &svFitPt_BBEC1_Down, "pt_sv_JetBBEC1_Down/F");
      TBranch *newBranchsvFitMass_BBEC1_Up = t->Branch("m_sv_JetBBEC1_Up", &svFitMass_BBEC1_Up, "m_sv_JetBBEC1_Up/F");
      TBranch *newBranchsvFitPt_BBEC1_Up = t->Branch("pt_sv_JetBBEC1_Up", &svFitPt_BBEC1_Up, "pt_sv_JetBBEC1_Up/F");

      TBranch *newBranchsvFitMass_BBEC1year_Down = t->Branch("m_sv_JetBBEC1year_Down", &svFitMass_BBEC1year_Down, "m_sv_JetBBEC1year_Down/F");
      TBranch *newBranchsvFitPt_BBEC1year_Down = t->Branch("pt_sv_JetBBEC1year_Down", &svFitPt_BBEC1year_Down, "pt_sv_JetBBEC1year_Down/F");
      TBranch *newBranchsvFitMass_BBEC1year_Up = t->Branch("m_sv_JetBBEC1year_Up", &svFitMass_BBEC1year_Up, "m_sv_JetBBEC1year_Up/F");
      TBranch *newBranchsvFitPt_BBEC1year_Up = t->Branch("pt_sv_JetBBEC1year_Up", &svFitPt_BBEC1year_Up, "pt_sv_JetBBEC1year_Up/F");

      TBranch *newBranchsvFitMass_EC2_Down = t->Branch("m_sv_JetEC2_Down", &svFitMass_EC2_Down, "m_sv_JetEC2_Down/F");
      TBranch *newBranchsvFitPt_EC2_Down = t->Branch("pt_sv_JetEC2_Down", &svFitPt_EC2_Down, "pt_sv_JetEC2_Down/F");
      TBranch *newBranchsvFitMass_EC2_Up = t->Branch("m_sv_JetEC2_Up", &svFitMass_EC2_Up, "m_sv_JetEC2_Up/F");
      TBranch *newBranchsvFitPt_EC2_Up = t->Branch("pt_sv_JetEC2_Up", &svFitPt_EC2_Up, "pt_sv_JetEC2_Up/F");

      TBranch *newBranchsvFitMass_EC2year_Down = t->Branch("m_sv_JetEC2year_Down", &svFitMass_EC2year_Down, "m_sv_JetEC2year_Down/F");
      TBranch *newBranchsvFitPt_EC2year_Down = t->Branch("pt_sv_JetEC2year_Down", &svFitPt_EC2year_Down, "pt_sv_JetEC2year_Down/F");
      TBranch *newBranchsvFitMass_EC2year_Up = t->Branch("m_sv_JetEC2year_Up", &svFitMass_EC2year_Up, "m_sv_JetEC2year_Up/F");
      TBranch *newBranchsvFitPt_EC2year_Up = t->Branch("pt_sv_JetEC2year_Up", &svFitPt_EC2year_Up, "pt_sv_JetEC2year_Up/F");

      TBranch *newBranchsvFitMass_FlavorQCD_Down = t->Branch("m_sv_JetFlavorQCD_Down", &svFitMass_FlavorQCD_Down, "m_sv_JetFlavorQCD_Down/F");
      TBranch *newBranchsvFitPt_FlavorQCD_Down = t->Branch("pt_sv_JetFlavorQCD_Down", &svFitPt_FlavorQCD_Down, "pt_sv_JetFlavorQCD_Down/F");
      TBranch *newBranchsvFitMass_FlavorQCD_Up = t->Branch("m_sv_JetFlavorQCD_Up", &svFitMass_FlavorQCD_Up, "m_sv_JetFlavorQCD_Up/F");
      TBranch *newBranchsvFitPt_FlavorQCD_Up = t->Branch("pt_sv_JetFlavorQCD_Up", &svFitPt_FlavorQCD_Up, "pt_sv_JetFlavorQCD_Up/F");

      TBranch *newBranchsvFitMass_HF_Down = t->Branch("m_sv_JetHF_Down", &svFitMass_HF_Down, "m_sv_JetHF_Down/F");
      TBranch *newBranchsvFitPt_HF_Down = t->Branch("pt_sv_JetHF_Down", &svFitPt_HF_Down, "pt_sv_JetHF_Down/F");
      TBranch *newBranchsvFitMass_HF_Up = t->Branch("m_sv_JetHF_Up", &svFitMass_HF_Up, "m_sv_JetHF_Up/F");
      TBranch *newBranchsvFitPt_HF_Up = t->Branch("pt_sv_JetHF_Up", &svFitPt_HF_Up, "pt_sv_JetHF_Up/F");

      TBranch *newBranchsvFitMass_HFyear_Down = t->Branch("m_sv_JetHFyear_Down", &svFitMass_HFyear_Down, "m_sv_JetHFyear_Down/F");
      TBranch *newBranchsvFitPt_HFyear_Down = t->Branch("pt_sv_JetHFyear_Down", &svFitPt_HFyear_Down, "pt_sv_JetHFyear_Down/F");
      TBranch *newBranchsvFitMass_HFyear_Up = t->Branch("m_sv_JetHFyear_Up", &svFitMass_HFyear_Up, "m_sv_JetHFyear_Up/F");
      TBranch *newBranchsvFitPt_HFyear_Up = t->Branch("pt_sv_JetHFyear_Up", &svFitPt_HFyear_Up, "pt_sv_JetHFyear_Up/F");

      TBranch *newBranchsvFitMass_RelBal_Down = t->Branch("m_sv_JetRelBal_Down", &svFitMass_RelBal_Down, "m_sv_JetRelBal_Down/F");
      TBranch *newBranchsvFitPt_RelBal_Down = t->Branch("pt_sv_JetRelBal_Down", &svFitPt_RelBal_Down, "pt_sv_JetRelBal_Down/F");
      TBranch *newBranchsvFitMass_RelBal_Up = t->Branch("m_sv_JetRelBal_Up", &svFitMass_RelBal_Up, "m_sv_JetRelBal_Up/F");
      TBranch *newBranchsvFitPt_RelBal_Up = t->Branch("pt_sv_JetRelBal_Up", &svFitPt_RelBal_Up, "pt_sv_JetRelBal_Up/F");

      TBranch *newBranchsvFitMass_RelSam_Down = t->Branch("m_sv_JetRelSam_Down", &svFitMass_RelSam_Down, "m_sv_JetRelSam_Down/F");
      TBranch *newBranchsvFitPt_RelSam_Down = t->Branch("pt_sv_JetRelSam_Down", &svFitPt_RelSam_Down, "pt_sv_JetRelSam_Down/F");
      TBranch *newBranchsvFitMass_RelSam_Up = t->Branch("m_sv_JetRelSam_Up", &svFitMass_RelSam_Up, "m_sv_JetRelSam_Up/F");
      TBranch *newBranchsvFitPt_RelSam_Up = t->Branch("pt_sv_JetRelSam_Up", &svFitPt_RelSam_Up, "pt_sv_JetRelSam_Up/F");

      TBranch *newBranchsvFitMass_Res_Down = t->Branch("m_sv_JetRes_Down", &svFitMass_Res_Down, "m_sv_JetRes_Down/F");
      TBranch *newBranchsvFitPt_Res_Down = t->Branch("pt_sv_JetRes_Down", &svFitPt_Res_Down, "pt_sv_JetRes_Down/F");
      TBranch *newBranchsvFitMass_Res_Up = t->Branch("m_sv_JetRes_Up", &svFitMass_Res_Up, "m_sv_JetRes_Up/F");
      TBranch *newBranchsvFitPt_Res_Up = t->Branch("pt_sv_JetRes_Up", &svFitPt_Res_Up, "pt_sv_JetRes_Up/F");

      TBranch *newBranchsvFitMass_UES_Down = t->Branch("m_sv_JetUES_Down", &svFitMass_UES_Down, "m_sv_JetUES_Down/F");
      TBranch *newBranchsvFitPt_UES_Down = t->Branch("pt_sv_JetUES_Down", &svFitPt_UES_Down, "pt_sv_JetUES_Down/F");
      TBranch *newBranchsvFitMass_UES_Up = t->Branch("m_sv_JetUES_Up", &svFitMass_UES_Up, "m_sv_JetUES_Up/F");
      TBranch *newBranchsvFitPt_UES_Up = t->Branch("pt_sv_JetUES_Up", &svFitPt_UES_Up, "pt_sv_JetUES_Up/F");

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

      // JES Uncertainties
      float met_reso_Up, met_reso_Down, met_resp_Up, met_resp_Down;
      float metphi_reso_Up, metphi_reso_Down, metphi_resp_Up, metphi_resp_Down;
      float met_JERUp, met_AbsoluteUp, met_AbsoluteyearUp, met_BBEC1Up, met_BBEC1yearUp, met_EC2Up, met_EC2yearUp, met_FlavorQCDUp, met_HFUp, met_HFyearUp,
            met_RelBalUp, met_RelSamUp, met_ResUp, met_UESUp, met_JERDown, met_AbsoluteDown, met_AbsoluteyearDown, met_BBEC1Down, met_BBEC1yearDown,
            met_EC2Down, met_EC2yearDown, met_FlavorQCDDown, met_HFDown, met_HFyearDown, met_RelBalDown, met_RelSamDown, met_ResDown,
            met_UESDown, metphi_JERUp, metphi_AbsoluteUp, metphi_AbsoluteyearUp, metphi_BBEC1Up, metphi_BBEC1yearUp, metphi_EC2Up, metphi_EC2yearUp,
            metphi_FlavorQCDUp, metphi_HFUp, metphi_HFyearUp, metphi_RelBalUp, metphi_RelSamUp, metphi_ResUp, metphi_UESUp, metphi_JERDown,
            metphi_AbsoluteDown, metphi_AbsoluteyearDown, metphi_BBEC1Down, metphi_BBEC1yearDown, metphi_EC2Down, metphi_EC2yearDown, metphi_FlavorQCDDown,
            metphi_HFDown, metphi_HFyearDown, metphi_RelBalDown, metphi_RelSamDown, metphi_ResDown, metphi_UESDown;

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

      t->SetBranchAddress("met_JERUp", &met_JERUp);
      t->SetBranchAddress("met_AbsoluteUp", &met_AbsoluteUp);
      t->SetBranchAddress("met_AbsoluteyearUp", &met_AbsoluteyearUp);
      t->SetBranchAddress("met_BBEC1Up", &met_BBEC1Up);
      t->SetBranchAddress("met_BBEC1yearUp", &met_BBEC1yearUp);
      t->SetBranchAddress("met_EC2Up", &met_EC2Up);
      t->SetBranchAddress("met_EC2yearUp", &met_EC2yearUp);
      t->SetBranchAddress("met_FlavorQCDUp", &met_FlavorQCDUp);
      t->SetBranchAddress("met_HFUp", &met_HFUp);
      t->SetBranchAddress("met_HFyearUp", &met_HFyearUp);
      t->SetBranchAddress("met_RelBalUp", &met_RelBalUp);
      t->SetBranchAddress("met_RelSamUp", &met_RelSamUp);
      t->SetBranchAddress("met_ResUp", &met_ResUp);
      t->SetBranchAddress("met_UESUp", &met_UESUp);
      t->SetBranchAddress("met_JERDown", &met_JERDown);
      t->SetBranchAddress("met_AbsoluteDown", &met_AbsoluteDown);
      t->SetBranchAddress("met_AbsoluteyearDown", &met_AbsoluteyearDown);
      t->SetBranchAddress("met_BBEC1Down", &met_BBEC1Down);
      t->SetBranchAddress("met_BBEC1yearDown", &met_BBEC1yearDown);
      t->SetBranchAddress("met_EC2Down", &met_EC2Down);
      t->SetBranchAddress("met_EC2yearDown", &met_EC2yearDown);
      t->SetBranchAddress("met_FlavorQCDDown", &met_FlavorQCDDown);
      t->SetBranchAddress("met_HFDown", &met_HFDown);
      t->SetBranchAddress("met_HFyearDown", &met_HFyearDown);
      t->SetBranchAddress("met_RelBalDown", &met_RelBalDown);
      t->SetBranchAddress("met_RelSamDown", &met_RelSamDown);
      t->SetBranchAddress("met_ResDown", &met_ResDown);
      t->SetBranchAddress("met_UESDown", &met_UESDown);
      t->SetBranchAddress("metphi_JERUp", &metphi_JERUp);
      t->SetBranchAddress("metphi_AbsoluteUp", &metphi_AbsoluteUp);
      t->SetBranchAddress("metphi_AbsoluteyearUp", &metphi_AbsoluteyearUp);
      t->SetBranchAddress("metphi_BBEC1Up", &metphi_BBEC1Up);
      t->SetBranchAddress("metphi_BBEC1yearUp", &metphi_BBEC1yearUp);
      t->SetBranchAddress("metphi_EC2Up", &metphi_EC2Up);
      t->SetBranchAddress("metphi_EC2yearUp", &metphi_EC2yearUp);
      t->SetBranchAddress("metphi_FlavorQCDUp", &metphi_FlavorQCDUp);
      t->SetBranchAddress("metphi_HFUp", &metphi_HFUp);
      t->SetBranchAddress("metphi_HFyearUp", &metphi_HFyearUp);
      t->SetBranchAddress("metphi_RelBalUp", &metphi_RelBalUp);
      t->SetBranchAddress("metphi_RelSamUp", &metphi_RelSamUp);
      t->SetBranchAddress("metphi_ResUp", &metphi_ResUp);
      t->SetBranchAddress("metphi_UESUp", &metphi_UESUp);
      t->SetBranchAddress("metphi_JERDown", &metphi_JERDown);
      t->SetBranchAddress("metphi_AbsoluteDown", &metphi_AbsoluteDown);
      t->SetBranchAddress("metphi_AbsoluteyearDown", &metphi_AbsoluteyearDown);
      t->SetBranchAddress("metphi_BBEC1Down", &metphi_BBEC1Down);
      t->SetBranchAddress("metphi_BBEC1yearDown", &metphi_BBEC1yearDown);
      t->SetBranchAddress("metphi_EC2Down", &metphi_EC2Down);
      t->SetBranchAddress("metphi_EC2yearDown", &metphi_EC2yearDown);
      t->SetBranchAddress("metphi_FlavorQCDDown", &metphi_FlavorQCDDown);
      t->SetBranchAddress("metphi_HFDown", &metphi_HFDown);
      t->SetBranchAddress("metphi_HFyearDown", &metphi_HFyearDown);
      t->SetBranchAddress("metphi_RelBalDown", &metphi_RelBalDown);
      t->SetBranchAddress("metphi_RelSamDown", &metphi_RelSamDown);
      t->SetBranchAddress("metphi_ResDown", &metphi_ResDown);
      t->SetBranchAddress("metphi_UESDown", &metphi_UESDown);

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

          metcorrJERDown_ex = met_JERDown*TMath::Cos(metphi_JERDown);
          metcorrJERUp_ex = met_JERUp*TMath::Cos(metphi_JERUp);
          metcorrAbsoluteDown_ex = met_AbsoluteDown*TMath::Cos(metphi_AbsoluteDown);
          metcorrAbsoluteUp_ex = met_AbsoluteUp*TMath::Cos(metphi_AbsoluteUp);
          metcorrAbsoluteyearDown_ex = met_AbsoluteyearDown*TMath::Cos(metphi_AbsoluteyearDown);
          metcorrAbsoluteyearUp_ex = met_AbsoluteyearUp*TMath::Cos(metphi_AbsoluteyearUp);
          metcorrBBEC1Down_ex = met_BBEC1Down*TMath::Cos(metphi_BBEC1Down);
          metcorrBBEC1Up_ex = met_BBEC1Up*TMath::Cos(metphi_BBEC1Up);
          metcorrBBEC1yearDown_ex = met_BBEC1yearDown*TMath::Cos(metphi_BBEC1yearDown);
          metcorrBBEC1yearUp_ex = met_BBEC1yearUp*TMath::Cos(metphi_BBEC1yearUp);
          metcorrEC2Down_ex = met_EC2Down*TMath::Cos(metphi_EC2Down);
          metcorrEC2Up_ex = met_EC2Up*TMath::Cos(metphi_EC2Up);
          metcorrEC2yearDown_ex = met_EC2yearDown*TMath::Cos(metphi_EC2yearDown);
          metcorrEC2yearUp_ex = met_EC2yearUp*TMath::Cos(metphi_EC2yearUp);
          metcorrFlavorQCDDown_ex = met_FlavorQCDDown*TMath::Cos(metphi_FlavorQCDDown);
          metcorrFlavorQCDUp_ex = met_FlavorQCDUp*TMath::Cos(metphi_FlavorQCDUp);
          metcorrHFDown_ex = met_HFDown*TMath::Cos(metphi_HFDown);
          metcorrHFUp_ex = met_HFUp*TMath::Cos(metphi_HFUp);
          metcorrHFyearDown_ex = met_HFyearDown*TMath::Cos(metphi_HFyearDown);
          metcorrHFyearUp_ex = met_HFyearUp*TMath::Cos(metphi_HFyearUp);
          metcorrRelBalDown_ex = met_RelBalDown*TMath::Cos(metphi_RelBalDown);
          metcorrRelBalUp_ex = met_RelBalUp*TMath::Cos(metphi_RelBalUp);
          metcorrRelSamDown_ex = met_RelSamDown*TMath::Cos(metphi_RelSamDown);
          metcorrRelSamUp_ex = met_RelSamUp*TMath::Cos(metphi_RelSamUp);
          metcorrResDown_ex = met_ResDown*TMath::Cos(metphi_ResDown);
          metcorrResUp_ex = met_ResUp*TMath::Cos(metphi_ResUp);
          metcorrUESDown_ex = met_UESDown*TMath::Cos(metphi_UESDown);
          metcorrUESUp_ex = met_UESUp*TMath::Cos(metphi_UESUp);

          metcorrJERDown_ey = met_JERDown*TMath::Sin(metphi_JERDown);
          metcorrJERUp_ey = met_JERUp*TMath::Sin(metphi_JERUp);
          metcorrAbsoluteDown_ey = met_AbsoluteDown*TMath::Sin(metphi_AbsoluteDown);
          metcorrAbsoluteUp_ey = met_AbsoluteUp*TMath::Sin(metphi_AbsoluteUp);
          metcorrAbsoluteyearDown_ey = met_AbsoluteyearDown*TMath::Sin(metphi_AbsoluteyearDown);
          metcorrAbsoluteyearUp_ey = met_AbsoluteyearUp*TMath::Sin(metphi_AbsoluteyearUp);
          metcorrBBEC1Down_ey = met_BBEC1Down*TMath::Sin(metphi_BBEC1Down);
          metcorrBBEC1Up_ey = met_BBEC1Up*TMath::Sin(metphi_BBEC1Up);
          metcorrBBEC1yearDown_ey = met_BBEC1yearDown*TMath::Sin(metphi_BBEC1yearDown);
          metcorrBBEC1yearUp_ey = met_BBEC1yearUp*TMath::Sin(metphi_BBEC1yearUp);
          metcorrEC2Down_ey = met_EC2Down*TMath::Sin(metphi_EC2Down);
          metcorrEC2Up_ey = met_EC2Up*TMath::Sin(metphi_EC2Up);
          metcorrEC2yearDown_ey = met_EC2yearDown*TMath::Sin(metphi_EC2yearDown);
          metcorrEC2yearUp_ey = met_EC2yearUp*TMath::Sin(metphi_EC2yearUp);
          metcorrFlavorQCDDown_ey = met_FlavorQCDDown*TMath::Sin(metphi_FlavorQCDDown);
          metcorrFlavorQCDUp_ey = met_FlavorQCDUp*TMath::Sin(metphi_FlavorQCDUp);
          metcorrHFDown_ey = met_HFDown*TMath::Sin(metphi_HFDown);
          metcorrHFUp_ey = met_HFUp*TMath::Sin(metphi_HFUp);
          metcorrHFyearDown_ey = met_HFyearDown*TMath::Sin(metphi_HFyearDown);
          metcorrHFyearUp_ey = met_HFyearUp*TMath::Sin(metphi_HFyearUp);
          metcorrRelBalDown_ey = met_RelBalDown*TMath::Sin(metphi_RelBalDown);
          metcorrRelBalUp_ey = met_RelBalUp*TMath::Sin(metphi_RelBalUp);
          metcorrRelSamDown_ey = met_RelSamDown*TMath::Sin(metphi_RelSamDown);
          metcorrRelSamUp_ey = met_RelSamUp*TMath::Sin(metphi_RelSamUp);
          metcorrResDown_ey = met_ResDown*TMath::Sin(metphi_ResDown);
          metcorrResUp_ey = met_ResUp*TMath::Sin(metphi_ResUp);
          metcorrUESDown_ey = met_UESDown*TMath::Sin(metphi_UESDown);
          metcorrUESUp_ey = met_UESUp*TMath::Sin(metphi_UESUp);

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

          runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0,
               svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          // MET systematics
          runSVFit(measuredTauLeptons, metcorrUESUp_ex, metcorrUESUp_ey, covMET, 0, svFitMass_UES_Up, svFitPt_UES_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrUESDown_ex, metcorrUESDown_ey, covMET, 0, svFitMass_UES_Down, svFitPt_UES_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJERUp_ex, metcorrJERUp_ey, covMET, 0, svFitMass_JER_Up, svFitPt_JER_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrJERDown_ex, metcorrJERDown_ey, covMET, 0, svFitMass_JER_Down, svFitPt_JER_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrAbsoluteUp_ex, metcorrAbsoluteUp_ey, covMET, 0, svFitMass_Absolute_Up, svFitPt_Absolute_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrAbsoluteDown_ex, metcorrAbsoluteDown_ey, covMET, 0, svFitMass_Absolute_Down, svFitPt_Absolute_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrAbsoluteyearUp_ex, metcorrAbsoluteyearUp_ey, covMET, 0, svFitMass_Absoluteyear_Up, svFitPt_Absoluteyear_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrAbsoluteyearDown_ex, metcorrAbsoluteyearDown_ey, covMET, 0, svFitMass_Absoluteyear_Down, svFitPt_Absoluteyear_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrBBEC1Up_ex, metcorrBBEC1Up_ey, covMET, 0, svFitMass_BBEC1_Up, svFitPt_BBEC1_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrBBEC1Down_ex, metcorrBBEC1Down_ey, covMET, 0, svFitMass_BBEC1_Down, svFitPt_BBEC1_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrBBEC1yearUp_ex, metcorrBBEC1yearUp_ey, covMET, 0, svFitMass_BBEC1year_Up, svFitPt_BBEC1year_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrBBEC1yearDown_ex, metcorrBBEC1yearDown_ey, covMET, 0, svFitMass_BBEC1year_Down, svFitPt_BBEC1year_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrEC2Up_ex, metcorrEC2Up_ey, covMET, 0, svFitMass_EC2_Up, svFitPt_EC2_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrEC2Down_ex, metcorrEC2Down_ey, covMET, 0, svFitMass_EC2_Down, svFitPt_EC2_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrEC2yearUp_ex, metcorrEC2yearUp_ey, covMET, 0, svFitMass_EC2year_Up, svFitPt_EC2year_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrEC2yearDown_ex, metcorrEC2yearDown_ey, covMET, 0, svFitMass_EC2year_Down, svFitPt_EC2year_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrFlavorQCDUp_ex, metcorrFlavorQCDUp_ey, covMET, 0, svFitMass_FlavorQCD_Up, svFitPt_FlavorQCD_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrFlavorQCDDown_ex, metcorrFlavorQCDDown_ey, covMET, 0, svFitMass_FlavorQCD_Down, svFitPt_FlavorQCD_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrHFUp_ex, metcorrHFUp_ey, covMET, 0, svFitMass_HF_Up, svFitPt_HF_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrHFDown_ex, metcorrHFDown_ey, covMET, 0, svFitMass_HF_Down, svFitPt_HF_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrHFyearUp_ex, metcorrHFyearUp_ey, covMET, 0, svFitMass_HFyear_Up, svFitPt_HFyear_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrHFyearDown_ex, metcorrHFyearDown_ey, covMET, 0, svFitMass_HFyear_Down, svFitPt_HFyear_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRelBalUp_ex, metcorrRelBalUp_ey, covMET, 0, svFitMass_RelBal_Up, svFitPt_RelBal_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRelBalDown_ex, metcorrRelBalDown_ey, covMET, 0, svFitMass_RelBal_Down, svFitPt_RelBal_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRelSamUp_ex, metcorrRelSamUp_ey, covMET, 0, svFitMass_RelSam_Up, svFitPt_RelSam_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrRelSamDown_ex, metcorrRelSamDown_ey, covMET, 0, svFitMass_RelSam_Down, svFitPt_RelSam_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrResUp_ex, metcorrResUp_ey, covMET, 0, svFitMass_Res_Up, svFitPt_Res_Up,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

          runSVFit(measuredTauLeptons, metcorrResDown_ex, metcorrResDown_ey, covMET, 0, svFitMass_Res_Down, svFitPt_Res_Down,
              svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);

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
            float ES_Up = 1 + tfes->getFES(decayMode2, eta2, gen_match_2, "up") + tfes->getTES(decayMode2, gen_match_2, "up");
            float ES_Down = 1 + tfes->getFES(decayMode2, eta2, gen_match_2, "down") + tfes->getTES(decayMode2, gen_match_2, "down");


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

              runSVFit(measuredTauUp, metcorr_ex, metcorr_ey, covMET, 0,
                svFitMass_MES_Up, svFitPt_MES_Up, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
              runSVFit(measuredTauDn, metcorr_ex, metcorr_ey, covMET, 0,
                svFitMass_MES_Down, svFitPt_MES_Down, svFitEta, svFitPhi, svFitMET, svFitTransverseMass, tau1, tau2);
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

            // tau DM0 shifted down
            if (gen_match_2 == 5 && decayMode2 == 0) {
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

    tau1_pt  = tau1.Pt();
    tau1_eta = tau1.Eta();
    tau1_phi = tau1.Phi();
    tau1_m   = tau1.M();
    tau2_pt  = tau2.Pt();
    tau2_eta = tau2.Eta();
    tau2_phi = tau2.Phi();
    tau2_m   = tau2.M();

    tauBranch1->Fill();
    tauBranch2->Fill();
    tauBranch3->Fill();
    tauBranch4->Fill();
    tauBranch5->Fill();
    tauBranch6->Fill();
    tauBranch7->Fill();
    tauBranch8->Fill();
    newBranch1->Fill();
    newBranch2->Fill();
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
    newMESBranch1->Fill();
    newMESBranch2->Fill();
    newMESBranch3->Fill();
    newMESBranch4->Fill();
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
    newBranchsvFitMass_JER_Down->Fill();
    newBranchsvFitPt_JER_Down->Fill();
    newBranchsvFitMass_JER_Up->Fill();
    newBranchsvFitPt_JER_Up->Fill();
    newBranchsvFitMass_Absolute_Down->Fill();
    newBranchsvFitPt_Absolute_Down->Fill();
    newBranchsvFitMass_Absolute_Up->Fill();
    newBranchsvFitPt_Absolute_Up->Fill();
    newBranchsvFitMass_Absoluteyear_Down->Fill();
    newBranchsvFitPt_Absoluteyear_Down->Fill();
    newBranchsvFitMass_Absoluteyear_Up->Fill();
    newBranchsvFitPt_Absoluteyear_Up->Fill();
    newBranchsvFitMass_BBEC1_Down->Fill();
    newBranchsvFitPt_BBEC1_Down->Fill();
    newBranchsvFitMass_BBEC1_Up->Fill();
    newBranchsvFitPt_BBEC1_Up->Fill();
    newBranchsvFitMass_BBEC1year_Down->Fill();
    newBranchsvFitPt_BBEC1year_Down->Fill();
    newBranchsvFitMass_BBEC1year_Up->Fill();
    newBranchsvFitPt_BBEC1year_Up->Fill();
    newBranchsvFitMass_EC2_Down->Fill();
    newBranchsvFitPt_EC2_Down->Fill();
    newBranchsvFitMass_EC2_Up->Fill();
    newBranchsvFitPt_EC2_Up->Fill();
    newBranchsvFitMass_EC2year_Down->Fill();
    newBranchsvFitPt_EC2year_Down->Fill();
    newBranchsvFitMass_EC2year_Up->Fill();
    newBranchsvFitPt_EC2year_Up->Fill();
    newBranchsvFitMass_FlavorQCD_Down->Fill();
    newBranchsvFitPt_FlavorQCD_Down->Fill();
    newBranchsvFitMass_FlavorQCD_Up->Fill();
    newBranchsvFitPt_FlavorQCD_Up->Fill();
    newBranchsvFitMass_HF_Down->Fill();
    newBranchsvFitPt_HF_Down->Fill();
    newBranchsvFitMass_HF_Up->Fill();
    newBranchsvFitPt_HF_Up->Fill();
    newBranchsvFitMass_HFyear_Down->Fill();
    newBranchsvFitPt_HFyear_Down->Fill();
    newBranchsvFitMass_HFyear_Up->Fill();
    newBranchsvFitPt_HFyear_Up->Fill();
    newBranchsvFitMass_RelBal_Down->Fill();
    newBranchsvFitPt_RelBal_Down->Fill();
    newBranchsvFitMass_RelBal_Up->Fill();
    newBranchsvFitPt_RelBal_Up->Fill();
    newBranchsvFitMass_RelSam_Down->Fill();
    newBranchsvFitPt_RelSam_Down->Fill();
    newBranchsvFitMass_RelSam_Up->Fill();
    newBranchsvFitPt_RelSam_Up->Fill();
    newBranchsvFitMass_Res_Down->Fill();
    newBranchsvFitPt_Res_Down->Fill();
    newBranchsvFitMass_Res_Up->Fill();
    newBranchsvFitPt_Res_Up->Fill();
    newBranchsvFitMass_UES_Down->Fill();
    newBranchsvFitPt_UES_Down->Fill();
    newBranchsvFitMass_UES_Up->Fill();
    newBranchsvFitPt_UES_Up->Fill();
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
