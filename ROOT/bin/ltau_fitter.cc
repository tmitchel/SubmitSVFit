#include <math.h>
#include <limits>
#include <string>
#include <vector>
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "../interface/vardump.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

// If doES       0 does not apply any ES shifts
//              1 applies ES shifts to TT channel, no effect on other channels
//
// If isWJets    0 no shift in number of jets used for recoil corrections
//              1 shifts njets + 1 for recoil corrections
//
// If metType    1 use mvamet
//        -1 use pf met

ClassicSVfit svfitAlgorithm;

int copyFiles(optutl::CommandLineParser parser, TFile *fOld, TFile *fNew);
void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES, int isWJets, int metType, double tesSize);
int CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source, optutl::CommandLineParser parser);
double tesUncertainties(unsigned int year, float decaymode);
double pt_shifted(float pt, double tesUnc, bool isDM, int updown);
double metcorr_shifted(double metcorr, float pt1, float phi1, bool isDM1, double tesUnc1, float pt2, float phi2, bool isDM2, double tesUnc2, int xory,
                       int updown);
void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> &measuredTauLeptons, double measuredMETx, double measuredMETy, TMatrixD &covMET,
              float num, float &svFitMass, float &svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass,
              TLorentzVector &tau1, TLorentzVector &tau2);
void four_vector(TLorentzVector p4, Float_t& pt, Float_t& eta, Float_t& phi, Float_t& mass);

int main(int argc, char *argv[]) {
    optutl::CommandLineParser parser("Sets Event Weights in the ntuple");
    parser.addOption("branch", optutl::CommandLineParser::kString);
    parser.addOption("newFile", optutl::CommandLineParser::kString);
    parser.addOption("inputFile", optutl::CommandLineParser::kString);
    parser.addOption("doES", optutl::CommandLineParser::kDouble);
    parser.addOption("isWJets", optutl::CommandLineParser::kDouble);
    parser.addOption("metType", optutl::CommandLineParser::kDouble);   // 1 = mvamet, -1 = pf met
    parser.addOption("tesSize", optutl::CommandLineParser::kDouble);  // Default TES = 1.2%
    parser.addOption("numEvents", optutl::CommandLineParser::kInteger);
    parser.parseArguments(argc, argv);

    std::cout << "EXTRA COMMANDS:"
              << "\n --- numEvents: " << parser.integerValue("numEvents") << "\n --- doES: " << parser.doubleValue("doES")
              << "\n --- isWJets: " << parser.doubleValue("isWJets") << "\n --- metType: " << parser.doubleValue("metType")
              << "\n --- tesSize: " << parser.doubleValue("tesSize") << std::endl;

    // Make sure a proper Met Type is chosen
    assert(parser.doubleValue("metType") == 1.0 || parser.doubleValue("metType") == -1.0);

    // No DiTauMass constraint
    svfitAlgorithm.setDiTauMassConstraint(-1.0);

    char TreeToUse[80] = "first";

    TFile *fProduce;  //= new TFile(parser.stringValue("newFile").c_str(),"UPDATE");

    TFile *f = new TFile(parser.stringValue("inputFile").c_str(), "READ");
    std::cout << "Creating new outputfile" << std::endl;
    std::string newFileName = parser.stringValue("newFile");

    fProduce = new TFile(newFileName.c_str(), "RECREATE");
    if (copyFiles(parser, f, fProduce) == 0) return -1;

    fProduce = new TFile(newFileName.c_str(), "UPDATE");
    std::cout << "listing the directories=================" << std::endl;
    fProduce->ls();
    readdir(fProduce, parser, TreeToUse, parser.doubleValue("doES"), parser.doubleValue("isWJets"), parser.doubleValue("metType"),
            parser.doubleValue("tesSize"));

    fProduce->Close();
    f->Close();
}

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES, int isWJets, int metType, double tesSize) {
    TLorentzVector tau1, tau2;
    // up systematics
    TLorentzVector tau1_Up, tau1_DM0_Up, tau1_DM1_Up, tau1_DM10_Up, tau1_UncMet_Up, tau1_ClusteredMet_Up;
    TLorentzVector tau1_JetEC2_Up, tau1_JetEta0to3_Up, tau1_JetEta0to5_Up, tau1_JetEta3to5_Up, tau1_JetRelativeBal_Up, tau1_JetRelativeSample_Up;
    TLorentzVector tau2_Up, tau2_DM0_Up, tau2_DM1_Up, tau2_DM10_Up, tau2_UncMet_Up, tau2_ClusteredMet_Up;
    TLorentzVector tau2_JetEC2_Up, tau2_JetEta0to3_Up, tau2_JetEta0to5_Up, tau2_JetEta3to5_Up, tau2_JetRelativeBal_Up, tau2_JetRelativeSample_Up;
    // down systematics
    TLorentzVector tau1_Down, tau1_DM0_Down, tau1_DM1_Down, tau1_DM10_Down, tau1_UncMet_Down, tau1_ClusteredMet_Down;
    TLorentzVector tau1_JetEC2_Down, tau1_JetEta0to3_Down, tau1_JetEta0to5_Down, tau1_JetEta3to5_Down, tau1_JetRelativeBal_Down,
        tau1_JetRelativeSample_Down;
    TLorentzVector tau2_Down, tau2_DM0_Down, tau2_DM1_Down, tau2_DM10_Down, tau2_UncMet_Down, tau2_ClusteredMet_Down;
    TLorentzVector tau2_JetEC2_Down, tau2_JetEta0to3_Down, tau2_JetEta0to5_Down, tau2_JetEta3to5_Down, tau2_JetRelativeBal_Down,
        tau2_JetRelativeSample_Down;

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
    while ((key = (TKey *)next())) {
        printf("Found key=%s \n", key->GetName());

        TObject *obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
            std::cout << "This is a directory, diving in!" << std::endl;
            // zero the processedNames vector, to allow processing trees with duplicate names in separate directories
            processedNames.clear();

            dir->cd(key->GetName());
            TDirectory *subdir = gDirectory;
            sprintf(TreeToUse,"%s",key->GetName());
            readdir(subdir, parser, TreeToUse, parser.doubleValue("doES"), parser.doubleValue("isWJets"), parser.doubleValue("metType"),
                    parser.doubleValue("tesSize"));

            dirsav->cd();
        } else if (obj->IsA()->InheritsFrom(TTree::Class())) {
            // check  if this tree was already processed
            std::vector<TString>::const_iterator it = find(processedNames.begin(), processedNames.end(), key->GetName());
            if (it != processedNames.end()) {
                std::cout << "This tree was already processed, skipping..." << std::endl;
                continue;
            }
            std::cout << "This is the tree! Start processing" << std::endl;
            processedNames.push_back(key->GetName());

            // Identify the process
            if (std::string(key->GetName()).find("tt") != std::string::npos) {
                decayType1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
                decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
                mass1 = 0.13957;
                mass2 = 0.13957;
                channel = "tt";
                std::cout << "Identified channel tt and using kappa = 5" << std::endl;
                svfitAlgorithm.addLogM_fixed(true, 5);
            } else if (std::string(key->GetName()).find("em") != std::string::npos) {
                std::cout << "EMu sample" << std::endl;
                decayType1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
                decayType2 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
                mass1 = 0.00051100;
                mass2 = 0.105658;
                channel = "em";
                std::cout << "Identified channel em and using kappa = 3" << std::endl;
                svfitAlgorithm.addLogM_fixed(true, 3);
            } else if (std::string(key->GetName()).find("et") != std::string::npos) {
                std::cout << "eleTauTree" << std::endl;
                decayType1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
                decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
                mass1 = 0.00051100;
                mass2 = 0;
                channel = "et";
                std::cout << "Identified channel et and using kappa = 4" << std::endl;
                svfitAlgorithm.addLogM_fixed(true, 4);
            } else if (std::string(key->GetName()).find("mt") != std::string::npos ||
                       std::string(key->GetName()).find("mutau") != std::string::npos) {
                std::cout << "muTauEvent" << std::endl;
                decayType1 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
                decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
                mass1 = 0.105658;
                mass2 = 0;
                channel = "mt";
                std::cout << "Identified channel mt and using kappa = 4" << std::endl;
                svfitAlgorithm.addLogM_fixed(true, 4);
            } else {
                std::cout << "Tree " << key->GetName() << " does not match ... Skipping!!" << std::endl;
                return;
            }

            TTree *t = (TTree *)obj;

            std::vector<std::pair<std::string, float*>> branch_info = {
                std::make_pair("m_sv", &svFitMass),
                std::make_pair("pt_sv", &svFitPt),
                std::make_pair("eta_sv", &svFitEta),
                std::make_pair("phi_sv", &svFitPhi),
                std::make_pair("met_sv", &svFitMET),
                std::make_pair("mt_sv", &svFitTransverseMass),
                std::make_pair("metcorr_ex", &metcorr_ex),
                std::make_pair("metcorr_ey", &metcorr_ey),
                std::make_pair("metcor", &metcor),
                std::make_pair("metcorphi", &metcorphi),
                std::make_pair("m_sv_Up", &svFitMass_Up),
                std::make_pair("pt_sv_Up", &svFitPt_Up),
                std::make_pair("eta_sv_Up", &svFitEta_Up),
                std::make_pair("phi_sv_Up", &svFitPhi_Up),
                std::make_pair("met_sv_Up", &svFitMET_Up),
                std::make_pair("mt_sv_Up", &svFitTransverseMass_Up),
                std::make_pair("m_sv_Down", &svFitMass_Down),
                std::make_pair("pt_sv_Down", &svFitPt_Down),
                std::make_pair("eta_sv_Down", &svFitEta_Down),
                std::make_pair("phi_sv_Down", &svFitPhi_Down),
                std::make_pair("met_sv_Down", &svFitMET_Down),
                std::make_pair("mt_sv_Down", &svFitTransverseMass_Down),
                std::make_pair("m_sv_DM0_Up", &svFitMass_DM0_Up),
                std::make_pair("pt_sv_DM0_Up", &svFitPt_DM0_Up),
                std::make_pair("eta_sv_DM0_Up", &svFitEta_DM0_Up),
                std::make_pair("phi_sv_DM0_Up", &svFitPhi_DM0_Up),
                std::make_pair("met_sv_DM0_Up", &svFitMET_DM0_Up),
                std::make_pair("mt_sv_DM0_Up", &svFitTransverseMass_DM0_Up),
                std::make_pair("m_sv_DM0_Down", &svFitMass_DM0_Down),
                std::make_pair("pt_sv_DM0_Down", &svFitPt_DM0_Down),
                std::make_pair("eta_sv_DM0_Down", &svFitEta_DM0_Down),
                std::make_pair("phi_sv_DM0_Down", &svFitPhi_DM0_Down),
                std::make_pair("met_sv_DM0_Down", &svFitMET_DM0_Down),
                std::make_pair("mt_sv_DM0_Down", &svFitTransverseMass_DM0_Down),
                std::make_pair("m_sv_DM1_Up", &svFitMass_DM1_Up),
                std::make_pair("pt_sv_DM1_Up", &svFitPt_DM1_Up),
                std::make_pair("eta_sv_DM1_Up", &svFitEta_DM1_Up),
                std::make_pair("phi_sv_DM1_Up", &svFitPhi_DM1_Up),
                std::make_pair("met_sv_DM1_Up", &svFitMET_DM1_Up),
                std::make_pair("mt_sv_DM1_Up", &svFitTransverseMass_DM1_Up),
                std::make_pair("m_sv_DM1_Down", &svFitMass_DM1_Down),
                std::make_pair("pt_sv_DM1_Down", &svFitPt_DM1_Down),
                std::make_pair("eta_sv_DM1_Down", &svFitEta_DM1_Down),
                std::make_pair("phi_sv_DM1_Down", &svFitPhi_DM1_Down),
                std::make_pair("met_sv_DM1_Down", &svFitMET_DM1_Down),
                std::make_pair("mt_sv_DM1_Down", &svFitTransverseMass_DM1_Down),
                std::make_pair("m_sv_DM10_Up", &svFitMass_DM10_Up),
                std::make_pair("pt_sv_DM10_Up", &svFitPt_DM10_Up),
                std::make_pair("eta_sv_DM10_Up", &svFitEta_DM10_Up),
                std::make_pair("phi_sv_DM10_Up", &svFitPhi_DM10_Up),
                std::make_pair("met_sv_DM10_Up", &svFitMET_DM10_Up),
                std::make_pair("mt_sv_DM10_Up", &svFitTransverseMass_DM10_Up),
                std::make_pair("m_sv_DM10_Down", &svFitMass_DM10_Down),
                std::make_pair("pt_sv_DM10_Down", &svFitPt_DM10_Down),
                std::make_pair("eta_sv_DM10_Down", &svFitEta_DM10_Down),
                std::make_pair("phi_sv_DM10_Down", &svFitPhi_DM10_Down),
                std::make_pair("met_sv_DM10_Down", &svFitMET_DM10_Down),
                std::make_pair("mt_sv_DM10_Down", &svFitTransverseMass_DM10_Down),
                std::make_pair("m_sv_UncMet_Up", &svFitMass_UncMet_Up),
                std::make_pair("pt_sv_UncMet_Up", &svFitPt_UncMet_Up),
                std::make_pair("eta_sv_UncMet_Up", &svFitEta_UncMet_Up),
                std::make_pair("phi_sv_UncMet_Up", &svFitPhi_UncMet_Up),
                std::make_pair("met_sv_UncMet_Up", &svFitMET_UncMet_Up),
                std::make_pair("mt_sv_UncMet_Up", &svFitTransverseMass_UncMet_Up),
                std::make_pair("m_sv_UncMet_Down", &svFitMass_UncMet_Down),
                std::make_pair("pt_sv_UncMet_Down", &svFitPt_UncMet_Down),
                std::make_pair("eta_sv_UncMet_Down", &svFitEta_UncMet_Down),
                std::make_pair("phi_sv_UncMet_Down", &svFitPhi_UncMet_Down),
                std::make_pair("met_sv_UncMet_Down", &svFitMET_UncMet_Down),
                std::make_pair("mt_sv_UncMet_Down", &svFitTransverseMass_UncMet_Down),
                std::make_pair("m_sv_ClusteredMet_Up", &svFitMass_ClusteredMet_Up),
                std::make_pair("pt_sv_ClusteredMet_Up", &svFitPt_ClusteredMet_Up),
                std::make_pair("eta_sv_ClusteredMet_Up", &svFitEta_ClusteredMet_Up),
                std::make_pair("phi_sv_ClusteredMet_Up", &svFitPhi_ClusteredMet_Up),
                std::make_pair("met_sv_ClusteredMet_Up", &svFitMET_ClusteredMet_Up),
                std::make_pair("mt_sv_ClusteredMet_Up", &svFitTransverseMass_ClusteredMet_Up),
                std::make_pair("m_sv_ClusteredMet_Down", &svFitMass_ClusteredMet_Down),
                std::make_pair("pt_sv_ClusteredMet_Down", &svFitPt_ClusteredMet_Down),
                std::make_pair("eta_sv_ClusteredMet_Down", &svFitEta_ClusteredMet_Down),
                std::make_pair("phi_sv_ClusteredMet_Down", &svFitPhi_ClusteredMet_Down),
                std::make_pair("met_sv_ClusteredMet_Down", &svFitMET_ClusteredMet_Down),
                std::make_pair("mt_sv_ClusteredMet_Down", &svFitTransverseMass_ClusteredMet_Down),
                std::make_pair("m_sv_JetEC2_Up", &svFitMass_JetEC2_Up),
                std::make_pair("pt_sv_JetEC2_Up", &svFitPt_JetEC2_Up),
                std::make_pair("eta_sv_JetEC2_Up", &svFitEta_JetEC2_Up),
                std::make_pair("phi_sv_JetEC2_Up", &svFitPhi_JetEC2_Up),
                std::make_pair("met_sv_JetEC2_Up", &svFitMET_JetEC2_Up),
                std::make_pair("mt_sv_JetEC2_Up", &svFitTransverseMass_JetEC2_Up),
                std::make_pair("m_sv_JetEC2_Down", &svFitMass_JetEC2_Down),
                std::make_pair("pt_sv_JetEC2_Down", &svFitPt_JetEC2_Down),
                std::make_pair("eta_sv_JetEC2_Down", &svFitEta_JetEC2_Down),
                std::make_pair("phi_sv_JetEC2_Down", &svFitPhi_JetEC2_Down),
                std::make_pair("met_sv_JetEC2_Down", &svFitMET_JetEC2_Down),
                std::make_pair("mt_sv_JetEC2_Down", &svFitTransverseMass_JetEC2_Down),
                std::make_pair("m_sv_JetEta0to3_Up", &svFitMass_JetEta0to3_Up),
                std::make_pair("pt_sv_JetEta0to3_Up", &svFitPt_JetEta0to3_Up),
                std::make_pair("eta_sv_JetEta0to3_Up", &svFitEta_JetEta0to3_Up),
                std::make_pair("phi_sv_JetEta0to3_Up", &svFitPhi_JetEta0to3_Up),
                std::make_pair("met_sv_JetEta0to3_Up", &svFitMET_JetEta0to3_Up),
                std::make_pair("mt_sv_JetEta0to3_Up", &svFitTransverseMass_JetEta0to3_Up),
                std::make_pair("m_sv_JetEta0to3_Down", &svFitMass_JetEta0to3_Down),
                std::make_pair("pt_sv_JetEta0to3_Down", &svFitPt_JetEta0to3_Down),
                std::make_pair("eta_sv_JetEta0to3_Down", &svFitEta_JetEta0to3_Down),
                std::make_pair("phi_sv_JetEta0to3_Down", &svFitPhi_JetEta0to3_Down),
                std::make_pair("met_sv_JetEta0to3_Down", &svFitMET_JetEta0to3_Down),
                std::make_pair("mt_sv_JetEta0to3_Down", &svFitTransverseMass_JetEta0to3_Down),
                std::make_pair("m_sv_JetEta0to5_Up", &svFitMass_JetEta0to5_Up),
                std::make_pair("pt_sv_JetEta0to5_Up", &svFitPt_JetEta0to5_Up),
                std::make_pair("eta_sv_JetEta0to5_Up", &svFitEta_JetEta0to5_Up),
                std::make_pair("phi_sv_JetEta0to5_Up", &svFitPhi_JetEta0to5_Up),
                std::make_pair("met_sv_JetEta0to5_Up", &svFitMET_JetEta0to5_Up),
                std::make_pair("mt_sv_JetEta0to5_Up", &svFitTransverseMass_JetEta0to5_Up),
                std::make_pair("m_sv_JetEta0to5_Down", &svFitMass_JetEta0to5_Down),
                std::make_pair("pt_sv_JetEta0to5_Down", &svFitPt_JetEta0to5_Down),
                std::make_pair("eta_sv_JetEta0to5_Down", &svFitEta_JetEta0to5_Down),
                std::make_pair("phi_sv_JetEta0to5_Down", &svFitPhi_JetEta0to5_Down),
                std::make_pair("met_sv_JetEta0to5_Down", &svFitMET_JetEta0to5_Down),
                std::make_pair("mt_sv_JetEta0to5_Down", &svFitTransverseMass_JetEta0to5_Down),
                std::make_pair("m_sv_JetEta3to5_Up", &svFitMass_JetEta3to5_Up),
                std::make_pair("pt_sv_JetEta3to5_Up", &svFitPt_JetEta3to5_Up),
                std::make_pair("eta_sv_JetEta3to5_Up", &svFitEta_JetEta3to5_Up),
                std::make_pair("phi_sv_JetEta3to5_Up", &svFitPhi_JetEta3to5_Up),
                std::make_pair("met_sv_JetEta3to5_Up", &svFitMET_JetEta3to5_Up),
                std::make_pair("mt_sv_JetEta3to5_Up", &svFitTransverseMass_JetEta3to5_Up),
                std::make_pair("m_sv_JetEta3to5_Down", &svFitMass_JetEta3to5_Down),
                std::make_pair("pt_sv_JetEta3to5_Down", &svFitPt_JetEta3to5_Down),
                std::make_pair("eta_sv_JetEta3to5_Down", &svFitEta_JetEta3to5_Down),
                std::make_pair("phi_sv_JetEta3to5_Down", &svFitPhi_JetEta3to5_Down),
                std::make_pair("met_sv_JetEta3to5_Down", &svFitMET_JetEta3to5_Down),
                std::make_pair("mt_sv_JetEta3to5_Down", &svFitTransverseMass_JetEta3to5_Down),
                std::make_pair("m_sv_JetRelativeBal_Up", &svFitMass_JetRelativeBal_Up),
                std::make_pair("pt_sv_JetRelativeBal_Up", &svFitPt_JetRelativeBal_Up),
                std::make_pair("eta_sv_JetRelativeBal_Up", &svFitEta_JetRelativeBal_Up),
                std::make_pair("phi_sv_JetRelativeBal_Up", &svFitPhi_JetRelativeBal_Up),
                std::make_pair("met_sv_JetRelativeBal_Up", &svFitMET_JetRelativeBal_Up),
                std::make_pair("mt_sv_JetRelativeBal_Up", &svFitTransverseMass_JetRelativeBal_Up),
                std::make_pair("m_sv_JetRelativeBal_Down", &svFitMass_JetRelativeBal_Down),
                std::make_pair("pt_sv_JetRelativeBal_Down", &svFitPt_JetRelativeBal_Down),
                std::make_pair("eta_sv_JetRelativeBal_Down", &svFitEta_JetRelativeBal_Down),
                std::make_pair("phi_sv_JetRelativeBal_Down", &svFitPhi_JetRelativeBal_Down),
                std::make_pair("met_sv_JetRelativeBal_Down", &svFitMET_JetRelativeBal_Down),
                std::make_pair("mt_sv_JetRelativeBal_Down", &svFitTransverseMass_JetRelativeBal_Down),
                std::make_pair("m_sv_JetRelativeSample_Up", &svFitMass_JetRelativeSample_Up),
                std::make_pair("pt_sv_JetRelativeSample_Up", &svFitPt_JetRelativeSample_Up),
                std::make_pair("eta_sv_JetRelativeSample_Up", &svFitEta_JetRelativeSample_Up),
                std::make_pair("phi_sv_JetRelativeSample_Up", &svFitPhi_JetRelativeSample_Up),
                std::make_pair("met_sv_JetRelativeSample_Up", &svFitMET_JetRelativeSample_Up),
                std::make_pair("mt_sv_JetRelativeSample_Up", &svFitTransverseMass_JetRelativeSample_Up),
                std::make_pair("m_sv_JetRelativeSample_Down", &svFitMass_JetRelativeSample_Down),
                std::make_pair("pt_sv_JetRelativeSample_Down", &svFitPt_JetRelativeSample_Down),
                std::make_pair("eta_sv_JetRelativeSample_Down", &svFitEta_JetRelativeSample_Down),
                std::make_pair("phi_sv_JetRelativeSample_Down", &svFitPhi_JetRelativeSample_Down),
                std::make_pair("met_sv_JetRelativeSample_Down", &svFitMET_JetRelativeSample_Down),
                std::make_pair("mt_sv_JetRelativeSample_Down", &svFitTransverseMass_JetRelativeSample_Down),
                std::make_pair("metcorClusteredDown", &metcorClusteredDown),
                std::make_pair("metcorphiClusteredDown", &metcorphiClusteredDown),
                std::make_pair("metcorClusteredUp", &metcorClusteredUp),
                std::make_pair("metcorphiClusteredUp", &metcorphiClusteredUp),
                std::make_pair("metcorUncDown", &metcorUncDown),
                std::make_pair("metcorphiUncDown", &metcorphiUncDown),
                std::make_pair("metcorUncUp", &metcorUncUp),
                std::make_pair("metcorphiUncUp", &metcorphiUncUp),
                std::make_pair("metcorJetEC2Down", &metcorJetEC2Down),
                std::make_pair("metcorphiJetEC2Down", &metcorphiJetEC2Down),
                std::make_pair("metcorJetEC2Up", &metcorJetEC2Up),
                std::make_pair("metcorphiJetEC2Up", &metcorphiJetEC2Up),
                std::make_pair("metcorJetEta0to3Down", &metcorJetEta0to3Down),
                std::make_pair("metcorphiJetEta0to3Down", &metcorphiJetEta0to3Down),
                std::make_pair("metcorJetEta0to3Up", &metcorJetEta0to3Up),
                std::make_pair("metcorphiJetEta0to3Up", &metcorphiJetEta0to3Up),
                std::make_pair("metcorJetEta0to5Down", &metcorJetEta0to5Down),
                std::make_pair("metcorphiJetEta0to5Down", &metcorphiJetEta0to5Down),
                std::make_pair("metcorJetEta0to5Up", &metcorJetEta0to5Up),
                std::make_pair("metcorphiJetEta0to5Up", &metcorphiJetEta0to5Up),
                std::make_pair("metcorJetEta3to5Down", &metcorJetEta3to5Down),
                std::make_pair("metcorphiJetEta3to5Down", &metcorphiJetEta3to5Down),
                std::make_pair("metcorJetEta3to5Up", &metcorJetEta3to5Up),
                std::make_pair("metcorphiJetEta3to5Up", &metcorphiJetEta3to5Up),
                std::make_pair("metcorJetRelativeBalDown", &metcorJetRelativeBalDown),
                std::make_pair("metcorphiJetRelativeBalDown", &metcorphiJetRelativeBalDown),
                std::make_pair("metcorJetRelativeBalUp", &metcorJetRelativeBalUp),
                std::make_pair("metcorphiJetRelativeBalUp", &metcorphiJetRelativeBalUp),
                std::make_pair("metcorJetRelativeSampleDown", &metcorJetRelativeSampleDown),
                std::make_pair("metcorphiJetRelativeSampleDown", &metcorphiJetRelativeSampleDown),
                std::make_pair("metcorJetRelativeSampleUp", &metcorJetRelativeSampleUp),
                std::make_pair("metcorphiJetRelativeSampleUp", &metcorphiJetRelativeSampleUp),
                std::make_pair("tau1_pt", &tau1_pt),
                std::make_pair("tau1_eta", &tau1_eta),
                std::make_pair("tau1_phi", &tau1_phi),
                std::make_pair("tau1_m", &tau1_m),
                std::make_pair("tau2_pt", &tau2_pt),
                std::make_pair("tau2_eta", &tau2_eta),
                std::make_pair("tau2_phi", &tau2_phi),
                std::make_pair("tau2_m", &tau2_m),
                std::make_pair("tau1_pt_Up", &tau1_pt_Up),
                std::make_pair("tau1_eta_Up", &tau1_eta_Up),
                std::make_pair("tau1_phi_Up", &tau1_phi_Up),
                std::make_pair("tau1_m_Up", &tau1_m_Up),
                std::make_pair("tau2_pt_Up", &tau2_pt_Up),
                std::make_pair("tau2_eta_Up", &tau2_eta_Up),
                std::make_pair("tau2_phi_Up", &tau2_phi_Up),
                std::make_pair("tau2_m_Up", &tau2_m_Up),
                std::make_pair("tau1_pt_Down", &tau1_pt_Down),
                std::make_pair("tau1_eta_Down", &tau1_eta_Down),
                std::make_pair("tau1_phi_Down", &tau1_phi_Down),
                std::make_pair("tau1_m_Down", &tau1_m_Down),
                std::make_pair("tau2_pt_Down", &tau2_pt_Down),
                std::make_pair("tau2_eta_Down", &tau2_eta_Down),
                std::make_pair("tau2_phi_Down", &tau2_phi_Down),
                std::make_pair("tau2_m_Down", &tau2_m_Down),
                std::make_pair("tau1_pt_DM0_Up", &tau1_pt_DM0_Up),
                std::make_pair("tau1_eta_DM0_Up", &tau1_eta_DM0_Up),
                std::make_pair("tau1_phi_DM0_Up", &tau1_phi_DM0_Up),
                std::make_pair("tau1_m_DM0_Up", &tau1_m_DM0_Up),
                std::make_pair("tau2_pt_DM0_Up", &tau2_pt_DM0_Up),
                std::make_pair("tau2_eta_DM0_Up", &tau2_eta_DM0_Up),
                std::make_pair("tau2_phi_DM0_Up", &tau2_phi_DM0_Up),
                std::make_pair("tau2_m_DM0_Up", &tau2_m_DM0_Up),
                std::make_pair("tau1_pt_DM0_Down", &tau1_pt_DM0_Down),
                std::make_pair("tau1_eta_DM0_Down", &tau1_eta_DM0_Down),
                std::make_pair("tau1_phi_DM0_Down", &tau1_phi_DM0_Down),
                std::make_pair("tau1_m_DM0_Down", &tau1_m_DM0_Down),
                std::make_pair("tau2_pt_DM0_Down", &tau2_pt_DM0_Down),
                std::make_pair("tau2_eta_DM0_Down", &tau2_eta_DM0_Down),
                std::make_pair("tau2_phi_DM0_Down", &tau2_phi_DM0_Down),
                std::make_pair("tau2_m_DM0_Down", &tau2_m_DM0_Down),
                std::make_pair("tau1_pt_DM1_Up", &tau1_pt_DM1_Up),
                std::make_pair("tau1_eta_DM1_Up", &tau1_eta_DM1_Up),
                std::make_pair("tau1_phi_DM1_Up", &tau1_phi_DM1_Up),
                std::make_pair("tau1_m_DM1_Up", &tau1_m_DM1_Up),
                std::make_pair("tau2_pt_DM1_Up", &tau2_pt_DM1_Up),
                std::make_pair("tau2_eta_DM1_Up", &tau2_eta_DM1_Up),
                std::make_pair("tau2_phi_DM1_Up", &tau2_phi_DM1_Up),
                std::make_pair("tau2_m_DM1_Up", &tau2_m_DM1_Up),
                std::make_pair("tau1_pt_DM1_Down", &tau1_pt_DM1_Down),
                std::make_pair("tau1_eta_DM1_Down", &tau1_eta_DM1_Down),
                std::make_pair("tau1_phi_DM1_Down", &tau1_phi_DM1_Down),
                std::make_pair("tau1_m_DM1_Down", &tau1_m_DM1_Down),
                std::make_pair("tau2_pt_DM1_Down", &tau2_pt_DM1_Down),
                std::make_pair("tau2_eta_DM1_Down", &tau2_eta_DM1_Down),
                std::make_pair("tau2_phi_DM1_Down", &tau2_phi_DM1_Down),
                std::make_pair("tau2_m_DM1_Down", &tau2_m_DM1_Down),
                std::make_pair("tau1_pt_DM10_Up", &tau1_pt_DM10_Up),
                std::make_pair("tau1_eta_DM10_Up", &tau1_eta_DM10_Up),
                std::make_pair("tau1_phi_DM10_Up", &tau1_phi_DM10_Up),
                std::make_pair("tau1_m_DM10_Up", &tau1_m_DM10_Up),
                std::make_pair("tau2_pt_DM10_Up", &tau2_pt_DM10_Up),
                std::make_pair("tau2_eta_DM10_Up", &tau2_eta_DM10_Up),
                std::make_pair("tau2_phi_DM10_Up", &tau2_phi_DM10_Up),
                std::make_pair("tau2_m_DM10_Up", &tau2_m_DM10_Up),
                std::make_pair("tau1_pt_DM10_Down", &tau1_pt_DM10_Down),
                std::make_pair("tau1_eta_DM10_Down", &tau1_eta_DM10_Down),
                std::make_pair("tau1_phi_DM10_Down", &tau1_phi_DM10_Down),
                std::make_pair("tau1_m_DM10_Down", &tau1_m_DM10_Down),
                std::make_pair("tau2_pt_DM10_Down", &tau2_pt_DM10_Down),
                std::make_pair("tau2_eta_DM10_Down", &tau2_eta_DM10_Down),
                std::make_pair("tau2_phi_DM10_Down", &tau2_phi_DM10_Down),
                std::make_pair("tau2_m_DM10_Down", &tau2_m_DM10_Down),
                std::make_pair("tau1_pt_UncMet_Up", &tau1_pt_UncMet_Up),
                std::make_pair("tau1_eta_UncMet_Up", &tau1_eta_UncMet_Up),
                std::make_pair("tau1_phi_UncMet_Up", &tau1_phi_UncMet_Up),
                std::make_pair("tau1_m_UncMet_Up", &tau1_m_UncMet_Up),
                std::make_pair("tau2_pt_UncMet_Up", &tau2_pt_UncMet_Up),
                std::make_pair("tau2_eta_UncMet_Up", &tau2_eta_UncMet_Up),
                std::make_pair("tau2_phi_UncMet_Up", &tau2_phi_UncMet_Up),
                std::make_pair("tau2_m_UncMet_Up", &tau2_m_UncMet_Up),
                std::make_pair("tau1_pt_UncMet_Down", &tau1_pt_UncMet_Down),
                std::make_pair("tau1_eta_UncMet_Down", &tau1_eta_UncMet_Down),
                std::make_pair("tau1_phi_UncMet_Down", &tau1_phi_UncMet_Down),
                std::make_pair("tau1_m_UncMet_Down", &tau1_m_UncMet_Down),
                std::make_pair("tau2_pt_UncMet_Down", &tau2_pt_UncMet_Down),
                std::make_pair("tau2_eta_UncMet_Down", &tau2_eta_UncMet_Down),
                std::make_pair("tau2_phi_UncMet_Down", &tau2_phi_UncMet_Down),
                std::make_pair("tau2_m_UncMet_Down", &tau2_m_UncMet_Down),
                std::make_pair("tau1_pt_ClusteredMet_Up", &tau1_pt_ClusteredMet_Up),
                std::make_pair("tau1_eta_ClusteredMet_Up", &tau1_eta_ClusteredMet_Up),
                std::make_pair("tau1_phi_ClusteredMet_Up", &tau1_phi_ClusteredMet_Up),
                std::make_pair("tau1_m_ClusteredMet_Up", &tau1_m_ClusteredMet_Up),
                std::make_pair("tau2_pt_ClusteredMet_Up", &tau2_pt_ClusteredMet_Up),
                std::make_pair("tau2_eta_ClusteredMet_Up", &tau2_eta_ClusteredMet_Up),
                std::make_pair("tau2_phi_ClusteredMet_Up", &tau2_phi_ClusteredMet_Up),
                std::make_pair("tau2_m_ClusteredMet_Up", &tau2_m_ClusteredMet_Up),
                std::make_pair("tau1_pt_ClusteredMet_Down", &tau1_pt_ClusteredMet_Down),
                std::make_pair("tau1_eta_ClusteredMet_Down", &tau1_eta_ClusteredMet_Down),
                std::make_pair("tau1_phi_ClusteredMet_Down", &tau1_phi_ClusteredMet_Down),
                std::make_pair("tau1_m_ClusteredMet_Down", &tau1_m_ClusteredMet_Down),
                std::make_pair("tau2_pt_ClusteredMet_Down", &tau2_pt_ClusteredMet_Down),
                std::make_pair("tau2_eta_ClusteredMet_Down", &tau2_eta_ClusteredMet_Down),
                std::make_pair("tau2_phi_ClusteredMet_Down", &tau2_phi_ClusteredMet_Down),
                std::make_pair("tau2_m_ClusteredMet_Down", &tau2_m_ClusteredMet_Down),
                std::make_pair("tau1_pt_JetEC2_Up", &tau1_pt_JetEC2_Up),
                std::make_pair("tau1_eta_JetEC2_Up", &tau1_eta_JetEC2_Up),
                std::make_pair("tau1_phi_JetEC2_Up", &tau1_phi_JetEC2_Up),
                std::make_pair("tau1_m_JetEC2_Up", &tau1_m_JetEC2_Up),
                std::make_pair("tau2_pt_JetEC2_Up", &tau2_pt_JetEC2_Up),
                std::make_pair("tau2_eta_JetEC2_Up", &tau2_eta_JetEC2_Up),
                std::make_pair("tau2_phi_JetEC2_Up", &tau2_phi_JetEC2_Up),
                std::make_pair("tau2_m_JetEC2_Up", &tau2_m_JetEC2_Up),
                std::make_pair("tau1_pt_JetEC2_Down", &tau1_pt_JetEC2_Down),
                std::make_pair("tau1_eta_JetEC2_Down", &tau1_eta_JetEC2_Down),
                std::make_pair("tau1_phi_JetEC2_Down", &tau1_phi_JetEC2_Down),
                std::make_pair("tau1_m_JetEC2_Down", &tau1_m_JetEC2_Down),
                std::make_pair("tau2_pt_JetEC2_Down", &tau2_pt_JetEC2_Down),
                std::make_pair("tau2_eta_JetEC2_Down", &tau2_eta_JetEC2_Down),
                std::make_pair("tau2_phi_JetEC2_Down", &tau2_phi_JetEC2_Down),
                std::make_pair("tau2_m_JetEC2_Down", &tau2_m_JetEC2_Down),
                std::make_pair("tau1_pt_JetEta0to3_Up", &tau1_pt_JetEta0to3_Up),
                std::make_pair("tau1_eta_JetEta0to3_Up", &tau1_eta_JetEta0to3_Up),
                std::make_pair("tau1_phi_JetEta0to3_Up", &tau1_phi_JetEta0to3_Up),
                std::make_pair("tau1_m_JetEta0to3_Up", &tau1_m_JetEta0to3_Up),
                std::make_pair("tau2_pt_JetEta0to3_Up", &tau2_pt_JetEta0to3_Up),
                std::make_pair("tau2_eta_JetEta0to3_Up", &tau2_eta_JetEta0to3_Up),
                std::make_pair("tau2_phi_JetEta0to3_Up", &tau2_phi_JetEta0to3_Up),
                std::make_pair("tau2_m_JetEta0to3_Up", &tau2_m_JetEta0to3_Up),
                std::make_pair("tau1_pt_JetEta0to3_Down", &tau1_pt_JetEta0to3_Down),
                std::make_pair("tau1_eta_JetEta0to3_Down", &tau1_eta_JetEta0to3_Down),
                std::make_pair("tau1_phi_JetEta0to3_Down", &tau1_phi_JetEta0to3_Down),
                std::make_pair("tau1_m_JetEta0to3_Down", &tau1_m_JetEta0to3_Down),
                std::make_pair("tau2_pt_JetEta0to3_Down", &tau2_pt_JetEta0to3_Down),
                std::make_pair("tau2_eta_JetEta0to3_Down", &tau2_eta_JetEta0to3_Down),
                std::make_pair("tau2_phi_JetEta0to3_Down", &tau2_phi_JetEta0to3_Down),
                std::make_pair("tau2_m_JetEta0to3_Down", &tau2_m_JetEta0to3_Down),
                std::make_pair("tau1_pt_JetEta0to5_Up", &tau1_pt_JetEta0to5_Up),
                std::make_pair("tau1_eta_JetEta0to5_Up", &tau1_eta_JetEta0to5_Up),
                std::make_pair("tau1_phi_JetEta0to5_Up", &tau1_phi_JetEta0to5_Up),
                std::make_pair("tau1_m_JetEta0to5_Up", &tau1_m_JetEta0to5_Up),
                std::make_pair("tau2_pt_JetEta0to5_Up", &tau2_pt_JetEta0to5_Up),
                std::make_pair("tau2_eta_JetEta0to5_Up", &tau2_eta_JetEta0to5_Up),
                std::make_pair("tau2_phi_JetEta0to5_Up", &tau2_phi_JetEta0to5_Up),
                std::make_pair("tau2_m_JetEta0to5_Up", &tau2_m_JetEta0to5_Up),
                std::make_pair("tau1_pt_JetEta0to5_Down", &tau1_pt_JetEta0to5_Down),
                std::make_pair("tau1_eta_JetEta0to5_Down", &tau1_eta_JetEta0to5_Down),
                std::make_pair("tau1_phi_JetEta0to5_Down", &tau1_phi_JetEta0to5_Down),
                std::make_pair("tau1_m_JetEta0to5_Down", &tau1_m_JetEta0to5_Down),
                std::make_pair("tau2_pt_JetEta0to5_Down", &tau2_pt_JetEta0to5_Down),
                std::make_pair("tau2_eta_JetEta0to5_Down", &tau2_eta_JetEta0to5_Down),
                std::make_pair("tau2_phi_JetEta0to5_Down", &tau2_phi_JetEta0to5_Down),
                std::make_pair("tau2_m_JetEta0to5_Down", &tau2_m_JetEta0to5_Down),
                std::make_pair("tau1_pt_JetEta3to5_Up", &tau1_pt_JetEta3to5_Up),
                std::make_pair("tau1_eta_JetEta3to5_Up", &tau1_eta_JetEta3to5_Up),
                std::make_pair("tau1_phi_JetEta3to5_Up", &tau1_phi_JetEta3to5_Up),
                std::make_pair("tau1_m_JetEta3to5_Up", &tau1_m_JetEta3to5_Up),
                std::make_pair("tau2_pt_JetEta3to5_Up", &tau2_pt_JetEta3to5_Up),
                std::make_pair("tau2_eta_JetEta3to5_Up", &tau2_eta_JetEta3to5_Up),
                std::make_pair("tau2_phi_JetEta3to5_Up", &tau2_phi_JetEta3to5_Up),
                std::make_pair("tau2_m_JetEta3to5_Up", &tau2_m_JetEta3to5_Up),
                std::make_pair("tau1_pt_JetEta3to5_Down", &tau1_pt_JetEta3to5_Down),
                std::make_pair("tau1_eta_JetEta3to5_Down", &tau1_eta_JetEta3to5_Down),
                std::make_pair("tau1_phi_JetEta3to5_Down", &tau1_phi_JetEta3to5_Down),
                std::make_pair("tau1_m_JetEta3to5_Down", &tau1_m_JetEta3to5_Down),
                std::make_pair("tau2_pt_JetEta3to5_Down", &tau2_pt_JetEta3to5_Down),
                std::make_pair("tau2_eta_JetEta3to5_Down", &tau2_eta_JetEta3to5_Down),
                std::make_pair("tau2_phi_JetEta3to5_Down", &tau2_phi_JetEta3to5_Down),
                std::make_pair("tau2_m_JetEta3to5_Down", &tau2_m_JetEta3to5_Down),
                std::make_pair("tau1_pt_JetRelativeBal_Up", &tau1_pt_JetRelativeBal_Up),
                std::make_pair("tau1_eta_JetRelativeBal_Up", &tau1_eta_JetRelativeBal_Up),
                std::make_pair("tau1_phi_JetRelativeBal_Up", &tau1_phi_JetRelativeBal_Up),
                std::make_pair("tau1_m_JetRelativeBal_Up", &tau1_m_JetRelativeBal_Up),
                std::make_pair("tau2_pt_JetRelativeBal_Up", &tau2_pt_JetRelativeBal_Up),
                std::make_pair("tau2_eta_JetRelativeBal_Up", &tau2_eta_JetRelativeBal_Up),
                std::make_pair("tau2_phi_JetRelativeBal_Up", &tau2_phi_JetRelativeBal_Up),
                std::make_pair("tau2_m_JetRelativeBal_Up", &tau2_m_JetRelativeBal_Up),
                std::make_pair("tau1_pt_JetRelativeBal_Down", &tau1_pt_JetRelativeBal_Down),
                std::make_pair("tau1_eta_JetRelativeBal_Down", &tau1_eta_JetRelativeBal_Down),
                std::make_pair("tau1_phi_JetRelativeBal_Down", &tau1_phi_JetRelativeBal_Down),
                std::make_pair("tau1_m_JetRelativeBal_Down", &tau1_m_JetRelativeBal_Down),
                std::make_pair("tau2_pt_JetRelativeBal_Down", &tau2_pt_JetRelativeBal_Down),
                std::make_pair("tau2_eta_JetRelativeBal_Down", &tau2_eta_JetRelativeBal_Down),
                std::make_pair("tau2_phi_JetRelativeBal_Down", &tau2_phi_JetRelativeBal_Down),
                std::make_pair("tau2_m_JetRelativeBal_Down", &tau2_m_JetRelativeBal_Down),
                std::make_pair("tau1_pt_JetRelativeSample_Up", &tau1_pt_JetRelativeSample_Up),
                std::make_pair("tau1_eta_JetRelativeSample_Up", &tau1_eta_JetRelativeSample_Up),
                std::make_pair("tau1_phi_JetRelativeSample_Up", &tau1_phi_JetRelativeSample_Up),
                std::make_pair("tau1_m_JetRelativeSample_Up", &tau1_m_JetRelativeSample_Up),
                std::make_pair("tau2_pt_JetRelativeSample_Up", &tau2_pt_JetRelativeSample_Up),
                std::make_pair("tau2_eta_JetRelativeSample_Up", &tau2_eta_JetRelativeSample_Up),
                std::make_pair("tau2_phi_JetRelativeSample_Up", &tau2_phi_JetRelativeSample_Up),
                std::make_pair("tau2_m_JetRelativeSample_Up", &tau2_m_JetRelativeSample_Up),
                std::make_pair("tau1_pt_JetRelativeSample_Down", &tau1_pt_JetRelativeSample_Down),
                std::make_pair("tau1_eta_JetRelativeSample_Down", &tau1_eta_JetRelativeSample_Down),
                std::make_pair("tau1_phi_JetRelativeSample_Down", &tau1_phi_JetRelativeSample_Down),
                std::make_pair("tau1_m_JetRelativeSample_Down", &tau1_m_JetRelativeSample_Down),
                std::make_pair("tau2_pt_JetRelativeSample_Down", &tau2_pt_JetRelativeSample_Down),
                std::make_pair("tau2_eta_JetRelativeSample_Down", &tau2_eta_JetRelativeSample_Down),
                std::make_pair("tau2_phi_JetRelativeSample_Down", &tau2_phi_JetRelativeSample_Down),
                std::make_pair("tau2_m_JetRelativeSample_Down", &tau2_m_JetRelativeSample_Down)
            };
            std::vector<TBranch*> fit_results;
            for (auto entry : branch_info) {
                fit_results.push_back(t->Branch(entry.first.c_str(), entry.second, (entry.first + "/F").c_str()));
            }
            std::cout << "That's a lot of tau 4-vector branches! N = " << fit_results.size() << std::endl;

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
            float decayMode = -999.;
            float decayMode2;
            float mvaCovMatrix00;
            float mvaCovMatrix10;
            float mvaCovMatrix01;
            float mvaCovMatrix11;
            float pfCovMatrix00;
            float pfCovMatrix10;
            float pfCovMatrix01;
            float pfCovMatrix11;
            // float mvamet_ex, // uncorrected mva met px (float)
            //  mvamet_ey, // uncorrected mva met py (float)

            int njets = -999.;  // number of jets (hadronic jet multiplicity) (int)

            // define MET
            float mvamet;
            float mvametphi;
            float pfmet;
            float pfmetphi;
            TLorentzVector TMet(0, 0, 0, 0);
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

            // ele/mu variables
            TBranch *pt1branch;

            t->SetBranchAddress("era", &era);
            t->SetBranchAddress("NtupleVer", &NtupleVer);
            t->SetBranchAddress("evt", &evt);
            t->SetBranchAddress("run", &run);
            t->SetBranchAddress("lumi", &lumi);
            t->SetBranchAddress("gen_match_1", &gen_match_1);
            t->SetBranchAddress("gen_match_2", &gen_match_2);
            if (channel == "tt") t->SetBranchAddress("t1_decayMode", &decayMode);
            if (channel == "tt") t->SetBranchAddress("t2_decayMode", &decayMode2);
            t->SetBranchAddress("pt_1", &pt1, &pt1branch);
            t->SetBranchAddress("eta_1", &eta1);
            t->SetBranchAddress("phi_1", &phi1);
            t->SetBranchAddress("pt_2", &pt2);
            t->SetBranchAddress("eta_2", &eta2);
            t->SetBranchAddress("phi_2", &phi2);
            t->SetBranchAddress("m_2", &m2);
            // t->SetBranchAddress("l1_decayMode",&decayMode);
            if (channel != "tt") t->SetBranchAddress("l2_decayMode", &decayMode2);
            t->SetBranchAddress("mvacov00", &mvaCovMatrix00);  // branch not stored in et/mt trees
            t->SetBranchAddress("mvacov01", &mvaCovMatrix01);  // branch not stored in et/mt trees
            t->SetBranchAddress("mvacov10", &mvaCovMatrix10);  // branch not stored in et/mt trees
            t->SetBranchAddress("mvacov11", &mvaCovMatrix11);  // branch not stored in et/mt trees
            t->SetBranchAddress("mvamet", &mvamet);            // branch not stored in et/mt trees
            t->SetBranchAddress("mvametphi", &mvametphi);      // branch not stored in et/mt trees
            t->SetBranchAddress("njets", &njets);
            t->SetBranchAddress("met", &pfmet);
            t->SetBranchAddress("metphi", &pfmetphi);
            // FOR PF MET ANALYSIS
            t->SetBranchAddress("metcov00", &pfCovMatrix00);
            t->SetBranchAddress("metcov01", &pfCovMatrix01);
            t->SetBranchAddress("metcov10", &pfCovMatrix10);
            t->SetBranchAddress("metcov11", &pfCovMatrix11);
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

            // double tesUp = 1.0 + tesSize;
            // double tesDown = 1.0 - tesSize;

            int nevents = t->GetEntries();
            if (parser.integerValue("numEvents") != -1) nevents = parser.integerValue("numEvents");
            for (Int_t i = 0; i < nevents; ++i) {
                t->GetEntry(i);

                // always PFMet
                TMet.SetPtEtaPhiM(pfmet, 0, pfmetphi, 0);
                metcorr_ex = pfmet * TMath::Cos(pfmetphi);
                metcorr_ey = pfmet * TMath::Sin(pfmetphi);
                // Shifted METs
                metcorrUncUp_ex = uncMetPtUp * TMath::Cos(uncMetPhiUp);
                metcorrUncUp_ey = uncMetPtUp * TMath::Sin(uncMetPhiUp);
                metcorrUncDown_ex = uncMetPtDown * TMath::Cos(uncMetPhiDown);
                metcorrUncDown_ey = uncMetPtDown * TMath::Sin(uncMetPhiDown);
                metcorrClusteredUp_ex = clusteredMetPtUp * TMath::Cos(clusteredMetPhiUp);
                metcorrClusteredUp_ey = clusteredMetPtUp * TMath::Sin(clusteredMetPhiUp);
                metcorrClusteredDown_ex = clusteredMetPtDown * TMath::Cos(clusteredMetPhiDown);
                metcorrClusteredDown_ey = clusteredMetPtDown * TMath::Sin(clusteredMetPhiDown);
                metcorrJetEC2Up_ex = met_JetEC2Up * TMath::Cos(metphi_JetEC2Up);
                metcorrJetEC2Up_ey = met_JetEC2Up * TMath::Sin(metphi_JetEC2Up);
                metcorrJetEC2Down_ex = met_JetEC2Down * TMath::Cos(metphi_JetEC2Down);
                metcorrJetEC2Down_ey = met_JetEC2Down * TMath::Sin(metphi_JetEC2Down);
                metcorrJetEta0to3Up_ex = met_JetEta0to3Up * TMath::Cos(metphi_JetEta0to3Up);
                metcorrJetEta0to3Up_ey = met_JetEta0to3Up * TMath::Sin(metphi_JetEta0to3Up);
                metcorrJetEta0to3Down_ex = met_JetEta0to3Down * TMath::Cos(metphi_JetEta0to3Down);
                metcorrJetEta0to3Down_ey = met_JetEta0to3Down * TMath::Sin(metphi_JetEta0to3Down);
                metcorrJetEta0to5Up_ex = met_JetEta0to5Up * TMath::Cos(metphi_JetEta0to5Up);
                metcorrJetEta0to5Up_ey = met_JetEta0to5Up * TMath::Sin(metphi_JetEta0to5Up);
                metcorrJetEta0to5Down_ex = met_JetEta0to5Down * TMath::Cos(metphi_JetEta0to5Down);
                metcorrJetEta0to5Down_ey = met_JetEta0to5Down * TMath::Sin(metphi_JetEta0to5Down);
                metcorrJetEta3to5Up_ex = met_JetEta3to5Up * TMath::Cos(metphi_JetEta3to5Up);
                metcorrJetEta3to5Up_ey = met_JetEta3to5Up * TMath::Sin(metphi_JetEta3to5Up);
                metcorrJetEta3to5Down_ex = met_JetEta3to5Down * TMath::Cos(metphi_JetEta3to5Down);
                metcorrJetEta3to5Down_ey = met_JetEta3to5Down * TMath::Sin(metphi_JetEta3to5Down);
                metcorrJetRelativeBalUp_ex = met_JetRelativeBalUp * TMath::Cos(metphi_JetRelativeBalUp);
                metcorrJetRelativeBalUp_ey = met_JetRelativeBalUp * TMath::Sin(metphi_JetRelativeBalUp);
                metcorrJetRelativeBalDown_ex = met_JetRelativeBalDown * TMath::Cos(metphi_JetRelativeBalDown);
                metcorrJetRelativeBalDown_ey = met_JetRelativeBalDown * TMath::Sin(metphi_JetRelativeBalDown);
                metcorrJetRelativeSampleUp_ex = met_JetRelativeSampleUp * TMath::Cos(metphi_JetRelativeSampleUp);
                metcorrJetRelativeSampleUp_ey = met_JetRelativeSampleUp * TMath::Sin(metphi_JetRelativeSampleUp);
                metcorrJetRelativeSampleDown_ex = met_JetRelativeSampleDown * TMath::Cos(metphi_JetRelativeSampleDown);
                metcorrJetRelativeSampleDown_ey = met_JetRelativeSampleDown * TMath::Sin(metphi_JetRelativeSampleDown);

                covMET[0][0] = pfCovMatrix00;
                covMET[1][0] = pfCovMatrix10;
                covMET[0][1] = pfCovMatrix01;
                covMET[1][1] = pfCovMatrix11;

                metcor = TMath::Sqrt(metcorr_ex * metcorr_ex + metcorr_ey * metcorr_ey);
                metcorphi = TMath::ATan2(metcorr_ey, metcorr_ex);
                std::cout << " - metcor " << metcor << " metcorphi " << metcorphi << std::endl;

                // Corrected MET values for saving
                metcorClusteredDown =
                    TMath::Sqrt(metcorrClusteredDown_ex * metcorrClusteredDown_ex + metcorrClusteredDown_ey * metcorrClusteredDown_ey);
                metcorphiClusteredDown = TMath::ATan2(metcorrClusteredDown_ey, metcorrClusteredDown_ex);

                metcorClusteredUp = TMath::Sqrt(metcorrClusteredUp_ex * metcorrClusteredUp_ex + metcorrClusteredUp_ey * metcorrClusteredUp_ey);
                metcorphiClusteredUp = TMath::ATan2(metcorrClusteredUp_ey, metcorrClusteredUp_ex);

                metcorUncDown = TMath::Sqrt(metcorrUncDown_ex * metcorrUncDown_ex + metcorrUncDown_ey * metcorrUncDown_ey);
                metcorphiUncDown = TMath::ATan2(metcorrUncDown_ey, metcorrUncDown_ex);

                metcorUncUp = TMath::Sqrt(metcorrUncUp_ex * metcorrUncUp_ex + metcorrUncUp_ey * metcorrUncUp_ey);
                metcorphiUncUp = TMath::ATan2(metcorrUncUp_ey, metcorrUncUp_ex);

                metcorJetEC2Down = TMath::Sqrt(metcorrJetEC2Down_ex * metcorrJetEC2Down_ex + metcorrJetEC2Down_ey * metcorrJetEC2Down_ey);
                metcorphiJetEC2Down = TMath::ATan2(metcorrJetEC2Down_ey, metcorrJetEC2Down_ex);

                metcorJetEC2Up = TMath::Sqrt(metcorrJetEC2Up_ex * metcorrJetEC2Up_ex + metcorrJetEC2Up_ey * metcorrJetEC2Up_ey);
                metcorphiJetEC2Up = TMath::ATan2(metcorrJetEC2Up_ey, metcorrJetEC2Up_ex);

                metcorJetEta0to3Down =
                    TMath::Sqrt(metcorrJetEta0to3Down_ex * metcorrJetEta0to3Down_ex + metcorrJetEta0to3Down_ey * metcorrJetEta0to3Down_ey);
                metcorphiJetEta0to3Down = TMath::ATan2(metcorrJetEta0to3Down_ey, metcorrJetEta0to3Down_ex);

                metcorJetEta0to3Up = TMath::Sqrt(metcorrJetEta0to3Up_ex * metcorrJetEta0to3Up_ex + metcorrJetEta0to3Up_ey * metcorrJetEta0to3Up_ey);
                metcorphiJetEta0to3Up = TMath::ATan2(metcorrJetEta0to3Up_ey, metcorrJetEta0to3Up_ex);

                metcorJetEta0to5Down =
                    TMath::Sqrt(metcorrJetEta0to5Down_ex * metcorrJetEta0to5Down_ex + metcorrJetEta0to5Down_ey * metcorrJetEta0to5Down_ey);
                metcorphiJetEta0to5Down = TMath::ATan2(metcorrJetEta0to5Down_ey, metcorrJetEta0to5Down_ex);

                metcorJetEta0to5Up = TMath::Sqrt(metcorrJetEta0to5Up_ex * metcorrJetEta0to5Up_ex + metcorrJetEta0to5Up_ey * metcorrJetEta0to5Up_ey);
                metcorphiJetEta0to5Up = TMath::ATan2(metcorrJetEta0to5Up_ey, metcorrJetEta0to5Up_ex);

                metcorJetEta3to5Down =
                    TMath::Sqrt(metcorrJetEta3to5Down_ex * metcorrJetEta3to5Down_ex + metcorrJetEta3to5Down_ey * metcorrJetEta3to5Down_ey);
                metcorphiJetEta3to5Down = TMath::ATan2(metcorrJetEta3to5Down_ey, metcorrJetEta3to5Down_ex);

                metcorJetEta3to5Up = TMath::Sqrt(metcorrJetEta3to5Up_ex * metcorrJetEta3to5Up_ex + metcorrJetEta3to5Up_ey * metcorrJetEta3to5Up_ey);
                metcorphiJetEta3to5Up = TMath::ATan2(metcorrJetEta3to5Up_ey, metcorrJetEta3to5Up_ex);

                metcorJetRelativeSampleDown = TMath::Sqrt(metcorrJetRelativeSampleDown_ex * metcorrJetRelativeSampleDown_ex +
                                                          metcorrJetRelativeSampleDown_ey * metcorrJetRelativeSampleDown_ey);
                metcorphiJetRelativeSampleDown = TMath::ATan2(metcorrJetRelativeSampleDown_ey, metcorrJetRelativeSampleDown_ex);

                metcorJetRelativeSampleUp = TMath::Sqrt(metcorrJetRelativeSampleUp_ex * metcorrJetRelativeSampleUp_ex +
                                                        metcorrJetRelativeSampleUp_ey * metcorrJetRelativeSampleUp_ey);
                metcorphiJetRelativeSampleUp = TMath::ATan2(metcorrJetRelativeSampleUp_ey, metcorrJetRelativeSampleUp_ex);

                metcorJetRelativeSampleDown = TMath::Sqrt(metcorrJetRelativeSampleDown_ex * metcorrJetRelativeSampleDown_ex +
                                                          metcorrJetRelativeSampleDown_ey * metcorrJetRelativeSampleDown_ey);
                metcorphiJetRelativeSampleDown = TMath::ATan2(metcorrJetRelativeSampleDown_ey, metcorrJetRelativeSampleDown_ex);

                metcorJetRelativeSampleUp = TMath::Sqrt(metcorrJetRelativeSampleUp_ex * metcorrJetRelativeSampleUp_ex +
                                                        metcorrJetRelativeSampleUp_ey * metcorrJetRelativeSampleUp_ey);
                metcorphiJetRelativeSampleUp = TMath::ATan2(metcorrJetRelativeSampleUp_ey, metcorrJetRelativeSampleUp_ex);

                if (channel == "mt" || channel == "et") {
                    mass2 = m2;
                    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons{
                        classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1, phi1, mass1),
                        classic_svFit::MeasuredTauLepton(decayType2, pt2, eta2, phi2, mass2, decayMode2)};

                    std::cout << "era: " << era << "evt: " << evt << " run: " << run << " lumi: " << lumi << " pt1 " << pt1 << " mass1 " << mass1
                              << " pt2: " << pt2 << " mass2: " << mass2 << std::endl;

                    runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET,
                             svFitTransverseMass, tau1, tau2);

                    // MET systematics
                    runSVFit(measuredTauLeptons, metcorrUncUp_ex, metcorrUncUp_ey, covMET, 0, svFitMass_UncMet_Up, svFitPt_UncMet_Up,
                             svFitEta_UncMet_Up, svFitPhi_UncMet_Up, svFitMET_UncMet_Up, svFitTransverseMass_UncMet_Up, tau1_UncMet_Up,
                             tau2_UncMet_Up);

                    runSVFit(measuredTauLeptons, metcorrUncDown_ex, metcorrUncDown_ey, covMET, 0, svFitMass_UncMet_Down, svFitPt_UncMet_Down,
                             svFitEta_UncMet_Down, svFitPhi_UncMet_Down, svFitMET_UncMet_Down, svFitTransverseMass_UncMet_Down, tau1_UncMet_Down,
                             tau2_UncMet_Down);

                    runSVFit(measuredTauLeptons, metcorrClusteredUp_ex, metcorrClusteredUp_ey, covMET, 0, svFitMass_ClusteredMet_Up,
                             svFitPt_ClusteredMet_Up, svFitEta_ClusteredMet_Up, svFitPhi_ClusteredMet_Up, svFitMET_ClusteredMet_Up,
                             svFitTransverseMass_ClusteredMet_Up, tau1_ClusteredMet_Up, tau2_ClusteredMet_Up);

                    runSVFit(measuredTauLeptons, metcorrClusteredDown_ex, metcorrClusteredDown_ey, covMET, 0, svFitMass_ClusteredMet_Down,
                             svFitPt_ClusteredMet_Down, svFitEta_ClusteredMet_Down, svFitPhi_ClusteredMet_Down, svFitMET_ClusteredMet_Down,
                             svFitTransverseMass_ClusteredMet_Down, tau1_ClusteredMet_Down, tau2_ClusteredMet_Down);

                    runSVFit(measuredTauLeptons, metcorrJetEC2Up_ex, metcorrJetEC2Up_ey, covMET, 0, svFitMass_JetEC2_Up, svFitPt_JetEC2_Up,
                             svFitEta_JetEC2_Up, svFitPhi_JetEC2_Up, svFitMET_JetEC2_Up, svFitTransverseMass_JetEC2_Up, tau1_JetEC2_Up,
                             tau2_JetEC2_Up);

                    runSVFit(measuredTauLeptons, metcorrJetEC2Down_ex, metcorrJetEC2Down_ey, covMET, 0, svFitMass_JetEC2_Down, svFitPt_JetEC2_Down,
                             svFitEta_JetEC2_Down, svFitPhi_JetEC2_Down, svFitMET_JetEC2_Down, svFitTransverseMass_JetEC2_Down, tau1_JetEC2_Down,
                             tau2_JetEC2_Down);

                    runSVFit(measuredTauLeptons, metcorrJetEta0to3Up_ex, metcorrJetEta0to3Up_ey, covMET, 0, svFitMass_JetEta0to3_Up,
                             svFitPt_JetEta0to3_Up, svFitEta_JetEta0to3_Up, svFitPhi_JetEta0to3_Up, svFitMET_JetEta0to3_Up,
                             svFitTransverseMass_JetEta0to3_Up, tau1_JetEta0to3_Up, tau2_JetEta0to3_Up);

                    runSVFit(measuredTauLeptons, metcorrJetEta0to3Down_ex, metcorrJetEta0to3Down_ey, covMET, 0, svFitMass_JetEta0to3_Down,
                             svFitPt_JetEta0to3_Down, svFitEta_JetEta0to3_Down, svFitPhi_JetEta0to3_Down, svFitMET_JetEta0to3_Down,
                             svFitTransverseMass_JetEta0to3_Down, tau1_JetEta0to3_Down, tau2_JetEta0to3_Down);

                    runSVFit(measuredTauLeptons, metcorrJetEta0to5Up_ex, metcorrJetEta0to5Up_ey, covMET, 0, svFitMass_JetEta0to5_Up,
                             svFitPt_JetEta0to5_Up, svFitEta_JetEta0to5_Up, svFitPhi_JetEta0to5_Up, svFitMET_JetEta0to5_Up,
                             svFitTransverseMass_JetEta0to5_Up, tau1_JetEta0to5_Up, tau2_JetEta0to5_Up);

                    runSVFit(measuredTauLeptons, metcorrJetEta0to5Down_ex, metcorrJetEta0to5Down_ey, covMET, 0, svFitMass_JetEta0to5_Down,
                             svFitPt_JetEta0to5_Down, svFitEta_JetEta0to5_Down, svFitPhi_JetEta0to5_Down, svFitMET_JetEta0to5_Down,
                             svFitTransverseMass_JetEta0to5_Down, tau1_JetEta0to5_Down, tau2_JetEta0to5_Down);

                    runSVFit(measuredTauLeptons, metcorrJetEta3to5Up_ex, metcorrJetEta3to5Up_ey, covMET, 0, svFitMass_JetEta3to5_Up,
                             svFitPt_JetEta3to5_Up, svFitEta_JetEta3to5_Up, svFitPhi_JetEta3to5_Up, svFitMET_JetEta3to5_Up,
                             svFitTransverseMass_JetEta3to5_Up, tau1_JetEta3to5_Up, tau2_JetEta3to5_Up);

                    runSVFit(measuredTauLeptons, metcorrJetEta3to5Down_ex, metcorrJetEta3to5Down_ey, covMET, 0, svFitMass_JetEta3to5_Down,
                             svFitPt_JetEta3to5_Down, svFitEta_JetEta3to5_Down, svFitPhi_JetEta3to5_Down, svFitMET_JetEta3to5_Down,
                             svFitTransverseMass_JetEta3to5_Down, tau1_JetEta3to5_Down, tau2_JetEta3to5_Down);

                    runSVFit(measuredTauLeptons, metcorrJetRelativeBalUp_ex, metcorrJetRelativeBalUp_ey, covMET, 0, svFitMass_JetRelativeBal_Up,
                             svFitPt_JetRelativeBal_Up, svFitEta_JetRelativeBal_Up, svFitPhi_JetRelativeBal_Up, svFitMET_JetRelativeBal_Up,
                             svFitTransverseMass_JetRelativeBal_Up, tau1_JetRelativeBal_Up, tau2_JetRelativeBal_Up);

                    runSVFit(measuredTauLeptons, metcorrJetRelativeBalDown_ex, metcorrJetRelativeBalDown_ey, covMET, 0, svFitMass_JetRelativeBal_Down,
                             svFitPt_JetRelativeBal_Down, svFitEta_JetRelativeBal_Down, svFitPhi_JetRelativeBal_Down, svFitMET_JetRelativeBal_Down,
                             svFitTransverseMass_JetRelativeBal_Down, tau1_JetRelativeBal_Down, tau2_JetRelativeBal_Down);

                    runSVFit(measuredTauLeptons, metcorrJetRelativeSampleUp_ex, metcorrJetRelativeSampleUp_ey, covMET, 0,
                             svFitMass_JetRelativeSample_Up, svFitPt_JetRelativeSample_Up, svFitEta_JetRelativeSample_Up,
                             svFitPhi_JetRelativeSample_Up, svFitMET_JetRelativeSample_Up, svFitTransverseMass_JetRelativeSample_Up,
                             tau1_JetRelativeSample_Up, tau2_JetRelativeSample_Up);

                    runSVFit(measuredTauLeptons, metcorrJetRelativeSampleDown_ex, metcorrJetRelativeSampleDown_ey, covMET, 0,
                             svFitMass_JetRelativeSample_Down, svFitPt_JetRelativeSample_Down, svFitEta_JetRelativeSample_Down,
                             svFitPhi_JetRelativeSample_Down, svFitMET_JetRelativeSample_Down, svFitTransverseMass_JetRelativeSample_Down,
                             tau1_JetRelativeSample_Down, tau2_JetRelativeSample_Down);

                    if (doES) {
                        // corrections only need to be done once
                        float ES_Up(1.), ES_Down(1.);  // shift TES
                        if (gen_match_2 == 5) {        // 0.6% uncertainty on hadronic tau
                            ES_Up = 1 + tesUncertainties(era, decayMode2);
                            ES_Down = 1 - tesUncertainties(era, decayMode2);
                        } else if (gen_match_2 < 5) {  // flat 3% on el/mu -> tau energy scale systematics
                            ES_Up = 1.03;
                            ES_Down = 0.97;
                        }

                        double pt_Up(pt2 * ES_Up), pt_Down(pt2 * ES_Down);  // shift tau pT by energy scale
                        double dx_Up(pt2 * TMath::Cos(phi2) * ((1. / ES_Up) - 1.)), dy_Up(pt2 * TMath::Sin(phi2) * ((1. / ES_Up) - 1.)),
                            dx_Down(pt2 * TMath::Cos(phi2) * ((1. / ES_Down) - 1.)), dy_Down(pt2 * TMath::Sin(phi2) * ((1. / ES_Down) - 1.));
                        double metcorr_ex_Up(metcorr_ex + dx_Up), metcorr_ey_Up(metcorr_ey + dy_Up), metcorr_ex_Down(metcorr_ex + dx_Down),
                            metcorr_ey_Down(metcorr_ey + dy_Down);

                        // leptons shifted up
                        std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUp{
                            classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1, phi1, mass1),
                            classic_svFit::MeasuredTauLepton(decayType2, pt_Up, eta2, phi2, mass2, decayMode2)};

                        // leptons shifted down
                        std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDown{
                            classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1, phi1, mass1),
                            classic_svFit::MeasuredTauLepton(decayType2, pt_Down, eta2, phi2, mass2, decayMode2)};

                        /////////////////////////////
                        // All upward shifts below //
                        /////////////////////////////

                        // all tau shift up
                        if (gen_match_2 < 6) {
                            runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_Up, svFitPt_Up, svFitEta_Up,
                                     svFitPhi_Up, svFitMET_Up, svFitTransverseMass_Up, tau1_Up, tau2_Up);
                        } else {
                            svFitMass_Up = svFitMass;
                            svFitPt_Up = svFitPt;
                            svFitEta_Up = svFitEta;
                            svFitPhi_Up = svFitPhi;
                            svFitMET_Up = svFitMET;
                            svFitTransverseMass_Up = svFitTransverseMass;
                            tau1_Up = tau1;
                            tau2_Up = tau2;
                        }

                        // tau DM0 shifted up
                        if (gen_match_2 == 5 && decayMode2 == 0) {
                            runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM0_Up, svFitPt_DM0_Up, svFitEta_DM0_Up,
                                     svFitPhi_DM0_Up, svFitMET_DM0_Up, svFitTransverseMass_DM0_Up, tau1_DM0_Up, tau2_DM0_Up);

                        } else {
                            svFitMass_DM0_Up = svFitMass;
                            svFitPt_DM0_Up = svFitPt;
                            svFitEta_DM0_Up = svFitEta;
                            svFitPhi_DM0_Up = svFitPhi;
                            svFitMET_DM0_Up = svFitMET;
                            svFitTransverseMass_DM0_Up = svFitTransverseMass;
                            tau1_DM0_Up = tau1;
                            tau2_DM0_Up = tau2;
                        }

                        // tau DM1 shifted up
                        if (gen_match_2 == 5 && decayMode2 == 1) {
                            runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM1_Up, svFitPt_DM1_Up, svFitEta_DM1_Up,
                                     svFitPhi_DM1_Up, svFitMET_DM1_Up, svFitTransverseMass_DM1_Up, tau1_DM1_Up, tau2_DM1_Up);
                        } else {
                            svFitMass_DM1_Up = svFitMass;
                            svFitPt_DM1_Up = svFitPt;
                            svFitEta_DM1_Up = svFitEta;
                            svFitPhi_DM1_Up = svFitPhi;
                            svFitMET_DM1_Up = svFitMET;
                            svFitTransverseMass_DM1_Up = svFitTransverseMass;
                            tau1_DM1_Up = tau1;
                            tau2_DM1_Up = tau2;
                        }

                        // tau DM10 shifted up
                        if (gen_match_2 == 5 && decayMode2 == 10) {
                            runSVFit(measuredTauLeptonsUp, metcorr_ex_Up, metcorr_ey_Up, covMET, 0, svFitMass_DM10_Up, svFitPt_DM10_Up,
                                     svFitEta_DM10_Up, svFitPhi_DM10_Up, svFitMET_DM10_Up, svFitTransverseMass_DM10_Up, tau1_DM10_Up, tau2_DM10_Up);
                        } else {
                            svFitMass_DM10_Up = svFitMass;
                            svFitPt_DM10_Up = svFitPt;
                            svFitEta_DM10_Up = svFitEta;
                            svFitPhi_DM10_Up = svFitPhi;
                            svFitMET_DM10_Up = svFitMET;
                            svFitTransverseMass_DM10_Up = svFitTransverseMass;
                            tau1_DM10_Up = tau1;
                            tau2_DM10_Up = tau2;
                        }

                        ///////////////////////////////
                        // All downward shifts below //
                        ///////////////////////////////

                        // all tau shift down
                        if (gen_match_2 < 6) {
                            runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_Down, svFitPt_Down, svFitEta_Down,
                                     svFitPhi_Down, svFitMET_Down, svFitTransverseMass_Down, tau1_Down, tau2_Down);
                        } else {
                            svFitMass_Down = svFitMass;
                            svFitPt_Down = svFitPt;
                            svFitEta_Down = svFitEta;
                            svFitPhi_Down = svFitPhi;
                            svFitMET_Down = svFitMET;
                            svFitTransverseMass_Down = svFitTransverseMass;
                            tau1_Down = tau1;
                            tau2_Down = tau2;
                        }

                        // tau DM0 shifted down
                        if (gen_match_2 == 5 && decayMode2 == 0) {
                            std::cout << "DM0 shift Down" << std::endl;
                            runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM0_Down, svFitPt_DM0_Down,
                                     svFitEta_DM0_Down, svFitPhi_DM0_Down, svFitMET_DM0_Down, svFitTransverseMass_DM0_Down, tau1_DM0_Down,
                                     tau2_DM0_Down);
                        } else {
                            svFitMass_DM0_Down = svFitMass;
                            svFitPt_DM0_Down = svFitPt;
                            svFitEta_DM0_Down = svFitEta;
                            svFitPhi_DM0_Down = svFitPhi;
                            svFitMET_DM0_Down = svFitMET;
                            svFitTransverseMass_DM0_Down = svFitTransverseMass;
                            tau1_DM0_Down = tau1;
                            tau2_DM0_Down = tau2;
                        }

                        // tau DM1 shifted down
                        if (gen_match_2 == 5 && decayMode2 == 1) {
                            std::cout << "DM1 shift down" << std::endl;
                            runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM1_Down, svFitPt_DM1_Down,
                                     svFitEta_DM1_Down, svFitPhi_DM1_Down, svFitMET_DM1_Down, svFitTransverseMass_DM1_Down, tau1_DM1_Down,
                                     tau2_DM1_Down);
                        } else {
                            svFitMass_DM1_Down = svFitMass;
                            svFitPt_DM1_Down = svFitPt;
                            svFitEta_DM1_Down = svFitEta;
                            svFitPhi_DM1_Down = svFitPhi;
                            svFitMET_DM1_Down = svFitMET;
                            svFitTransverseMass_DM1_Down = svFitTransverseMass;
                            tau1_DM1_Down = tau1;
                            tau2_DM1_Down = tau2;
                        }

                        // tau DM10 shifted down
                        if (gen_match_2 == 5 && decayMode2 == 10) {
                            std::cout << "DM10 shift down" << std::endl;
                            runSVFit(measuredTauLeptonsDown, metcorr_ex_Down, metcorr_ey_Down, covMET, 0, svFitMass_DM10_Down, svFitPt_DM10_Down,
                                     svFitEta_DM10_Down, svFitPhi_DM10_Down, svFitMET_DM10_Down, svFitTransverseMass_DM10_Down, tau1_DM10_Down,
                                     tau2_DM10_Down);
                        } else {
                            svFitMass_DM10_Down = svFitMass;
                            svFitPt_DM10_Down = svFitPt;
                            svFitEta_DM10_Down = svFitEta;
                            svFitPhi_DM10_Down = svFitPhi;
                            svFitMET_DM10_Down = svFitMET;
                            svFitTransverseMass_DM10_Down = svFitTransverseMass;
                            tau1_DM10_Down = tau1;
                            tau2_DM10_Down = tau2;
                        }
                    }  // end doES
                } else {
                    svFitMass = -100;
                    svFitPt = -100;
                    svFitEta = -100;
                    svFitPhi = -100;
                    svFitMET = -100;
                    svFitTransverseMass = -100;
                    tau1.SetPtEtaPhiM(0, 0, 0, 0);
                    tau2.SetPtEtaPhiM(0, 0, 0, 0);
                }
                // fill the tau 4-vector parameters for branch-filling
                four_vector(tau1, tau1_pt, tau1_eta, tau1_phi, tau1_m);
                four_vector(tau2, tau2_pt, tau2_eta, tau2_phi, tau2_m);

                four_vector(tau1_Up, tau1_pt_Up, tau1_eta_Up, tau1_phi_Up, tau1_m_Up);
                four_vector(tau2_Up, tau2_pt_Up, tau2_eta_Up, tau2_phi_Up, tau2_m_Up);
                four_vector(tau1_Down, tau1_pt_Down, tau1_eta_Down, tau1_phi_Down, tau1_m_Down);
                four_vector(tau2_Down, tau2_pt_Down, tau2_eta_Down, tau2_phi_Down, tau2_m_Down);

                four_vector(tau1_DM0_Up, tau1_pt_DM0_Up, tau1_eta_DM0_Up, tau1_phi_DM0_Up, tau1_m_DM0_Up);
                four_vector(tau2_DM0_Up, tau2_pt_DM0_Up, tau2_eta_DM0_Up, tau2_phi_DM0_Up, tau2_m_DM0_Up);
                four_vector(tau1_DM0_Down, tau1_pt_DM0_Down, tau1_eta_DM0_Down, tau1_phi_DM0_Down, tau1_m_DM0_Down);
                four_vector(tau2_DM0_Down, tau2_pt_DM0_Down, tau2_eta_DM0_Down, tau2_phi_DM0_Down, tau2_m_DM0_Down);

                four_vector(tau1_DM1_Up, tau1_pt_DM1_Up, tau1_eta_DM1_Up, tau1_phi_DM1_Up, tau1_m_DM1_Up);
                four_vector(tau2_DM1_Up, tau2_pt_DM1_Up, tau2_eta_DM1_Up, tau2_phi_DM1_Up, tau2_m_DM1_Up);
                four_vector(tau1_DM1_Down, tau1_pt_DM1_Down, tau1_eta_DM1_Down, tau1_phi_DM1_Down, tau1_m_DM1_Down);
                four_vector(tau2_DM1_Down, tau2_pt_DM1_Down, tau2_eta_DM1_Down, tau2_phi_DM1_Down, tau2_m_DM1_Down);

                four_vector(tau1_DM10_Up, tau1_pt_DM10_Up, tau1_eta_DM10_Up, tau1_phi_DM10_Up, tau1_m_DM10_Up);
                four_vector(tau2_DM10_Up, tau2_pt_DM10_Up, tau2_eta_DM10_Up, tau2_phi_DM10_Up, tau2_m_DM10_Up);
                four_vector(tau1_DM10_Down, tau1_pt_DM10_Down, tau1_eta_DM10_Down, tau1_phi_DM10_Down, tau1_m_DM10_Down);
                four_vector(tau2_DM10_Down, tau2_pt_DM10_Down, tau2_eta_DM10_Down, tau2_phi_DM10_Down, tau2_m_DM10_Down);

                four_vector(tau1_UncMet_Up, tau1_pt_UncMet_Up, tau1_eta_UncMet_Up, tau1_phi_UncMet_Up, tau1_m_UncMet_Up);
                four_vector(tau2_UncMet_Up, tau2_pt_UncMet_Up, tau2_eta_UncMet_Up, tau2_phi_UncMet_Up, tau2_m_UncMet_Up);
                four_vector(tau1_UncMet_Down, tau1_pt_UncMet_Down, tau1_eta_UncMet_Down, tau1_phi_UncMet_Down, tau1_m_UncMet_Down);
                four_vector(tau2_UncMet_Down, tau2_pt_UncMet_Down, tau2_eta_UncMet_Down, tau2_phi_UncMet_Down, tau2_m_UncMet_Down);

                four_vector(tau1_ClusteredMet_Up, tau1_pt_ClusteredMet_Up, tau1_eta_ClusteredMet_Up, tau1_phi_ClusteredMet_Up, tau1_m_ClusteredMet_Up);
                four_vector(tau2_ClusteredMet_Up, tau2_pt_ClusteredMet_Up, tau2_eta_ClusteredMet_Up, tau2_phi_ClusteredMet_Up, tau2_m_ClusteredMet_Up);
                four_vector(tau1_ClusteredMet_Down, tau1_pt_ClusteredMet_Down, tau1_eta_ClusteredMet_Down, tau1_phi_ClusteredMet_Down, tau1_m_ClusteredMet_Down);
                four_vector(tau2_ClusteredMet_Down, tau2_pt_ClusteredMet_Down, tau2_eta_ClusteredMet_Down, tau2_phi_ClusteredMet_Down, tau2_m_ClusteredMet_Down);

                four_vector(tau1_JetEC2_Up, tau1_pt_JetEC2_Up, tau1_eta_JetEC2_Up, tau1_phi_JetEC2_Up, tau1_m_JetEC2_Up);
                four_vector(tau2_JetEC2_Up, tau2_pt_JetEC2_Up, tau2_eta_JetEC2_Up, tau2_phi_JetEC2_Up, tau2_m_JetEC2_Up);
                four_vector(tau1_JetEC2_Down, tau1_pt_JetEC2_Down, tau1_eta_JetEC2_Down, tau1_phi_JetEC2_Down, tau1_m_JetEC2_Down);
                four_vector(tau2_JetEC2_Down, tau2_pt_JetEC2_Down, tau2_eta_JetEC2_Down, tau2_phi_JetEC2_Down, tau2_m_JetEC2_Down);

                four_vector(tau1_JetEta0to3_Up, tau1_pt_JetEta0to3_Up, tau1_eta_JetEta0to3_Up, tau1_phi_JetEta0to3_Up, tau1_m_JetEta0to3_Up);
                four_vector(tau2_JetEta0to3_Up, tau2_pt_JetEta0to3_Up, tau2_eta_JetEta0to3_Up, tau2_phi_JetEta0to3_Up, tau2_m_JetEta0to3_Up);
                four_vector(tau1_JetEta0to3_Down, tau1_pt_JetEta0to3_Down, tau1_eta_JetEta0to3_Down, tau1_phi_JetEta0to3_Down, tau1_m_JetEta0to3_Down);
                four_vector(tau2_JetEta0to3_Down, tau2_pt_JetEta0to3_Down, tau2_eta_JetEta0to3_Down, tau2_phi_JetEta0to3_Down, tau2_m_JetEta0to3_Down);

                four_vector(tau1_JetEta0to5_Up, tau1_pt_JetEta0to5_Up, tau1_eta_JetEta0to5_Up, tau1_phi_JetEta0to5_Up, tau1_m_JetEta0to5_Up);
                four_vector(tau2_JetEta0to5_Up, tau2_pt_JetEta0to5_Up, tau2_eta_JetEta0to5_Up, tau2_phi_JetEta0to5_Up, tau2_m_JetEta0to5_Up);
                four_vector(tau1_JetEta0to5_Down, tau1_pt_JetEta0to5_Down, tau1_eta_JetEta0to5_Down, tau1_phi_JetEta0to5_Down, tau1_m_JetEta0to5_Down);
                four_vector(tau2_JetEta0to5_Down, tau2_pt_JetEta0to5_Down, tau2_eta_JetEta0to5_Down, tau2_phi_JetEta0to5_Down, tau2_m_JetEta0to5_Down);

                four_vector(tau1_JetEta3to5_Up, tau1_pt_JetEta3to5_Up, tau1_eta_JetEta3to5_Up, tau1_phi_JetEta3to5_Up, tau1_m_JetEta3to5_Up);
                four_vector(tau2_JetEta3to5_Up, tau2_pt_JetEta3to5_Up, tau2_eta_JetEta3to5_Up, tau2_phi_JetEta3to5_Up, tau2_m_JetEta3to5_Up);
                four_vector(tau1_JetEta3to5_Down, tau1_pt_JetEta3to5_Down, tau1_eta_JetEta3to5_Down, tau1_phi_JetEta3to5_Down, tau1_m_JetEta3to5_Down);
                four_vector(tau2_JetEta3to5_Down, tau2_pt_JetEta3to5_Down, tau2_eta_JetEta3to5_Down, tau2_phi_JetEta3to5_Down, tau2_m_JetEta3to5_Down);

                four_vector(tau1_JetRelativeBal_Up, tau1_pt_JetRelativeBal_Up, tau1_eta_JetRelativeBal_Up, tau1_phi_JetRelativeBal_Up, tau1_m_JetRelativeBal_Up);
                four_vector(tau2_JetRelativeBal_Up, tau2_pt_JetRelativeBal_Up, tau2_eta_JetRelativeBal_Up, tau2_phi_JetRelativeBal_Up, tau2_m_JetRelativeBal_Up);
                four_vector(tau1_JetRelativeBal_Down, tau1_pt_JetRelativeBal_Down, tau1_eta_JetRelativeBal_Down, tau1_phi_JetRelativeBal_Down, tau1_m_JetRelativeBal_Down);
                four_vector(tau2_JetRelativeBal_Down, tau2_pt_JetRelativeBal_Down, tau2_eta_JetRelativeBal_Down, tau2_phi_JetRelativeBal_Down, tau2_m_JetRelativeBal_Down);

                four_vector(tau1_JetRelativeSample_Up, tau1_pt_JetRelativeSample_Up, tau1_eta_JetRelativeSample_Up, tau1_phi_JetRelativeSample_Up, tau1_m_JetRelativeSample_Up);
                four_vector(tau2_JetRelativeSample_Up, tau2_pt_JetRelativeSample_Up, tau2_eta_JetRelativeSample_Up, tau2_phi_JetRelativeSample_Up, tau2_m_JetRelativeSample_Up);
                four_vector(tau1_JetRelativeSample_Down, tau1_pt_JetRelativeSample_Down, tau1_eta_JetRelativeSample_Down, tau1_phi_JetRelativeSample_Down, tau1_m_JetRelativeSample_Down);
                four_vector(tau2_JetRelativeSample_Down, tau2_pt_JetRelativeSample_Down, tau2_eta_JetRelativeSample_Down, tau2_phi_JetRelativeSample_Down, tau2_m_JetRelativeSample_Down);

                std::cout << "\n\n" << std::endl;
                // std::cout << "\n\nex: " << metcorr_ex << "   ey: " << metcorr_ey <<  " phi: " << metcorphi<<"\n"<<std::endl;
                for (auto branch : fit_results) {
                    branch->Fill();
                }
            }
            dir->cd();
            t->Write("", TObject::kOverwrite);
            delete t;
        }  // if the iterator of the key is a TTree
    }
}

void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> &measuredTauLeptons, double measuredMETx, double measuredMETy, TMatrixD &covMET,
              float num, float &svFitMass, float &svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass,
              TLorentzVector &tau1, TLorentzVector &tau2) {
    svfitAlgorithm.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
    if (svfitAlgorithm.isValidSolution()) {
        svFitMass = static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->getMass();
        svFitPt = static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->getPt();
        svFitEta = static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->getEta();
        svFitPhi = static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->getPhi();
        svFitTransverseMass = static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->getTransverseMass();

        classic_svFit::HistogramAdapterTau *h_tau1 =
            static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->tau1();
        classic_svFit::HistogramAdapterTau *h_tau2 =
            static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->tau2();
        const classic_svFit::LorentzVector tau1_p4 = h_tau1->getP4();
        const classic_svFit::LorentzVector tau2_p4 = h_tau2->getP4();
        tau1.SetPtEtaPhiM(tau1_p4.Pt(), tau1_p4.Eta(), tau1_p4.Phi(), tau1_p4.M());
        tau2.SetPtEtaPhiM(tau2_p4.Pt(), tau2_p4.Eta(), tau2_p4.Phi(), tau2_p4.M());

        const classic_svFit::LorentzVector ditau = static_cast<classic_svFit::HistogramAdapterDiTau *>(svfitAlgorithm.getHistogramAdapter())->getP4();
        svFitMET = (ditau - (measuredTauLeptons[0].p4() + measuredTauLeptons[1].p4())).Pt();

        // TLorentzVector testTau1, testTau2;
        // testTau1.SetPtEtaPhiM(measuredTauLeptons[0].p4().Pt(), measuredTauLeptons[0].p4().Eta(), measuredTauLeptons[0].p4().Phi(),
        // measuredTauLeptons[0].p4().M()); testTau2.SetPtEtaPhiM(measuredTauLeptons[1].p4().Pt(), measuredTauLeptons[1].p4().Eta(),
        // measuredTauLeptons[1].p4().Phi(), measuredTauLeptons[1].p4().M());

        // TLorentzVector ditau;
        // ditau.SetPtEtaPhiM(svFitPt, svFitEta, svFitPhi, svFitMass);
    }
}

// Thank you Renee Brun :)
void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
    // copy all objects and subdirs of directory source as a subdir of the current directory
    TDirectory *savdir = gDirectory;
    TDirectory *adir = savdir;
    if (source->GetName() != parser.stringValue("inputFile")) {
        adir = savdir->mkdir(source->GetName());
        std::cout << "Source name is not outputfile name" << std::endl;
        adir->cd();
    } else {
        // adir = savdir->mkdir("input");
        adir->cd();
    }

    // loop on all entries of this directory
    TKey *key;
    TIter nextkey(source->GetListOfKeys());
    while ((key = (TKey *)nextkey())) {
        const char *classname = key->GetClassName();
        TClass *cl = gROOT->GetClass(classname);
        if (!cl) continue;
        if (cl->InheritsFrom(TDirectory::Class())) {
            source->cd(key->GetName());
            TDirectory *subdir = gDirectory;
            adir->cd();
            CopyDir(subdir, parser);
            adir->cd();
        } else if (cl->InheritsFrom(TTree::Class())) {
            TTree *T = (TTree *)source->Get(key->GetName());
            adir->cd();
            TTree *newT = T->CloneTree(-1);
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
    // Copy all objects and subdirs of file fname as a subdir of the current directory
    TDirectory *target = gDirectory;
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
        printf("Cannot copy file: %s\n", fname);
        target->cd();
        return 0;
    }
    target->cd();
    CopyDir(f, parser);
    delete f;
    target->cd();
    return 1;
}
int copyFiles(optutl::CommandLineParser parser, TFile *fOld, TFile *fNew) {
    // prepare files to be copied
    if (gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
        gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
    }

    fNew->cd();
    if (CopyFile(parser.stringValue("inputFile").c_str(), parser) == 0) return 0;
    fNew->ls();
    fNew->Close();
    return 1;
}

double tesUncertainties(unsigned int year, float decaymode) {
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingLegacyRun2#Tau_energy_scale_uncertainty
    double tesSize = -1000;
    if (year == 2016) {
        if (decaymode == 0)
            tesSize = 0.010;
        else if (decaymode == 1)
            tesSize = 0.009;
        else if (decaymode == 10)
            tesSize = 0.011;
    }
    if (year == 2017) {
        if (decaymode == 0)
            tesSize = 0.008;
        else if (decaymode == 1)
            tesSize = 0.008;
        else if (decaymode == 10)
            tesSize = 0.009;
    }
    if (year == 2018) {
        if (decaymode == 0)
            tesSize = 0.011;
        else if (decaymode == 1)
            tesSize = 0.008;
        else if (decaymode == 10)
            tesSize = 0.009;
    }

    return tesSize;
}

void four_vector(TLorentzVector p4, Float_t& pt, Float_t& eta, Float_t& phi, Float_t& mass) {
    pt = p4.Pt();
    eta = p4.Eta();
    phi = p4.Phi();
    mass = p4.M();
}