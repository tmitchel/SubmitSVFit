#include <map>
#include <string>

class scenario {
 private:
  // Input Objects
  Float_t pt_1, pt_2, met, metphi;

 public:
  scenario (TTree*, std::string);
  virtual ~scenario () {};
  Float_t get_pt_1()          { return pt_1;        };
  Float_t get_pt_2()          { return pt_2;        };
  Float_t get_met()           { return met;          };
  Float_t get_metphi()        { return metphi;       };

  // Output Objects
  TLorentzVector tau1, tau2;
  Float_t svFitMass, svFitPt, svFitEta, svFitPhi, svFitTransverseMass;
  Float_t svFitMET, metcorr_ex, metcorr_ey, metcor, metcorphi;
  Float_t tau1_pt, tau1_eta, tau1_phi, tau1_m;
  Float_t tau2_pt, tau2_eta, tau2_phi, tau2_m;
};

std::string Replace(std::string &str, const std::string& from, const std::string& to){
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();
  }
  return str;
}      


// read branches from namu
scenario::scenario(TTree* namu, std::string unc) {
  // set nominal as a defalt
  std::string pt_1_unc = "pt_1";
  std::string pt_2_unc = "pt_2";
  std::string met_unc = "met";
  std::string metphi_unc = "metphi";
  std::cout << "-------------- " << unc << " --------------" << std::endl;

  // if branches for unc is availble, replace to them
  if (unc=="nominal") unc = "";
  else unc = "_"+unc;
  // JES Branch names are different 
  std::string uncReplaced = unc;
  if (unc.find("Jet") != std::string::npos || unc.find("JER") != std::string::npos) {
    if (unc.find("_Up") != std::string::npos)  Replace(uncReplaced,"_Up","Up");
    else if (unc.find("Down") != std::string::npos) Replace(uncReplaced,"_Down","Down");
  }
  // UncMet
  

  TIter next(namu->GetListOfBranches());
  TBranch *branch;
  while ((branch = (TBranch*)next())) {
    if ((pt_1_unc+uncReplaced)==branch->GetName())        pt_1_unc=branch->GetName();
    else if ((pt_2_unc+uncReplaced)==branch->GetName())   pt_2_unc=branch->GetName();
    else if ((met_unc+uncReplaced)==branch->GetName())    met_unc=branch->GetName();
    else if ((metphi_unc+uncReplaced)==branch->GetName()) metphi_unc=branch->GetName();
  }
  // print used tau and jet kinematics
  std::cout << std::endl << " >> Searching " << uncReplaced << std::endl << std::endl;
  std::cout << pt_1_unc << std::endl;
  std::cout << pt_2_unc << std::endl;
  std::cout << met_unc << std::endl;
  std::cout << metphi_unc << std::endl;
  std::cout << "-------------- used branches --------------" << std::endl << std::endl;
  // pick up branche
  namu -> SetBranchAddress( pt_1_unc.c_str()         , &pt_1       );
  namu -> SetBranchAddress( pt_2_unc.c_str()         , &pt_2       );
  namu -> SetBranchAddress( met_unc.c_str()          , &met         );
  namu -> SetBranchAddress( metphi_unc.c_str()       , &metphi      );

  // initializing outputs
  svFitMass = -999.0;
  svFitPt = -999.0;
  svFitEta = -999.0;
  svFitPhi = -999.0;
  svFitTransverseMass = -999.0;
  svFitMET = -999.0;
  metcorr_ex = -999.0;
  metcorr_ey = -999.0;
  metcor = -999.0;
  metcorphi = -999.0;
  tau1_pt = -999.0;
  tau1_eta = -999.0;
  tau1_phi = -999.0;
  tau1_m = -999.0;
  tau2_pt = -999.0;
  tau2_eta = -999.0;
  tau2_phi = -999.0;
  tau2_m = -999.0;
}



