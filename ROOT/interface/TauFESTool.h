// Copyright Tyler Mitchell [2020]

#ifndef ROOT_INTERFACE_TAUFESTOOL_H_
#define ROOT_INTERFACE_TAUFESTOOL_H_

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TGraph.h"

typedef std::vector<Float_t> fvector;
std::vector<Float_t> percentage_to_weight(std::vector<Float_t>);
std::vector<Float_t> percentage_to_decimal(std::vector<Float_t>);

class TauFESTool {
   public:
    TauFESTool(std::string, std::string, std::string, bool);
    Float_t getFES(Float_t, Float_t, Float_t, std::string);
    Float_t getTES(Float_t, Float_t, Bool_t);

   private:
    std::unordered_map<std::string, fvector> fes_sfs;
    std::unordered_map<std::string, fvector> tes_sfs;
    std::unordered_map<int, int> dm_map;
};

TauFESTool::TauFESTool(std::string year, std::string id, std::string path, bool isEmbed)
    : fes_sfs{{"nominal", fvector{}}, {"up", fvector{}}, {"down", fvector{}}},
      tes_sfs{{"nominal", fvector{}}, {"syst", fvector{}}},
      dm_map{{0, 0}, {1, 1}, {10, 2}, {11, 3}} {
    std::string filepath = "$CMSSW_BASE/src/" + path + "TauFES_eta-dm_" + id + "_" + year + ".root";
    auto tes_file = new TFile(filepath.c_str(), "READ");
    TGraph* tes_graph = reinterpret_cast<TGraph*>(tes_file->Get("fes"));

    // i = 0 : barrel dm 0
    // i = 1 : barrel dm 1
    // i = 2 : endcap dm 0
    // i = 3 : endcap dm 1
    for (auto i = 0; i < 4; i++) {
        fes_sfs["nominal"].push_back(tes_graph->GetY()[i]);
        fes_sfs["up"].push_back(tes_graph->GetErrorYhigh(i) - tes_graph->GetY()[i]);
        fes_sfs["down"].push_back(tes_graph->GetErrorYlow(i) - tes_graph->GetY()[i]);
    }

    // make_pair(energy scale, systematic shift)
    // order: dm0, dm1, dm10, dm11
    std::vector<Float_t> shifts;
    if (isEmbed) {
        // embedded only
        Float_t dm0_embed(0.975), dm0_embed_syst(0.008);
        Float_t dm1_embed(0.975 * 1.051), dm1_embed_syst(.016);
        Float_t dm10_embed(0.975 * 0.975 * 0.975), dm10_embed_syst(.024);
        tes_sfs["nominal"] = {dm0_embed, dm1_embed, dm10_embed, 1.};
        tes_sfs["syst"] = {dm0_embed_syst, dm1_embed_syst, dm10_embed_syst, 0.};
    } else {
        std::vector<Float_t> nominal;
        if (year == "2016Legacy") {
            tes_sfs["nominal"] = percentage_to_weight({-0.6, -0.5, 0., -0.1});
            tes_sfs["syst"] = percentage_to_decimal({1., 0.9, 1.1, 1.});
        } else if (year == "2017ReReco") {
            tes_sfs["nominal"] = percentage_to_weight({0.7, -0.2, 0.1, -0.1});
            tes_sfs["syst"] = percentage_to_decimal({0.8, 0.8, 0.9, 1.});
        } else if (year == "2018ReReco") {
            tes_sfs["nominal"] = percentage_to_weight({-1.3, -0.5, -1.2, -0.1});
            tes_sfs["syst"] = percentage_to_decimal({1.1, 0.9, 0.8, 1.});
        }
    }
}

Float_t TauFESTool::getTES(Float_t dm, Float_t gen_match, Bool_t syst = false) {
    if (gen_match != 5) {
        return syst ? 0. : 1.;  // 0 for syst, 1 for nominal
    } else if (dm == 5 || dm == 6) {
      return syst ? 0. : 1.;  // these will be removed anyways
    }

    auto index = dm_map.at(dm);
    std::string key = "nominal";
    if (syst) {
        key = "syst";
    }

    return tes_sfs.at(key).at(index);
}

Float_t TauFESTool::getFES(Float_t dm, Float_t eta, Float_t gen_match, std::string syst = "") {
    // SF is 1 if not an electron
    if (gen_match != 1 && gen_match != 3) {
        return syst == "" || syst == "nomnial" ? 1. : 0.;
    } else if (dm != 0 && dm != 1) {
      return syst == "" || syst == "nomnial" ? 1. : 0.;  // these will be removed anyways
    }

    // handle default case
    if (syst == "") {
        syst = "nominal";
    }

    // handle barrel vs endcap
    int index = 0;
    if (eta > 1.5) {
        index = 2;
    }

    // decay mode 0 or 1
    index += dm;
    return fes_sfs[syst].at(index);
}

std::vector<Float_t> percentage_to_weight(std::vector<Float_t> input) {
    std::vector<Float_t> output;
    for (auto in : input) {
        output.push_back(1 + (0.01 * in));
    }
    return output;
}

std::vector<Float_t> percentage_to_decimal(std::vector<Float_t> input) {
    std::vector<Float_t> output;
    for (auto in : input) {
        output.push_back(0.01 * in);
    }
    return output;
}

#endif  // ROOT_INTERFACE_TAUFESTOOL_H_
