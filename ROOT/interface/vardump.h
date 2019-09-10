float svFitMass = -10;
float svFitPt = -10;
float svFitEta = -10;
float svFitPhi = -10;
float svFitMET = -10;
float svFitTransverseMass = -10;

float metcorr_ex = -10;               // corrected met px (float)
float metcorr_ey = -10;               // corrected met py (float)
float metcorrUncUp_ex = -10;          // corrUncUpected met px (float)
float metcorrUncUp_ey = -10;          // corrUncUpected met py (float)
float metcorrUncDown_ex = -10;        // corrUncDownected met px (float)
float metcorrUncDown_ey = -10;        // corrUncDownected met py (float)
float metcorrClusteredUp_ex = -10;    // corrClusteredUpected met px (float)
float metcorrClusteredUp_ey = -10;    // corrClusteredUpected met py (float)
float metcorrClusteredDown_ex = -10;  // corrClusteredDownected met px (float)
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
float metcor = -10;     // corrected metcor
float metcorphi = -10;  // corrected metcorphi
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
float tau1_pt = -10;
float tau1_eta = -10;
float tau1_phi = -10;
float tau1_m = -10;
float tau2_pt = -10;
float tau2_eta = -10;
float tau2_phi = -10;
float tau2_m = -10;
// up (whatever it is) ============
float tau1_pt_Up = -10;
float tau1_eta_Up = -10;
float tau1_phi_Up = -10;
float tau1_m_Up = -10;
float tau2_pt_Up = -10;
float tau2_eta_Up = -10;
float tau2_phi_Up = -10;
float tau2_m_Up = -10;
// down
float tau1_pt_Down = -10;
float tau1_eta_Down = -10;
float tau1_phi_Down = -10;
float tau1_m_Down = -10;
float tau2_pt_Down = -10;
float tau2_eta_Down = -10;
float tau2_phi_Down = -10;
float tau2_m_Down = -10;
// up DM0 =========================
float tau1_pt_DM0_Up = -10;
float tau1_eta_DM0_Up = -10;
float tau1_phi_DM0_Up = -10;
float tau1_m_DM0_Up = -10;
float tau2_pt_DM0_Up = -10;
float tau2_eta_DM0_Up = -10;
float tau2_phi_DM0_Up = -10;
float tau2_m_DM0_Up = -10;
// down
float tau1_pt_DM0_Down = -10;
float tau1_eta_DM0_Down = -10;
float tau1_phi_DM0_Down = -10;
float tau1_m_DM0_Down = -10;
float tau2_pt_DM0_Down = -10;
float tau2_eta_DM0_Down = -10;
float tau2_phi_DM0_Down = -10;
float tau2_m_DM0_Down = -10;
// up DM1 =========================
float tau1_pt_DM1_Up = -10;
float tau1_eta_DM1_Up = -10;
float tau1_phi_DM1_Up = -10;
float tau1_m_DM1_Up = -10;
float tau2_pt_DM1_Up = -10;
float tau2_eta_DM1_Up = -10;
float tau2_phi_DM1_Up = -10;
float tau2_m_DM1_Up = -10;
// down
float tau1_pt_DM1_Down = -10;
float tau1_eta_DM1_Down = -10;
float tau1_phi_DM1_Down = -10;
float tau1_m_DM1_Down = -10;
float tau2_pt_DM1_Down = -10;
float tau2_eta_DM1_Down = -10;
float tau2_phi_DM1_Down = -10;
float tau2_m_DM1_Down = -10;
// up DM10 =========================
float tau1_pt_DM10_Up = -10;
float tau1_eta_DM10_Up = -10;
float tau1_phi_DM10_Up = -10;
float tau1_m_DM10_Up = -10;
float tau2_pt_DM10_Up = -10;
float tau2_eta_DM10_Up = -10;
float tau2_phi_DM10_Up = -10;
float tau2_m_DM10_Up = -10;
// down
float tau1_pt_DM10_Down = -10;
float tau1_eta_DM10_Down = -10;
float tau1_phi_DM10_Down = -10;
float tau1_m_DM10_Down = -10;
float tau2_pt_DM10_Down = -10;
float tau2_eta_DM10_Down = -10;
float tau2_phi_DM10_Down = -10;
float tau2_m_DM10_Down = -10;
// up UncMet =========================
float tau1_pt_UncMet_Up = -10;
float tau1_eta_UncMet_Up = -10;
float tau1_phi_UncMet_Up = -10;
float tau1_m_UncMet_Up = -10;
float tau2_pt_UncMet_Up = -10;
float tau2_eta_UncMet_Up = -10;
float tau2_phi_UncMet_Up = -10;
float tau2_m_UncMet_Up = -10;
// down
float tau1_pt_UncMet_Down = -10;
float tau1_eta_UncMet_Down = -10;
float tau1_phi_UncMet_Down = -10;
float tau1_m_UncMet_Down = -10;
float tau2_pt_UncMet_Down = -10;
float tau2_eta_UncMet_Down = -10;
float tau2_phi_UncMet_Down = -10;
float tau2_m_UncMet_Down = -10;
// up ClusteredMet =========================
float tau1_pt_ClusteredMet_Up = -10;
float tau1_eta_ClusteredMet_Up = -10;
float tau1_phi_ClusteredMet_Up = -10;
float tau1_m_ClusteredMet_Up = -10;
float tau2_pt_ClusteredMet_Up = -10;
float tau2_eta_ClusteredMet_Up = -10;
float tau2_phi_ClusteredMet_Up = -10;
float tau2_m_ClusteredMet_Up = -10;
// down
float tau1_pt_ClusteredMet_Down = -10;
float tau1_eta_ClusteredMet_Down = -10;
float tau1_phi_ClusteredMet_Down = -10;
float tau1_m_ClusteredMet_Down = -10;
float tau2_pt_ClusteredMet_Down = -10;
float tau2_eta_ClusteredMet_Down = -10;
float tau2_phi_ClusteredMet_Down = -10;
float tau2_m_ClusteredMet_Down = -10;
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