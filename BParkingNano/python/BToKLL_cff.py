import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

electronPairsForKee = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.5 && userFloat("unBiased") >= 3'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.00'
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-4'),
)

BToKee = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('electronPairsForKee'),
    leptonTransientTracks = electronPairsForKee.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    kaonSelection = cms.string(''),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.00 '
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.005 '
        '&& userFloat("fitted_cos_theta_2D") >= 0 '
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    ),
    doJpsiConstr = cms.bool(True),    
)

muonPairsForKmumu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.0'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 4.5 '
                                 '&& mass() > 1.5 && charge() == 0 && userFloat("lep_deltaR") > 0.0'),
    postVtxSelection = electronPairsForKee.postVtxSelection,
)

BToKmumu = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('muonPairsForKmumu'),
    leptonTransientTracks = muonPairsForKmumu.transientTracksSrc,
    kaons = BToKee.kaons,
    kaonsTransientTracks = BToKee.kaonsTransientTracks,
    beamSpot = cms.InputTag("offlineBeamSpot"),
    kaonSelection = cms.string(''),
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.0'
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.005 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    ),
    doJpsiConstr = cms.bool(True),    
)

BToKeeTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKee"),
    cut = cms.string(""),
    name = cms.string("BToKEE"),
    doc = cms.string("BToKEE Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
        kIdx = uint('k_idx'),
        minDR = ufloat('min_dr'),
        maxDR = ufloat('max_dr'),
        # fit and vtx info
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work
        mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        fit_l1_pt = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        fit_l2_pt = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        fit_k_pt = ufloat('fitted_k_pt'),
        fit_k_eta = ufloat('fitted_k_eta'),
        fit_k_phi = ufloat('fitted_k_phi'),
    )
)

BToKmumuTable = BToKeeTable.clone(
    src = cms.InputTag("BToKmumu"),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable")
)


CountBToKee = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToKee")
)    
CountBToKmumu = CountBToKee.clone(
    minNumber = cms.uint32(1),
    src = cms.InputTag("BToKmumu")
)


BToKMuMuSequence = cms.Sequence(
    (muonPairsForKmumu * BToKmumu)
)
BToKEESequence = cms.Sequence(
    (electronPairsForKee * BToKee)
)

BToKLLSequence = cms.Sequence(
    (electronPairsForKee * BToKee) +
    (muonPairsForKmumu * BToKmumu)
)
BToKLLTables = cms.Sequence(BToKeeTable + BToKmumuTable)

