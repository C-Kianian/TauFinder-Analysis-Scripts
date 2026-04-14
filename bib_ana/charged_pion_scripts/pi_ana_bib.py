# Comes from ethanmar/MuColl-TauStudy/analysis/pfo_matching.py
# Run this code after the reco step on a pion gun
from pyLCIO import IOIMPL
import ROOT
from ROOT import TH1F, TFile, TCanvas, gPad
import math
from argparse import ArgumentParser
import os
from array import array

ROOT.gStyle.SetOptFit(111)
ROOT.gStyle.SetOptStat("nemruo") #for uf/of info

##################
# CHANGE BINS HERE
rebins = list(range(0, 1001, 50))
n_rebins = len(rebins) - 1
##################

# Args
parser = ArgumentParser()
parser.add_argument('-c', '--charge', type=str, required=True,
                    help='Charge state of the pions ("plus", "minus", or "both")') # Input plus, minus, or both for this arg
parser.add_argument('-i', '--inputFile', type=str, default='output_reco.slcio')
parser.add_argument('-o', '--outputFile', type=str, default='pi_bib_ana.root')
args = parser.parse_args()

# get charge pion we are interested in
charge = str(args.charge).lower()
allowed_pdgs = {
    'plus': {211},
    'minus': {-211},
    'both': {211, -211}
}
latex = {
    'plus': "Simulated #pi^{+} Gun",
    'minus': "Simulated #pi^{-} Gun",
    'both': "Simulated #pi^{+} and #pi^{-} Guns"
}
if charge not in allowed_pdgs: raise ValueError("Invalid charge state. Must be 'plus', 'minus', or 'both'.")
selected_pdgs = allowed_pdgs[charge]

# Hist setup
PT_MIN, PT_MAX, PT_BINS = 0.0, 1000.0, 160

ETA_MIN, ETA_MAX, ETA_BINS = 0.0, 2.5, 20

THETA_BINS = 20

PHI_BINS = 24

EFF_MIN, EFF_MAX = 0.0, 1.2

# Hists
hists = []

def book(h):
    h.SetDirectory(0)
    hists.append(h)
    return h

# All reco pions
fAllPt = book(TH1F('reco_pfo_pt', 'All Reco Charged Pion p_{T};p_{T}_reco [GeV/c];Entries', PT_BINS, PT_MIN, PT_MAX))
fAllEta = book(TH1F('reco_pfo_eta', 'All Reco Charged Pion |#eta|;|#eta_reco|;Entries', ETA_BINS, ETA_MIN, ETA_MAX))
fAllTheta = book(TH1F('reco_pfo_theta', 'All Reco Charged Pion #theta;#theta_reco [rad];Entries', THETA_BINS, 0, math.pi))
fAllPhi = book(TH1F('reco_pfo_phi', 'All Reco Charged Pion #phi;#phi_reco;Entries', PHI_BINS, -math.pi, math.pi))
fAllE = book(TH1F('reco_pfo_e', 'All Reco Charged Pion E;E_reco [GeV];Entries', PT_BINS, PT_MIN, PT_MAX))

# MC pion per event
fMCPt = book(TH1F('mc_pion_pt', 'MC Charged Pion p_{T};p_{T}_reco;Entries', PT_BINS, PT_MIN, PT_MAX))
fMCEta = book(TH1F('mc_pion_eta', 'MC Charged Pion |#eta|;|#eta_reco|;Entries', ETA_BINS, ETA_MIN, ETA_MAX))
fMCTheta = book(TH1F('mc_pion_theta', 'MC Charged Pion #theta;#theta_reco [rad];Entries', THETA_BINS, 0, math.pi))
fMCPhi = book(TH1F('mc_pion_phi', 'MC Charged Pion #phi;#phi_reco;Entries', PHI_BINS, -math.pi, math.pi))
fMCE = book(TH1F('mc_pion_e', 'MC Charged Pion E;E_reco [GeV];Entries', PT_BINS, PT_MIN, PT_MAX))

# Best reco pion matched to MC pion
fMatchedPt = book(TH1F('mc_matched_pt', 'Matched Best Reco Charged Pion MC p_{T};p_{T}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
fMatchedEta = book(TH1F('mc_matched_eta', 'Matched Best Reco Charged Pion MC |#eta|;|#eta_true|;Entries', ETA_BINS, ETA_MIN, ETA_MAX))
fMatchedTheta = book(TH1F('mc_matched_theta', 'Matched Best Charged Reco Pion MC #theta;#theta_true [rad];Entries', THETA_BINS, 0, math.pi))
fMatchedPhi = book(TH1F('mc_matched_phi', 'Matched Best Reco Charged Pion MC #phi;#phi_true;Entries', PHI_BINS, -math.pi, math.pi))
fMatchedE = book(TH1F('mc_matched_e', 'Matched Best Reco Charged Pion MC E;E_true [GeV];Entries', PT_BINS, PT_MIN, PT_MAX))

# Best reco pion matched to MC pion stricter DR
fMatchedPtStrict = book(TH1F('mc_matched_pt_strict', '(Strict) Matched Best Reco Charged Pion MC p_{T};p_{T}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
fMatchedEtaStrict = book(TH1F('mc_matched_eta_strict', '(Strict) Matched Best Reco Charged Pion MC |#eta|;|#eta_true|;Entries', ETA_BINS, ETA_MIN, ETA_MAX))
fMatchedThetaStrict = book(TH1F('mc_matched_theta_strict', '(Strict) Matched Best Charged Reco Pion MC #theta;#theta_true [rad];Entries', THETA_BINS, 0, math.pi))
fMatchedPhiStrict = book(TH1F('mc_matched_phi_strict', '(Strict) Matched Best Reco Charged Pion MC #phi;#phi_true;Entries', PHI_BINS, -math.pi, math.pi))
fMatchedEStrict = book(TH1F('mc_matched_e_strict', '(Strict) Matched Best Reco Charged Pion MC E;E_true [GeV];Entries', PT_BINS, PT_MIN, PT_MAX))

# Residuals
fResidualPt = book(TH1F('residual_pt', 'Residual Transverse Momentum;True - Reco p_{T} [GeV/c];Entries', 400, -100, 100))
fResidualPt.SetMaximum(1000)

fResidualE = book(TH1F('residual_e', 'Residual Energy;True - Reco E [GeV];Entries', 400, -100, 100))
fResidualE.SetMaximum(1000)

# Strict residuals
fResidualEStrict = book(TH1F('residual_e_strict', '(Strict) Residual Energy;True - Reco E [GeV];Entries', 400, -100, 100))
fResidualEStrict.SetMaximum(1000)

fResidualPtStrict = book(TH1F('residual_pt_strict', '(Strict) Residual Transverse Momentum;True - Reco p_{T} [GeV/c];Entries', 400, -100, 100))
fResidualPtStrict.SetMaximum(1000)

# Resolution
fResPt = TH1F('resolution_pt', 'p_{T} Resolution;(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
fResPt.SetMaximum(1000)

fResE = TH1F('resolution_e', 'E Resolution;(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
fResE.SetMaximum(1000)

# Strict resolution
fResEStrict = TH1F('resolution_e_strict', '(Strict) E Resolution;(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
fResEStrict.SetMaximum(1000)

fResPtStrict = TH1F('resolution_pt_strict', '(Strict) p_{T} Resolution;(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
fResPtStrict.SetMaximum(1000)

# Regional pT histograms
regions = ['barrel', 'centbarrel', 'transition', 'endcap']
regional_max = [800, 600, 400, 200]
fAllPtReg = {}
fMatchedPtReg = {}
fMatchedPtRegStrict = {}
fResPtReg = {}
fResEReg = {}
fResPtRegStrict = {}
fResERegStrict = {}

for r in regions:
    # Matched info
    fAllPtReg[r] = book(TH1F(f'mc_pion_pt_{r}', f'MC Charged Pion p_{{T}} ({r})', PT_BINS, PT_MIN, PT_MAX))
    fMatchedPtReg[r] = book(TH1F(f'mc_matched_pt_{r}',f'Matched MC Charged Pion p_{{T}} ({r})',PT_BINS, PT_MIN, PT_MAX))
    fMatchedPtRegStrict[r] = book(TH1F(f'mc_matched_pt_{r}_strict',f'(Strict) Matched MC Charged Pion p_{{T}} ({r})',PT_BINS, PT_MIN, PT_MAX))

    # Regional resolutions
    fResPtReg[r] = TH1F(f'resolution_pt_{r}', 'p_{T} Resolution ('+r+');(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
    fResPtReg[r].SetMaximum(regional_max[regions.index(r)])
    fResEReg[r] = TH1F(f'resolution_e_{r}', 'Energy Resolution (' + r + ');(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
    fResEReg[r].SetMaximum(regional_max[regions.index(r)])
    fResPtRegStrict[r] = TH1F(f'resolution_pt_{r}_strict', '(Strict) p_{T} Resolution ('+r+');(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
    fResPtRegStrict[r].SetMaximum(regional_max[regions.index(r)])
    fResERegStrict[r] = TH1F(f'resolution_e_{r}_strict', '(Strict) Energy Resolution ('+r+');(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
    fResERegStrict[r].SetMaximum(regional_max[regions.index(r)])


# Get theta regions and eta
def eta(theta):
    # Added to prevent 0 crashes
    if theta < 1e-10:
        theta = 1e-10
    elif (math.pi - theta) < 1e-10:
        theta = math.pi - 1e-10

    return -math.log(math.tan(theta / 2.0))

def theta_region(theta):
    regs = []
    if 0.70 < theta < 2.45:
        regs.append('barrel')
    if 1.0 < theta < 2.0:
        regs.append('centbarrel')
    if (0.577 < theta < 1.0) or (2.0 < theta < 2.56):
        regs.append('transition')
    if theta < 0.577 or theta > 2.56:
        regs.append('endcap')
    if len(regs) > 0: return regs
    return None

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2*math.pi
    while dphi < -math.pi:
        dphi += 2*math.pi
    return dphi


# Get input files
to_process = []
if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)

# Event loop
for file in to_process:

    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:

        pfos = event.getCollection('PandoraPFOs')
        mcs  = event.getCollection('MCParticle')

        # Best MC charged pion
        best_mc = None
        best_mc_pt = -1.0

        for mc in mcs:
            if mc.getPDG() not in selected_pdgs: continue

            mcMom = mc.getMomentum()
            mcPt = math.hypot(mcMom[0], mcMom[1])

            if mcPt > best_mc_pt:
                best_mc_pt = mcPt
                best_mc = mc

        if best_mc is None:
            continue

        # MC kinematics
        mcPDG = best_mc.getPDG()
        mcMom = best_mc.getMomentum()
        mcPt = best_mc_pt
        mcTheta = math.acos(mcMom[2] / math.sqrt(mcPt**2 + mcMom[2]**2))
        mcEta = eta(mcTheta)
        mcPhi = math.atan2(mcMom[1], mcMom[0])
        mcE = best_mc.getEnergy()

        regs = theta_region(mcTheta)

        # MC info, denominators for eff plots
        fMCE.Fill(mcE)
        fMCPt.Fill(mcPt)
        fMCEta.Fill(abs(mcEta))
        fMCTheta.Fill(mcTheta)
        fMCPhi.Fill(mcPhi)

        if regs:
            for reg in regs:
                fAllPtReg[reg].Fill(mcPt)

        # Best reconstructed charged pion and PFOs
        best_reco_pi = None
        best_reco_pi_pt = -1.0

        for pfo in pfos:
            if pfo.getType() != mcPDG: continue

            # Pion kinematics
            recoPiMom = pfo.getMomentum()
            recoPiPt = math.hypot(recoPiMom[0], recoPiMom[1])
            pfoE = pfo.getEnergy()
            recoPiTheta = math.acos(recoPiMom[2] / math.sqrt(recoPiPt ** 2 + recoPiMom[2] ** 2))
            recoPiEta = eta(recoPiTheta)
            recoPiPhi = math.atan2(recoPiMom[1], recoPiMom[0])

            # Fill all pion histograms
            fAllPt.Fill(recoPiPt)
            fAllEta.Fill(abs(recoPiEta))
            fAllTheta.Fill(recoPiTheta)
            fAllPhi.Fill(recoPiPhi)
            fAllE.Fill(pfoE)

            if recoPiPt > best_reco_pi_pt:
                best_reco_pi_pt = recoPiPt
                best_reco_pi = pfo

        if best_reco_pi is None:
            continue

        # Charged Pi kinematics
        recoPiMom = best_reco_pi.getMomentum()
        recoPiPt = best_reco_pi_pt
        recoPiTheta = math.acos(recoPiMom[2] / math.sqrt(recoPiPt ** 2 + recoPiMom[2] ** 2))
        recoPiEta = eta(recoPiTheta)
        recoPiPhi = math.atan2(recoPiMom[1], recoPiMom[0])
        recoPiE = best_reco_pi.getEnergy()

        # dR matching
        dphi = delta_phi(mcPhi, recoPiPhi)
        dR = math.sqrt(dphi*dphi + (mcEta - recoPiEta)**2)

        if dR < 0.25:
            fMatchedPt.Fill(mcPt)
            fMatchedEta.Fill(abs(mcEta))
            fMatchedTheta.Fill(mcTheta)
            fMatchedPhi.Fill(mcPhi)
            fMatchedE.Fill(mcE)
            fResidualPt.Fill(mcPt - recoPiPt)
            fResidualE.Fill(mcE - recoPiE)
            fResPt.Fill((mcPt - recoPiPt) / mcPt)
            fResE.Fill((mcE - recoPiE) / mcE)

            if regs:
                for reg in regs:
                    fMatchedPtReg[reg].Fill(mcPt)
                    fResPtReg[reg].Fill((mcPt - recoPiPt) / mcPt)
                    fResEReg[reg].Fill((mcE - recoPiE) / mcE)

        if dR < 0.1:
            fMatchedPtStrict.Fill(mcPt)
            fMatchedEtaStrict.Fill(abs(mcEta))
            fMatchedThetaStrict.Fill(mcTheta)
            fMatchedPhiStrict.Fill(mcPhi)
            fMatchedEStrict.Fill(mcE)
            fResidualPtStrict.Fill(mcPt - recoPiPt)
            fResidualEStrict.Fill(mcE - recoPiE)
            fResPtStrict.Fill((mcPt - recoPiPt) / mcPt)
            fResEStrict.Fill((mcE - recoPiE) / mcE)

            if regs:
                for reg in regs:
                    fMatchedPtRegStrict[reg].Fill(mcPt)
                    fResPtRegStrict[reg].Fill((mcPt - recoPiPt) / mcPt)
                    fResERegStrict[reg].Fill((mcE - recoPiE) / mcE)

    reader.close()

# Eff plots
def make_eff(num, den, name, title, xtitle, rebin=False):
    if rebin:
        rebin_arr = array('d', rebins) # make rebin array into cpp array

        num_rebin = num.Rebin(n_rebins, num.GetName()+"_rebin", rebin_arr)
        den_rebin = den.Rebin(n_rebins, den.GetName()+"_rebin", rebin_arr)

        eff = num_rebin.Clone(name)
        eff.Divide(num_rebin, den_rebin, 1, 1, 'B')
    else:
        eff = num.Clone(name)
        eff.Divide(num, den, 1, 1, 'B')


    eff.SetLineWidth(2)
    eff.SetTitle(title)
    eff.GetXaxis().SetTitle(xtitle)
    eff.GetYaxis().SetTitle('#epsilon')
    eff.GetYaxis().SetRangeUser(EFF_MIN, EFF_MAX)
    eff.SetStats(0)
    book(eff)

# Function to fit gaussians to resolutions
def fit_gauss(hist):
    hist.Fit("gaus", "Q", "", -0.1, 0.1)
    book(hist)

# Eff plots
make_eff(fMatchedPt, fMCPt, 'pt_eff', 'Charged Pion Efficiency vs p_{T}', 'True p_{T} [GeV/c]', rebin=True)
make_eff(fMatchedEta, fMCEta, 'eta_eff', 'Charged Pion Efficiency vs |#eta|', 'True |#eta|')
make_eff(fMatchedTheta, fMCTheta, 'theta_eff', 'Charged Pion Efficiency vs #theta', 'True #theta')
make_eff(fMatchedPhi, fMCPhi, 'phi_eff', 'Charged Pion Efficiency vs #phi', 'True #phi')
make_eff(fMatchedE, fMCE, 'e_eff', 'Charged Pion Efficiency vs E', 'True E [GeV]')

# Strict eff plots
make_eff(fMatchedPtStrict, fMCPt, 'pt_eff_strict', '(Strict) Charged Pion Efficiency vs p_{T}', 'True p_{T} [GeV/c]', rebin=True)
make_eff(fMatchedEtaStrict, fMCEta, 'eta_eff_strict', '(Strict) Charged Pion Efficiency vs |#eta|', 'True |#eta|')
make_eff(fMatchedThetaStrict, fMCTheta, 'theta_eff_strict', '(Strict) Charged Pion Efficiency vs #theta', 'True #theta')
make_eff(fMatchedPhiStrict, fMCPhi, 'phi_eff_strict', '(Strict) Charged Pion Efficiency vs #phi', 'True #phi')
make_eff(fMatchedEStrict, fMCE, 'e_eff_strict', '(Strict) Charged Pion Efficiency vs E', 'True E [GeV]')

# Fit guassians
fit_gauss(fResPt)
fit_gauss(fResE)
fit_gauss(fResPtStrict)
fit_gauss(fResEStrict)

# Regional hists
for r in regions:
    # Regional eff
    make_eff(fMatchedPtReg[r], fAllPtReg[r],f'pt_{r}_eff',f'Charged Pion Efficiency vs p_{{T}} ({r})','True p_{T} [GeV/c]', rebin=True)
    make_eff(fMatchedPtRegStrict[r],fAllPtReg[r],f'pt_{r}_eff_strict',f'(Strict) Charged Pion Efficiency vs p_{{T}} ({r})','True p_{T} [GeV/c]', rebin=True)

    # Regional resoltuions
    fit_gauss(fResEReg[r])
    fit_gauss(fResPtReg[r])
    fit_gauss(fResPtRegStrict[r])
    fit_gauss(fResERegStrict[r])

# Root output
output = TFile(args.outputFile, 'RECREATE')
for h in hists:
    h.Write()
output.Close()

# png plots
c = TCanvas() # Canvas

text = ROOT.TLatex() # Text settings
text.SetNDC()
text.SetTextFont(42)
text.SetTextSize(0.04)
text.SetTextAlign(13)

for h in hists:
    c.Clear() # Clear canvas
    h.Draw()
    gPad.Update()

    text.DrawLatex(0.12, 0.93, "#bf{#it{MAIA Detector Concept}}")
    text.DrawLatex(0.12, 0.88, latex[charge])
    c.Update()

    c.SaveAs(h.GetName() + '.png')

