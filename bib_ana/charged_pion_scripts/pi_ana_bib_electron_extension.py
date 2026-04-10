# Extended from ethanmar/MuColl-TauStudy/analysis/pfo_matching.py
# MC pion matched to best reconstructed charged pion OR electron
# Run this code after the reco step on a pion gun
from pyLCIO import IOIMPL
import ROOT
from ROOT import TH1F, TFile, TCanvas, gPad
import math
from argparse import ArgumentParser
import os

ROOT.gStyle.SetOptFit(111)
ROOT.gStyle.SetOptStat("nemruo") #for uf/of info

# Args
parser = ArgumentParser()
parser.add_argument('-i', '--inputFile', type=str, default='output_reco.slcio')
parser.add_argument('-o', '--outputFile', type=str, default='piEle_bib_ana.root')
args = parser.parse_args()

# Hist setup
PT_MIN, PT_MAX, PT_BINS = 0.0, 1000.0, 20

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
fAllPt = book(TH1F('reco_piEle_pt', 'All Reco Charged Pion or Electron p_{T};p_{T}^{reco} [GeV];Entries', PT_BINS, PT_MIN, PT_MAX))
fAllEta = book(TH1F('reco_piEle_eta', 'All Reco Charged Pion or Electron |#eta|;|#eta^{reco}|;Entries', ETA_BINS, ETA_MIN, ETA_MAX))
fAllTheta = book(TH1F('reco_piEle_theta', 'All Reco Charged Pion or Electron #theta;#theta^{reco};Entries', THETA_BINS, 0, math.pi))
fAllPhi = book(TH1F('reco_piEle_phi', 'All Reco Charged Pion or Electron #phi;#phi^{reco};Entries', PHI_BINS, -math.pi, math.pi))
fAllE = book(TH1F('reco_piEle_e', 'All Reco Charged Pion or Electron E;E^{reco} [GeV];Entries',PT_BINS, PT_MIN, PT_MAX))

# MC pion per event
fMCPt = book(TH1F('mc_pion_pt', 'MC Charged Pion p_{T};p_{T}^{true};Entries',PT_BINS, PT_MIN, PT_MAX))
fMCEta = book(TH1F('mc_pion_eta', 'MC Charged Pion |#eta|;|#eta^{true}|;Entries', ETA_BINS, ETA_MIN, ETA_MAX))
fMCTheta = book(TH1F('mc_pion_theta', 'MC Charged Pion #theta;#theta^{true};Entries', THETA_BINS, 0, math.pi))
fMCPhi = book(TH1F('mc_pion_phi', 'MC Charged Pion #phi;#phi^{true};Entries',PHI_BINS, -math.pi, math.pi))
fMCE = book(TH1F('mc_pion_e', 'MC Charged Pion E;E^{true};Entries',PT_BINS, PT_MIN, PT_MAX))

# Best reco pion or e matched to MC pion
fMatchedPt = book(TH1F('mc_matched_pt_piEle','Matched MC Pion to Reco #pi or e p_{T};p_{T}^{true};Entries', PT_BINS, PT_MIN, PT_MAX))
fMatchedEta = book(TH1F('mc_matched_eta_piEle', 'Matched MC Pion to Reco #pi or e |#eta|;|#eta_true|;Entries',ETA_BINS, ETA_MIN, ETA_MAX))
fMatchedTheta = book(TH1F('mc_matched_theta_piEle', 'Matched MC Pion to Reco #pi or e #theta;#theta_true [rad];Entries', THETA_BINS, 0, math.pi))
fMatchedPhi = book(TH1F('mc_matched_phi_piEle', 'Matched MC Pion to Reco #pi or e #phi;#phi_true;Entries', PHI_BINS, -math.pi, math.pi))
fMatchedE = book(TH1F('mc_matched_e_piEle', 'Matched MC Pion to Reco #pi or e E;E_true [GeV];Entries', PT_BINS, PT_MIN, PT_MAX))

# Best reco pion or e matched to MC pion stricter DR
fMatchedPtStrict = book(TH1F('mc_matched_pt_piEle_strict', '(Strict) Matched MC Pion to Reco #pi or e p_{T};p_{T}^{true};Entries',PT_BINS, PT_MIN, PT_MAX))
fMatchedEtaStrict = book(TH1F('mc_matched_eta_piEle_strict', '(Strict) Matched MC Pion to Reco #pi or e |#eta|;|#eta_true|;Entries',ETA_BINS, ETA_MIN, ETA_MAX))
fMatchedThetaStrict = book(TH1F('mc_matched_theta_piEle_strict', '(Strict) Matched MC Pion to Reco #pi or e #theta;#theta_true [rad];Entries', THETA_BINS, 0, math.pi))
fMatchedPhiStrict = book(TH1F('mc_matched_phi_piEle_strict', '(Strict) Matched MC Pion to Reco #pi or e #phi;#phi_true;Entries', PHI_BINS, -math.pi, math.pi))
fMatchedEStrict = book(TH1F('mc_matched_e_piEle_strict', '(Strict) Matched MC Pion to Reco #pi or e E;E_true [GeV];Entries', PT_BINS, PT_MIN, PT_MAX))

# Residual
fResidualPt = book(TH1F('residual_pt_piEle', 'Residual p_{T};True - Reco p_{T} [GeV/c];Entries', 400, -100, 100))
fResidualPt.SetMaximum(1000)

fResidualE = book(TH1F('residual_e_piEle', 'Residual E;True - Reco E [GeV];Entries', 400, -100, 100))
fResidualE.SetMaximum(1000)

# Strict residual
fResidualPtStrict = book(TH1F('residual_pt_piEle_strict', 'Residual p_{T};True - Reco p_{T} [GeV/c];Entries', 400, -100, 100))
fResidualPtStrict.SetMaximum(1000)

fResidualEStrict = book(TH1F('residual_e_piEle_strict', 'Residual E;True - Reco E [GeV];Entries', 400, -100, 100))
fResidualEStrict.SetMaximum(1000)

# Resolution
fResPt = TH1F('resolution_pt_piEle', 'p_{T} Resolution;(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
fResPt.SetMaximum(1000)

fResE = TH1F('resolution_e_piEle', 'Energy Resolution;(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
fResE.SetMaximum(1000)

# Strict resolution
fResPtStrict = TH1F('resolution_pt_piEle_strict', 'p_{T} Resolution;(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
fResPtStrict.SetMaximum(1000)

fResEStrict = TH1F('resolution_e_piEle_strict', 'Energy Resolution;(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
fResEStrict.SetMaximum(1000)

# Truth E vs residual
#fEvsResidual = book(TH1F('e_vs_residual_piEle', 'E vs Residual;True E;True - Reco E', 100, 0, 100))
#fEvsResidualStrict = book(TH1F('e_vs_residual_piEle_strict', 'E vs Residual;True E;True - Reco E', 100, 0, 100))

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
    fAllPtReg[r] = book(TH1F(f'mc_pion_pt_{r}', f'MC Pion p_{{T}} ({r})', PT_BINS, PT_MIN, PT_MAX))
    fMatchedPtReg[r] = book(TH1F(f'mc_matched_pt_piEle_{r}', f'Matched MC Pion to Reco #pi or e p_{{T}} ({r})', PT_BINS, PT_MIN, PT_MAX))
    fMatchedPtRegStrict[r] = book(TH1F(f'mc_matched_pt_piEle_{r}_strict', f'(Strict) Matched MC Pion to Reco #pi or e p_{{T}} ({r})', PT_BINS, PT_MIN, PT_MAX))

    # Regional resolutions
    fResPtReg[r] = TH1F(f'resolution_pt_piEle_{r}', 'p_{T} Resolution ('+r+');(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
    fResPtReg[r].SetMaximum(regional_max[regions.index(r)])
    fResEReg[r] = TH1F(f'resolution_e_piEle_{r}', 'Energy Resolution (' + r + ');(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
    fResEReg[r].SetMaximum(regional_max[regions.index(r)])
    fResPtRegStrict[r] = TH1F(f'resolution_pt_piEle_{r}_strict', '(Strict) p_{T} Resolution ('+r+');(True - Reco p_{T})/(True p_{T});Entries', 100, -0.1, 0.1)
    fResPtRegStrict[r].SetMaximum(regional_max[regions.index(r)])
    fResERegStrict[r] = TH1F(f'resolution_e_piEle_{r}_strict', '(Strict) Energy Resolution ('+r+');(True - Reco E)/(True E);Entries', 100, -0.1, 0.1)
    fResERegStrict[r].SetMaximum(regional_max[regions.index(r)])

# Get theta regions and eta
def eta(theta):
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
    return regs if regs else None

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

        # Best MC pion
        best_mc = None
        best_mc_pt = -1

        for mc in mcs:
            if abs(mc.getPDG()) != 211:
                continue

            mcMom = mc.getMomentum()
            mcPt = math.hypot(mcMom[0], mcMom[1])

            if mcPt > best_mc_pt:
                best_mc_pt = mcPt
                best_mc = mc

        if best_mc is None:
            continue

        mcMom = best_mc.getMomentum()
        mcPt = best_mc_pt
        mcTheta = math.acos(mcMom[2] / math.sqrt(mcPt**2 + mcMom[2]**2))
        mcEta = eta(mcTheta)
        mcPhi = math.atan2(mcMom[1], mcMom[0])
        mcE = best_mc.getEnergy()

        regs = theta_region(mcTheta)

        fMCPt.Fill(mcPt)
        fMCEta.Fill(abs(mcEta))
        fMCTheta.Fill(mcTheta)
        fMCPhi.Fill(mcPhi)
        fMCE.Fill(mcE)

        if regs:
            for reg in regs:
                fAllPtReg[reg].Fill(mcPt)

        # Best reco pi or e
        best_reco = None
        best_reco_pt = -1

        for pfo in pfos:
            if abs(pfo.getType()) not in (211, 11):
                continue

            mom = pfo.getMomentum()
            pt = math.hypot(mom[0], mom[1])
            energy = pfo.getEnergy()

            theta = math.acos(mom[2] / math.sqrt(pt**2 + mom[2]**2))
            eta_val = eta(theta)
            phi = math.atan2(mom[1], mom[0])

            fAllPt.Fill(pt)
            fAllEta.Fill(abs(eta_val))
            fAllTheta.Fill(theta)
            fAllPhi.Fill(phi)
            fAllE.Fill(energy)

            if pt > best_reco_pt:
                best_reco_pt = pt
                best_reco = pfo

        if best_reco is None:
            continue

        recoMom = best_reco.getMomentum()
        recoPt = best_reco_pt
        recoE = best_reco.getEnergy()

        recoTheta = math.acos(
            recoMom[2] / math.sqrt(recoPt**2 + recoMom[2]**2)
        )
        recoEta = eta(recoTheta)
        recoPhi = math.atan2(recoMom[1], recoMom[0])

        # ----- ΔR matching -----
        dphi = delta_phi(mcPhi, recoPhi)
        dR = math.sqrt(dphi*dphi + (mcEta - recoEta)**2)

        if dR < 0.25:

            fMatchedPt.Fill(mcPt)
            fMatchedEta.Fill(abs(mcEta))
            fMatchedTheta.Fill(mcTheta)
            fMatchedPhi.Fill(mcPhi)
            fMatchedE.Fill(mcE)

            fResidualPt.Fill(mcPt - recoPt)
            fResidualE.Fill(mcE - recoE)

            if mcPt > 0:
                fResPt.Fill((mcPt - recoPt)/mcPt)
            if mcE > 0:
                fResE.Fill((mcE - recoE)/mcE)

            if regs:
                for reg in regs:
                    fMatchedPtReg[reg].Fill(mcPt)
                    fResPtReg[reg].Fill((mcPt - recoPt)/mcPt)
                    fResEReg[reg].Fill((mcE - recoE) / mcE)

        if dR < 0.1:
            fMatchedPtStrict.Fill(mcPt)
            fMatchedEtaStrict.Fill(abs(mcEta))
            fMatchedThetaStrict.Fill(mcTheta)
            fMatchedPhiStrict.Fill(mcPhi)
            fMatchedEStrict.Fill(mcE)

            fResidualPtStrict.Fill(mcPt - recoPt)
            fResidualEStrict.Fill(mcE - recoE)

            if mcPt > 0:
                fResPtStrict.Fill((mcPt - recoPt)/mcPt)
            if mcE > 0:
                fResEStrict.Fill((mcE - recoE)/mcE)

            if regs:
                for reg in regs:
                    fMatchedPtRegStrict[reg].Fill(mcPt)
                    fResPtRegStrict[reg].Fill((mcPt - recoPt)/mcPt)
                    fResERegStrict[reg].Fill((mcE - recoE)/recoE)

    reader.close()

# Eff plots
def make_eff(num, den, name, title, xtitle):
    eff = num.Clone(name)
    eff.Divide(num, den, 1, 1, 'B')
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
make_eff(fMatchedPt, fMCPt, 'pt_eff_piEle', 'Charged #pi (+e) Efficiency vs p_{T}', 'True p_{T} [GeV]')
make_eff(fMatchedEta, fMCEta, 'eta_eff_piEle', 'Charged Pion Efficiency vs |#eta|', 'True |#eta|')
make_eff(fMatchedTheta, fMCTheta, 'theta_eff_piEle', 'Charged Pion Efficiency vs #theta', '#theta [rad]')
make_eff(fMatchedPhi, fMCPhi, 'phi_eff_piEle', 'Charged Pion Efficiency vs #phi', '#phi [rad]')
make_eff(fMatchedE, fMCE, 'e_eff_piEle', 'Charged Pion Efficiency vs E', 'E [GeV]')

# Strict eff plots
make_eff(fMatchedPtStrict, fMCPt, 'pt_eff_piEle_strict', '(Strict) Charged #pi (+e) Efficiency vs p_{T}', 'True p_{T} [GeV]')
make_eff(fMatchedEtaStrict, fMCEta, 'eta_eff_piEle_strict', '(Strict) Charged Pion Efficiency vs |#eta|', 'True |#eta|')
make_eff(fMatchedThetaStrict, fMCTheta, 'theta_eff_piEle_strict', '(Strict) Charged Pion Efficiency vs #theta', '#theta [rad]')
make_eff(fMatchedPhiStrict, fMCPhi, 'phi_eff_piEle_strict', '(Strict) Charged Pion Efficiency vs #phi', '#phi [rad]')
make_eff(fMatchedEStrict, fMCE, 'e_eff_piEle_strict', '(Strict) Charged Pion Efficiency vs E', 'E [GeV]')

# Fit guassians
fit_gauss(fResPt)
fit_gauss(fResE)
fit_gauss(fResPtStrict)
fit_gauss(fResEStrict)

for r in regions:
    # Regional eff
    make_eff(fMatchedPtReg[r], fAllPtReg[r],f'pt_{r}_eff_piEle',f'Charged #pi (+e) Efficiency ({r})','True p_{T} [GeV]')
    make_eff(fMatchedPtRegStrict[r], fAllPtReg[r],f'pt_{r}_eff_piEle_strict',f'(Strict) Charged #pi (+e) Efficiency ({r})','True p_{T} [GeV]')

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
for h in hists:
    c = TCanvas()
    h.Draw()
    gPad.Update()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.04)
    text.SetTextAlign(13)
    text.DrawLatex(0.12, 0.93, "#bf{#it{MAIA Detector Concept}}")
    text.DrawLatex(0.12, 0.88, "Simulated #pi^{+} Gun")
    c.Update()

    c.SaveAs(h.GetName() + '.png')
    del c

