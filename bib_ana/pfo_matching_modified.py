from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
import os

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------
parser = ArgumentParser()

parser.add_argument(
    '-i', '--inputFile',
    help='--inputFile output_reco.slcio',
    type=str,
    default='output_reco.slcio'
)

parser.add_argument(
    '-o', '--outputFile',
    help='--outputFile pfo_matching.root',
    type=str,
    default='pfo_matching.root'
)

args = parser.parse_args()

# --------------------------------------------------
# Histogram configuration (FIXED, NOT DYNAMIC)
# --------------------------------------------------
PT_MIN  = 0.0
PT_MAX  = 1000.0   # 1 TeV
PT_BINS = 40

ETA_MIN  = 0.0
ETA_MAX  = 2.5
ETA_BINS = 20

EFF_MIN = 0.0
EFF_MAX = 1.0

# --------------------------------------------------
# Histograms
# --------------------------------------------------
hists = []

fMatchedPt = TH1F(
    'matched_pt',
    'Matched PFO True Pt',
    PT_BINS,
    PT_MIN,
    PT_MAX
)
hists.append(fMatchedPt)

fMatchedEta = TH1F(
    'matched_eta',
    'Matched PFO True Eta',
    ETA_BINS,
    ETA_MIN,
    ETA_MAX
)
hists.append(fMatchedEta)

fResidualPt = TH1F(
    'res_pt',
    'Residual Transverse Momentum',
    20,
    -100,
    100
)
fResidualPt.SetXTitle('True - Reco Pt [GeV/c]')
hists.append(fResidualPt)

fAllPt = TH1F(
    'all_pt',
    'All True Pt',
    PT_BINS,
    PT_MIN,
    PT_MAX
)
hists.append(fAllPt)

fAllEta = TH1F(
    'all_eta',
    'All True Eta',
    ETA_BINS,
    ETA_MIN,
    ETA_MAX
)
hists.append(fAllEta)

for hist in hists:
    hist.SetDirectory(0)

# --------------------------------------------------
# Utility functions
# --------------------------------------------------
def eta(theta):
    return -math.log(math.tan(theta / 2.0))

# --------------------------------------------------
# Collect input files
# --------------------------------------------------
to_process = []

if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)

# --------------------------------------------------
# Event loop
# --------------------------------------------------
for file in to_process:

    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for ievt, event in enumerate(reader):

        pfos = event.getCollection('PandoraPFOs')
        mcs  = event.getCollection('MCParticle')

        # ----------------------------------------------
        # Find highest-pT charged MC pion
        # ----------------------------------------------
        best_mc = None
        best_mc_pt = -1.0

        for mc in mcs:
            if abs(mc.getPDG()) != 211:
                continue

            mcMom = mc.getMomentum()
            mcPt = math.sqrt(mcMom[0]**2 + mcMom[1]**2)

            if mcPt > best_mc_pt:
                best_mc_pt = mcPt
                best_mc = mc

        if best_mc is None:
            continue

        # MC kinematics
        mcMom = best_mc.getMomentum()
        mcPt = best_mc_pt
        mcTheta = math.acos(mcMom[2] / math.sqrt(mcPt**2 + mcMom[2]**2))
        mcEta = eta(mcTheta)
        mcPhi = math.atan2(mcMom[1], mcMom[0])

        # Fill denominators
        fAllPt.Fill(mcPt)
        fAllEta.Fill(abs(mcEta))

        # ----------------------------------------------
        # Find highest-pT charged PFO pion
        # ----------------------------------------------
        best_pfo = None
        best_pfo_pt = -1.0

        for pfo in pfos:
            if abs(pfo.getType()) != 211:
                continue

            pfoMom = pfo.getMomentum()
            pfoPt = math.sqrt(pfoMom[0]**2 + pfoMom[1]**2)

            if pfoPt > best_pfo_pt:
                best_pfo_pt = pfoPt
                best_pfo = pfo

        if best_pfo is None:
            continue

        # PFO kinematics
        pfoMom = best_pfo.getMomentum()
        pfoPt = best_pfo_pt
        pfoTheta = math.acos(pfoMom[2] / math.sqrt(pfoPt**2 + pfoMom[2]**2))
        pfoEta = eta(pfoTheta)
        pfoPhi = math.atan2(pfoMom[1], pfoMom[0])

        # ----------------------------------------------
        # ΔR matching
        # ----------------------------------------------
        dR = math.sqrt((mcPhi - pfoPhi)**2 + (mcEta - pfoEta)**2)

        if dR < 0.25:
            fMatchedPt.Fill(mcPt)
            fMatchedEta.Fill(abs(mcEta))
            fResidualPt.Fill(mcPt - pfoPt)

    reader.close()

# --------------------------------------------------
# Efficiency plots (FIXED AXES)
# --------------------------------------------------
fPiPtEff = fMatchedPt.Clone('pt_eff')
fPiPtEff.Divide(fPiPtEff, fAllPt, 1, 1, 'B')
fPiPtEff.SetLineColor(6)
fPiPtEff.SetLineWidth(2)
fPiPtEff.SetTitle('Charged Pion Efficiency vs Pt')
fPiPtEff.GetXaxis().SetTitle('True pT [GeV/c]')
fPiPtEff.GetYaxis().SetTitle('#epsilon')
fPiPtEff.GetYaxis().SetRangeUser(EFF_MIN, EFF_MAX)
fPiPtEff.GetXaxis().SetRangeUser(PT_MIN, PT_MAX)
fPiPtEff.SetStats(0)
hists.append(fPiPtEff)

fPiEtaEff = fMatchedEta.Clone('eta_eff')
fPiEtaEff.Divide(fPiEtaEff, fAllEta, 1, 1, 'B')
fPiEtaEff.SetLineColor(7)
fPiEtaEff.SetLineWidth(2)
fPiEtaEff.SetTitle('Charged Pion Efficiency vs Eta')
fPiEtaEff.GetXaxis().SetTitle('True #eta')
fPiEtaEff.GetYaxis().SetTitle('#epsilon')
fPiEtaEff.GetYaxis().SetRangeUser(EFF_MIN, EFF_MAX)
fPiEtaEff.GetXaxis().SetRangeUser(ETA_MIN, ETA_MAX)
fPiEtaEff.SetStats(0)
hists.append(fPiEtaEff)

# --------------------------------------------------
# Output ROOT file
# --------------------------------------------------
output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
    hist.Write()
output_file.Close()

# --------------------------------------------------
# Save plots (non-dynamic axes)
# --------------------------------------------------
for hist in hists:
    canvas = TCanvas()
    hist.Draw('HIST')
    canvas.SaveAs(hist.GetName() + '.png')

