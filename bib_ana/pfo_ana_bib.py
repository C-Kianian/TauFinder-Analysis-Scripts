# Comes from ethanmar/MuColl-TauStudy/analysis/pfo_matching.py
# Run this code after the reco step on a pion gun
from pyLCIO import IOIMPL
from ROOT import TH1F, TFile, TCanvas, TLegend, gPad
import math
from argparse import ArgumentParser
import os
import numpy as np

# args
parser = ArgumentParser()
parser.add_argument('-i', '--inputFile', type=str, default='output_reco.slcio')
parser.add_argument('-o', '--outputFile', type=str, default='pfo_ana_bib.root')
args = parser.parse_args()

# function to get angle between two particles
def angle_between(p1, p2):
    # np.array copies the C++ data immediately
    v1 = np.array(p1)
    v2 = np.array(p2)

    mag1 = np.linalg.norm(v1)
    mag2 = np.linalg.norm(v2)

    if mag1 == 0 or mag2 == 0:
        return float("inf")

    dot = np.dot(v1, v2)
    cos_theta = np.clip(dot / (mag1 * mag2), -1.0, 1.0)
    return np.arccos(cos_theta)


# histogram bookkeeping
hists = []
def book(hist):
    hist.SetDirectory(0)
    hists.append(hist)
    return hist

# theta hist
hPfoTheta = book(TH1F(
    "pfo_theta",
    "Number of Reconstructed PFOs vs #theta;#theta [rad];Counts",
    50, 0, math.pi
))

# PDG bookkeeping
pdg_counts = {} # pdg counts
pdg_close_counts = {} # close pdg counts
hEnergy = {} # pdg energy histogram
hEnergyClose = {} # close pdg energy histogram


ENERGY_BINS = 10000
ENERGY_MAX = 10000

total_pfos = 0
total_close_pfos = 0

# simple charged pion and electron counters
counters = {
    "pi_no_e": 0,
    "e_no_pi": 0,
    "pi_and_e": 0,
    "highest_pt_charged_pi": 0,
    "highest_pt_charged_e": 0,
    "any_pfo": 0
}

# get input files
to_process = []
if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)

# event loop
for file in to_process:

    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:

        MCParticles = event.getCollection('MCParticle')
        if not MCParticles:
            continue
        MCPis = [mc for mc in MCParticles if mc.getPDG() == 211] # get only mc pions
        if not MCPis:
            continue

         # get highest pt pion
        best_mc_pt = -1.0
        MCPi = None

        for mc in MCPis:
            # Compute pt
            MC_mom = mc.getMomentum()
            px, py, pz = MC_mom[0], MC_mom[1], MC_mom[2]
            pt = math.hypot(px, py)

            # Check if this is the highest pt so far
            if pt > best_mc_pt:
                best_mc_pt = pt
                MCPi = mc

        # pfos
        pfos = event.getCollection('PandoraPFOs')
        if len(pfos) == 0:
            continue

        pfos_close = []
        if MCPi is not None:
            mcMom = MCPi.getMomentum()
            mc_px, mc_py, mc_pz = mcMom[0], mcMom[1], mcMom[2]
            if mcMom:  # Only proceed if MCPi momentum is valid
                for pfo in pfos:
                    pfoMom = pfo.getMomentum()
                    pfo_px, pfo_py, pfo_pz = pfoMom[0], pfoMom[1], pfoMom[2]
                    if pfoMom:  # Only proceed if PFO momentum is valid
                        try:
                            ang = angle_between((mc_px, mc_py, mc_pz), (pfo_px, pfo_py, pfo_pz))
                            if ang < 0.1:
                                pfos_close.append(pfo)
                        except Exception as e:
                            print("Skipping PFO due to error:", e)


        total_pfos += len(pfos)
        total_close_pfos += len(pfos_close)

        counters["any_pfo"] += 1

        pi_pfos = []
        e_pfos = []

        highest_pt = 0
        highest_pt_type = None

        for pfo in pfos:

            mom = pfo.getMomentum()
            px, py, pz = mom[0], mom[1], mom[2]
            pt = math.hypot(px, py)
            theta = math.atan2(pt, pz)

            pdg = abs(pfo.getType())
            energy = pfo.getEnergy()

            # Theta
            hPfoTheta.Fill(theta)

            # PDG counts
            pdg_counts[pdg] = pdg_counts.get(pdg, 0) + 1
            pdg_close_counts[pdg] = pdg_close_counts.get(pdg, 0) + (1 if pfo in pfos_close else 0)

            # energy hists
            if pdg not in hEnergy:
                hEnergy[pdg] = book(
                    TH1F(f"energy_pdg_{pdg}",
                         f"PFO Energy PDG {pdg};Energy [GeV];Counts",
                         ENERGY_BINS, 0, ENERGY_MAX)
                )
            if pdg not in hEnergyClose and pfo in pfos_close:
                hEnergyClose[pdg] = book(
                    TH1F(f"energy_close_pdg_{pdg}",
                         f"Close PFO Energy PDG {pdg};Energy [GeV];Counts",
                         ENERGY_BINS, 0, ENERGY_MAX)
                )

            hEnergy[pdg].Fill(energy)
            if pfo in pfos_close: hEnergyClose[pdg].Fill(energy)

            # pi vs electron events
            if pdg == 211:
                pi_pfos.append(pt)
                if pt > highest_pt:
                    highest_pt = pt
                    highest_pt_type = pdg
            elif pdg == 11:
                e_pfos.append(pt)
                if pt > highest_pt:
                    highest_pt = pt
                    highest_pt_type = pdg

        # event categories
        if pi_pfos and not e_pfos:
            counters["pi_no_e"] += 1
        elif e_pfos and not pi_pfos:
            counters["e_no_pi"] += 1
        elif pi_pfos and e_pfos:
            counters["pi_and_e"] += 1

        if highest_pt_type == 211:
            counters["highest_pt_charged_pi"] += 1
        elif highest_pt_type == 11:
            counters["highest_pt_charged_e"] += 1

    reader.close()

# PDG count bar histogram
pdgs_sorted = sorted(pdg_counts.keys(), reverse=True)
pdgs_close_sorted = sorted(pdg_close_counts.keys(), reverse=True)

hPdgCounts = book(
    TH1F("pdg_counts",
         "PFO Counts by PDG;PDG ID;Counts",
         len(pdgs_sorted), 0, len(pdgs_sorted))
)
hPdgCountsClose = book(
    TH1F("pdg_counts_close",
         "Close PFO Counts by PDG;PDG ID;Counts",
         len(pdgs_close_sorted), 0, len(pdgs_close_sorted))
)

hPdgCounts.GetYaxis().SetTitle("Counts")
hPdgCounts.GetXaxis().SetTitle("PDG ID")
hPdgCountsClose.GetYaxis().SetTitle("Counts")
hPdgCountsClose.GetXaxis().SetTitle("PDG ID")

for i, pdg in enumerate(pdgs_sorted, start=1):
    hPdgCounts.SetBinContent(i, pdg_counts[pdg])
    hPdgCounts.GetXaxis().SetBinLabel(i, str(pdg))

for i, pdg in enumerate(pdgs_close_sorted, start=1):
    hPdgCountsClose.SetBinContent(i, pdg_close_counts[pdg])
    hPdgCountsClose.GetXaxis().SetBinLabel(i, str(pdg))

# Overlay ALL PDG energy spectra
c_energy = TCanvas("energy_overlay", "PFO Energy Overlay", 900, 700)
c_energy.SetLogy(True)

legend = TLegend(0.62, 0.60, 0.88, 0.88)
legend.SetBorderSize(0)

color = 1
first = True

for pdg in pdgs_sorted:
    h = hEnergy[pdg]
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetStats(0)
    h.SetMinimum(0.5)
    h.GetXaxis().SetTitle("Energy [GeV]")
    h.GetYaxis().SetTitle("Counts")

    legend.AddEntry(h, f"PDG {pdg}", "l")

    if first:
        h.Draw()
        h.GetXaxis().SetTitle("Energy [GeV]")
        h.GetYaxis().SetTitle("Counts")
        first = False
    else:
        h.Draw("SAME")

    color += 1

legend.Draw()
c_energy.SaveAs("energy_overlay_all_pdgs_log.png")

# Overlay close PDG energy spectra
c_energy_close = TCanvas("energy_overlay_close", "Close PFO Energy Overlay", 900, 700)
c_energy_close.SetLogy(True)

legend = TLegend(0.62, 0.60, 0.88, 0.88)
legend.SetBorderSize(0)

color = 1
first = True

for pdg in pdgs_close_sorted:
    h = hEnergyClose[pdg]
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetStats(0)
    h.SetMinimum(0.5)
    h.GetXaxis().SetTitle("Energy [GeV]")
    h.GetYaxis().SetTitle("Counts")

    legend.AddEntry(h, f"PDG {pdg}", "l")

    if first:
        h.Draw()
        h.GetXaxis().SetTitle("Energy [GeV]")
        h.GetYaxis().SetTitle("Counts")
        first = False
    else:
        h.Draw("SAME")

    color += 1

legend.Draw()
c_energy_close.SaveAs("energy_overlay_close_pdgs_log.png")

# Save PDG bar chart
c_pdg = TCanvas("pdg_counts", "PDG Counts", 800, 600)
c_pdg.SetLogy(True)
hPdgCounts.SetFillColor(38)
hPdgCounts.SetStats(0)
hPdgCounts.Draw("BAR")
c_pdg.SaveAs("pdg_counts.png")

# Save close PDG bar chart
c_pdg_close = TCanvas("pdg_counts_close", "PDG Counts", 800, 600)
c_pdg_close.SetLogy(True)
hPdgCountsClose.SetFillColor(38)
hPdgCountsClose.SetStats(0)
hPdgCountsClose.Draw("BAR")
c_pdg_close.SaveAs("pdg_counts_close.png")

# normalize PDG bar chart
hPdgCountsNorm = hPdgCounts.Clone("hPdgCountsNorm")
hPdgCountsNorm.Scale(1.0/total_pfos)
hPdgCountsNorm.GetYaxis().SetRangeUser(0, 1.1)
c_pdg_norm = TCanvas("pdg_counts_norm", "PDG Counts Normalized", 800, 600)
hPdgCountsNorm.Draw("BAR")
c_pdg_norm.SaveAs("pdg_counts_norm.png")

# normalize close PDG bar chart
hPdgCountsCloseNorm = hPdgCountsClose.Clone("hPdgCountsCloseNorm")
hPdgCountsCloseNorm.Scale(1.0/total_close_pfos)
hPdgCountsCloseNorm.GetYaxis().SetRangeUser(0, 1.1)
c_pdg_norm = TCanvas("pdg_counts_close_norm", "Close PDG Counts Normalized", 800, 600)
hPdgCountsCloseNorm.Draw("BAR")
c_pdg_norm.SaveAs("pdg_counts_close_norm.png")

# Save theta plot
c_theta = TCanvas("theta", "PFO Theta", 800, 600)
hPfoTheta.Draw()
c_theta.SaveAs("pfo_theta.png")

# Output ROOT file
output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
    hist.Write()
output_file.Close()

# Print counters
print("\n--- Charged pion → electron mis-ID study ---")
for k, v in counters.items():
    print(f"{k}: {v}")

