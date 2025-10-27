from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas, TLegend
import os
from argparse import ArgumentParser
import matplotlib.pyplot as plt

from tau_mc_link import getLinkedMCTau, getDecayMode, getNRecoQPis

# Command line arguments
parser = ArgumentParser()
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')
parser.add_argument('--outputFile', type=str, default='reco_tau_prong_combinations.root')
args = parser.parse_args()

# Discover input files
to_process = []
if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)

hist_list = []
counter_table = {}
all_prongs = set()
file_prong_dicts = {}

for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    n_reco_3p_events = 0
    n_reco_matched_3p = 0
    prong_dict = {}

    for event in reader:
        # Get collections
        reco_taus = event.getCollection('RecoTaus')
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauPFOLink')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        prong_list = []
        counted = False

        for reco_tau in reco_taus:
            mcTau = getLinkedMCTau(reco_tau, relationNavigatorTau, relationNavigatorRecoMC)
            decayMode = getDecayMode(mcTau)

            if decayMode == 4 and not counted:
                n_reco_3p_events += 1
                n_reco_matched_3p += 1
                counted = True
                prong_list.append(getNRecoQPis(reco_tau))
            elif decayMode == 4:
                prong_list.append(getNRecoQPis(reco_tau))
                n_reco_matched_3p += 1

        if prong_list:
            prong_tuple = tuple(sorted(prong_list))
            prong_dict[prong_tuple] = prong_dict.get(prong_tuple, 0) + 1

    # Store results
    file_prong_dicts[file] = {
        "prong_dict": prong_dict,
        "n_reco_3p_events": n_reco_3p_events,
        "n_reco_matched_3p": n_reco_matched_3p
    }
    all_prongs.update(prong_dict.keys())

# All prong combinations from all files
all_prongs = sorted(all_prongs)
bin_labels = {prong: i+1 for i, prong in enumerate(all_prongs)}

for file, results in file_prong_dicts.items():
    prong_dict = results["prong_dict"]

    hProngCombinations = TH1F(
        f'{os.path.basename(file)}',
        'Different Prong Combinations',
        len(all_prongs), 0, len(all_prongs)
    )

    # set bin labels consistently
    for prong, bin_id in bin_labels.items():
        label = "-".join(map(str, prong))
        hProngCombinations.GetXaxis().SetBinLabel(bin_id, label)

    # fill contents for this file
    for prong_combination, n_occur in prong_dict.items():
        bin_id = bin_labels[prong_combination]
        hProngCombinations.SetBinContent(bin_id, n_occur)

    hist_list.append(hProngCombinations)
    counter_table[file] = {
        "n_reco_3p_events": results["n_reco_3p_events"],
        "n_reco_matched_3p": results["n_reco_matched_3p"]
    }

# Save everything in one ROOT file
outFile = TFile(args.outputFile, "RECREATE")
for h in hist_list:
    h.Write()
outFile.Close()

# Output hists
i = 0
names = ["Static (0.1 rad) Merge", "No Merge"]
for h in hist_list:
    c = TCanvas("c", h.GetName(), 800, 600)
    c.SetLogy(True)

    # Clean style
    h.GetYaxis().SetRangeUser(0.1, 1500)
    h.SetTitle(f"{names[i]} Prong Combinations")
    h.SetLineColor(i+1)
    h.SetLineWidth(2)
    h.SetFillColor(0)
    h.SetStats(0)
    h.GetXaxis().SetTitle(f"Reco Tau Prongs")
    h.GetYaxis().SetTitle("Entries")

    # Draw and save
    h.Draw("HIST")
    c.SaveAs(f"{names[i]}.png")  # save as PNG
    i += 1
    del c

# Overlay histograms
if len(hist_list) > 1:
    c = TCanvas("c", "Overlay", 800, 600)
    c.SetLogy(True)
    legend = TLegend(0.15, 0.75, 0.55, 0.9)

    for i, h in enumerate(hist_list):
        h.SetLineColor(i+1)
        h.SetLineWidth(2)
        h.SetStats(0)
        if i == 0:
            h.SetTitle("Different Prong Combinations")
            h.GetXaxis().SetTitle("Reco Tau Prongs")
            h.GetYaxis().SetTitle("Entries")
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")
        legend.AddEntry(h, names[i], "l")

    legend.Draw()
    c.SaveAs("overlay_histograms.png")

# Build table data
headers = ["Name", "# Reco 3P Events", "# Reco $\\tau^{-}$ Matched to MC 3P"]

file_keys = list(counter_table.keys())
rows = []

for i in range(len(names)):
    name = names[i]
    # get counts from corresponding file
    counts_dict = counter_table[file_keys[i]]
    rows.append([name, counts_dict["n_reco_3p_events"], counts_dict["n_reco_matched_3p"]])


# Plot table
fig, ax = plt.subplots()
ax.axis('off')
table = ax.table(cellText=rows, colLabels=headers, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.2)

plt.savefig("counter_table.png", dpi=200, bbox_inches="tight")

