from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import os
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import math

from tau_mc_link import getLinkedMCTau, getDecayMode, getNRecoQPis

# Command line arguments
parser = ArgumentParser()
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')
parser.add_argument('--outputFile', type=str, default='info_1p-2p_reco_tau_combo.root')
args = parser.parse_args()

# read in file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(args.inputFile)

# ROOT histograms
hist_pt_1p = TH1F("pt_1prong", "1-prong tau pT; p_{T}_{reco} [GeV]; Events", 50, 0, 100)
hist_pt_2p = TH1F("pt_2prong", "2-prong tau pT; p_{T}_{reco} [GeV]; Events", 50, 0, 100)

hist_theta_1p = TH1F("theta_1prong", "1-prong tau theta; #theta_{reco} [rad]; Events", 50, 0, 3.2)
hist_theta_2p = TH1F("theta_2prong", "2-prong tau theta; #theta_{reco} [rad]; Events", 50, 0, 3.2)

hist_phi_1p = TH1F("phi_1prong", "1-prong tau phi; #phi_{reco} [rad]; Events", 50, 0, 3.2)
hist_phi_2p = TH1F("phi_2prong", "2-prong tau phi; #phi_{reco} [rad]; Events", 50, 0, 3.2)

# Table storage
event_table = []
event_index = 0

for event in reader:
    try:
        reco_taus = event.getCollection('RecoTaus')
        tauRecoLink = event.getCollection('TauPFOLink')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        # Skip events with fewer than 2 reco taus
        #if reco_taus.getNumberOfElements() < 2:
        #    continue

        # Count how many are from a 3-prong MC tau
        n_3p_mc = 0
        for reco_tau in reco_taus:
            mcTau = getLinkedMCTau(reco_tau, relationNavigatorTau, relationNavigatorRecoMC)
            if mcTau and getDecayMode(mcTau) == 4:
                n_3p_mc += 1

        if n_3p_mc < 1: # make sure both reco taus are from a 3-prong MC tau
            continue

        # Try to find a 1p and a 2p tau in the reco collection
        one_prong = None
        two_prong = None
        total_q_pis = 0

        for reco_tau in reco_taus:
            n_q_pis = getNRecoQPis(reco_tau)
            if n_q_pis == 1 and one_prong is None:
                one_prong = reco_tau
                total_q_pis += 1
            elif n_q_pis == 2 and two_prong is None:
                two_prong = reco_tau
                total_q_pis += 2

        if one_prong and two_prong and total_q_pis == 3:
            # 1p info
            p1 = one_prong.getMomentum()
            px1, py1, pz1 = p1[0], p1[1], p1[2]
            pt_1p = math.sqrt(px1**2 + py1**2)
            theta_1p = math.acos(pz1 / math.sqrt(px1**2 + py1**2 + pz1**2))
            phi_1p = abs(math.atan2(py1, px1))

            # 2p info
            p2 = two_prong.getMomentum()
            px2, py2, pz2 = p2[0], p2[1], p2[2]
            pt_2p = math.sqrt(px2**2 + py2**2)
            theta_2p = math.acos(pz2 / math.sqrt(px2**2 + py2**2 + pz2**2))
            phi_2p = abs(math.atan2(py2, px2))

            # Fill ROOT histograms
            hist_pt_1p.Fill(pt_1p)
            hist_pt_2p.Fill(pt_2p)
            hist_theta_1p.Fill(theta_1p)
            hist_theta_2p.Fill(theta_2p)
            hist_phi_1p.Fill(phi_1p)
            hist_phi_2p.Fill(phi_2p)

            # Add to table
            event_table.append(f"{event_index:<8} | {pt_1p:7.2f} | {pt_2p:7.2f} | {theta_1p:8.3f} | {theta_2p:8.3f} | {phi_1p:8.3f} | {phi_2p:8.3f}")
            event_index += 1

    except Exception as e:
        print(f"Error in event {event_index}: {e}")
        continue

# Save ROOT histograms to .root file
output_file = TFile(args.outputFile, "RECREATE")
hist_pt_1p.Write()
hist_pt_2p.Write()
hist_theta_1p.Write()
hist_theta_2p.Write()
hist_phi_1p.Write()
hist_phi_2p.Write()
output_file.Close()

# Save ROOT histograms as PNGs
canvas = TCanvas("c", "c", 800, 600)

def save_hist_as_png(hist, filename):
    canvas.Clear()
    hist.Draw()
    canvas.SaveAs(filename)

save_hist_as_png(hist_pt_1p, "pt_1prong.png")
save_hist_as_png(hist_pt_2p, "pt_2prong.png")
save_hist_as_png(hist_theta_1p, "theta_1prong.png")
save_hist_as_png(hist_theta_2p, "theta_2prong.png")
save_hist_as_png(hist_phi_1p, "phi_1prong.png")
save_hist_as_png(hist_phi_2p, "phi_2prong.png")

# Create matplotlib table
header = ["Event #", "1p reco pt", "2p reco pt", "1p reco theta", "2p reco theta", "1p reco phi", "2p reco phi"]
table_data = [row.split('|') for row in event_table]
table_data = [[col.strip() for col in row] for row in table_data]

fig, ax = plt.subplots(figsize=(10, len(table_data) * 0.4 + 1))
ax.axis('off')

table = plt.table(cellText=table_data,
                  colLabels=header,
                  loc='center',
                  cellLoc='center')

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.2)

plt.savefig("tau_reco_event_table.png", bbox_inches='tight')
plt.close()

