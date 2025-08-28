from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import ROOT
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np

from tau_mc_link import getLinkedMCTau, getVisibleProperties, getDecayMode, getNRecoNeutralPis, getNRecoQPis

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

args = parser.parse_args()

# Check if input file is a directory or a single file
to_process = []

decay_modes = {
    0: '1P0N',
    1: '1P + N',
    2: '1P + N',
    3: '1P + N',
    5: '3P + N'
}

hists_dict = {}

for decay_num, decay_mode in decay_modes.items():
    if decay_num == 2 or decay_num == 3: continue

    # Reco Photons by Reco Pt
    hRecoPhotonRecoPt = TH1F(f"{decay_mode}_pt_reco_photons", f"{decay_mode} Reco Photons", 20, 0, 320)
    hRecoPhotonRecoPt.GetXaxis().SetTitle('p_{T} [GeV/c]')

    # truth photon pt from true pi0
    hTruePhotonsFromPi0Pt = TH1F(f"{decay_mode}_pt_true_photons_from_true_pi0", f"{decay_mode} True Photons From True Pi0", 20, 0, 320)
    hTruePhotonsFromPi0Pt.GetXaxis().SetTitle('p_{T} [GeV/c]')

    # truth photons pt not from true pi0
    hTruePhotonsNotFromPi0Pt = TH1F(f"{decay_mode}_pt_true_photons_not_from_pi0", f"{decay_mode} True Photons Not From Pi0", 20, 0, 320)
    hTruePhotonsNotFromPi0Pt.GetXaxis().SetTitle('p_{T} [GeV/c]')

    # add to dictionary
    hists_dict[decay_num] = {
        'reco_gamma_pt': hRecoPhotonRecoPt,
        'true_pi0_photon_pt' : hTruePhotonsFromPi0Pt,
        'true_not_pi0_photon_pt' : hTruePhotonsNotFromPi0Pt
    }

bins = ['1P0N', '1P + N', '3P0N']
n_bins = len(bins)
mode_to_bin_dict = {
    0: 1,  # '1P0N'
    1: 2,  # '1P + N' (also used for 2 and 3)
    5: 3   # '3P0N'
}

def style_general_hists(h, bin_names):
    h.GetXaxis().SetTitle('Decay Mode')
    h.GetYaxis().SetTitle('# of Entries')
    for i, label in enumerate(bin_names, 1):
        h.GetXaxis().SetBinLabel(i, label)


#general counting hists
# true not p0 photons
hNTruePhotonNotFromPi0 = TH1F(f"n_true_photons_not_from_pi0", f"N True Photons Direct Daughters of Tau", n_bins, 0, n_bins)
style_general_hists(hNTruePhotonNotFromPi0, bins)
# true pi0 photons
hNTruePhotonFromPi0 = TH1F(f"n_true_photons_from_pi0", f"N True Photons From Pi0", n_bins, 0, n_bins)
style_general_hists(hNTruePhotonFromPi0, bins)
# true pi_0
hNPi0True = TH1F(f"n_true_pi0", f"N True Pi0s", n_bins, 0, n_bins)
style_general_hists(hNPi0True, bins)
# reco photons
hNPhotonReco = TH1F(f"n_reco_photons", f"N Reconstructed Photons", n_bins, 0, n_bins)
style_general_hists(hNPhotonReco, bins)

if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)

for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for ievt, event in enumerate(reader):

        reco_taus = event.getCollection('TauRec_PFO')
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauRecLink_PFO')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        mc_taus = [t for t in mcParticles if abs(t.getPDG()) == 15]
        reco_gamma_list = [pfo for pfo in pfos if pfo.getType() == 22]
        linked_gamma_list = []

        # Loop through MC taus
        for mc_tau in mc_taus:

            # Get tau decay mode
            decayMode = getDecayMode(mc_tau)

            if decayMode in decay_modes:
                if decayMode == 2 or decayMode == 3: decayMode = 1

                #truth pi0
                tau_daughters = mc_tau.getDaughters()
                true_pi0s = [d for d in tau_daughters if abs(d.getPDG()) == 111]
                hNPi0True.Fill(mode_to_bin_dict[decayMode], len(true_pi0s))

                #truth photons from pi0
                true_pi0_gammas = [] # get mc photons from mc pi0
                for mc_pi0 in true_pi0s:
                    for d in mc_pi0.getDaughters():
                        if abs(d.getPDG()) == 22:
                            true_pi0_gammas.append(d)
                            p = d.getMomentum() # get pt info
                            px, py, pz = p[0], p[1], p[2]
                            pt = math.sqrt(px**2 + py**2)
                            hists_dict[decayMode]['true_pi0_photon_pt'].Fill(pt)
                hNTruePhotonFromPi0.Fill(mode_to_bin_dict[decayMode], len(true_pi0_gammas))

                true_photons = [d for d in tau_daughters if abs(d.getPDG()) == 22]
                true_gammas_not_from_pi0 = [g for g in true_photons if g not in true_pi0_gammas] # all photons in tau not from pi0, only direct daughters
                hNTruePhotonNotFromPi0.Fill(mode_to_bin_dict[decayMode], len(true_gammas_not_from_pi0))

                for g in true_gammas_not_from_pi0:
                    p = g.getMomentum() # get pt info
                    px, py, pz = p[0], p[1], p[2]
                    pt = math.sqrt(px**2 + py**2)
                    hists_dict[decayMode]['true_not_pi0_photon_pt'].Fill(pt)

            # Find reco taus linked to this MC tau
            linked_reco_taus = relationNavigatorTau.getRelatedToObjects(mc_tau)

            # For each linked reco tau, find linked PFOs and count photons
            n_reco_photons = 0

            for reco_tau in linked_reco_taus:
                # Find PFOs linked to this reco tau
                linked_pfos = relationNavigatorTau.getRelatedToObjects(reco_tau)
                for pfo in linked_pfos:
                    if pfo.getType() == 22:  # photon
                        n_reco_photons += 1

                        # Fill reco photon pt histogram
                        p = pfo.getMomentum()
                        px, py, pz = p[0], p[1], p[2]
                        pt = math.sqrt(px**2 + py**2)
                        hists_dict[decayMode]['reco_gamma_pt'].Fill(pt)

            # Fill general reco photon count histogram
            if decayMode in mode_to_bin_dict:
                hNPhotonReco.Fill(mode_to_bin_dict[decayMode], n_reco_photons)


    reader.close()

    # overlay the general hists
    # Style the histograms so they are distinguishable
    hNTruePhotonNotFromPi0.SetLineColor(ROOT.kRed)
    hNTruePhotonNotFromPi0.SetLineWidth(2)

    hNTruePhotonFromPi0.SetLineColor(ROOT.kBlue)
    hNTruePhotonFromPi0.SetLineWidth(2)

    hNPi0True.SetLineColor(ROOT.kGreen+2)
    hNPi0True.SetLineWidth(2)

    hNPhotonReco.SetLineColor(ROOT.kBlack)
    hNPhotonReco.SetLineWidth(2)

    # Create canvas
    c_general = ROOT.TCanvas("c_general", "General Photon Counts", 800, 600)

    # Draw first histogram, then overlay others
    hNTruePhotonNotFromPi0.Draw()
    hNTruePhotonFromPi0.Draw("SAME")
    hNPi0True.Draw("SAME")
    hNPhotonReco.Draw("SAME")

    # Create and configure legend
    legend = ROOT.TLegend(0.65, 0.7, 0.9, 0.9)
    legend.AddEntry(hNTruePhotonNotFromPi0, "True Photons Not from Pi0", "l")
    legend.AddEntry(hNTruePhotonFromPi0, "True Photons from Pi0", "l")
    legend.AddEntry(hNPi0True, "True Pi0s", "l")
    legend.AddEntry(hNPhotonReco, "Reconstructed Photons", "l")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    c_general.Update()
    c_general.SaveAs("general_photon_counts_overlay.png")


    #overlay decay specific hists
    for decay_num, decay_mode in decay_modes.items():
        if decay_num == 2 or decay_num == 3:
            continue

        h_reco = hists_dict[decay_num]['reco_gamma_pt']
        h_true_pi0 = hists_dict[decay_num]['true_pi0_photon_pt']
        h_true_not_pi0 = hists_dict[decay_num]['true_not_pi0_photon_pt']

        # Set axis titles and ranges
        h_reco.GetXaxis().SetTitle('p_{T} [GeV/c]')
        h_reco.GetXaxis().SetRangeUser(0, 320)
        h_reco.SetLineColor(ROOT.kBlue)
        h_reco.SetLineWidth(2)

        h_true_pi0.SetLineColor(ROOT.kRed)
        h_true_pi0.SetLineWidth(2)

        h_true_not_pi0.SetLineColor(ROOT.kGreen+2)
        h_true_not_pi0.SetLineWidth(2)

        # Draw
        c = ROOT.TCanvas(f"c_{decay_mode}", f"Photon pT Overlay {decay_mode}", 800, 600)
        h_reco.Draw()
        h_true_pi0.Draw("SAME")
        h_true_not_pi0.Draw("SAME")

        legend = ROOT.TLegend(0.65, 0.7, 0.9, 0.9)
        legend.AddEntry(h_reco, "Reco Photons", "l")
        legend.AddEntry(h_true_pi0, "True Photons from Pi0", "l")
        legend.AddEntry(h_true_not_pi0, "True Photons Not from Pi0", "l")
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()

        c.Update()
        c.SaveAs(f"{decay_mode}_photon_pt_overlay.png")
        c.Close()


    output_file = TFile("output_histograms.root", "RECREATE")

    # Write general histograms
    hNTruePhotonNotFromPi0.Write()
    hNTruePhotonFromPi0.Write()
    hNPi0True.Write()
    hNPhotonReco.Write()

    # Write decay mode specific histograms
    for decay_num in hists_dict:
        for hist in hists_dict[decay_num].values():
            hist.Write()

    output_file.Close()
