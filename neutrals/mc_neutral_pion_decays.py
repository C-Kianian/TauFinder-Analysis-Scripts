from pyLCIO import IOIMPL, EVENT, UTIL
import ROOT
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt

from tau_mc_link import getDecayMode

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='pi_0_decay_ana.root')

args = parser.parse_args()

# Check if input file is a directory or a single file
to_process = []

if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)


#calculate the angle between two photons
def calc_angle(gamma_list):
    g1 = gamma_list[0]
    g2 = gamma_list[1]

    p1 = g1.getMomentum()
    p2 = g2.getMomentum()

    dot = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]
    norm1 = math.sqrt(p1[0]**2 + p1[1]**2 + p1[2]**2)
    norm2 = math.sqrt(p2[0]**2 + p2[1]**2 + p2[2]**2)

    cos_angle = dot / (norm1 * norm2)
    cos_angle = min(1.0, max(-1.0, cos_angle))
    angle = math.acos(cos_angle)
    return angle


decay_modes = {
    1: '1P + Neutrals',
    2: '1P + Neutrals',
    3: '1P + Neutrals',
    5: '3P + Neutrals'
}

# Initialize histograms
hists_dict = {}

for decay_num, decay_mode in decay_modes.items():
    if decay_num == 2 or decay_num == 3: continue
    # true angle between photons
    hPiPhotonAngleTrue = TH1F(f"{decay_mode}_true_pi_photons_angle", f"{decay_mode} True Photon Angle", 50, 0, 0.5)
    hPiPhotonAngleTrue.GetXaxis().SetTitle('#theta [rad]')

    # reco angle between reco photons
    hPhotonAngleReco = TH1F(f"{decay_mode}_reco_photons_min_angle", f"{decay_mode} Reconstructed Minimum Photon Angle", 50, 0, 0.5)
    hPhotonAngleReco.GetXaxis().SetTitle('#theta_{min} [rad]')

    # add to dictionary
    hists_dict[decay_num] = {
        "true": hPiPhotonAngleTrue,
        'reco_gamma' : hPhotonAngleReco
    }

# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')

        taus = [t for t in mcParticles if abs(t.getPDG()) == 15] # get MC taus

        # Loop through MC taus
        for tau in taus:

            # Get tau decay mode
            decayMode = getDecayMode(tau)

            if decayMode in decay_modes:
                if decayMode == 2 or decayMode == 3: decayMode = 1

                #get reco photons, check all possible pairings
                photons = [pfo for pfo in pfos if abs(pfo.getType()) == 22]
                if len(photons) >= 2:
                    angle = 3
                    for i in range(len(photons)):
                        for j in range(i + 1, len(photons)):
                            new_angle = calc_angle([photons[i], photons[j]])
                            angle = min(angle, new_angle) # get min separation angle, likely photons that decay from same pi

                        hists_dict[decayMode]["reco_gamma"].Fill(angle) # add separation angle

                #truth
                daughters = tau.getDaughters()
                pis_0_true = [d for d in daughters if abs(d.getPDG()) == 111]

                if len(pis_0_true) > 0:
                    for pi_true in pis_0_true:
                        # check MC photon
                        daughters = pi_true.getDaughters()
                        gammas = [d for d in daughters if abs(d.getPDG()) == 22]
                        if len(gammas) != 2:
                            continue

                        pt_pi = math.sqrt(pi_true.getMomentum()[0]**2+pi_true.getMomentum()[1]**2)
                        angle = calc_angle(gammas)
                        hists_dict[decayMode]["true"].Fill(angle)

    # Close file
    reader.close()

for decay_num, decay_mode in decay_modes.items():
    if decay_num == 2 or decay_num == 3: continue

    hMC = hists_dict[decay_num]['true']
    hRecoPhoton = hists_dict[decay_num]['reco_gamma']
    hRecoPhoton.GetXaxis().SetRangeUser(0.0, 0.5)

    hMC.SetTitle(f"{decay_mode} Photon Theta")
    hMC.GetXaxis().SetTitle('#theta [rad]')
    hMC.GetXaxis().SetRangeUser(0.0, 0.5)


    hRecoPhoton.SetLineColor(ROOT.kBlue)
    hMC.SetLineColor(ROOT.kRed)

    hRecoPhoton.SetLineWidth(3)
    hMC.SetLineWidth(3)

    hRecoPhoton.SetStats(0)
    hMC.SetStats(0)

    c = ROOT.TCanvas('c', 'overlay', 800, 600)

    hMC.Draw()
    hRecoPhoton.Draw('SAME')

    c.Update()

    legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.90)
    legend.AddEntry(hRecoPhoton, 'Reco Photon #theta_{min}', 'lp')
    legend.AddEntry(hMC, 'MC Photon #theta', 'lp')
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()
    filename = f"{decay_mode}_photon_angles.png"
    c.SaveAs(filename)

    del c

# Write to output file
output_file = TFile(args.outputFile, 'RECREATE')
for decay_num, decay_mode in decay_modes.items():
    if decay_num == 2 or decay_num == 3: continue
    hists_dict[decay_num]['true'].Write()
    hists_dict[decay_num]['reco_gamma'].Write()
output_file.Close()
