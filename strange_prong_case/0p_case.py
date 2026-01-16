from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas, TLegend
import os
from argparse import ArgumentParser
import matplotlib.pyplot as plt

from tau_mc_link import getLinkedMCTau, getDecayMode, getNRecoQPis

# Command line arguments
parser = ArgumentParser()
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')
#parser.add_argument('--outputFile', type=str, default='reco_tau_prong_combinations.root')
args = parser.parse_args()

# Discover input files
to_process = []
if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)


def getNChargedParticles(recoTau): # counter for electrons, muons, pions
    recoTauParticles = recoTau.getParticles()
    nChargedParticles = 0
    for particle in recoTauParticles:
        if abs(particle.getType()) == 13 or abs(particle.getType()) == 11 or abs(particle.getType()) == 211:
            nChargedParticles += 1
    return nChargedParticles

prong_decay_modes = [0, 1, 2, 3, 4, 5]
total_muons = 0
total_electrons = 0
total_events_with_at_least_one_electron_or_muon = 0
total_0p_events = 0
total_0p_linked_1p = 0
n_1_charged_particles_in_linked_1p = 0
total_0p_linked_3p = 0
n_3_charged_particles_in_linked_3p = 0


for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:
        # Get collections
        reco_taus = event.getCollection('RecoTaus')
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauPFOLink')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        for reco_tau in reco_taus:

            # Get linked MC tau
            mcTau = getLinkedMCTau(reco_tau, relationNavigatorTau, relationNavigatorRecoMC)

            # Get linked MC tau decay mode
            decayMode = getDecayMode(mcTau)

            if decayMode in prong_decay_modes: # only study the hadronic decays
                if getNRecoQPis(reco_tau) == 0: # study only 0p taus
                    total_0p_events += 1
                    if decayMode == 4 or decayMode == 5:
                        total_0p_linked_3p += 1
                        if getNChargedParticles(reco_tau) == 3: n_3_charged_particles_in_linked_3p += 1
                    else:
                        total_0p_linked_1p += 1
                        if getNChargedParticles(reco_tau) == 1: n_1_charged_particles_in_linked_1p += 1
                    counted = False # only count one muons or electron per tau
                    reco_tau_daughters = reco_tau.getParticles()
                    for daughter in reco_tau_daughters:
                        if abs(daughter.getType()) == 13:
                            total_muons += 1
                            if not counted:
                                counted = True
                                total_events_with_at_least_one_electron_or_muon += 1
                        elif abs(daughter.getType()) == 11:
                            total_electrons += 1
                            if not counted:
                                counted = True
                                total_events_with_at_least_one_electron_or_muon += 1



    reader.close()

print("# of 0p taus: " + str(total_0p_events))
print("# of muons in 0p taus: " + str(total_muons))
print("# of electrons in 0p taus: " + str(total_electrons))
print("# of 0p taus with one electron or muon: " + str(total_events_with_at_least_one_electron_or_muon))
print("# of 0p taus linked to 1p MC tau: " + str(total_0p_linked_1p))
print("# of 1 charged particle reco taus linked to 1p tau: " + str(n_1_charged_particles_in_linked_1p))
print("# of 0p taus linked to 3p MC tau: " + str(total_0p_linked_3p))
print("# of 3 charged particle reco taus linked to 3p tau: " + str(n_3_charged_particles_in_linked_3p))

