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

prong_decay_modes = [0, 1, 2, 3, 4, 5]

def getNChargedParticles(recoTau): # counter for electrons, muons, pions
    recoTauParticles = recoTau.getParticles()
    nChargedParticles = 0
    for particle in recoTauParticles:
        if abs(particle.getType()) == 13 or abs(particle.getType()) == 11 or abs(particle.getType()) == 211:
            nChargedParticles += 1
    return nChargedParticles

for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    n_reco_hadronic_events = 0
    total_muons_or_electrons_in_hadronic_events = 0
    n_reco_hadronic_events_with_muon_or_electron = 0
    n_3_charged_particle_events_linked_3p = 0
    n_4_charged_particle_events_including_electron_or_muon = 0

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
            mcTau = getLinkedMCTau(reco_tau, relationNavigatorTau, relationNavigatorRecoMC)
            decayMode = getDecayMode(mcTau)

            if decayMode in prong_decay_modes:
                per_event_total_electrons_or_muons = 0
                n_reco_hadronic_events += 1
                nQPis = getNRecoQPis(reco_tau)
                charged_particles = getNChargedParticles(reco_tau)
                if 0 < nQPis < 4: # only look at events that will have both a reco pion and electron
                    reco_tau_daughters = reco_tau.getParticles()
                    counted = False # only count one muons or electron per tau
                    if charged_particles != nQPis:
                        per_event_total_electrons_or_muons = abs(nQPis - charged_particles)
                        total_muons_or_electrons_in_hadronic_events += abs(nQPis - charged_particles)
                        if not counted:
                            counted = True
                            n_reco_hadronic_events_with_muon_or_electron += 1
                        if charged_particles == 3 and decayMode == 4: n_3_charged_particle_events_linked_3p += 1
                        if charged_particles == 4: n_4_charged_particle_events_including_electron_or_muon += 1

        reader.close()

    print("# of hadronic reco taus (with > 0 pions): " + str(n_reco_hadronic_events))
    print("# of muons or electrons in hadronic reco tau: " + str(total_muons_or_electrons_in_hadronic_events))
    print("# of hadronic reco taus with at least one muon or electron: " + str(n_reco_hadronic_events_with_muon_or_electron))
    print("# of 3 charged particle reco taus linked to 3-prong MC tau: " + str(n_3_charged_particle_events_linked_3p))
    print("# of 4 charged particle reco taus with one electron or muon: " + str(n_4_charged_particle_events_including_electron_or_muon))

