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

for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    n_reco_4p_events = 0
    n_reco_4p_matched_3p = 0
    total_muons_or_electrons = 0
    n_reco_4p_with_electron_or_muon = 0
    non_3p_reco_4p_events = []

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

            if getNRecoQPis(reco_tau) == 4:
                n_reco_4p_events += 1
                if decayMode == 4 or decayMode == 5: n_reco_4p_matched_3p += 1
                else: non_3p_reco_4p_events.append(decayMode)

    reader.close()


    print("# of 4p reco taus: " + str(n_reco_4p_events))
    print("# of 4p reco taus matched to 3p MC tau: " + str(n_reco_4p_matched_3p))
    print("Non-3p decay modes: " + str(non_3p_reco_4p_events))



