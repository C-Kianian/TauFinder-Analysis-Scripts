from argparse import ArgumentParser
from tau_mc_link import getDecayMode, getLinkedMCTau, getNRecoNeutralPis, getNRecoQPis
import os
from pyLCIO import IOIMPL, EVENT, UTIL
import matplotlib.pyplot as plt

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

args = parser.parse_args()

# Check if input file is a directory or a single file
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

    num_true_taus_linked = 0 #mc decay
    num_reco_taus = 0
    num_reco_taus_linked = 0
    num_3_prong_reco_taus = 0
    true_num_3_prong_reco = 0
    num_prongs_per_reco_tau_linked_to_3_prong = []

    for ievt, event in enumerate(reader):

        reco_taus = event.getCollection('TauRec_PFO')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauRecLink_PFO')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        reco_tau_in_event_links_to_mc_tau = False # check if we link more than one tau per event
        reco_tau_in_event_links_to_3_prong = False # check if we link more than one tau per 3-prong event

        # Instantiate relation navigators to parse tauReco and RecoMC links
        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        num_reco_taus += len(reco_taus)

        for reco_tau in reco_taus:
            mcTauLinked = getLinkedMCTau(reco_tau, relationNavigatorTau, relationNavigatorRecoMC)
            decayMode = getDecayMode(mcTauLinked)

            if mcTauLinked is not None:
                num_reco_taus_linked += 1
            if mcTauLinked is not None and reco_tau_in_event_links_to_mc_tau is False:
                reco_tau_in_event_links_to_mc_tau = True
                num_true_taus_linked += 1
            if decayMode == 4 and reco_tau_in_event_links_to_3_prong is False:
                reco_tau_in_event_links_to_3_prong = True
                true_num_3_prong_reco += 1
            if decayMode == 4:
                num_prongs_per_reco_tau_linked_to_3_prong.append(getNRecoQPis(reco_tau))
                num_3_prong_reco_taus += 1


    print(f"Total number of MC taus matched: {num_true_taus_linked}")
    print(f"Total number of reco taus: {num_reco_taus}")
    print(f"Total number of reco taus matched: {num_reco_taus_linked}")
    print(f"Total number of 3-prong reco taus: {num_3_prong_reco_taus}")
    print(f"Total number of 3-prong reco taus matched to MC 3 prong: {true_num_3_prong_reco}")
    print(f"Average number of prongs per reco tau matched to 3-prong: {sum(num_prongs_per_reco_tau_linked_to_3_prong)/len(num_prongs_per_reco_tau_linked_to_3_prong)}")

    plt.hist(num_prongs_per_reco_tau_linked_to_3_prong, bins=[0,1,2,3,4,5,6], color='magenta', alpha=0.7)
    plt.xlabel(r'# of reco prongs per matched reco 3-prong $\tau^-$', fontsize=12, fontweight='bold')
    plt.ylabel(r'# of reco $\tau^-$', fontsize=12, fontweight='bold')
    plt.title(r'# of reco 3-prongs $\tau^-$ vs. # of reco prongs', fontsize=14, fontweight='bold')
    plt.savefig('reco_prongs_per_reco_linked_3_prong.png')

    reader.close()
    del reader

