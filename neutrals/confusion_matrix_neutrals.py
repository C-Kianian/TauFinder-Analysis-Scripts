from argparse import ArgumentParser
from tau_mc_link import getDecayMode, getLinkedMCTau, getNRecoNeutralPis, getNRecoQPis
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import os
from pyLCIO import IOIMPL, EVENT, UTIL
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

uw_red = '#c5050c'
light_base = '#fff0f0'
uw_cmap = LinearSegmentedColormap.from_list("uw_cmap", ['white', uw_red])


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

    y_true = [] #mc decay
    y_pred = [] #reco classified decay

    for ievt, event in enumerate(reader):

        taus = event.getCollection('TauRec_PFO')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauRecLink_PFO')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        decay_dict = {0: '1P0N', 1:'1P + Ns', 2:'1P + Ns', 3:'1P + Ns', 4: '3P0N'} #, 5:'Other', 6:'Other', 7:'Other'}

        mc_taus = [mcParticle for mcParticle in mcParticles if mcParticle.getPDG() == 15]

        for mc_tau in mc_taus:
            # Determine true label
            true_mode = getDecayMode(mc_tau)
            if true_mode == 2 or true_mode == 3: true_mode = 1
            true_label = decay_dict.get(true_mode, "Other")

            # Try to find a reco tau matched to this MC tau
            matched_tau = None
            for tau in taus:
                mcTauLinked = getLinkedMCTau(tau, relationNavigatorTau, relationNavigatorRecoMC)
                if mcTauLinked == mc_tau:
                    matched_tau = tau
                    break

            if matched_tau is not None:
                nRecoNeutralPis = getNRecoNeutralPis(matched_tau) # count reco photons
                nQPis = getNRecoQPis(matched_tau) # count reco charge pis

                if nQPis == 1 and nRecoNeutralPis == 0:
                    pred_label = '1P0N'
                elif nQPis == 1 and nRecoNeutralPis > 0:
                    pred_label = '1P + Ns'
                #elif nQPis == 2 and nRecoNeutralPis == 0:
                #    pred_label = '2P'
                #elif nQPis == 2 and nRecoNeutralPis > 0:
                #    pred_label = '2P + Neutrals'
                elif nQPis == 3:
                    pred_label = '3P0N'
                else:
                    pred_label = 'Other'
            else:
                pred_label = 'Other'

            y_true.append(true_label)
            y_pred.append(pred_label)

    reader.close()

    labels = ["1P0N", "1P + Ns", "3P0N", "Other"]
    conf_matrix = confusion_matrix(y_true=y_true, y_pred=y_pred, labels=labels, normalize='true')
    display = ConfusionMatrixDisplay(confusion_matrix=conf_matrix, display_labels=labels)
    display.plot(cmap=uw_cmap)
    display.ax_.set_xlabel("Reconstructed Decay Mode")
    display.ax_.set_ylabel("True Decay Mode")
    display.ax_.set_xticklabels(["1P0N", "1P + Ns", "3P0N", "Not Matched"])
    plt.title("Tau Decay Mode Classification")
    plt.savefig('confusion_matrix.png')
    print(conf_matrix)
    plt.close()
    del plt
