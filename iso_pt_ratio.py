from ROOT import TFile
from argparse import ArgumentParser


# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='Taus_loose.root')

args = parser.parse_args()

file = TFile.Open(args.inputFile)

tree = file.Get('anatree')

if not tree:
    print('TTree not found in file!')
    exit()


tau_pt = []
tau_isoE = []
tau_ratios = []
num_reco_taus = 0

for i in range(tree.GetEntries()):
    tree.GetEntry(i)

    for j in range(tree.ntau):
        pt = tree.t_pt[j]
        tau_pt.append(pt)
        isoE = tree.t_isoE[j]
        tau_isoE.append(isoE)
        tau_ratios.append(float(isoE) / float(pt))
        num_reco_taus += 1

ratio_dict = {0.05: 0, 0.1: 0, 0.15: 0, 0.2: 0}

for ratio in tau_ratios:
    for key in ratio_dict:
        if ratio < key:
            ratio_dict[key] += 1

for key in ratio_dict:
    print(f'Number of reco Taus accepted by IsoE/Pt < {key}: {ratio_dict[key]}')

print(f'Total number of reco Taus: {num_reco_taus}')
