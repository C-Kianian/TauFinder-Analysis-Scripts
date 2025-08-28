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
tau_invM = []
event_num = []
pt_by_5 = [0] * 60
nrej_isoE_tot = 0
n_isoE_more_5 = 0
n_invM_more_2 = 0
n_Pt_less_5 = 0
n_Pt_more_5_less_10 = 0


for i in range(tree.GetEntries()): # get info for the tau candidate parameters
    tree.GetEntry(i)
    nrej_isoE_tot += tree.nrej_isoE

    for j in range(tree.ntau):
        tau_pt.append(tree.t_pt[j])
        tau_isoE.append(tree.t_isoE[j])
        tau_invM.append(tree.t_minv[j])
        event_num.append(tree.event_num[j])

for isoE in tau_isoE:
    if isoE > 5:
        n_isoE_more_5 += 1

n_candidates = len(event_num)

print("")
print("--------------------------Isolation Energy------------------------------")
print("Number of Tau candidates rejected with isolation energy greater than 5: " + str(n_isoE_more_5))
print("Percent of Tau candidates rejected with isolation energy greater than 5: " + str(n_isoE_more_5/n_candidates * 100))
print("")


for invM in tau_invM:
    if invM > 2:
        n_invM_more_2 += 1

print("")
print("--------------------------Invariant Mass--------------------------------")
print("Number of Tau candidates rejected with invariant mass greater than 2: " + str(n_invM_more_2))
print("Percent of Tau candidates rejected with invariant mass greater than 2: " + str(round(n_invM_more_2/n_candidates * 100 , 2)) + "%")
print("")

for pt in tau_pt:
    if pt <= 5:
        n_Pt_less_5 += 1
    for i in range(60):
        if 5 * i < pt <= 5 + 5 * i: # for each bin of 5 pt find num of tau candidates
            pt_by_5[i] += 1

print("---------------------------Transverse Momentum--------------------------")
print("Number of Tau candidates with transverse momentum less than 5: " + str(n_Pt_less_5))
print("Percent of Tau candidates with transverse momentum less than 5: " + str(round(n_Pt_less_5/n_candidates * 100, 2)) + "%")
print("")

print("-------------------------Transverse Momentum By 5------------------------")

for i in range(60):
    print("Number of Tau candidates with transverse momentum more than " + str(5 * i) + " less than " + str(5 * i + 5) + ": " + str(pt_by_5[i]))
    print("Percent of Tau candidates with transverse momentum more than " + str(5 * i) + " less than " + str(5 * i + 5) + ": " + str(round(pt_by_5[i]/n_candidates * 100, 2))+ "%")
    print("")


