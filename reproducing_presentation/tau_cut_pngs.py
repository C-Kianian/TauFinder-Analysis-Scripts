from ROOT import TFile
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

uw_red = '#c5050c'
light_base = '#ffe6e6'
uw_cmap = LinearSegmentedColormap.from_list("uw_cmap", [light_base, uw_red])

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='Taus_default.root')

args = parser.parse_args()

file = TFile.Open(args.inputFile)

tree = file.Get('anatree')

if not tree:
    print('TTree not found in file!')
    exit()

tau_pt = []
tau_isoE = []
tau_invM = []
isoE_100 = []
event_num = []
event_num_isoE_100 = []
nrej_isoE_tot = 0
n_isoE_more_5 = 0
n_invM_more_2 = 0


for i in range(tree.GetEntries()):
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

for invM in tau_invM:
    if invM > 2:
        n_invM_more_2 += 1

print(f'Number of taus that fail isoE: {nrej_isoE_tot}')
print(f'Fraction of taus with isoE > 5 GeV: {n_isoE_more_5/len(tau_isoE)}')
print(f'Fraction of taus with invM > 2 GeV/c^2: {n_invM_more_2/len(tau_invM)}')

for i in range(len(event_num)):
    if tau_isoE[i] > 100:
        event_num_isoE_100.append(event_num[i])
        isoE_100.append(tau_isoE[i])

print(f'Events with isoE > 100 GeV:')
for i in range(len(event_num_isoE_100)):
    print(f'Event number: {event_num_isoE_100[i]}     Isolation Energy: {isoE_100[i]}')


#isoE
plt.hexbin(
    tau_pt,
    tau_isoE,
    gridsize=50,
    extent=[0, 300, 0, 10],  # Match your axis limits
    mincnt=1,  # Ignore empty bins
    cmap=uw_cmap
)
plt.colorbar(label=r'# of reco $\tau^-$')

#plt.scatter(tau_pt, tau_isoE, color='magenta', s=10, alpha=0.2)

#plt.scatter(tau_pt, tau_isoE, color='magenta', s=10)

#plt.hist2d(tau_pt, tau_isoE, bins=100, cmap=uw_cmap)
#plt.colorbar(label='Count')

#default cut
plt.axhline(y=5, linestyle='--', color='black', label='5 GeV')

#scale
plt.xlim(0,300)
plt.ylim(0,10)

#axis labels
plt.xlabel(r'$p_T$ [GeV/c]', fontsize=12, fontweight='bold')
plt.ylabel(r'$E_{iso}$ [GeV]', fontsize=12, fontweight='bold')
plt.title(r'$\tau^-$ Isolation Energy vs Transverse Momentum', fontsize=15, fontweight='bold')
plt.legend()

#bold
plt.tick_params(axis='both', labelsize=10, width=1.5)
for label in plt.gca().get_xticklabels() + plt.gca().get_yticklabels():
    label.set_fontweight('bold')
plt.savefig('hex_pt_iso_E.png')


#plt.show()
plt.clf()

#invM
plt.hexbin(
    tau_pt,
    tau_invM,
    gridsize=50,
    extent=[0, 300, 0, 3],  # Match your axis limits
    mincnt=1,  # Ignore empty bins
    cmap=uw_cmap
)
plt.colorbar(label=r'# of reco $\tau^-$')

#plt.scatter(tau_pt, tau_invM, color='magenta', s=10, alpha=0.2)

#plt.scatter(tau_pt, tau_invM, color='magenta', s=10)

##Convert to numpy arrays if not already
#tau_pt = np.array(tau_pt)
#tau_invM = np.array(tau_invM)
#
##Filter out NaN or inf values
#mask = np.isfinite(tau_pt) & np.isfinite(tau_invM)
#
#plt.hist2d(tau_pt[mask], tau_invM[mask], bins=100, cmap=uw_cmap)
#plt.colorbar(label='Count')

#default cut
plt.axhline(y=2, linestyle='--', color='black', label='2 $GeV/c^2$')

#axis titles
plt.xlabel(r'$p_T$ [GeV/c]', fontsize=12, fontweight='bold')
plt.ylabel(r'$M_{inv}$ [$GeV/c^2$]', fontsize=12, fontweight='bold')
plt.title(r'$\tau^-$ Invariant Mass vs Transverse Momentum', fontsize=15, fontweight='bold')

#scale
plt.xlim(0,300)
plt.ylim(0,3)

plt.legend()

#bold
plt.tick_params(axis='both', labelsize=10, width=1.5)
for label in plt.gca().get_xticklabels() + plt.gca().get_yticklabels():
    label.set_fontweight('bold')
plt.savefig('hex_pt_inv_M.png')

#plt.show()


file.Close()
