import math
import os
from array import array
from argparse import ArgumentParser
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TFile, TTree

# Args
parser = ArgumentParser(description="Convert LCIO Reco Taus to ROOT Tree for BDT training")
parser.add_argument('--input', type=str, required=True, help='TauFinder .slcio file')
parser.add_argument('--output', type=str, default='tau_bdt_training.root', help='Name output .root file')
parser.add_argument('--isBackground', action='store_true', help='Set if processing background (BIB/Neutrino Gun), sets isSignal=0')
args = parser.parse_args()

# Pion info
def calculate_delta_r(eta1, phi1, eta2, phi2):
    d_phi = phi1 - phi2
    while d_phi > math.pi: d_phi -= 2*math.pi
    while d_phi < -math.pi: d_phi += 2*math.pi
    return math.sqrt((eta1 - eta2)**2 + d_phi**2)

def get_kinematics(p):
    px, py, pz = p[0], p[1], p[2]
    pt = math.sqrt(px**2 + py**2)
    mag = math.sqrt(px**2 + py**2 + pz**2)

    # Theta (polar angle)
    theta = math.acos(pz/mag) if mag > 0 else 0
    # Pseudorapidity (eta)
    eta = -math.log(math.tan(theta / 2.0)) if (0 < theta < math.pi) else 0
    phi = math.atan2(py, px)

    return pt, eta, phi, theta

def compute_dynamic_cone(pt): # Following default TauFinder.cc implementation
    if pt < 10.0:
        return 0.6
    elif 10.0 <= pt <= 120.0:
        return max(6.0 / pt, 0.05)
    else:
        return 0.05

# Root tree
out_file = TFile(args.output, "RECREATE")
tree = TTree("tau_tree", "Tau BDT Training Data")

# Define Buffers
b_pt = array('f', [0.])
b_eta = array('f', [0.])
b_phi = array('f', [0.])
b_theta = array('f', [0.])
b_mass_inv = array('f', [0.])
b_nTracks = array('i', [0])
b_nCharged = array('i', [0])
b_nNeutral = array('i', [0])
b_isoE = array('f', [0.])
b_isSignal = array('i', [0])
b_daughter_types = array('i', [0] * 10)

# Branches
tree.Branch("pt", b_pt, "pt/F")
tree.Branch("theta", b_theta, "theta/F")
tree.Branch("eta", b_eta, "eta/F")
tree.Branch("phi", b_phi, "phi/F")
tree.Branch("invMass", b_mass_inv, "mass/F")
tree.Branch("nTracks", b_nTracks, "nTracks/I")
tree.Branch("nCharged", b_nCharged, "nCharged/I")
tree.Branch("nNeutral", b_nNeutral, "nNeutral/I")
tree.Branch("isoE", b_isoE, "isoE/F")
tree.Branch("isSignal", b_isSignal, "isSignal/I")
tree.Branch("daughterTypes", b_daughter_types, "daughterTypes[10]/I")

# Read slcio file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
try:
    reader.open(args.input)
except Exception as e:
    print(f"Error opening file: {e}")
    exit(1)

print(f"Starting on: {args.input}")

num_true_taus = 0 # For debugging, should be ~18-19k

# Event loop
for i, event in enumerate(reader):
    if i % 1000 == 0: print(f"Processing Event {i}...")

    try:
        reco_taus = event.getCollection("RecoTaus")
        pfos = event.getCollection("PandoraPFOs")
        mc_particles = event.getCollection("MCParticle")
    except Exception:
        print(f"Missing collection in event: {i}")
        continue

    # MC taus for truth matching
    true_taus = []
    if not args.isBackground:
        for mcp in mc_particles:
            if abs(mcp.getPDG()) == 15: #and mcp.getGeneratorStatus() == 1:
                _, m_eta, m_phi, _ = get_kinematics(mcp.getMomentum())
                num_true_taus += 1
                true_taus.append((m_eta, m_phi))

    for tau in reco_taus:
        t_pt, t_eta, t_phi, t_theta = get_kinematics(tau.getMomentum())

        # Get daughter's tracks, inv mass, types (neutral/charged/pdg)
        daughters = tau.getParticles()
        n_tracks = 0
        n_charged = 0
        n_neutral = 0
        sum_e, sum_px, sum_py, sum_pz = 0.0, 0.0, 0.0, 0.0
        daughter_types = []

        for d in daughters:
            # 4 vector sum
            sum_e += d.getEnergy()
            dp = d.getMomentum()
            sum_px += dp[0]; sum_py += dp[1]; sum_pz += dp[2]

            # Type/PDG counting
            daughter_types.append(d.getType())
            if d.getCharge() != 0:
                n_charged += 1
                if d.getTracks().size() > 0:
                    n_tracks += 1
            else:
                n_neutral += 1

        # Calculate inv mass
        p2 = sum_px**2 + sum_py**2 + sum_pz**2
        m2 = sum_e**2 - p2
        calc_mass = math.sqrt(m2) if m2 > 0 else 0.0

        # Dynamic isolation cone
        inner_cone = compute_dynamic_cone(t_pt)
        outer_cone = inner_cone + 0.2
        iso_e = 0.0
        for pfo in pfos:
            p_mom = pfo.getMomentum()
            p_pt = math.sqrt(p_mom[0]**2 + p_mom[1]**2)
            if p_pt < 0.2: continue
            _, p_eta, p_phi, _ = get_kinematics(p_mom)
            dr = calculate_delta_r(t_eta, t_phi, p_eta, p_phi)
            if inner_cone < dr < outer_cone:
                iso_e += pfo.getEnergy()

        # Fill buffers
        b_pt[0], b_theta[0], b_eta[0], b_phi[0] = t_pt, t_theta, t_eta, t_phi
        b_mass_inv[0] = calc_mass
        b_nTracks[0] = n_tracks
        b_nCharged[0] = n_charged
        b_nNeutral[0] = n_neutral
        b_isoE[0] = iso_e

        # Clear buffer
        for j in range(10): 
            b_daughter_types[j] = 0

        # Add daughter types
        for j, dtype in enumerate(daughter_types):
            if j < 10:
                b_daughter_types[j] = int(dtype)
            else:
                break # TauFinder limits to 10

        # Truth matching
        matched = False
        for m_eta, m_phi in true_taus:
            if calculate_delta_r(t_eta, t_phi, m_eta, m_phi) < 0.1:
                matched = True; break
        b_isSignal[0] = 1 if (matched and not args.isBackground) else 0

        tree.Fill()

reader.close()
out_file.Write()
out_file.Close()
print(f"Processed {num_true_taus} true taus")
print(f"Tree saved successfully to {args.output}")

