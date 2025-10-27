#taken from https://github.com/ethanmar/MuColl-TauStudy, added a Pi0 counting method
from pyLCIO import IOIMPL, EVENT, UTIL

# Get mc tau associated with reco tau
def getLinkedMCTau(recoTau, relationNavigatorTau, relationNavigatorRecoMC):

    mcTau = None
    
    # Get linked reco tau daughters
    linkedRecoTauDaughters = relationNavigatorTau.getRelatedToObjects(recoTau)

    # Loop over linked reco tau daughters
    for linkedRecoTauDaughter in linkedRecoTauDaughters:

        # Get linked mc tau daughters
        linkedMCTauDaughters = relationNavigatorRecoMC.getRelatedToObjects(linkedRecoTauDaughter)

        # Loop over linked mc tau daughters
        for linkedMCTauDaughter in linkedMCTauDaughters:


            # Check whether linkedMCTauDaughter parent is really a tau

            # Initialize dummy MCParticle
            dummy = linkedMCTauDaughter

            # Initialize parent MCParticle which will be set to tau if linkedMCTauDaughter has tau parent
            parent = linkedMCTauDaughter

            # Get number of parents
            nParents = len(linkedMCTauDaughter.getParents())

            # Loop over parents until we find a tau or we run out of parents (no tau parents)
            while (nParents != 0):
                dummy = parent.getParents()[0]
                nParents = len(dummy.getParents())
                parent = dummy
                if (abs(parent.getPDG()) == 15):
                    nParents = 0
            if (abs(parent.getPDG()) == 15):
                mcTau = parent

    return mcTau

# Get visible properties of MC tau (or any other particle)
def getVisibleProperties(mcTau):

    # Get MC tau daughters
    daughters = mcTau.getDaughters()

    # Instantiate visible energy and momentum to zero
    E_vis = 0
    px_vis = 0
    py_vis = 0
    pz_vis = 0
    n_daughters_vis = 0
    vis_daughter_types = []

    # Loop through daughters
    for daughter in daughters:

        # Get PDG
        pdg = abs(daughter.getPDG())

        # Ignore neutrinos
        if (pdg == 12 or pdg == 14 or pdg == 16):
            continue

        # Sum visible observables
        E_vis += daughter.getEnergy()
        px_vis += daughter.getMomentum()[0]
        py_vis += daughter.getMomentum()[1]
        pz_vis += daughter.getMomentum()[2]
        n_daughters_vis += 1
        vis_daughter_types.append(pdg)
        
    return E_vis, px_vis, py_vis, pz_vis, n_daughters_vis, vis_daughter_types

# Get decat mode of MC tau
def getDecayMode(mcTau):

    '''
    Decay modes for tau+:
    0: -211, 16, (pi-, nu-tau) (10.82%)
    1: -211, 16, 111 (pi-, pi0, nu-tau) (25.49%)
    2: -211, 16, 111, 111 (pi-, pi0x2, nu-tau) (9.26%)
    3: -211, 16, 111, 111, 111 (pi-, pi0x3, nu-tau) (1.04%)
    4: -211, -211, 16, 211 (3-prong plus nu-tau) (8.99%)
    5: -211, -211, 16, 111, 211 (3-prong plus pi0 and nu-tau) (2.74%)
    6: -12, 11, 16 (nu-e-bar, e-, nu-tau) (17.82%)
    7: -14, 13, 16 (nu-mu-bar, mu-, nu-tau) (17.39%) 
    '''
    # Initialize list to store daughter pdgs
    daughter_pdgs = []

    n_daughters = len(mcTau.getDaughters())
    
    # Loop over daughters and store pdgs
    for daughter in mcTau.getDaughters():
        daughter_pdgs.append(daughter.getPDG())

    if (mcTau.getPDG() == 15):

        # Sort daughter pdgs from lowest to highest
        daughter_pdgs.sort()


        # Hadronic decays
        if (daughter_pdgs[0] == -211):
            if (daughter_pdgs[1] == -211):
                # 3-prong
                if (daughter_pdgs[3] == 211):
                    # 3-prong no neutrals
                    return 4
                elif (daughter_pdgs[3] == 111):
                    # 3-prong with neutral
                    return 5
            elif (daughter_pdgs[1] == 16):
                # 1-prong
                if (n_daughters == 2):
                    return 0
                elif (daughter_pdgs[2] == 111):
                    if (n_daughters == 3):
                        return 1
                    elif (n_daughters == 4):
                        return 2
                    elif (n_daughters == 5):
                        return 3
        elif (daughter_pdgs[0] == -12):
            return 6
        elif (daughter_pdgs[0] == -14):
            return 7

# Get number of reco charged pions
def getNRecoQPis(recoTau):

    # Initialize variable to keep track of number of reco charged pions
    nRecoQPis = 0
    
    # Get reco daughters
    recoDaughters = recoTau.getParticles()

    # Loop through reco daughters
    for recoDaughter in recoDaughters:

        # Get type
        type_ = abs(recoDaughter.getType())

        # Check if charged pion
        if (type_ == 211):

            # Increment number of charged pions
            nRecoQPis += 1

    return nRecoQPis

# Get number of reco charged pions (only difference between this and https://github.com/ethanmar/MuColl-TauStudy)
def getNRecoNeutralPis(recoTau):

    # Initialize variable to keep track of number of reco photons or psuedo pi0s
    nRecoPi0s = 0

    # Get reco daughters
    recoDaughters = recoTau.getParticles()

    # Loop through reco daughters
    for recoDaughter in recoDaughters:

        # Get type
        type_ = abs(recoDaughter.getType())

        # Check if photon
        if (type_ == 22):

            # Increment number of photons
            nRecoPi0s += 1

    return nRecoPi0s
