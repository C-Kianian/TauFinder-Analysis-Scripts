from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np

from tau_mc_link import getLinkedMCTau, getVisibleProperties, getDecayMode, getNRecoNeutralPis, getNRecoQPis

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='tau_neutral_ana.root')

args = parser.parse_args()

# Initialize histograms
general_hists = []
hists_dict = {1: {}, 5:{}}

# Decay modes
decay_modes = {1: '1P + N',
               5: '3P + N'}
modes = [1,2,3,5]

for decay_num, decay_mode in decay_modes.items(): # decay specific histograms
    name = decay_mode.replace(' ','_').lower()
    #linked mc tau properties
    hLinkedMCTauPt = TH1F(f'linked_{name}_mc_tau_pt', f'Linked {decay_mode} MC Tau Pt', 20, 0, 320)
    hists_dict[decay_num]['mc_tau_pt'] = hLinkedMCTauPt

    hLinkedMCTauTheta = TH1F(f'linked_{name}_mc_tau_theta', f'Linked {decay_mode} MC Tau Theta', 20, 0, math.pi)
    hists_dict[decay_num]['mc_tau_theta'] = hLinkedMCTauTheta

    hLinkedMCTauPhi = TH1F(f'linked_{name}_mc_tau_phi', f'Linked {decay_mode} MC Tau Phi', 20, 0, math.pi)
    hists_dict[decay_num]['mc_tau_phi'] = hLinkedMCTauPhi

    #linked reco tau properties
    hLinkedMCTauPtReco = TH1F(f'linked_{name}_reco_mc_tau_pt', f'Linked Reco-{decay_mode} MC Tau Pt', 20, 0, 320)
    hists_dict[decay_num]['linked_correct_reco_mc_tau_pt'] = hLinkedMCTauPtReco

    hLinkedMCTauThetaReco = TH1F(f'linked_{name}_reco_mc_tau_theta', f'Linked Reco-{decay_mode} MC Tau Theta', 20, 0, math.pi)
    hists_dict[decay_num]['linked_correct_reco_mc_tau_theta'] = hLinkedMCTauThetaReco

    hLinkedMCTauPhiReco = TH1F(f'linked_{name}_reco_mc_tau_phi', f'Linked Reco-{decay_mode} MC Tau Phi', 20, 0, math.pi)
    hists_dict[decay_num]['linked_correct_reco_mc_tau_phi'] = hLinkedMCTauPhiReco

    # reconstructed tau properties
    hTauEReco = TH1F(f'{name}_tau_energy', f'Reconstructed {decay_mode} Tau Energy', 500, 0, 1500)
    hists_dict[decay_num]['reco_tau_energy'] = hTauEReco

    hTauPtReco = TH1F(f'{name}_tau_pt', f'Reconstructed {decay_mode} Tau Pt', 20, 0, 320)
    hists_dict[decay_num]['reco_tau_pt'] = hTauPtReco

    hTauThetaReco = TH1F(f'{name}_tau_theta', f'Reconstructed {decay_mode} Tau Theta', 20, 0, math.pi)
    hists_dict[decay_num]['reco_tau_theta'] = hTauThetaReco

    hTauPhiReco = TH1F(f'{name}_tau_phi', f'Reconstructed {decay_mode} Tau Phi', 20, 0, math.pi)
    hists_dict[decay_num]['reco_tau_phi'] = hTauPhiReco

    #reco daughters
    hTauNDaughtersReco = TH1F(f'n_{name}_tau_daughters', f'Number of Reconstructed {decay_mode} Tau Daughters', 10, 0, 10)
    hists_dict[decay_num]['n_linked_correct_reco_tau_daughters'] = hTauNDaughtersReco

    hTauDaughterTypeReco = TH1F(f'{name}_tau_daughter_types', f'Reconstructed {decay_mode} Tau Daughter Types', 2200, 0, 2200)
    hists_dict[decay_num]['reco_tau_daughter_types'] = hTauDaughterTypeReco

    #mc daughters
    hTauNDaughtersMC = TH1F(f'n_{name}_true_tau_daughters', f'Number of Reconstructed {decay_mode} Tau Daughters', 10, 0, 10)
    hists_dict[decay_num]['n_mc_tau_daughters'] = hTauNDaughtersMC

    hTauDaughterTypeMC = TH1F(f'{name}_true_tau_daughter_types', f'Reconstructed {decay_mode} Tau Daughter Types', 2200, 0, 2200)
    hists_dict[decay_num]['mc_tau_daughter_types'] = hTauDaughterTypeMC

    #mc true visible and total properties
    hTauVisPtTrue = TH1F(f'{name}_true_vis_pT', f'True {decay_mode} Tau Visible Pt', 20, 0, 320)
    hists_dict[decay_num]['true_vis_pT'] = hTauVisPtTrue

    hTauVisThetaTrue = TH1F(f'{name}_true_vis_theta', f'True {decay_mode} Tau Visible Theta', 20, 0, math.pi)
    hists_dict[decay_num]['true_vis_theta'] = hTauVisThetaTrue

    hTauVisPhiTrue = TH1F(f'{name}_true_vis_phi', f'True {decay_mode} Tau Visible Phi', 20, 0, math.pi)
    hists_dict[decay_num]['true_vis_phi'] = hTauVisPhiTrue

    hTauTotPtTrue = TH1F(f'{name}_true_tot_pT', f'True {decay_mode} Tau Total Pt', 20, 0, 320)
    hists_dict[decay_num]['true_tot_pT'] = hTauTotPtTrue

    hTauTotThetaTrue = TH1F(f'{name}_true_tot_theta', f'True {decay_mode} Tau Total Theta', 20, 0, math.pi)
    hists_dict[decay_num]['true_tot_theta'] = hTauTotThetaTrue

    hTauTotPhiTrue = TH1F(f'{name}_true_tot_phi', f'True {decay_mode} Tau Total Phi', 20, 0, math.pi)
    hists_dict[decay_num]['true_tot_phi'] = hTauTotPhiTrue

    # neutral pions
    hPi0PtTrue = TH1F(f'{name}_pi_0_true_pT', f'True {decay_mode} Neutral Pion Pt', 20, 0, 320)
    hists_dict[decay_num]['pi_0_true_pT'] = hPi0PtTrue

    hPi0PhiTrue = TH1F(f'{name}_pi_0_true_phi', f'True {decay_mode} Neutral Pion Phi', 20, 0, math.pi)
    hists_dict[decay_num]['pi_0_true_phi'] = hPi0PhiTrue

    hPi0ThetaTrue = TH1F(f'{name}_pi_0_true_theta', f'True {decay_mode} Neutral Pion Theta', 20, 0, math.pi)
    hists_dict[decay_num]['pi_0_true_theta'] = hPi0ThetaTrue

    hPi0PtMatched = TH1F(f'{name}_pi_0_matched_pT', f'Matched {decay_mode} Neutral Pion Pt', 20, 0, 320)
    hists_dict[decay_num]['pi_0_matched_pT'] = hPi0PtMatched

    hPi0PhiMatched = TH1F(f'{name}_pi_0_matched_phi', f'Matched {decay_mode} Neutral Pion Phi', 20, 0, math.pi)
    hists_dict[decay_num]['pi_0_matched_phi'] = hPi0PhiMatched

    hPi0ThetaMatched = TH1F(f'{name}_pi_0_matched_theta', f'Matched {decay_mode} Neutral Pion Theta', 20, 0, math.pi)
    hists_dict[decay_num]['pi_0_matched_theta'] = hPi0ThetaMatched

    # charged pions
    hPiPtTrue = TH1F(f'{name}_pi_true_pT', f'True {decay_mode} Pion Pt', 20, 0, 320)
    hists_dict[decay_num]['pi_true_pT'] = hPiPtTrue

    hPiPhiTrue = TH1F(f'{name}_pi_true_phi', f'True {decay_mode} Pion Phi', 20, 0, math.pi)
    hists_dict[decay_num]['pi_true_phi'] = hPiPhiTrue

    hPiThetaTrue = TH1F(f'{name}_pi_true_theta', f'True {decay_mode} Pion Theta', 20, 0, math.pi)
    hists_dict[decay_num]['pi_true_theta'] = hPiThetaTrue

    hPiPtMatched = TH1F(f'{name}_pi_matched_pT', f'Matched {decay_mode} Pion Pt', 20, 0, 320)
    hists_dict[decay_num]['pi_matched_pT'] = hPiPtMatched

    hPiPhiMatched = TH1F(f'{name}_pi_matched_phi', f'Matched {decay_mode} Pion Phi', 20, 0, math.pi)
    hists_dict[decay_num]['pi_matched_phi'] = hPiPhiMatched

    hPiThetaMatched = TH1F(f'{name}_pi_matched_theta', f'Matched {decay_mode} Pion Theta', 20, 0, math.pi)
    hists_dict[decay_num]['pi_matched_theta'] = hPiThetaMatched


# general Reco
hTauE = TH1F('tau_energy', 'Reconstructed Tau Energy', 500, 0, 1500)
general_hists.append(hTauE)

hTauPt = TH1F('tau_pt', 'Reconstructed Tau Pt', 20, 0, 320)
general_hists.append(hTauPt)

hTauTheta = TH1F('tau_theta', 'Reconstructed Tau Theta', 20, 0, math.pi)
general_hists.append(hTauTheta)

hTauPhi = TH1F('tau_phi', 'Reconstructed Tau Phi', 20, 0, math.pi)
general_hists.append(hTauPhi)

hTauNDaughters = TH1F('n_tau_daughters', 'Number of Reconstructed Tau Daughters', 10, 0, 10)
general_hists.append(hTauNDaughters)

hTauDaughterType = TH1F('tau_daughter_types', 'Reconstructed Tau Daughter Types', 1100, 0, 2200)
general_hists.append(hTauDaughterType)

hLinkedMCTauPt = TH1F('linked_mc_tau_pt', 'Linked MC Tau Pt', 20, 0, 320)
general_hists.append(hLinkedMCTauPt)

hLinkedMCTauTheta = TH1F('linked_mc_tau_theta', 'Linked MC Tau Theta', 20, 0, math.pi)
general_hists.append(hLinkedMCTauTheta)

hLinkedMCTauPhi = TH1F('linked_mc_tau_phi', 'Linked MC Tau Phi', 20, 0, math.pi)
general_hists.append(hLinkedMCTauPhi)

# general mc propterties
hTauVisETrue = TH1F('tau_true_energy', 'True Tau Visible Energy', 500, 0, 1500)
general_hists.append(hTauVisETrue)

hTauVisPtTrue = TH1F('tau_true_pT', 'True Tau Visible Pt', 20, 0, 320)
general_hists.append(hTauVisPtTrue)

hTauVisThetaTrue = TH1F('tau_true_theta', 'True Tau Visible Theta', 20, 0, math.pi)
general_hists.append(hTauVisThetaTrue)

hTauVisPhiTrue = TH1F('tau_phi_true', 'True Tau Visible Phi', 20, 0, math.pi)
general_hists.append(hTauVisPhiTrue)

hTauNVisDaughtersTrue = TH1F('n_mc_visible_tau_daughters', 'Number of True Visible Tau Daughters', 10, 0, 10)
general_hists.append(hTauNVisDaughtersTrue)

hTauVisDaughterTypeTrue = TH1F('visible_mc_tau_daughter_types', 'True Visible Tau Daughter Types', 300, 0, 300)
general_hists.append(hTauVisDaughterTypeTrue)


# Check if input file is a directory or a single file
to_process = []

if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)


# Keep track of one-pion + neutrals reco daughter types
pion_types = {1: {}, 5: {}}

# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        reco_taus = event.getCollection('TauRec_PFO')
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauRecLink_PFO')
        recoMCLink = event.getCollection('RecoMCTruthLink')
        reco_photons = [pfo for pfo in pfos if abs(pfo.getType()) == 22]
        reco_pis = [pfo for pfo in pfos if abs(pfo.getType()) == 211]

        # Instantiate relation navigators to parse tauReco and RecoMC links
        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        # Loop through tau PFOs
        for reco_tau in reco_taus:

            # Instantiate and calculate tau observables
            E = reco_tau.getEnergy()
            p = reco_tau.getMomentum()
            px, py, pz = p[0], p[1], p[2]
            pt = math.sqrt(px**2 + py**2)
            theta = math.acos(p[2]/(math.sqrt(pt**2+p[2]**2)))
            phi = math.acos(p[0]/pt)
            reco_tau_daughters = reco_tau.getParticles()

            # Fill general reco tau hists
            hTauE.Fill(E)
            hTauPt.Fill(pt)
            hTauTheta.Fill(theta)
            hTauPhi.Fill(phi)
            hTauNDaughters.Fill(len(reco_tau_daughters))

            # Loop over reco daughters and fill reco daughter hists
            for daughter in reco_tau_daughters:
                type_ = abs(daughter.getType())
                hTauDaughterType.Fill(type_)

            # Get linked MC tau
            mcTau = getLinkedMCTau(reco_tau, relationNavigatorTau, relationNavigatorRecoMC)

            # Get visible properties of linked MC tau
            E_vis, px_vis, py_vis, pz_vis, n_daughters_vis, vis_daughter_types = getVisibleProperties(mcTau)
            pt_vis = math.sqrt(px_vis**2 + py_vis**2)
            theta_vis = math.acos(pz_vis/(math.sqrt(pt_vis**2 + pz_vis**2)))
            phi_vis = math.acos(px_vis/pt_vis)

            # Fill linked MC tau hists
            hLinkedMCTauPt.Fill(pt_vis)
            hLinkedMCTauTheta.Fill(theta_vis)
            hLinkedMCTauPhi.Fill(phi_vis)

            # Get linked MC tau decay mode
            decayMode = getDecayMode(mcTau)

            # Fill neutral hists
            if decayMode in modes:
                #fill mc info for all reco taus based on decay mode, general reco info
                if decayMode == 2 or decayMode == 3: decayMode = 1 # combine these modes
                hists_dict[decayMode]['mc_tau_pt'].Fill(pt_vis)
                hists_dict[decayMode]['mc_tau_theta'].Fill(theta_vis)
                hists_dict[decayMode]['mc_tau_phi'].Fill(phi_vis)
                hists_dict[decayMode]['n_mc_tau_daughters'].Fill(n_daughters_vis)
                for vis_daughter_type in vis_daughter_types:
                    hTauNVisDaughtersTrue.Fill(vis_daughter_type)

                # Get number of reco neutral pions in reco tau daughters
                nRecoNeutralPis = getNRecoNeutralPis(reco_tau)
                nQRecopis = getNRecoQPis(reco_tau)

                # correctly reconstructed tau, ie. links to proper type and satisfies the criteria for that type
                if nRecoNeutralPis > 0 and ((decayMode == 1 and nQRecopis == 1) or (decayMode == 5 and nQRecopis == 3)):
                    # Fill truth info hists for linked reco-neutral that pass criteria
                    hists_dict[decayMode]['linked_correct_reco_mc_tau_pt'].Fill(pt_vis)
                    hists_dict[decayMode]['linked_correct_reco_mc_tau_theta'].Fill(theta_vis)
                    hists_dict[decayMode]['linked_correct_reco_mc_tau_phi'].Fill(phi_vis)
                    hists_dict[decayMode]['n_linked_correct_reco_tau_daughters'].Fill(nRecoNeutralPis)

                    # Fill reco info hists with reco neutral info
                    hists_dict[decayMode]['reco_tau_energy'].Fill(E)
                    hists_dict[decayMode]['reco_tau_pt'].Fill(pt)
                    hists_dict[decayMode]['reco_tau_theta'].Fill(theta)
                    hists_dict[decayMode]['reco_tau_phi'].Fill(phi)
                    hTauNDaughters.Fill(len(reco_tau_daughters))
                    for daughter in reco_tau_daughters: # get daughters
                        type_ = abs(daughter.getType())
                        hists_dict[decayMode]['reco_tau_daughter_types'].Fill(type_)
                        if str(type_) in pion_types[decayMode]: pion_types[decayMode][str(type_)] += 1
                        else: pion_types[decayMode][str(type_)] = 1

        # Loop through MC particles
        for mcParticle in mcParticles:

            # Get MC taus
            pdg = abs(mcParticle.getPDG())
            if pdg == 15:

                # Get visible properties
                E_vis, px_vis, py_vis, pz_vis, n_daughters_vis, vis_daughter_types = getVisibleProperties(mcParticle)
                pt_vis = math.sqrt(px_vis**2 + py_vis**2)
                theta_vis = math.acos(pz_vis/(math.sqrt(pt_vis**2 + pz_vis**2)))
                phi_vis = math.acos(px_vis/pt_vis)

                # Fill visible mc tau hists
                hTauVisETrue.Fill(E_vis)
                hTauVisPtTrue.Fill(pt_vis)
                hTauVisThetaTrue.Fill(theta_vis)
                hTauVisPhiTrue.Fill(phi_vis)
                hTauNVisDaughtersTrue.Fill(n_daughters_vis)
                for vis_daughter_type in vis_daughter_types:
                    hTauVisDaughterTypeTrue.Fill(vis_daughter_type)

                # Get total properties
                tot_mom = mcParticle.getMomentum()
                tot_px, tot_py, tot_pz  = tot_mom[0], tot_mom[1], tot_mom[2]
                tot_pt = math.sqrt(tot_px**2 + tot_py**2)
                tot_theta = math.acos(tot_pz/(math.sqrt(tot_pt**2 + tot_pz**2)))
                tot_phi = math.acos(tot_px/tot_pt)

                # Get mc particle decay mode
                decayMode = getDecayMode(mcParticle)

                if decayMode in modes:
                    if decayMode == 2 or decayMode == 3: decayMode = 1
                    # Fill mc visible tau hists
                    hists_dict[decayMode]['true_vis_pT'].Fill(pt_vis)
                    hists_dict[decayMode]['true_vis_theta'].Fill(theta_vis)
                    hists_dict[decayMode]['true_vis_phi'].Fill(phi_vis)
                    hists_dict[decayMode]['true_tot_pT'].Fill(tot_pt)
                    hists_dict[decayMode]['true_tot_theta'].Fill(tot_theta)
                    hists_dict[decayMode]['true_tot_phi'].Fill(tot_phi)

                    mc_pis_0 = [daughter for daughter in mcParticle.getDaughters() if abs(daughter.getPDG()) == 111]
                    mc_pis = [daughter for daughter in mcParticle.getDaughters() if abs(daughter.getPDG()) == 211]

                    for mc_pi_0 in mc_pis_0: # neutral pi info
                        # Fill true histogram using mc tau info
                        hists_dict[decayMode]['pi_0_true_pT'].Fill(pt_vis)
                        hists_dict[decayMode]['pi_0_true_phi'].Fill(phi_vis)
                        hists_dict[decayMode]['pi_0_true_theta'].Fill(theta_vis)
                    for mc_pi in mc_pis: #charged pi
                        hists_dict[decayMode]['pi_true_pT'].Fill(pt_vis)
                        hists_dict[decayMode]['pi_true_phi'].Fill(phi_vis)
                        hists_dict[decayMode]['pi_true_theta'].Fill(theta_vis)

                    linked_mc_pi_0 = []
                    linked_mc_pis = []
                    for reco_photon in reco_photons:
                        linkedMCParticles = relationNavigatorRecoMC.getRelatedToObjects(reco_photon)
                        for linkedMCParticle in linkedMCParticles:
                            if abs(linkedMCParticle.getPDG()) == 111 and linkedMCParticle not in linked_mc_pi_0:
                                # if reco photon matches a mc pion and mc pion has not been matched before
                                linked_mc_pi_0.append(linkedMCParticle)
                                hists_dict[decayMode]['pi_0_matched_pT'].Fill(pt_vis)
                                hists_dict[decayMode]['pi_0_matched_phi'].Fill(phi_vis)
                                hists_dict[decayMode]['pi_0_matched_theta'].Fill(theta_vis)

                    for reco_pi in reco_pis: # get reco charged pi info
                        linkedMCParticles = relationNavigatorRecoMC.getRelatedToObjects(reco_pi)
                        for linkedMCParticle in linkedMCParticles:
                            if abs(linkedMCParticle.getPDG()) == 211 and linkedMCParticle not in linked_mc_pis:
                                linked_mc_pis.append(linkedMCParticle)

                    if (len(linked_mc_pis) == 1 and decayMode == 1) or (len(linked_mc_pis) == 3 and decayMode == 5):
                        hists_dict[decayMode]['pi_matched_pT'].Fill(pt_vis)
                        hists_dict[decayMode]['pi_matched_phi'].Fill(phi_vis)
                        hists_dict[decayMode]['pi_matched_theta'].Fill(theta_vis)

    # Close file
    reader.close()

# function to style a hist based on specifications and add to hist list
def style_hist(hist, color, width, title, x_label, y_label):
    hist.SetLineColor(color)
    hist.SetLineWidth(width)
    hist.SetTitle(title)
    if "Efficiency" in title: hist.GetYaxis().SetRangeUser(0.0, 1.2)
    hist.GetXaxis().SetTitle(x_label)
    hist.GetYaxis().SetTitle(y_label)
    hist.SetStats(0)
    general_hists.append(hist)

for decay_num, decay_mode in decay_modes.items():
    name = decay_mode.replace(' ','_').lower()
    #general reco efficiency hists
    hPtEff = hists_dict[decay_num]["mc_tau_pt"].Clone(f'{name}_pt_eff')
    hPtEff.Divide(hPtEff, hists_dict[decay_num]['true_vis_pT'], 1, 1, 'B')
    style_hist(hPtEff, 6, 2, f'{decay_mode} Reco Efficiency vs Pt', 'True Visible P_{t} [GeV/c]', '#epsilon')

    hThetaEff = hists_dict[decay_num]["mc_tau_theta"].Clone(f'{name}_theta_eff')
    hThetaEff.Divide(hThetaEff, hists_dict[decay_num]['true_vis_theta'], 1, 1, 'B')
    style_hist(hThetaEff, 7, 2, f'{decay_mode} Reco Efficiency vs Theta', 'True Visible #theta [rad]', '#epsilon')

    hPhiEff = hists_dict[decay_num]["mc_tau_phi"].Clone(f'{name}_phi_eff')
    hPhiEff.Divide(hPhiEff, hists_dict[decay_num]['true_vis_phi'], 1, 1, 'B')
    style_hist(hPhiEff, 418, 2, f'{decay_mode} Reco Efficiency vs Phi', 'True Visible #phi [rad]', '#epsilon')

    #reco reco efficiency hists
    hPtEffReco = hists_dict[decay_num]["linked_correct_reco_mc_tau_pt"].Clone(f'{name}_pt_eff_reco')
    hPtEffReco.Divide(hPtEffReco, hists_dict[decay_num]['true_vis_pT'], 1, 1, 'B')
    style_hist(hPtEffReco, 6, 2, f'Reco-{decay_mode} Reco Efficiency vs Pt', 'True Visible P_{t} [GeV/c]', '#epsilon')

    hThetaEffReco = hists_dict[decay_num]["linked_correct_reco_mc_tau_theta"].Clone(f'{name}_theta_eff_reco')
    hThetaEffReco.Divide(hThetaEffReco, hists_dict[decay_num]['true_vis_theta'], 1, 1, 'B')
    style_hist(hThetaEffReco, 7, 2, f'Reco-{decay_mode} Reco Efficiency vs Theta', 'True Visible #theta [rad]', '#epsilon')

    hPhiEffReco = hists_dict[decay_num]["linked_correct_reco_mc_tau_phi"].Clone(f'{name}_phi_eff_reco')
    hPhiEffReco.Divide(hPhiEffReco, hists_dict[decay_num]['true_vis_phi'], 1, 1, 'B')
    style_hist(hPhiEffReco, 418, 2, f'Reco-{decay_mode} Reco Efficiency vs Phi', 'True Visible #phi [rad]', '#epsilon')

    # fake eff hists
    hPtEffFake = hPtEff.Clone(f'fake_{name}_pt_eff')
    hPtEffFake.Add(hPtEffReco, -1)
    style_hist(hPtEffFake, 6, 2, f'Fake-{decay_mode} Efficiency vs Pt', 'True Visible P_{t} [GeV/c]', '#epsilon')

    hThetaEffFake = hThetaEff.Clone(f'fake_{name}_theta_eff')
    hThetaEffFake.Add(hThetaEffReco, -1)
    style_hist(hThetaEffFake, 7, 2, f'Fake-{decay_mode} Reco Efficiency vs Theta', 'True Visible #theta [rad]', '#epsilon')

    hPhiEffFake = hPhiEff.Clone(f'fake_{name}_phi_eff')
    hPhiEffFake.Add(hPhiEffReco, -1)
    style_hist(hPhiEffFake, 418, 2, f'Fake-{decay_mode} Reco Efficiency vs Phi', 'True Visible #phi [rad]', '#epsilon')

    # neutral pion eff
    hPi0PtEff = hists_dict[decay_num]["pi_0_matched_pT"].Clone(f'{name}_pi_0_pt_eff')
    hPi0PtEff.Divide(hPi0PtEff, hists_dict[decay_num]['pi_0_true_pT'], 1, 1, 'B')
    style_hist(hPi0PtEff, 6, 2, f'{decay_mode} Neutral Pion Reco Efficiency vs Pt', 'True Visible P_{t} [GeV/c]', '#epsilon')

    hPi0ThetaEff = hists_dict[decay_num]["pi_0_matched_theta"].Clone(f'{name}_pi_0_theta_eff')
    hPi0ThetaEff.Divide(hPi0ThetaEff, hists_dict[decay_num]['pi_0_true_theta'], 1, 1, 'B')
    style_hist(hPi0ThetaEff, 7, 2, f'{decay_mode} Neutral Pion Reco Efficiency vs Theta', 'True Visible #theta [rad]', '#epsilon')

    hPi0PhiEff = hists_dict[decay_num]["pi_0_matched_phi"].Clone(f'{name}_pi_0_phi_eff')
    hPi0PhiEff.Divide(hPi0PhiEff, hists_dict[decay_num]['pi_0_true_phi'], 1, 1, 'B')
    style_hist(hPi0PhiEff, 418, 2, f'{decay_mode} Neutral Pion Reco Efficiency vs Phi', 'True Visible #phi [rad]', '#epsilon')

    # charge pion eff
    hPiPtEff = hists_dict[decay_num]["pi_matched_pT"].Clone(f'{name}_pi_pt_eff')
    hPiPtEff.Divide(hPiPtEff, hists_dict[decay_num]['pi_true_pT'], 1, 1, 'B')
    style_hist(hPiPtEff, 6, 2, f'{decay_mode} Pion Reco Efficiency vs Pt', 'True Visible P_{t} [GeV/c]', '#epsilon')

    hPiThetaEff = hists_dict[decay_num]["pi_matched_theta"].Clone(f'{name}_pi_theta_eff')
    hPiThetaEff.Divide(hPiThetaEff, hists_dict[decay_num]['pi_true_theta'], 1, 1, 'B')
    style_hist(hPiThetaEff, 7, 2, f'{decay_mode} Pion Reco Efficiency vs Theta', 'True Visible #theta [rad]', '#epsilon')

    hPiPhiEff = hists_dict[decay_num]["pi_matched_phi"].Clone(f'{name}_pi_phi_eff')
    hPiPhiEff.Divide(hPiPhiEff, hists_dict[decay_num]['pi_true_phi'], 1, 1, 'B')
    style_hist(hPiPhiEff, 418, 2, f'{decay_mode} Pion Reco Efficiency vs Phi', 'True Visible #phi [rad]', '#epsilon')

# Write to output file
output_file = TFile(args.outputFile, 'RECREATE')
for decay_num in hists_dict.keys():
    for key, value in hists_dict[decay_num].items():
        value.Write()
for hist in general_hists:
    hist.Write()

output_file.Close()

# Draw hists and save as PNG
for hist in general_hists:
    filename = hist.GetTitle() + '.png'
    canvas = TCanvas()
    hist.Draw()
    canvas.SaveAs(filename)

for decay_num in hists_dict.keys():
    for hist in hists_dict[decay_num].values():
        filename = hist.GetTitle() + '.png'
        canvas = TCanvas()
        hist.Draw()
        canvas.SaveAs(filename)

#output the daughter types of the pions
for decay_key in pion_types.keys():
    mode = decay_modes[decay_key]
    name = mode.replace(' ', '_')
    with open(f'{name}_daughter_types.txt', 'w') as file:
        file.write(f'{mode} Daughter Types\n')
        for key, value in pion_types[decay_key].items():
            print(f'Type: {key}, Total Number: {value}\n', file=file)
