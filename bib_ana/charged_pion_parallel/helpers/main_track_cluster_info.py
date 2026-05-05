#!/usr/bin/env python3

import glob
import math
import numpy as np
import ROOT
from ROOT import TH1F, TFile
import pyLCIO
from pyLCIO import IOIMPL
from .geometry import eta, theta_region, delta_phi
#from .track_truth_match import build_rel_nav, system_to_relname

ROOT.gStyle.SetOptFit(111)
ROOT.gStyle.SetOptStat("nemruo") #for uf/of info

# helper functions
def get_theta(px, py, pz):
    pt = math.sqrt(px**2 + py**2)
    return math.atan2(pt, pz)

def get_phi(px, py):
    return math.atan2(py, px)

def getDr(MCP, track):
    trk_omega, trk_tan_lambda, trk_phi = (
        track.getOmega(), #signed curvature of track in 1/mm
        track.getTanLambda(), #"dip angle" in r-z at reference point
        track.getPhi(),
    )
    trk_theta = (math.pi / 2) - math.atan(trk_tan_lambda)
    mcx, mcy, mcz = (MCP.getMomentum()[0], MCP.getMomentum()[1], MCP.getMomentum()[2])
    mc_theta = get_theta(mcx, mcy, mcz)
    mc_phi = get_phi(mcx,mcy)

    dtheta = mc_theta - trk_theta
    dphi = mc_phi - trk_phi

    if dphi > math.pi:
        dphi = math.fabs(dphi - 2 * math.pi)

    return math.sqrt(dtheta * dtheta + dphi * dphi)

BFIELD = 5.0 # Taking 5 T for MAIA and 3.57 T for MuColl_v1.
FACTOR = 3e-4 # conversion factor to take T calculation to GeV pT

# main processing function
# TODO: Clean some of the below parameters
# they should be input args in some sense

hit_collection_names = [
    "VBTrackerHitsConed",
    "VETrackerHitsConed",
    "IBTrackerHitsConed",
    "IETrackerHitsConed",
    "OBTrackerHitsConed",
    "OETrackerHitsConed",
    "VBTrackerHitsRelationsConed",
    "VETrackerHitsRelationsConed",
    "IBTrackerHitsRelationsConed",
    "IETrackerHitsRelationsConed",
    "OBTrackerHitsRelationsConed",
    "OETrackerHitsRelationsConed",
    "VertexBarrelCollectionConed",
    "VertexEndcapCollectionConed",
    "InnerTrackerBarrelCollectionConed",
    "InnerTrackerEndcapCollectionConed",
    "OuterTrackerBarrelCollectionConed",
    "OuterTrackerEndcapCollectionConed"
]

# hit_collection_names = [
#     "VBTrackerHits",
#     "VETrackerHits",
#     "IBTrackerHits",
#     "IETrackerHits",
#     "OBTrackerHits",
#     "OETrackerHits",
#     "VBTrackerHitsRelations",
#     "VETrackerHitsRelations",
#     "IBTrackerHitsRelations",
#     "IETrackerHitsRelations",
#     "OBTrackerHitsRelations",
#     "OETrackerHitsRelations",
#     "VertexBarrelCollection",
#     "VertexEndcapCollection",
#     "InnerTrackerBarrelCollection",
#     "InnerTrackerEndcapCollection",
#     "OuterTrackerBarrelCollection",
#     "OuterTrackerEndcapCollection"
# ]

# copying this in from the tracking script although I don't think it works
hit_collection_mask = {key:True for key in hit_collection_names}

allowed_pdgs = {
    'plus': {211},
    'minus': {-211},
    'both': {211, -211},
    'none': {211, -211}
}

# hard-coding a positively charged pion
charge = "none"
ignore_charge = True

# Hist setup
PT_MIN, PT_MAX, PT_BINS = 0.0, 1000.0, 160

ETA_MIN, ETA_MAX, ETA_BINS = 0.0, 2.5, 20

THETA_BINS = 20

PHI_BINS = 24

EFF_MIN, EFF_MAX = 0.0, 1.2

def book(h):
    h.SetDirectory(0)
    return h

# Regional pT histograms
regions = ['barrel', 'centbarrel', 'transition', 'endcap']
regional_maxs = {
    'plus': [800, 600, 400, 200],
    'minus': [800, 600, 400, 200],
    'both': [1600, 1200, 800, 400],
    'none': [800, 600, 400, 200]
}

regional_max = regional_maxs[charge]

def process_set(pattern, max_events):

    hists = {}

    # booking a ton of TH1Fs
    for region in regions:
        # initializing ROOT TH1Fs for filling
        # Histograms for counting number of MCPs
        hists[f"fMCPt_{region}"] = book(TH1F(f'mc_pion_pt_{region}', f'MC Charged Pion p_{{T}} ({region})', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fMCTheta_{region}"] = book(TH1F(f'mc_pion_theta_{region}', f'MC Charged Pion #theta ({region});#theta_reco [rad];Entries', THETA_BINS, 0, np.pi))

        # Histograms for counting tracking efficiency
        hists[f"fTrackPt_{region}"] = book(TH1F(f'matched_track_pt_{region}', f'Matched Track p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fTrackTheta_{region}"] = book(TH1F(f'matched_track_theta_{region}', f'Matched Track #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

        ####### NEW HISTOGRAMS FOR TRACK-CLUSTER MATCHING FAILURE ########
        # Histograms for counting all/matched/failed track-cluster matching events
        # pfos
        hists[f"fTrkClsPFOs_all_{region}"] = book(TH1F(f'trk_cls_all_PFOs_{region}', f'# Charged PFOs per All (Matched + Failed) Track-Cluster Event ({region});Number of Charged PFOs;Entries', 5, 0, 5))
        hists[f"fTrkClsPFOs_match_{region}"] = book(TH1F(f'trk_cls_match_PFOs_{region}', f'# Charged PFOs per  Matched Track-Cluster Event ({region});Number of Charged PFOs;Entries', 5, 0, 5))
        hists[f"fTrkClsPFOs_fail_{region}"] = book(TH1F(f'trk_cls_fail_PFOs_{region}', f'# Charged PFOs per Failed Track-Cluster Event ({region});Number of Charged PFOs;Entries', 5, 0, 5))
        # tracks
        hists[f"fTrkClsTracks_all_{region}"] = book(TH1F(f'trk_cls_all_tracks_{region}', f'# Charged Tracks per All (Matched + Failed) Track-Cluster Event ({region});Number of Charged PFOs;Entries', 5, 0, 5))
        hists[f"fTrkClsTracks_match_{region}"] = book(TH1F(f'trk_cls_match_track_{region}', f'# Charged Tracks per Matched Track-Cluster Event({region});Number of Charged PFOs;Entries', 5, 0, 5))
        hists[f"fTrkClsTracks_fail_{region}"] = book(TH1F(f'trk_cls_fail_track_{region}', f'# Charged Tracks per Failed Track-Cluster Event({region});Number of Charged PFOs;Entries', 5, 0, 5))

        ####### NEW HISTOGRAMS FOR TRACK-CLUSTER MATCHING QUALITY ########
        # Histograms for checking track-cluster matched event quality (normalized and unnormalized)
        hists[f"fTrkClsQualNorm_{region}"] = book(TH1F(f'trk_cls_match_qual_norm_{region}', f'Normalized Matched Track-Cluster p_{{T}}_{{track}} - E_{{cluster}} ({region});p_{{T}}_{{track}} - E_{{cluster}} / p_{{T}}_{{track}};Entries', 100, -0.3, 0.3))
        hists[f"fTrkClsQual_{region}"] = book(TH1F(f'trk_cls_match_qual_{region}', f'Matched Track-Cluster p_{{T}}_{{track}} - E_{{cluster}} ({region});p_{{T}}_{{track}}-E_{{cluster}};Entries', 500, -50, 50))
        # hists for cluster counts
        hists[f"fTrkClsNumCls_{region}"] = book(TH1F(f'trk_cls_match_num_cls_{region}', f'# of Clusters per Matched Track-Cluster Event ({region});# of Clusters;Entries', 5, 0, 5))

        # Histograms for counting track-cluster matching efficiency
        hists[f"fTrkClsPt_{region}"] = book(TH1F(f'trk_cls_match_pt_{region}', f'Matched Track-Cluster p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fTrkClsTheta_{region}"] = book(TH1F(f'trk_cls_match_theta_{region}', f'Matched Track-Cluster p_{{T}} #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

        # Histograms for counting reco charged pions
        hists[f"fMatchedPt_{region}"] = book(TH1F(f'mc_matched_pt_{region}', f'Matched Best Reco Charged Pion MC p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fMatchedTheta_{region}"] = book(TH1F(f'mc_matched_theta_{region}', f'Matched Best Charged Reco Pion MC #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

    files = sorted(glob.glob(pattern))
    selected_pdgs = allowed_pdgs[charge]
    event_count = 0

    (f"Found {len(files)} files")

    for fname in files:

        if event_count >= max_events:
            break

        reader = IOIMPL.LCFactory.getInstance().createLCReader()
        print("fname: ", fname)
        collection_names = ['MCParticle', 'PandoraPFOs', 'SelectedTracks']#, 'MCParticle_SelectedTracks']# + hit_collection_names
        #reader.setReadCollectionNames(collection_names)
        reader.open(fname)

        evt = reader.readNextEvent()
        event_count += 1
        print("Event count: ", event_count)

        if "MCParticle" not in evt.getCollectionNames():
            # there seems to be an issue with some non-BIB files. The total number is small, so skipping them should be fine.
            print("Event seems bugged! Skipping.")
            continue

        mcs = evt.getCollection('MCParticle')
        # Best MC charged pion
        best_mc = mcs[0]
        if best_mc is None:
            continue
        mcPDG = best_mc.getPDG()
        mcMom = best_mc.getMomentum()
        mcPt = math.hypot(mcMom[0], mcMom[1])
        mcTheta = math.acos(mcMom[2] / math.sqrt(mcPt**2 + mcMom[2]**2))
        mcEta = eta(mcTheta)
        mcPhi = math.atan2(mcMom[1], mcMom[0])
        mcE = best_mc.getEnergy()
        # find what region the MCP is in
        regs = theta_region(mcTheta)

        # fill histograms according to region
        if regs:
            for reg in regs:
                hists[f"fMCPt_{reg}"].Fill(mcPt)
                hists[f"fMCTheta_{reg}"].Fill(mcTheta)

        tracks = evt.getCollection('SelectedTracks')
        relationsContainer = evt.getCollection('MCParticle_SelectedTracks')
        relation = pyLCIO.UTIL.LCRelationNavigator(relationsContainer)
        related_tracks = relation.getRelatedToObjects(best_mc)
        print("number of relation tracks: ", len(related_tracks))
        # auto-continue if there are no tracks or truth-matched tracks in the event
        if len(tracks) == 0 or len(related_tracks) == 0: continue
        # get num charge tracks:
        charged_tracks = [t for t in tracks if abs(t.getOmega()) > 1e-12]
        num_charge_tracks = len(charged_tracks)


        #### Below is for proper track truth matching. Currently commented out. See comments for why. ####

        # # build relation between hit collections and sub-detector
        # rel_nav = build_rel_nav(evt)

        # hit_collections = []
        # for hname in hit_collection_names:
        #     if(not hit_collection_mask[hname]):
        #         print("I should never hit this, right??")
        #         continue
        #     try:
        #         hit_collections.append(evt.getCollection(hname))
        #     except:
        #         hit_collection_mask[hname] = False
        #         print('\tDid not find hit collection: {}. Disabling...'.format(hname))
        #         pass

        # apparently SelectedTracks don't have any hits associated
        # this is a known bug
        # the tracking results assume that these tracks have high hit purity
        # therefore they will pass the truth-matching requirement
        # therefore take any event with a SelectedTrack as passing tracking requirements
        # for track in tracks:
        #     print("Looping over track")
        #     print("Track has number of hits:", len(track.getTrackerHits()))
        #     print("track omega: ", track.getOmega())
        #     truth_matched_hits = 0
        #     for hit in track.getTrackerHits():
        #         position = hit.getPosition()
        #         print("hit pos:", position)
        #         encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
        #         decoder = pyLCIO.UTIL.BitField64(encoding)
        #         cellID = int(hit.getCellID0())
        #         decoder.setValue(cellID)
        #         detector = decoder["system"].value()
        #         layer = decoder['layer'].value()
        #         if detector == 1 or detector == 2:
        #             LC_pixel_nhit += 1
        #         if detector == 3 or detector == 4:
        #             LC_inner_nhit += 1
        #         if detector == 5 or detector == 6:
        #             LC_outer_nhit += 1
        #         print("detector:", detector)
        #         print("len(getrel objects):", len(rel_nav[system_to_relname[detector]].getRelatedToObjects(hit)))
        #         for sim_hit in rel_nav[system_to_relname[detector]].getRelatedToObjects(hit):
        #             print(type(sim_hit))
        #             print(sim_hit)
        #             mcp_true = sim_hit.getMCParticle()
        #             if mcp_true and abs(mcp_true.getPDG()) == 211:
        #                 truth_matched_hits += 1
        #         truth_hit_ratio = truth_matched_hits / len(track.getTrackerHits())
        #         print("truth_hit_ratio: ", truth_hit_ratio)

        #### end track truth matching section ####


        # this satisfies our tracking efficiency requirements
        # fill tracking efficiency plots
        if regs:
            for reg in regs:
                hists[f"fTrackPt_{reg}"].Fill(mcPt)
                hists[f"fTrackTheta_{reg}"].Fill(mcTheta)

        # initialize reco pis, to be found
        best_reco_charged = None
        best_reco_charged_pt = -1.0

        pfos = evt.getCollection('PandoraPFOs')
        charge_counts = 0

        for pfo in pfos:
            # allowing all charged particles, to filter by charged pions later
            if abs(pfo.getType()) != abs(mcPDG) and abs(pfo.getType()) != 11 and abs(pfo.getType()) != 13: continue # no charge matching case

            charge_counts += 1
            # Pion kinematics
            recoChargedMomDefault = pfo.getMomentum()
            recoChargedPtDefault = math.hypot(recoChargedMomDefault[0], recoChargedMomDefault[1])
            recoChargedTheta = math.acos(recoChargedMomDefault[2] / math.sqrt(recoChargedPtDefault ** 2 + recoChargedMomDefault[2] ** 2))
            recoChargedEta = eta(recoChargedTheta)
            recoChargedPhi = math.atan2(recoChargedMomDefault[1], recoChargedMomDefault[0])

            # dR matching
            dphi = delta_phi(mcPhi, recoChargedPhi)
            dR = math.sqrt(dphi*dphi + (mcEta - recoChargedEta)**2)

            if dR < 0.1:
                if recoChargedPtDefault > best_reco_charged_pt:
                    best_reco_charged_pt = recoChargedPtDefault
                    best_reco_charged = pfo

        if best_reco_charged is None: # check if no default match was found, if so, skip this event
            if regs: # fill hists with num of charged pfos for track-cluster all/failed events
                for reg in regs:
                    # charged pfos
                    hists[f"fTrkClsPFOs_all_{reg}"].Fill(charge_counts) # all
                    hists[f"fTrkClsPFOs_fail_{reg}"].Fill(charge_counts) # fail
                    # charged tracks
                    hists[f"fTrkClsTracks_all_{reg}"].Fill(num_charge_tracks) # all
                    hists[f"fTrkClsTracks_fail_{reg}"].Fill(num_charge_tracks) # fail
            continue

        # fill histograms according to region
        if regs:
            for reg in regs:
                # charged pfos
                hists[f"fTrkClsPFOs_all_{reg}"].Fill(charge_counts) # all
                hists[f"fTrkClsPFOs_match_{reg}"].Fill(charge_counts) # match
                # charged tracks
                hists[f"fTrkClsTracks_all_{reg}"].Fill(num_charge_tracks) # all
                hists[f"fTrkClsTracks_match_{reg}"].Fill(num_charge_tracks) # match

                # track-cluster matching quality
                tracks = best_reco_charged.getTracks()
                best_reco_charged_track_pt = 0.0
                if tracks is None or len(tracks) == 0: continue
                if len(tracks) > 1: print("PFO has more than one track!")
                if len(tracks) != 0:
                    dRTracks = []
                    for i in range(0,len(tracks)):
                        dRTracks.append(getDr(best_mc, tracks[i]))
                    # find closest track to truth pion
                    myTrack = tracks[dRTracks.index(min(dRTracks))]
                    best_reco_charged_track_pt = (BFIELD * FACTOR) / abs(myTrack.getOmega())

                #total_px = 0.0
                #total_py = 0.0
                #for trk in tracks:
                #    omega = trk.getOmega()
                #    if abs(omega) < 1e-12:
                #        continue
                #
                #    pt = (BFIELD * FACTOR) / abs(omega)
                #
                #    phi = trk.getPhi()
                #
                #    px = pt * math.cos(phi)
                #    py = pt * math.sin(phi)
                #
                #    total_px += px
                #    total_py += py
                #
                #best_reco_charged_track_pt = math.hypot(total_px, total_py)

                best_reco_cluster_e = 0.0
                clusters = best_reco_charged.getClusters()
                if clusters is None or len(clusters) == 0: continue
                for cluster in clusters:
                    best_reco_cluster_e += cluster.getEnergy()

                hists[f"fTrkClsNumCls_{reg}"].Fill(len(clusters))

                trkCls = best_reco_charged_track_pt - best_reco_cluster_e
                #ratio = best_reco_cluster_e / best_reco_charged_track_pt
                hists[f"fTrkClsQualNorm_{reg}"].Fill(trkCls/best_reco_cluster_e)
                hists[f"fTrkClsQual_{reg}"].Fill(trkCls)

                hists[f"fTrkClsPt_{reg}"].Fill(mcPt)
                hists[f"fTrkClsTheta_{reg}"].Fill(mcTheta)
                # now add to the charged pion ID histogram
                if abs(best_reco_charged.getType()) == abs(mcPDG):
                    hists[f"fMatchedPt_{reg}"].Fill(mcPt)
                    hists[f"fMatchedTheta_{reg}"].Fill(mcTheta)

        del evt
        del mcs
        del best_mc
        del tracks
        del relationsContainer
        del related_tracks
        del pfos
        del best_reco_charged

    return hists

