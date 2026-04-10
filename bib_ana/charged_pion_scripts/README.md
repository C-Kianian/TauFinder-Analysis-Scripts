Instructions on running the BIB charged pion related scripts:

Get local copy of the scripts:
git clone https://github.com/C-Kianian/TauFinder-Analysis-Scripts.git

cd TauFinder-Analysis-Scripts/bib_ana/charged_pion_scripts

Enter the container and bind to the shared directory: 
apptainer run \
  -B /ospool/uc-shared/project/futurecolliders/:/host/futurecolliders \
  /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.11-amd64

Run a script:
python ./pi_ana_bib.py --inputFile=/host/futurecolliders/gpenn/pions/<path to either pi+/pi-, bib/non bib samples directories> --outputFile=pi_bib_ana.root --charge=<charge of the pion samples, either plus, minus or both if the sample is mixed>

If you encounter a segfault/memory error running over all samples at once, there are three options:

1. Run over less samples, ~5k
2. Exit the container and copy an already merged file, example: cp /home/kianian/future_collider/pions/digi_reco/pi+_digi_reco/bib_results/merged_bib_pi+_digi_reco_9.7k.slcio 
3. Write or ask an LLM for a script to batch merge 10k files with the lcio_merge_events command, if you try to merge 10k files all at once you will encounter the same error as above
