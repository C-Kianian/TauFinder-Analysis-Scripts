import ROOT
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('--inputFile', type=str, default='./tau_ana_default.root')

args = parser.parse_args()

root_file_default = ROOT.TFile(args.inputFile, 'READ')

#overlaying pt hists

#nums = [1,3]
decay_modes = {1: '1P + Neutrals',
               5: '3P + Neutrals'}
#for num in nums:
for decay_num, decay_mode in decay_modes.items():
    name = decay_mode.replace(' ','_').lower()

    h_vis = root_file_default.Get(f'{name}_true_vis_pT')
    h_tot = root_file_default.Get(f'{name}_true_tot_pT')
    #h_vis = root_file_default.Get(f'{num}p_true_vis_pT')
    #h_tot = root_file_default.Get(f'{num}p_true_tot_pT')

    h_vis.SetTitle("True " + decay_mode + " Tau p_{T}")
    h_vis.GetXaxis().SetTitle('p_{T} [GeV/c]')
    #if num == 1: h_vis.GetYaxis().SetRangeUser(0, 300)
    #else: h_vis.GetYaxis().SetRangeUser(0, 200)

    h_vis.SetLineColor(ROOT.kBlue)
    h_tot.SetLineColor(ROOT.kRed)

    h_vis.SetLineWidth(3)
    h_tot.SetLineWidth(3)

    h_vis.SetStats(0)
    h_tot.SetStats(0)

    c = ROOT.TCanvas('c', 'overlay', 800, 600)

    h_vis.Draw()
    h_tot.Draw('SAME')

    c.Update()

    legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.90)
    legend.AddEntry(h_vis, 'Visible p_{T}', 'lp')
    legend.AddEntry(h_tot, 'Total p_{T}', 'lp')
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()
    filename = "True " + decay_mode + " Tau p_{T}" + ".png"
    c.SaveAs(filename)

    del c


#overlaying eff hists

names = ['pt', 'phi', 'theta']
variables = ['p_{T}', '#phi', '#theta']
units = ['GeV/c', 'rad', 'rad']

for i in range(3):
    name1 = decay_modes[1].replace(' ','_').lower()
    name3 = decay_modes[5].replace(' ','_').lower()

    h_one = root_file_default.Get(f'{name1}_{names[i]}_eff')
    h_three = root_file_default.Get(f'{name3}_{names[i]}_eff')

    #h_one = root_file_default.Get('reco_1p_' + names[i] + '_eff')
    #h_three = root_file_default.Get('reco_3p_' + names[i] + '_eff')

    h_one.SetTitle('Benchmark Reconstruction Efficiencies vs ' + variables[i])
    h_one.GetXaxis().SetTitle('True Visible #tau^{-} ' + variables[i] + ' [' + units[i] + ']')
    h_one.GetYaxis().SetRangeUser(0.0, 1.0)

    h_one.SetLineColor(ROOT.kRed)
    h_three.SetLineColor(ROOT.kOrange)

    h_one.SetLineWidth(2)
    h_three.SetLineWidth(2)

    h_one.SetMarkerStyle(ROOT.kDot)
    h_three.SetMarkerStyle(ROOT.kDot)

    h_one.SetMarkerColor(ROOT.kRed)
    h_three.SetMarkerColor(ROOT.kOrange)

    h_one.SetMarkerSize(1.5)
    h_three.SetMarkerSize(1.5)

    h_one.SetStats(0)
    h_three.SetStats(0)

    c = ROOT.TCanvas('c', 'overlay', 800, 600)

    h_one.Draw()
    h_three.Draw('SAME')

    c.Update()

    legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.90)
    legend.AddEntry(h_one, 'One-Prong', 'lp')
    legend.AddEntry(h_three, 'Three-Prong', 'lp')
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.Draw()

    c.Update()
    filename = "Benchmark Tau Reconstruction Efficiencies " + variables[i] + '.png'
    c.SaveAs(filename)

    del c


