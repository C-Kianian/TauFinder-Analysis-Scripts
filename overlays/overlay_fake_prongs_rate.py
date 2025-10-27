import ROOT
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('--inputFileDefault', type=str, default='./tau_ana_default.root')
parser.add_argument('--inputFileLoose', type=str, default='./tau_ana_loose.root')
parser.add_argument('--inputFileDynamic', type=str, default='./tau_ana_dynamic.root')

args = parser.parse_args()

root_file_default = ROOT.TFile(args.inputFileDefault, 'READ')
root_file_loose = ROOT.TFile(args.inputFileLoose, 'READ')
root_file_dynamic = ROOT.TFile(args.inputFileDynamic, 'READ')



nums = ['1', '3']
names = ['pt', 'phi', 'theta']
variables = ['p_{T}', '#phi', '#theta']
units = ['GeV/c', 'rad', 'rad']

decay_modes = {1: '1P + N',
               5: '3P + N'}

for num in nums:
    for i in range(3):

        #3P0N, 1P0N overlay
        hDefault = root_file_default.Get(f'fake_{num}p_{names[i]}')
        hLoose = root_file_loose.Get(f'fake_{num}p_{names[i]}')
        hDynamic = root_file_dynamic.Get(f'fake_{num}p_{names[i]}')


        #hDefault.SetTitle("Fake " + decay_mode + ' Rate Efficiencies vs ' + variables[i])
        hDefault.SetTitle(f"{num}-Prong Migration Rates vs {variables[i]}")
        hDefault.GetXaxis().SetTitle('True Visible #tau^{-} ' + variables[i] + ' [' + units[i] + ']')
        hDefault.GetYaxis().SetRangeUser(0.0, 1.2)

        hDefault.SetLineColor(ROOT.kRed)
        hLoose.SetLineColor(ROOT.kBlue)
        hDynamic.SetLineColor(ROOT.kViolet)
        #hPion0.SetLineColor(ROOT.kViolet)

        hDefault.SetLineWidth(2)
        hLoose.SetLineWidth(2)
        #hPion0.SetLineWidth(2)
        hDynamic.SetLineWidth(2)


        hDefault.SetMarkerStyle(ROOT.kFullDotLarge)
        hLoose.SetMarkerStyle(ROOT.kMultiply)
        hDynamic.SetMarkerStyle(ROOT.kCircle)
        #hPion0.SetMarkerStyle(ROOT.kCircle)

        hDefault.SetMarkerColor(ROOT.kRed)
        hLoose.SetMarkerColor(ROOT.kBlue)
        hDynamic.SetMarkerColor(ROOT.kViolet)
        #hPion0.SetMarkerColor(ROOT.kViolet)

        # Bold axis fonts
        hDefault.GetXaxis().SetTitleFont(62)
        hDefault.GetXaxis().SetLabelFont(62)
        hDefault.GetYaxis().SetTitleFont(62)
        hDefault.GetYaxis().SetLabelFont(62)

        hDefault.SetMarkerSize(1.5)
        hLoose.SetMarkerSize(1.5)

        hDefault.SetStats(0)
        hLoose.SetStats(0)
        #hPion0.SetStats(0)
        hDynamic.SetStats(0)

        c = ROOT.TCanvas('c', 'overlay', 800, 600)

        hDefault.Draw()
        hLoose.Draw('SAME')
        #hPion0.Draw('SAME')
        hDynamic.Draw('SAME')

        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(42)
        text.SetTextSize(0.04)
        text.SetTextAlign(13)
        text.DrawLatex(0.12, 0.85, "#bf{#it{MAIA Detector Concept}}")
        text.DrawLatex(0.12, 0.80, "Simulated #tau^{-} Gun (no BIB)")

        c.Update()

        legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.90)
        legend.AddEntry(hDynamic, 'Shrinking Cone', 'lp')
        legend.AddEntry(hLoose, 'Loose', 'lp')
        legend.AddEntry(hDefault, 'Default', 'lp')
        legend.SetBorderSize(1)
        legend.SetTextSize(0.028)
        legend.SetFillStyle(0)
        legend.Draw()


        c.Update()

        filename = hDefault.GetTitle() + '.png'
        c.SaveAs(filename)

        del c
