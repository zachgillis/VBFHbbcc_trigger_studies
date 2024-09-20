## Determine trigger efficiencies and create plots, VBF H -> bb/cc analysis
## Zach Gillis


## IMPORTING FILES & PARAMETER SETTINGS

# import libraries
import ROOT
import numpy as np

# create dictionaries to store TTree objects & histograms
data22 = {'mu': {1: 'data22_mu'}, 'vbf': {1: 'data22_vbf'}}
data23 = {'mu': {1: 'data23_mu'},'vbf': {1: 'data23_vbf'}}
mc23a = {'mu': {1: 'mc23a_mu'},'vbf': {1: 'mc23a_vbf'}}
mc23d = {'mu': {1: 'mc23d_mu'},'vbf': {1: 'mc23d_vbf'}}

# import all files
file1 = ROOT.TFile.Open('Ntuples2/version_1/run3_data22_Zmumu_MUONtrigger_20240815_v1.0/data-MVATree/data22.root', 'READ')
data22['vbf'][0] = file1.Get('Nominal')

file2 = ROOT.TFile.Open('Ntuples2/version_1/run3_data22_Zmumu_VBFtrig_20240813_v1.0/data-MVATree/Zmumu_data22.root', 'READ')
data22['mu'][0] = file2.Get('Nominal')

file3 = ROOT.TFile.Open('Ntuples2/version_1/run3_data23_Zmumu_MUONtrigger_20240815_v1.0/data-MVATree/data23.root', 'READ')
data23['vbf'][0] = file3.Get('Nominal')

file4 = ROOT.TFile.Open('Ntuples2/version_1/run3_data23_Zmumu_VBFtrig_20240813_v1.0/data-MVATree/data23.root', 'READ')
data23['mu'][0] = file4.Get('Nominal')

file5 = ROOT.TFile.Open('Ntuples2/version_1/run3_MC23a_Zmumu_MUONtrigger_20240815_v1.0/data-MVATree/data-CxAOD.root', 'READ')
mc23a['mu'][0] = file5.Get('Nominal')

file6 = ROOT.TFile.Open('Ntuples2/version_1/run3_MC23a_Zmumu_VBFtrigger_20240815_v1.0/data-MVATree/data-CxAOD.root', 'READ')
mc23a['vbf'][0] = file6.Get('Nominal')

file7 = ROOT.TFile.Open('Ntuples2/version_1/run3_MC23d_Zmumu_MUONtrigger_20240815_v1.0/data-MVATree/data-CxAOD.root', 'READ')
mc23d['mu'][0] = file7.Get('Nominal')

file8 = ROOT.TFile.Open('Ntuples2/version_1/run3_MC23d_Zmumu_VBFtrigger_20240815_v1.0/data-MVATree/data-CxAOD.root', 'READ')
mc23d['vbf'][0] = file8.Get('Nominal')

# save directory for plots
directory = 'trig_eff_0823'

# parameters for which to determine trigger efficiency
params = ['pTJ1', 'pTJ2', 'etaJ1', 'etaJ2', 'mJJ', 'dPhiJJ', 'dEtaJJ']
files = [data22, data23, mc23a, mc23d]
files_str = ['data22', 'data23', 'mc23a', 'mc23d']

# for pT, m histograms, MeV -> GeV
factor = {
    'pTJ1': 1e3,
    'pTJ2': 1e3,
    'etaJ1': 1,
    'etaJ2': 1,
    'mJJ': 1e3,
    'dPhiJJ': 1,
    'dEtaJJ': 1,
}

# histogram binning
bins = {
    'pTJ1': [80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220],
    'pTJ2': [80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220],
    'etaJ1': [-5, -4, -3.5, -3.25, -3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1, 0, 1, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 5],
    'etaJ2': [-5, -4, -3.5, -3.25, -3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1, 0, 1, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 5],
    'mJJ': [1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2750, 3000, 4000],
    'dPhiJJ': [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
    'dEtaJJ': [4, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.25, 6.75],
}

# unify bin datatype
for key, value in bins.items():
    new = np.array(value).astype(np.float64)
    bins[key] = new
    
# x-axis histogram label
x_label = {
    'pTJ1': '#it{p}_{T, J1} [GeV]',
    'pTJ2': '#it{p}_{T, J2} [GeV]',
    'etaJ1': '#it{#eta}_{J1}',
    'etaJ2': '#it{#eta}_{J2}',
    'mJJ': '#it{m}_{JJ} [GeV]',
    'dPhiJJ': '#Delta#it{#phi}_{JJ}',
    'dEtaJJ': '#Delta#it{#eta}_{JJ}',
}

# create dictionaries to store tefficiency objects
teff = {
    'data22': {},
    'data23': {},
    'mc23a': {},
    'mc23d': {}
}

# MC eventWeight off by 0.1
fill_fact = {
    'mc23a': 0.1,
    'mc23d': 0.1,
    'data22': 1.0,
    'data23': 1.0,
}


## PLOTTING FUNCTIONS

# function to plot ratio plots for two TGraph objects
def tgraph_ratio_plot(graph1c, graph2c, graph1_label, graph2_label, xlabel, ylabel, rlabel, text, save):
    graph1 = graph1c.Clone("graph1_clone")
    graph2 = graph2c.Clone("graph2_clone")

    c = ROOT.TCanvas("c", "c", 800, 800)
    
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.Draw()
    pad1.cd()
    
    graph1.SetLineColor(2)
    graph1.SetLineWidth(2)
    graph1.SetMarkerStyle(25)
    graph1.SetMarkerColor(2)
    
    graph2.SetMinimum(0.5)
    graph2.SetMaximum(1.1)
    
    graph2.SetLineColor(1)
    graph2.SetLineWidth(2)
    graph2.SetMarkerStyle(8)
    graph2.SetMarkerColor(1)
    
    graph2.Draw("AP")
    graph1.Draw("P SAME")
    
    legend = ROOT.TLegend(0.72, 0.78, 0.92, 0.88)
    legend.SetTextSize(0.035)
    legend.AddEntry(graph1, graph1_label, "lp")
    legend.AddEntry(graph2, graph2_label, "lp")
    legend.SetBorderSize(0)
    legend.Draw()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.055)
    latex.DrawLatex(0.2, 0.83, "#it{ATLAS} #bf{Internal}")
    latex.SetTextSize(0.035)
    latex.DrawLatex(0.2, 0.78, f'#bf{{{text}}}')
    
    title_size = 0.1
    label_size = 0.1
    
    graph2.GetYaxis().SetTitle(ylabel)
    graph2.GetYaxis().SetTitleSize(0.5 * title_size)
    graph2.GetYaxis().SetLabelSize(0.5 * label_size)
    graph2.GetXaxis().SetLabelSize(0)
    graph2.GetXaxis().SetTitleSize(0)
    
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.35)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()
    
    ratio_graph = ROOT.TGraphAsymmErrors(graph1.GetN())
    ratio_graph.SetTitle("")
    
    for i in range(graph1.GetN()):
        x1, y1 = np.array([0.0]), np.array([0.0])
        x2, y2 = np.array([0.0]), np.array([0.0])

        graph1.GetPoint(i, x1, y1)
        graph2.GetPoint(i, x2, y2)

        x1, y1 = float(x1), float(y1)
        x2, y2 = float(x2), float(y2)

        if y1 != 0:
            ratio_y = y2 / y1
            ratio_graph.SetPoint(i, x1, ratio_y)

            ey1_low = graph1.GetErrorYlow(i)
            ey1_high = graph1.GetErrorYhigh(i)
            ey2_low = graph2.GetErrorYlow(i)
            ey2_high = graph2.GetErrorYhigh(i)

            ey_low = ratio_y * np.sqrt((ey2_low/y2)**2 + (ey1_low/y1)**2)
            ey_high = ratio_y * np.sqrt((ey2_high/y2)**2 + (ey1_high/y1)**2)

            ex1_low = graph1.GetErrorXlow(i)
            ex1_high = graph1.GetErrorXhigh(i)
            ex2_low = graph2.GetErrorXlow(i)
            ex2_high = graph2.GetErrorXhigh(i)

            ratio_graph.SetPointError(i, ex2_low, ex2_high, ey_low, ey_high)
        else:
            ratio_graph.SetPoint(i, x1, 0)
            ratio_graph.SetPointError(i, 0, 0, 0, 0)

    
    ratio_graph.SetLineColor(ROOT.kBlack)
    ratio_graph.SetLineWidth(2)
    ratio_graph.SetMarkerStyle(22)
    ratio_graph.SetMarkerColor(ROOT.kBlack)
    
    ratio_graph.GetXaxis().SetTitle(xlabel)
    ratio_graph.GetYaxis().SetTitle(rlabel)
    ratio_graph.GetYaxis().SetNdivisions(505)
    
    ratio_graph.GetXaxis().SetLabelSize(label_size)
    ratio_graph.GetXaxis().SetTitleSize(title_size)
    ratio_graph.GetYaxis().SetLabelSize(label_size)
    ratio_graph.GetYaxis().SetTitleSize(title_size)
    
    ratio_graph.GetYaxis().SetTitleOffset(0.5)
    
    ratio_graph.SetMinimum(0.5)
    ratio_graph.SetMaximum(1.5)
    
    ratio_graph.Draw("AP")

    xmin = ratio_graph.GetXaxis().GetXmin()
    xmax = ratio_graph.GetXaxis().GetXmax()

    line = ROOT.TLine(xmin, 1, xmax, 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)
    line.Draw()
    
    pad1.cd()
    graph1.GetXaxis().SetLimits(xmin, xmax)
    graph2.GetXaxis().SetLimits(xmin, xmax)
    
    pad2.cd()
    ratio_graph.GetXaxis().SetLimits(xmin, xmax)
    
    c.SaveAs(save)

# function to plot ratio plots for two TH1D objects
def histogram_ratio_plot(hist1c, hist2c, hist1_label, hist2_label, xlabel, ylabel, rlabel, text, save):
    hist1 = hist1c.Clone("hist1_clone")
    hist2 = hist2c.Clone("hist2_clone")

    c = ROOT.TCanvas("c", "c", 800, 800)
    
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.Draw()
    pad1.cd()
    
    hist1.SetStats(0)
    hist2.SetStats(0)
    
    hist1.SetLineColor(ROOT.kRed+1)
    hist1.SetLineWidth(2)
    hist1.SetMarkerStyle(20)
    hist1.SetMarkerColor(ROOT.kRed+1)
    
    hist2.SetLineColor(ROOT.kBlue+1)
    hist2.SetLineWidth(2)
    hist2.SetMarkerStyle(21)
    hist2.SetMarkerColor(ROOT.kBlue+1)
    
    y_max = max(hist1.GetMaximum(), hist2.GetMaximum()) * 1.2
    hist1.SetMaximum(y_max)
    
    hist1.Draw("HIST")
    hist2.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.7, 0.78, 0.92, 0.88)
    legend.SetTextSize(0.035)
    legend.AddEntry(hist1, hist1_label, "l")
    legend.AddEntry(hist2, hist2_label, "l")
    legend.SetBorderSize(0)
    legend.Draw()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.055)
    latex.DrawLatex(0.2, 0.83, "#it{ATLAS} #bf{Internal}")
    latex.SetTextSize(0.035)
    latex.DrawLatex(0.2, 0.78, f'#bf{{{text}}}')
    
    title_size = 0.1
    label_size = 0.1
    
    hist1.GetYaxis().SetTitle(ylabel)
    hist1.GetYaxis().SetTitleSize(0.5 * title_size)
    hist1.GetYaxis().SetLabelSize(0.5 * label_size)
    hist1.GetXaxis().SetLabelSize(0)
    hist1.GetXaxis().SetTitleSize(0)
    
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.35)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()
    
    ratio_hist = hist1.Clone("ratio_hist")
    ratio_hist.SetTitle("")
    
    ratio_hist.Divide(hist2)
    
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetLineWidth(2)
    ratio_hist.SetMarkerStyle(22)
    ratio_hist.SetMarkerColor(ROOT.kBlack)
    
    ratio_hist.GetXaxis().SetTitle(xlabel)
    ratio_hist.GetYaxis().SetTitle(rlabel)
    ratio_hist.GetYaxis().SetNdivisions(505)
    
    ratio_hist.GetXaxis().SetLabelSize(label_size)
    ratio_hist.GetXaxis().SetTitleSize(title_size)
    ratio_hist.GetYaxis().SetLabelSize(label_size)
    ratio_hist.GetYaxis().SetTitleSize(title_size)
    
    ratio_hist.GetYaxis().SetTitleOffset(0.5)
    
    ratio_hist.SetMinimum(0.5)
    ratio_hist.SetMaximum(1.5)
    
    ratio_hist.Draw("HIST")
    
    line = ROOT.TLine(ratio_hist.GetXaxis().GetXmin(), 1, ratio_hist.GetXaxis().GetXmax(), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)
    line.Draw()
    
    pad1.cd()
    x_min = min(hist1.GetXaxis().GetXmin(), hist2.GetXaxis().GetXmin())
    x_max = max(hist1.GetXaxis().GetXmax(), hist2.GetXaxis().GetXmax())
    
    hist1.GetXaxis().SetLimits(x_min, x_max)
    hist2.GetXaxis().SetLimits(x_min, x_max)
    
    pad2.cd()
    ratio_hist.GetXaxis().SetLimits(x_min, x_max)
    
    c.SaveAs(save)

# function to plot ratio plots for four TH1D objects (for graphs to compare two data taking periods)
def tgraph_ratio_plot_4(graph1c, graph2c, graph3c, graph4c, graph1_label, graph2_label, graph3_label, graph4_label, xlabel, ylabel, rlabel1, rlabel2, text, save):
    graph1 = graph1c.Clone("graph1_clone")
    graph2 = graph2c.Clone("graph2_clone")
    graph3 = graph3c.Clone("graph3_clone")
    graph4 = graph4c.Clone("graph4_clone")

    c = ROOT.TCanvas("c", "c", 800, 800)
    
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.Draw()
    pad1.cd()
    
    graph1.SetLineColor(2)
    graph1.SetLineWidth(2)
    graph1.SetMarkerStyle(25)
    graph1.SetMarkerColor(2)
    
    graph2.SetLineColor(1)
    graph2.SetLineWidth(2)
    graph2.SetMarkerStyle(24)
    graph2.SetMarkerColor(1)

    graph3.SetLineColor(ROOT.kRed)
    graph3.SetLineWidth(2)
    graph3.SetMarkerStyle(21)
    graph3.SetMarkerColor(ROOT.kRed)
    
    graph4.SetLineColor(1)
    graph4.SetLineWidth(2)
    graph4.SetMarkerStyle(8)
    graph4.SetMarkerColor(1)
    
    graph2.SetMinimum(0.5)
    graph2.SetMaximum(1.1)
    
    graph2.Draw("APX")
    graph1.Draw("PX SAME")
    graph3.Draw("PX SAME")
    graph4.Draw("PX SAME")
    
    legend = ROOT.TLegend(0.6, 0.77, 0.92, 0.88)
    legend.SetTextSize(0.035)
    legend.SetNColumns(2)
    legend.AddEntry(graph1, graph1_label, "p")
    legend.AddEntry(graph2, graph2_label, "p")
    legend.AddEntry(graph3, graph3_label, "p")
    legend.AddEntry(graph4, graph4_label, "p")
    legend.SetBorderSize(0)
    legend.Draw()
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.055)
    latex.DrawLatex(0.2, 0.83, "#it{ATLAS} #bf{Internal}")
    latex.SetTextSize(0.035)
    latex.DrawLatex(0.2, 0.78, f'#bf{{{text}}}')
    
    title_size = 0.1
    label_size = 0.1
    
    graph2.GetYaxis().SetTitle(ylabel)
    graph2.GetYaxis().SetTitleSize(0.5 * title_size)
    graph2.GetYaxis().SetLabelSize(0.5 * label_size)
    graph2.GetXaxis().SetLabelSize(0)
    graph2.GetXaxis().SetTitleSize(0)
    
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.35)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()
    
    ratio_graph1 = ROOT.TGraphAsymmErrors(graph1.GetN())
    ratio_graph2 = ROOT.TGraphAsymmErrors(graph3.GetN())
    ratio_graph1.SetTitle("")
    ratio_graph2.SetTitle("")
    
    for i in range(graph1.GetN()):
        x1, y1 = np.array([0.0]), np.array([0.0])
        x2, y2 = np.array([0.0]), np.array([0.0])
        x3, y3 = np.array([0.0]), np.array([0.0])
        x4, y4 = np.array([0.0]), np.array([0.0])

        graph1.GetPoint(i, x1, y1)
        graph2.GetPoint(i, x2, y2)
        graph3.GetPoint(i, x3, y3)
        graph4.GetPoint(i, x4, y4)

        x1, y1 = float(x1), float(y1)
        x2, y2 = float(x2), float(y2)
        x3, y3 = float(x3), float(y3)
        x4, y4 = float(x4), float(y4)

        if y1 != 0:
            ratio_y1 = y2 / y1
            ratio_graph1.SetPoint(i, x1, ratio_y1)

            ey1_low = graph1.GetErrorYlow(i)
            ey1_high = graph1.GetErrorYhigh(i)
            ey2_low = graph2.GetErrorYlow(i)
            ey2_high = graph2.GetErrorYhigh(i)

            ey_low1 = ratio_y1 * np.sqrt((ey2_low/y2)**2 + (ey1_low/y1)**2)
            ey_high1 = ratio_y1 * np.sqrt((ey2_high/y2)**2 + (ey1_high/y1)**2)

            ratio_graph1.SetPointError(i, graph1.GetErrorXlow(i), graph1.GetErrorXhigh(i), ey_low1, ey_high1)
        else:
            ratio_graph1.SetPoint(i, x1, 0)
            ratio_graph1.SetPointError(i, 0, 0, 0, 0)

        if y3 != 0:
            ratio_y2 = y4 / y3
            ratio_graph2.SetPoint(i, x3, ratio_y2)

            ey3_low = graph3.GetErrorYlow(i)
            ey3_high = graph3.GetErrorYhigh(i)
            ey4_low = graph4.GetErrorYlow(i)
            ey4_high = graph4.GetErrorYhigh(i)

            ey_low2 = ratio_y2 * np.sqrt((ey4_low/y4)**2 + (ey3_low/y3)**2)
            ey_high2 = ratio_y2 * np.sqrt((ey4_high/y4)**2 + (ey3_high/y3)**2)

            ratio_graph2.SetPointError(i, graph3.GetErrorXlow(i), graph3.GetErrorXhigh(i), ey_low2, ey_high2)
        else:
            ratio_graph2.SetPoint(i, x3, 0)
            ratio_graph2.SetPointError(i, 0, 0, 0, 0)
    
    ratio_graph1.SetLineColor(ROOT.kBlue)
    ratio_graph1.SetLineWidth(2)
    ratio_graph1.SetMarkerStyle(22)
    ratio_graph1.SetMarkerColor(ROOT.kBlue)
    
    ratio_graph2.SetLineColor(ROOT.kGreen)
    ratio_graph2.SetLineWidth(2)
    ratio_graph2.SetMarkerStyle(23)
    ratio_graph2.SetMarkerColor(ROOT.kGreen)
    
    ratio_graph1.GetXaxis().SetTitle(xlabel)
    ratio_graph1.GetYaxis().SetTitle(rlabel1)
    ratio_graph1.GetYaxis().SetNdivisions(505)
    
    ratio_graph1.GetXaxis().SetLabelSize(label_size)
    ratio_graph1.GetXaxis().SetTitleSize(title_size)
    ratio_graph1.GetYaxis().SetLabelSize(label_size)
    ratio_graph1.GetYaxis().SetTitleSize(title_size)
    
    ratio_graph1.GetYaxis().SetTitleOffset(0.5)
    
    ratio_graph1.SetMinimum(0.5)
    ratio_graph1.SetMaximum(1.5)
    
    ratio_graph1.Draw("APX")

    ratio_graph2.GetXaxis().SetLabelSize(label_size)
    ratio_graph2.GetXaxis().SetTitleSize(title_size)
    ratio_graph2.GetYaxis().SetLabelSize(label_size)
    ratio_graph2.GetYaxis().SetTitleSize(title_size)
    
    ratio_graph2.GetYaxis().SetTitleOffset(0.5)
    
    ratio_graph2.Draw("PX SAME")

    ratio_legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    ratio_legend.SetTextSize(0.08)
    ratio_legend.AddEntry(ratio_graph1, f'{graph2_label} / {graph1_label}', "p")
    ratio_legend.AddEntry(ratio_graph2, f'{graph4_label} / {graph3_label}', "p")
    ratio_legend.SetBorderSize(0)
    ratio_legend.Draw()

    xmin = ratio_graph1.GetXaxis().GetXmin()
    xmax = ratio_graph1.GetXaxis().GetXmax()

    line1 = ROOT.TLine(xmin, 1, xmax, 1)
    line1.SetLineColor(ROOT.kBlack)
    line1.SetLineStyle(2)
    line1.Draw()

    pad1.cd()
    graph1.GetXaxis().SetLimits(xmin, xmax)
    graph2.GetXaxis().SetLimits(xmin, xmax)
    graph3.GetXaxis().SetLimits(xmin, xmax)
    graph4.GetXaxis().SetLimits(xmin, xmax)
    
    pad2.cd()
    ratio_graph1.GetXaxis().SetLimits(xmin, xmax)
    
    c.SaveAs(save)


## FILLING HISTOGRAMS / CREATING TEFFICIENCY OBJECTS

# dictionary to store VBF events indexed by eventNumber
events_vbf = {}
        
for i, file in enumerate(files_str):
    # create TEfficiency objects, 
    for param in params:
        teff[file][param] = ROOT.TEfficiency('eff', f';{x_label[param]};Trigger Efficiency', len(bins[param])-1, np.array(bins[param], dtype='double'))
        files[i]['mu'][param] = ROOT.TH1D(f'{file[1]}, {param}', f';{x_label[param]};Events / bin width', len(bins[param])-1, bins[param])
        files[i]['vbf'][param] = ROOT.TH1D(f'{file[1]}, {param}', f';{x_label[param]};Events / bin width', len(bins[param])-1, bins[param])
        
    # VBF events indexed by eventNumber
    events_vbf[file] = {event_vbf.eventNumber: event_vbf for event_vbf in files[i]['vbf'][0]}

    # if event from muon TTree in VBF TTree, (1) remove event from VBF event dictionary (for speed), (2) fill TEfficiency & VBF histogram. If not, add 'False' to TEfficiency histogram
    for event_mu in files[i]['mu'][0]:
        if event_mu.eventNumber in events_vbf[file]:
            events_vbf[file].pop(event_mu.eventNumber)
            for param in params:
                teff[file][param].Fill(True, getattr(event_mu, param) / factor[param], getattr(event_mu, 'eventWeight', 1) * fill_fact[file])
                files[i]['vbf'][param].Fill(getattr(event_mu, param) / factor[param], getattr(event_mu, 'eventWeight', 1) * fill_fact[file])
        else:
            for param in params:
                teff[file][param].Fill(False, getattr(event_mu, param) / factor[param], getattr(event_mu, 'eventWeight', 1) * fill_fact[file])
        for param in params:
            files[i]['mu'][param].Fill(getattr(event_mu, param) / factor[param], getattr(event_mu, 'eventWeight', 1) * fill_fact[file])
        
    # scale histograms by bin width
    for param in params:
        files[i]['mu'][param].Scale(1.0, 'width')
        files[i]['vbf'][param].Scale(1.0, 'width')

## CREATING PLOTS (change file paths)

# MC histogram ratio plots (for MC23a/MC23d histogram comparison )
for trig in ['vbf', 'mu']:
    for param in params:
        histogram_ratio_plot(
            mc23a[trig][param], 
            mc23d[trig][param], 
            'MC23a', 
            'MC23d', 
            x_label[param], 
            'Events / bin width', 
            'Ratio', 
            'Z#rightarrow #mu#mu CR', 
            f'{directory}/mc_ratio{trig}_{param}.pdf'
        )

# VBF / Mu histogram plots
for param in params:
    for i, file in enumerate(files):
        histogram_ratio_plot(
            file['vbf'][param], 
            file['mu'][param], 
            'VBF & muon triggers', 
            'Muon trigger', 
            x_label[param], 
            'Events / bin width', 
            'VBF Trig. Eff.', 
            f'{files_str[i]}; Z#rightarrow #mu#mu CR', 
            f'{directory}/vbf_mu_hist_{files_str[i]}_{param}.pdf'
        )

# trigger efficiency ratio plots
for param in params:
    tgraph_ratio_plot(
        teff['mc23a'][param].CreateGraph(), 
        teff['data22'][param].CreateGraph(), 
        'MC23a', 
        'Data22', 
        x_label[param], 
        'Trigger Efficiency', 
        'Data / MC', 
        'Z#rightarrow #mu#mu CR', 
        f'{directory}/22_trig_ratio_{param}.pdf'
    )
    tgraph_ratio_plot(
        teff['mc23d'][param].CreateGraph(), 
        teff['data23'][param].CreateGraph(), 
        'MC23d', 
        'Data23', 
        x_label[param], 
        'Trigger Efficiency', 
        'Data / MC', 
        'Z#rightarrow #mu#mu CR', 
        f'{directory}/23_trig_ratio_{param}.pdf'
    )

# trigger efficiency ratio plots, ALL
for param in params:
    tgraph_ratio_plot_4(
        teff['mc23a'][param].CreateGraph(), 
        teff['data22'][param].CreateGraph(), 
        teff['mc23d'][param].CreateGraph(), 
        teff['data23'][param].CreateGraph(), 
        'MC23a', 
        'Data22', 
        'MC23d', 
        'Data23', 
        x_label[param], 
        'Trigger Efficiency', 
        'Data / MC', 
        'Data / MC', 
        'Z#rightarrow #mu#mu CR', 
        f'{directory}/comb_trig_ratio_{param}.pdf'
    )