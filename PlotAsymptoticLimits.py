import os,json
import ROOT
import numpy as np
ROOT.gROOT.SetBatch(True)

if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--dirs",    dest="dirs",     default='prod_231129_M900')
    parser.add_option("--feat",    dest="feat",     default='ZZKinFit_mass')
    parser.add_option("--ver",     dest="ver",      default='prod_231129')
    parser.add_option("--ch",      dest="ch",       default='combination')
    (options, args) = parser.parse_args()

    dirs = options.dirs.split(',') if ',' in options.dirs else [options.dirs]
    feat = options.feat
    ver = options.ver
    ch = options.ch

    maindir = os.getcwd() + '/ResLimits'
    os.system('mkdir -p ' + maindir)

    mass_dict = []
    
    mass  = []
    exp   = []
    m1s_t = []
    p1s_t = []
    m2s_t = []
    p2s_t = []

    for dir in dirs:

        limit_file = maindir + f'/{dir}/{ch}/{feat}/limits.json'
        with open(limit_file, 'r') as json_file:
            mass_dict = json.load(json_file)
        m = list(mass_dict.keys())[0]
        mass.append(m)
        exp.append(mass_dict[m]['exp'])
        m1s_t.append(mass_dict[m]['m1s_t'])
        p1s_t.append(mass_dict[m]['p1s_t'])
        m2s_t.append(mass_dict[m]['m2s_t'])
        p2s_t.append(mass_dict[m]['p2s_t'])

    mass  = np.array(mass, dtype=float)
    exp   = np.array(exp, dtype=float)
    m1s_t = np.array(m1s_t, dtype=float)
    p1s_t = np.array(p1s_t, dtype=float)
    m2s_t = np.array(m2s_t, dtype=float)
    p2s_t = np.array(p2s_t, dtype=float)  
        
    G_exp = ROOT.TGraph()
    G_sig1 = ROOT.TGraphAsymmErrors()
    G_sig2 = ROOT.TGraphAsymmErrors()

    ipt = 0
    for m, e, m1_t, p1_t, m2_t, p2_t in zip(mass, exp, m1s_t, p1s_t, m2s_t, p2s_t):

        p2 = p2_t - e
        p1 = p1_t - e
        m2 = e - m2_t
        m1 = e - m1_t

        G_exp.SetPoint(ipt, m, e)
        G_sig1.SetPoint(ipt, m, e)
        G_sig1.SetPointError(ipt, 0, 0, m1, p1)
        G_sig2.SetPoint(ipt, m, e)
        G_sig2.SetPointError(ipt, 0, 0, m2, p2)

        ipt += 1

    G_exp.SetMarkerStyle(24)
    G_exp.SetMarkerColor(4)
    G_exp.SetMarkerSize(0.8)
    G_exp.SetLineColor(ROOT.kBlack)
    G_exp.SetLineWidth(3)
    G_exp.SetLineStyle(2)
    G_exp.SetFillColor(0)

    G_sig1.SetMarkerStyle(0)
    G_sig1.SetMarkerColor(3)
    G_sig1.SetFillColor(ROOT.kGreen+1)
    G_sig1.SetLineColor(ROOT.kGreen+1)
    G_sig1.SetFillStyle(1001)
    
    G_sig2.SetMarkerStyle(0)
    G_sig2.SetMarkerColor(5)
    G_sig2.SetFillColor(ROOT.kOrange)
    G_sig2.SetLineColor(ROOT.kOrange)
    G_sig2.SetFillStyle(1001)

    canvas = ROOT.TCanvas('canvas', 'canvas', 650, 500)
    canvas.SetFrameLineWidth(3)
    canvas.SetBottomMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetLeftMargin(0.15)
    canvas.SetGridx()
    canvas.SetGridy()

    # Outside frame
    frame_bounds = 170, 1030
    hframe = ROOT.TH1F('hframe', '',
                       100, frame_bounds[0], frame_bounds[1])
    hframe.SetMinimum(0.1)

    hframe.GetYaxis().SetTitleSize(0.047)
    hframe.GetXaxis().SetTitleSize(0.055)
    hframe.GetYaxis().SetLabelSize(0.045)
    hframe.GetXaxis().SetLabelSize(0.045)
    hframe.GetXaxis().SetLabelOffset(0.012)
    hframe.GetYaxis().SetTitleOffset(1.2)
    hframe.GetXaxis().SetTitleOffset(1.1)
    hframe.GetYaxis().SetTitle("95% CL on #sigma #times #bf{#it{#Beta}}(X#rightarrowZZ#rightarrow bb#tau#tau) [pb]")
    hframe.GetXaxis().SetTitle("m_{X} [GeV]")
    hframe.SetStats(0)
    ROOT.gPad.SetTicky()
    hframe.Draw()

    legend = ROOT.TLegend(0.7,0.7,0.9,0.89)
    legend.SetBorderSize(0)
    legend.AddEntry(G_sig1, '68% exp.', 'f')
    legend.AddEntry(G_sig2, '95% exp.', 'f')
    legend.AddEntry(G_exp, 'Expected', 'l')
    G_sig2.Draw('3same')
    G_sig1.Draw('3same')
    G_exp.Draw('Lsame')
    
    canvas.Update()
    canvas.SetLogy()

    legend.Draw()

    save_name = maindir + f'/Limits_'+feat+'_'+ver+'_'+ch
    canvas.SaveAs(save_name+'.png')
    canvas.SaveAs(save_name+'.pdf')
