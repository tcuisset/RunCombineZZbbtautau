import os,sys
import ROOT
import numpy as np
import concurrent.futures
ROOT.gROOT.SetBatch(True)

#######################################################################
######################### SCRIPT BODY #################################
#######################################################################

def CanvasCreator(dims, margins=0.11):
    if len(dims) == 2: canvas = ROOT.TCanvas( "c", "c", dims[0], dims[1] )
    elif len(dims) == 4: canvas = ROOT.TCanvas( "c", "c", dims[0], dims[1], dims[2], dims[3] )
    else: sys.exit("[ERROR] dims argument (1) either list of len 2 or 4")

    canvOptions = {"SetTitle": "", "SetGrid": True}
    for opt in canvOptions.keys():
        getattr(canvas, opt)(canvOptions[opt])

    ROOT.gPad.SetRightMargin(margins)
    ROOT.gPad.SetLeftMargin(margins)
    ROOT.gPad.SetBottomMargin(margins)
    ROOT.gPad.SetTopMargin(margins)
    ROOT.gPad.SetFrameLineWidth(3)

    canvas.Update()
    return canvas

'''
python3 run_combine.py --ver ul_2016_ZZ_v12,ul_2016_HIPM_ZZ_v12,ul_2017_ZZ_v12,ul_2018_ZZ_v12 \
    --cat cat_ZZ_elliptical_cut_90_etau,cat_ZZ_elliptical_cut_90_mutau,cat_ZZ_elliptical_cut_90_tautau \
    --feat dnn_ZZbbtt_kl_1 --prd prod_240328
'''

if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run",     dest="run",      default=True)
    parser.add_option("--ver",     dest="ver",      default='')
    parser.add_option("--cat",     dest="cat",      default='')
    parser.add_option("--prd",     dest="prd",      default='')
    parser.add_option("--feat",    dest="feat",     default='dnn_ZZbbtt_kl_1')
    parser.add_option("--grp",     dest="grp",      default='datacard_zz')
    parser.add_option("--singleThread", action="store_false", help="Don't run in parallel, disable for debugging")
    (options, args) = parser.parse_args()

    if ',' in options.ver:
        versions = options.ver.split(',')
    else:
        versions = [options.ver]
    
    if ',' in options.cat:
        categories = options.cat.split(',')
    else:
        categories = [options.cat]

    if ',' in options.feat:
        features = options.feat.split(',')
    else:
        features = [options.feat]

    prd = options.prd
    grp = options.grp

    if os.environ["USER"] == 'evernazza':
        basedir = '/data_CMS/cms/' + os.environ["USER"][1:] + '/cmt/CreateDatacards/'
    else:
        basedir = '/data_CMS/cms/' + os.environ["USER"] + '/cmt/CreateDatacards/'

    maindir = os.getcwd() 


    def run_limit(feature, version, category):
        if 'etau' in category:      ch = 'etau'
        if 'mutau' in category:     ch = 'mutau'
        if 'tautau' in category:    ch = 'tautau'
        odir = maindir + f'/NonRes/{version}/{prd}/{feature}'
        datadir = basedir + f'/{version}/{category}/{prd}'
        datafile = datadir + f'{feature}_{grp}_{ch}_os_iso.txt'
        run = int(options.run) == 1

        ch_dir = odir + f'/{ch}'
        os.system('mkdir -p ' + ch_dir)

        print(" ### INFO: Move to data directory ", datadir)
        if run: os.chdir(datadir)
        # os.system('ls') #DEBUG

        print(" ### INFO: Create workspace")
        cmd = f'text2workspace.py {feature}_{grp}_{ch}_os_iso.txt -o {ch_dir}/model.root'
        if run: os.system(cmd)

        print(" ### INFO: Run Delta Log Likelihood Scan")
        if "ZZ" in version:
            r_range = "--rMin 0 --rMax 2"
        elif "ZbbHtt" in version or "ZttHbb" in version:
            if "ZHKinFit_mass" in feature:
                r_range = "--rMin -100 --rMax 100"
            else: # DNN
                r_range = "--rMin -20 --rMax 25"
        else:
            raise ValueError("COuld not determine ZZ or ZH analysis")
        cmd = f'combine -M MultiDimFit {ch_dir}/model.root --algo=grid --points 100 {r_range} --preFitValue 1 --expectSignal 1 -t -1'
        if run: os.chdir(ch_dir)
        if run: os.system(cmd)

        LS_file = f'{ch_dir}/higgsCombineTest.MultiDimFit.mH120.root'
        f = ROOT.TFile(LS_file)
        limit = f.Get("limit")

        to_draw = ROOT.TString("2*deltaNLL:r")
        n = limit.Draw( to_draw.Data(), "", "l")

        x = np.ndarray((n), 'd', limit.GetV2())[1:]
        y = np.ndarray((n), 'd', limit.GetV1())[1:]

        graphScan = ROOT.TGraph(x.size,x,y)

        graphScan.SetTitle("")
        graphScan.SetLineWidth(3)
        graphScan.SetLineColor(ROOT.kRed)
        graphScan.GetYaxis().SetTitle("-2 #Delta LL")
        graphScan.GetXaxis().SetTitle('#mu')

        x_min = graphScan.GetXaxis().GetXmin()
        x_max = graphScan.GetXaxis().GetXmax()

        o_sigma = ROOT.TLine(x_min, 1, x_max, 1)
        o_sigma.SetLineStyle(7)
        o_sigma.SetLineWidth(2)
        o_sigma.SetLineColor(ROOT.kGray+2)
        t_sigma = ROOT.TLine(x_min, 3.84, x_max, 3.84)
        t_sigma.SetLineStyle(7)
        t_sigma.SetLineWidth(2)
        t_sigma.SetLineColor(ROOT.kGray+2)

        c = CanvasCreator([800,800], margins=0.11)
        c.cd()
        graphScan.Draw()

        t1 = ROOT.TLatex(0.11, 0.91, "#scale[1.5]{CMS} Private work ("+ch+")")
        t1.SetTextSize(0.03)
        t1.SetNDC(True)
        t1.Draw("SAME")

        t2 = ROOT.TLatex(0.6, 0.91, "2018, (13 TeV) 59.7 fb^{-1}")
        t2.SetTextSize(0.03)
        t2.SetNDC(True)
        t2.Draw("SAME")

        o_sigma.Draw("SAME")
        t_sigma.Draw("SAME")
        c.Update()

        OS = ROOT.TLatex()
        OS.SetTextFont(42)
        OS.SetTextSize(0.03)
        OS.DrawLatex( x[0]*1.1, 1+0.1, '68% C.L.' )
        TS = ROOT.TLatex()
        TS.SetTextFont(42)
        TS.SetTextSize(0.03)
        TS.DrawLatex( x[0]*1.1, 3.84+0.1, '95% C.L.' )

        c.SaveAs(f"{ch_dir}/DeltaNLL.png")
        c.SaveAs(f"{ch_dir}/DeltaNLL.pdf")
        c.Close()

        print(" ### INFO: Run significance extraction")
        cmd = f'combine -M Significance {feature}_{grp}_{ch}_os_iso.txt -t -1 --expectSignal=1 --pvalue &> {ch_dir}/PValue.log'
        os.chdir(datadir)
        if run: os.system(cmd)
        if run: os.system(f'mv {datadir}/higgsCombineTest.Significance.mH120.root {ch_dir}/higgsCombineTest.Significance.mH120.pvalue.root')

        LS_file = f'{ch_dir}/higgsCombineTest.Significance.mH120.pvalue.root'
        f = ROOT.TFile(LS_file)
        limit = f.Get("limit")
        limit.GetEntry(0)
        a = limit.limit

        cmd = f'combine -M Significance {feature}_{grp}_{ch}_os_iso.txt -t -1 --expectSignal=1 &> {ch_dir}/Significance.log'
        os.chdir(datadir)
        if run: os.system(cmd)
        if run: os.system(f'mv {datadir}/higgsCombineTest.Significance.mH120.root {ch_dir}/higgsCombineTest.Significance.mH120.significance.root')

        LS_file = f'{ch_dir}/higgsCombineTest.Significance.mH120.significance.root'
        f = ROOT.TFile(LS_file)
        limit = f.Get("limit")
        limit.GetEntry(0)
        b = limit.limit

        print(" ### INFO: Results for", ch)
        print(" ### p-value     = ", a)
        print(" ### significane = ", b)

        if run: os.chdir(ch_dir)

    if options.singleThread:
        for feature in features:
            for version in versions:
                for category in categories:
                    run_limit(feature, version, category)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
            for feature in features:
                for version in versions:
                    for category in categories:
                        exc.submit(run_limit, feature, version, category)
 ### INFO: Results for etau
 ### p-value     =  0.0005750532017919205
 ### significane =  3.2509733438179635

 ### INFO: Results for mutau
 ### p-value     =  1.0108180003358628e-09
 ### significane =  5.9960589968567835

 ### INFO: Results for tautau
 ### p-value     =  6.192409029815178e-06
 ### significane =  4.370701189865533

#################### New

 ### INFO: Results for etau
 ### p-value     =  9.918416846514268e-07
 ### significane =  4.755079503776725
# r :    +1.000   -0.212/+0.215 (68%)

 ### INFO: Results for mutau
 ### p-value     =  4.072216424367219e-14
 ### significane =  7.467965274049631
# r :    +1.000   -0.136/+0.137 (68%)

 ### INFO: Results for tautau
 ### p-value     =  9.954731305062541e-10
 ### significane =  5.998543959380103
# r :    +1.000   -0.173/+0.178 (68%)

# For the combination:
# combineCards.py Name1=dnn_new_zzbbtt_kl_1_etau_os_iso.txt.txt Name2=dnn_new_zzbbtt_kl_1_mutau_os_iso.txt Name3=dnn_new_zzbbtt_kl_1_tautau_os_iso.txt  > dnn_comb.txt
# text2workspace.py dnn_comb.txt -o model.root
# combine -M MultiDimFit model.root --algo=singles --rMin 0 --rMax 2 --preFitValue 1 --expectSignal 1 -t -1
# combine -M MultiDimFit model.root --algo=grid --points 100 --rMin 0 --rMax 2 --preFitValue 1 --expectSignal 1 -t -1
# combine -M Significance dnn_comb.txt -t -1 --expectSignal=1

 ### INFO: Results for combination
 ### p-value     =  2.53842e-26
 ### significane =  10.5501
# r :    +1.000   -0.097/+0.098 (68%)