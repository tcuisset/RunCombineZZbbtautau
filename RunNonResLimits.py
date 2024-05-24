import os,sys
import ROOT
import numpy as np
import concurrent.futures
ROOT.gROOT.SetBatch(True)

import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.CMS)

'''
python3 RunNonResLimits.py --ver ul_2018_ZZ_v12 \
    --cat cat_ZZ_elliptical_cut_90_resolved_1b,cat_ZZ_elliptical_cut_90_resolved_2b,cat_ZZ_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZZbbtt_kl_1 --prd prod_240523 --grp datacard_zz
'''

#######################################################################
######################### SCRIPT BODY #################################
#######################################################################

if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run",     dest="run",      default=True)
    parser.add_option("--ver",     dest="ver",      default='')
    parser.add_option("--cat",     dest="cat",      default='')
    parser.add_option("--prd",     dest="prd",      default='')
    parser.add_option("--feat",    dest="feat",     default='dnn_ZZbbtt_kl_1')
    parser.add_option("--grp",     dest="grp",      default='datacard_zz')
    parser.add_option("--channels", dest="channels", default="etau,mutau,tautau")
    parser.add_option("--singleThread", action="store_true", help="Don't run in parallel, disable for debugging")
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

    if ',' in options.channels:
        channels = options.channels.split(',')
    else:
        channels = [options.channels]

    prd = options.prd
    grp = options.grp

    if os.environ["USER"] == 'evernazza':
        basedir = '/data_CMS/cms/' + os.environ["USER"][1:] + '/cmt/CreateDatacards/'
    else:
        basedir = '/data_CMS/cms/' + os.environ["USER"] + '/cmt/CreateDatacards/'

    maindir = os.getcwd() 

    cat_name_dict = {
        "cat_ZZ_elliptical_cut_90_resolved_1b": "",
        "cat_ZZ_elliptical_cut_90_resolved_2b": "",
        "cat_ZZ_elliptical_cut_90_boosted_noPNet": "",
        "ZbbHtt_elliptical_cut_90_resolved_1b": "",
        "ZbbHtt_elliptical_cut_90_resolved_1b": "",
        "ZbbHtt_elliptical_cut_90_resolved_2b": "",
        "ZbbHtt_elliptical_cut_90_boosted": "",
    }

    def run_limit(feature, version, category, ch):

        if "boosted" in category:        cat_name = r"Boosted"
        elif "resolved_1b" in category:  cat_name = r"Res 1b"
        elif "resolved_2b" in category:  cat_name = r"Res 2b"
        else:                            cat_name = "category"

        if ch == "etau":                 ch_name = "$\\tau_{e}\\tau_{h}$"
        elif ch == "mutau":              ch_name = "$\\tau_{\\mu}\\tau_{h}$"
        elif ch == "tautau":             ch_name = "$\\tau_{h}\\tau_{h}$"
        else:                            ch_name = ch

        odir = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}'
        datadir = basedir + f'/{version}/{category}/{prd}'
        datafile = datadir + f'{feature}_{grp}_{ch}_os_iso.txt'
        run = int(options.run) == 1

        ch_dir = odir + f'/{ch}'
        os.system('mkdir -p ' + ch_dir)

        print(" ### INFO: Create workspace")
        cmd = f'cd {datadir} && text2workspace.py {feature}_{grp}_{ch}_os_iso.txt -o {ch_dir}/model.root &>{ch_dir}/text2workspace.log'
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
        cmd = f'cd {ch_dir} && combine -M MultiDimFit {ch_dir}/model.root --algo=grid --points 100 {r_range} --preFitValue 1 --expectSignal 1 -t -1 &>multiDimFit.log'
        if run: os.system(cmd)

        LS_file = f'{ch_dir}/higgsCombineTest.MultiDimFit.mH120.root'
        f = ROOT.TFile(LS_file)
        limit = f.Get("limit")

        to_draw = ROOT.TString("2*deltaNLL:r")
        n = limit.Draw( to_draw.Data(), "", "l")

        x = np.ndarray((n), 'd', limit.GetV2())[1:]
        y = np.ndarray((n), 'd', limit.GetV1())[1:]

        plt.figure(figsize=(8, 8))
        plt.plot(x, y, label='Data', color='red', linewidth=3)
        plt.axhline(y=1, color='gray', linestyle='--', linewidth=2)
        plt.axhline(y=3.84, color='gray', linestyle='--', linewidth=2)
        plt.text(x[0] + 0.05, 1 + 0.1, '68% C.L.', fontsize=18)
        plt.text(x[0] + 0.05, 3.84 + 0.1, '95% C.L.', fontsize=18)
        plt.text(0.97, 0.97, cat_name, ha="right", va="top", transform=plt.gca().transAxes, color="black")
        plt.text(0.97, 0.92, ch_name, ha="right", va="top", transform=plt.gca().transAxes, color="black")
        mplhep.cms.label(data=False)

        plt.xlabel('$\\mu$')
        plt.ylabel('-2 $\\Delta LL$')
        plt.title("")
        plt.xlim(x[0], x[-1])
        plt.ylim(0,5)

        plt.savefig(f"{ch_dir}/DeltaNLL.png")
        plt.savefig(f"{ch_dir}/DeltaNLL.pdf")

        if options.singleThread:

            # Creates problems when parallelizing, but it's only to print out the significance so we can drop it
            print(" ### INFO: Run significance extraction")
            cmd = f'cd {datadir} && combine -M Significance {feature}_{grp}_{ch}_os_iso.txt -t -1 --expectSignal=1 --pvalue &> {ch_dir}/PValue.log'
            if run: os.system(cmd)
            if run: os.system(f'mv {datadir}/higgsCombineTest.Significance.mH120.root {ch_dir}/higgsCombineTest.Significance.mH120.pvalue.root')

            LS_file = f'{ch_dir}/higgsCombineTest.Significance.mH120.pvalue.root'
            f_pv = ROOT.TFile(LS_file)
            limit = f.Get("limit")
            limit.GetEntry(0)
            a = limit.limit
            f_pv.Close()

            cmd = f'cd {datadir} && combine -M Significance {feature}_{grp}_{ch}_os_iso.txt -t -1 --expectSignal=1 &> {ch_dir}/Significance.log'
            if run: os.system(cmd)
            if run: os.system(f'mv {datadir}/higgsCombineTest.Significance.mH120.root {ch_dir}/higgsCombineTest.Significance.mH120.significance.root')

            LS_file = f'{ch_dir}/higgsCombineTest.Significance.mH120.significance.root'
            f_sig = ROOT.TFile(LS_file)
            limit = f.Get("limit")
            limit.GetEntry(0)
            b = limit.limit
            f_sig.Close()

            print(" ### INFO: Results for", ch)
            print(" ### p-value     = ", a)
            print(" ### significane = ", b)

        if run: os.chdir(ch_dir)

    if options.singleThread:
        for feature in features:
            for version in versions:
                for category in categories:
                    for channel in channels:
                        run_limit(feature, version, category, channel)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
            futures = []
            for feature in features:
                for version in versions:
                    for category in categories:
                        for channel in channels:
                            futures.append(exc.submit(run_limit, feature, version, category, channel))
            for res in concurrent.futures.as_completed(futures):
                res.result()
