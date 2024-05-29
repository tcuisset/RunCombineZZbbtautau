import os,sys
import numpy as np
import concurrent.futures

import uproot
import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.CMS)

import warnings, logging
warnings.filterwarnings("ignore", message=".*Type 3 font.*")
logging.getLogger('matplotlib').setLevel(logging.ERROR)

# IF YOU ONLY NEED TO CHANGE THE PLOTTING STYLE, RUN WITH THE OPTION: --plot_only

'''
ulimit -s unlimited

python3 RunNonResLimits.py --ver ul_2016_ZZ_v12,ul_2016_HIPM_ZZ_v12,ul_2017_ZZ_v12,ul_2018_ZZ_v12 \
    --cat cat_ZZ_elliptical_cut_90_resolved_1b,cat_ZZ_elliptical_cut_90_resolved_2b,cat_ZZ_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZZbbtt_kl_1 --prd prod_240523 --grp datacard_zz --move_eos --user_eos evernazz

python3 RunAsymptoticLimits.py --ver ul_2016_HIPM_ZbbHtt_v12,ul_2016_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_resolved_1b,cat_ZbbHtt_elliptical_cut_90_resolved_2b,cat_ZbbHtt_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZbbHtt_kl_1 --featureDependsOnMass --prd prod_... --grp datacard_zbbhtt \
    --move_eos --user_eos cuisset

python3 RunAsymptoticLimits.py --ver ul_2016_HIPM_ZttHbb_v12,ul_2016_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_resolved_1b,cat_ZttHbb_elliptical_cut_90_resolved_2b,cat_ZttHbb_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZttHbb_kl_1 --featureDependsOnMass --prd prod_... --grp datacard_ztthbb \
    --move_eos --user_eos cuisset
'''

#######################################################################
######################### SCRIPT BODY #################################
#######################################################################

if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run",          dest="run",          default=True)
    parser.add_option("--ver",          dest="ver",          default='')
    parser.add_option("--cat",          dest="cat",          default='')
    parser.add_option("--prd",          dest="prd",          default='')
    parser.add_option("--feat",         dest="feat",         default='dnn_ZZbbtt_kl_1')
    parser.add_option("--grp",          dest="grp",          default='datacard_zz')
    parser.add_option("--channels",     dest="channels",     default="etau,mutau,tautau")
    parser.add_option("--singleThread", action="store_true", help="Don't run in parallel, disable for debugging")
    parser.add_option("--run_one",      dest="run_one",      default=True,             help='Run each channel or not')
    parser.add_option("--run_ch",       dest="run_ch",       default=True,             help='Combine channels or not')
    parser.add_option("--run_cat",      dest="run_cat",      default=True,             help='Combine categories or not')
    parser.add_option("--run_year",     dest="run_year",     default=True,             help='Combine years or not')
    parser.add_option("--plot_only",    dest="plot_only",    default=False,            action='store_true')
    parser.add_option("--move_eos",     dest="move_eos",     default=False,            action='store_true')
    parser.add_option("--user_eos",     dest="user_eos",     default='evernazz',       help='User Name for lxplus account')
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
    run = int(options.run) == 1
    run_one = options.run_one == True
    run_ch = options.run_ch == True
    run_cat = options.run_cat == True
    run_year = options.run_year == True

    if os.environ["USER"] == 'evernazza':
        cmtdir = '/data_CMS/cms/' + os.environ["USER"][1:] + '/cmt/CreateDatacards/'
    else:
        cmtdir = '/data_CMS/cms/' + os.environ["USER"] + '/cmt/CreateDatacards/'

    maindir = os.getcwd() 

    if "ZZ" in options.ver:
        o_name = 'ZZbbtt'; fancy_name = '$ZZ_{bb\tau\tau}$'
        r_range_single = r_range_single_KinFit = r_range_comb = r_range_comb_KinFit = "--rMin 0 --rMax 2"
        r_range_setPR_comb = r_range_setPR_KinFit = "--setParameterRanges r=0,2"
    else: 
        if "ZbbHtt" in options.ver:
            o_name = 'ZbbHtt'; fancy_name = '$Z_{bb}H_{\tau\tau}$'
        elif "ZttHbb" in options.ver:
            o_name = 'ZttHbb'; fancy_name = '$Z_{\tau\tau}H_{bb}$'
        r_range_single = "--rMin -20 --rMax 25"
        r_range_comb = "--rMin -10 --rMax 15"
        r_range_setPR_comb = "--setParameterRanges r=-10,15"
        r_range_single_KinFit = "--rMin -100 --rMax 100"
        r_range_comb_KinFit = "--rMin -50 --rMax 50"
        r_range_setPR_KinFit = "--setParameterRanges r=-50,50"

    dict_ch_name = {"etau": "$\\tau_{e}\\tau_{h}$", "mutau": "$\\tau_{\\mu}\\tau_{h}$", "tautau": "$\\tau_{h}\\tau_{h}$"}

    def SetStyle(fig, x, line1="", line2="", max=5):
        plt.axhline(y=1, color='gray', linestyle='--', linewidth=2)
        plt.axhline(y=3.84, color='gray', linestyle='--', linewidth=2)
        plt.text(x[0] + 0.05, 1 + 0.1, '68% C.L.', fontsize=18)
        plt.text(x[0] + 0.05, 3.84 + 0.1, '95% C.L.', fontsize=18)
        plt.text(0.03, 0.97, line1, ha="left", va="top", transform=plt.gca().transAxes, color="black", bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
        plt.text(0.03, 0.92, line2, ha="left", va="top", transform=plt.gca().transAxes, color="black", bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
        mplhep.cms.label(data=False)
        plt.xlabel('$\\mu$')
        plt.ylabel('-2 $\\Delta LL$')
        plt.title("")
        plt.xlim(x[0], x[-1])
        plt.ylim(0,max)
        plt.grid()

    def WriteResults (fig, x, y, x_stat, y_stat, sig_file, sig=True, round = 2):
        central = x[np.argmin(y)]
        interval_1sigma = x[np.where(y < 1)]
        min_1sigma = np.abs(min(interval_1sigma)-central)
        max_1sigma = np.abs(max(interval_1sigma)-central)
        interval_1sigma_stat = x_stat[np.where(y_stat < 1)]
        min_1sigma_stat = np.abs(min(interval_1sigma_stat)-central)
        max_1sigma_stat = np.abs(max(interval_1sigma_stat)-central)
        r = central
        up = max_1sigma
        down = min_1sigma
        up_stat = max_1sigma_stat
        down_stat = min_1sigma_stat
        up_syst = np.sqrt(max_1sigma**2 - max_1sigma_stat**2)
        down_syst = np.sqrt(min_1sigma**2 - min_1sigma_stat**2)
        sig = GetLimit(sig_file)

        text = fr"$\mu = 1.00^{{+{up_syst:.{round}f}}}_{{-{down_syst:.{round}f}}}(syst)^{{+{up_stat:.{round}f}}}_{{-{down_stat:.{round}f}}}(stat)$"
        if sig: text += f"\nSignificance = {sig:.{round}f}$\sigma$"
        plt.text(0.03, 0.91, text, ha='left', va='top', transform=plt.gca().transAxes, fontsize='small',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
        sig = GetLimit(sig_file)
        # plt.text(0.03, 0.85, f"Significance = {sig:.{round}f}$\sigma$", ha='left', va='top', transform=plt.gca().transAxes, fontsize='small',
        #     bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
    
    def GetDeltaLL(LS_file):
        file = uproot.open(LS_file)
        limit = file["limit"]
        r = limit.array("r")
        deltaNLL = 2 * limit.array("deltaNLL")
        sorted_indices = np.argsort(r)
        r_sorted = r[sorted_indices]
        deltaNLL_sorted = deltaNLL[sorted_indices]
        return r_sorted[1:], deltaNLL_sorted[1:]

    def GetLimit(LS_file):
        file = uproot.open(LS_file)
        limit_tree = file["limit"]
        limit_value = limit_tree["limit"].array()[0]
        return limit_value

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    def run_single_limit(feature, version, category, ch):

        if "boosted" in category:        cat_name = r"Boosted"
        elif "resolved_1b" in category:  cat_name = r"Res 1b"
        elif "resolved_2b" in category:  cat_name = r"Res 2b"
        else:                            cat_name = "category"

        if "KinFit" in feature: r_range = r_range_single_KinFit
        else:                   r_range = r_range_single

        odir = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}'
        datadir = cmtdir + f'/{version}/{category}/{prd}'
        datafile = datadir + f'/{feature}_{grp}_{ch}_os_iso.txt'

        ch_dir = odir + f'/{ch}'
        os.system('mkdir -p ' + ch_dir)

        print(" ### INFO: Create workspace")
        cmd = f'cd {datadir} && text2workspace.py {datafile} -o {ch_dir}/model.root &>{ch_dir}/text2workspace.log'
        if run: os.system(cmd)

        print(" ### INFO: Run Delta Log Likelihood Scan")
        cmd = f'cd {ch_dir} && combine -M MultiDimFit {ch_dir}/model.root --algo=grid --points 100 {r_range} --preFitValue 1 --expectSignal 1 -t -1 &>multiDimFit.log'
        if run: os.system(cmd)

        LS_file = f'{ch_dir}/higgsCombineTest.MultiDimFit.mH120.root'
        x, y = GetDeltaLL(LS_file)

        fig = plt.figure(figsize=(10, 10))
        plt.plot(x, y, label='Data', color='red', linewidth=3)
        SetStyle(fig, x, cat_name, dict_ch_name[ch])
        ver_short = version.split("ul_")[1].split("_Z")[0] ; cat_short = category.split("_cut_90_")[1]
        plt.savefig(f"{ch_dir}/DeltaNLL_{ver_short}_{cat_short}_{ch}.png")
        plt.savefig(f"{ch_dir}/DeltaNLL_{ver_short}_{cat_short}_{ch}.pdf")

        if options.singleThread:

            # Creates problems when parallelizing, but it's only to print out the significance so we can drop it
            print(" ### INFO: Run significance extraction")
            cmd = f'cd {datadir} && combine -M Significance {datafile} -t -1 --expectSignal=1 --pvalue &> {ch_dir}/PValue.log'
            if run: os.system(cmd)
            if run: os.system(f'mv {datadir}/higgsCombineTest.Significance.mH120.root {ch_dir}/higgsCombineTest.Significance.mH120.pvalue.root')
            LS_file = f'{ch_dir}/higgsCombineTest.Significance.mH120.pvalue.root'
            a = GetLimit(LS_file)

            cmd = f'cd {datadir} && combine -M Significance {datafile} -t -1 --expectSignal=1 &> {ch_dir}/Significance_{ver_short}_{cat_short}_{ch}.log'
            if run: os.system(cmd)
            if run: os.system(f'mv {datadir}/higgsCombineTest.Significance.mH120.root {ch_dir}/higgsCombineTest.Significance.mH120.significance.root')
            LS_file = f'{ch_dir}/higgsCombineTest.Significance.mH120.significance.root'
            b = GetLimit(LS_file)

            print(" ### INFO: Results for", feature, version, category, ch)
            print(" ### p-value     = ", a)
            print(" ### significane = ", b, "\n")

        if run: os.chdir(ch_dir)

    if run_one:

        print("\n ### INFO: Run Single Limits\n")

        if options.singleThread:
            for feature in features:
                for version in versions:
                    for category in categories:
                        for channel in channels:
                            run_single_limit(feature, version, category, channel)
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                futures = []
                for feature in features:
                    for version in versions:
                        for category in categories:
                            for channel in channels:
                                futures.append(exc.submit(run_single_limit, feature, version, category, channel))
                for res in concurrent.futures.as_completed(futures):
                    res.result()

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    def run_comb_channels(feature, version, category):

        combdir = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch'
        print(" ### INFO: Saving combination in ", combdir)
        if run: os.system('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for ch in channels:
            ch_file = cmtdir + f'/{version}/{category}/{prd}/{feature}_{grp}_{ch}_os_iso.txt'
            ch_root = cmtdir + f'/{version}/{category}/{prd}/{feature}_{grp}_{ch}_os_iso.root'
            os.system(f'cp {ch_file} {combdir}/{version}_{category}_{feature}_{grp}_{ch}_os_iso.txt')
            os.system(f'cp {ch_root} {combdir}')
            if os.path.exists(ch_file):
                cmd += f' {ch}={version}_{category}_{feature}_{grp}_{ch}_os_iso.txt'
        cmd += f' > {version}_{feature}_{category}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

        if "KinFit" in feature: r_range = r_range_comb_KinFit ; r_range_setPR = r_range_setPR_KinFit
        else:                   r_range = r_range_comb ; r_range_setPR = r_range_setPR_comb
    
        cmd = f'text2workspace.py {version}_{feature}_{category}_os_iso.txt -o model.root'
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit model.root --algo=singles {r_range} --preFitValue 1 --expectSignal 1 -t -1'
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit model.root --algo=grid --points 100 {r_range} --preFitValue 1 --expectSignal 1 -t -1'
        if run: os.system(cmd)
        cmd = f'combine -M Significance {version}_{feature}_{category}_os_iso.txt -t -1 --expectSignal=1 &> Significance_{version}_{feature}_{category}.log'
        if run: os.system(cmd)

        cmd = f'combine -M MultiDimFit model.root -m 125 -n .bestfit.with_syst {r_range_setPR} --saveWorkspace'\
            ' --preFitValue 1 --expectSignal 1 -t -1 &>MultiDimFit.log'
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit higgsCombine.bestfit.with_syst.MultiDimFit.mH125.root {r_range_setPR} '\
            '--saveWorkspace --preFitValue 1 --expectSignal 1 -t -1 -n .scan.with_syst.statonly_correct --algo grid '\
            '--points 100 --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances &>MultiDimFit_statOnly.log'
        if run: os.system(cmd)

    if run_ch:

        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Combination of channels \n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                    for version in versions:
                        for category in categories:
                            run_comb_channels(feature, version, category)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for version in versions:
                            for category in categories:
                                futures.append(exc.submit(run_comb_channels, feature, version, category))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

        ############################################################################
        print("\n ### INFO: Plot Combination of channels \n")
        ############################################################################

        for feature in features:
            for version in versions:
                for category in categories:

                    if "boosted" in category:        cat_name = r"Boosted"
                    elif "resolved_1b" in category:  cat_name = r"Res 1b"
                    elif "resolved_2b" in category:  cat_name = r"Res 2b"
                    else:                            cat_name = category

                    fig = plt.figure(figsize=(10, 10))
                    cmap = plt.get_cmap('tab10')
                    for i, ch in enumerate(channels):
                        LS_file = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/{ch}/higgsCombineTest.MultiDimFit.mH120.root'
                        x, y = GetDeltaLL(LS_file)
                        plt.plot(x, y, label=dict_ch_name[ch], linewidth=3, color=cmap(i))
                    LS_file = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/higgsCombineTest.MultiDimFit.mH120.root'
                    x, y = GetDeltaLL(LS_file)
                    plt.plot(x, y, label='Combination', linewidth=3, color=cmap(i+1))
                    LS_file = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/higgsCombine.scan.with_syst.statonly_correct.MultiDimFit.mH120.root'
                    x_stat, y_stat = GetDeltaLL(LS_file)
                    plt.plot(x_stat, y_stat, label='Stat-only', linewidth=3, linestyle='--', color=cmap(i+1))
                    plt.legend(loc='upper right', fontsize=18, frameon=True)
                    SetStyle(fig, x, cat_name, "", 8)
                    WriteResults(fig, x, y, x_stat, y_stat, maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/higgsCombineTest.Significance.mH120.root')
                    ver_short = version.split("ul_")[1].split("_Z")[0] ; cat_short = category.split("_cut_90_")[1]
                    plt.savefig(maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/DeltaNLL_{ver_short}_{cat_short}.png')
                    plt.savefig(maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/DeltaNLL_{ver_short}_{cat_short}.pdf')

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    def run_comb_categories(feature, version):

        combdir = maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat'
        print(" ### INFO: Saving combination in ", combdir)
        if run: os.system('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for category in categories:
            cat_file = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/{version}_{feature}_{category}_os_iso.txt'
            cat_short = category.split("_cut_90_")[1]
            cmd += f' {cat_short}={cat_file}'
        cmd += f' > {version}_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

        if "KinFit" in feature: r_range = r_range_comb_KinFit ; r_range_setPR = r_range_setPR_KinFit
        else:                   r_range = r_range_comb ; r_range_setPR = r_range_setPR_comb
    
        cmd = f'text2workspace.py {version}_{feature}_os_iso.txt -o model.root'
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit model.root --algo=singles {r_range} --preFitValue 1 --expectSignal 1 -t -1'
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit model.root --algo=grid --points 100 {r_range} --preFitValue 1 --expectSignal 1 -t -1'
        if run: os.system(cmd)
        cmd = f'combine -M Significance {version}_{feature}_os_iso.txt -t -1 --expectSignal=1 &> Significance_{version}_{feature}.log'
        if run: os.system(cmd)

        cmd = f'combine -M MultiDimFit model.root -m 125 -n .bestfit.with_syst {r_range_setPR} --saveWorkspace'\
            ' --preFitValue 1 --expectSignal 1 -t -1 &>MultiDimFit.log'
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit higgsCombine.bestfit.with_syst.MultiDimFit.mH125.root {r_range_setPR} '\
            '--saveWorkspace --preFitValue 1 --expectSignal 1 -t -1 -n .scan.with_syst.statonly_correct --algo grid '\
            '--points 100 --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances &>MultiDimFit_statOnly.log'
        if run: os.system(cmd)

    if run_cat:

        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Combination of categories \n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                    for version in versions:
                        run_comb_categories(feature, version)

            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for version in versions:
                            futures.append(exc.submit(run_comb_categories, feature, version))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

        ############################################################################
        print("\n ### INFO: Plot Combination of categories \n")
        ############################################################################

        for feature in features:
            for version in versions:
                fig = plt.figure(figsize=(10, 10))
                cmap = plt.get_cmap('tab10')
                for i, category in enumerate(categories):
                    LS_file = maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/higgsCombineTest.MultiDimFit.mH120.root'
                    x, y = GetDeltaLL(LS_file)
                    if "boosted" in category:        cat_name = r"Boosted"
                    elif "resolved_1b" in category:  cat_name = r"Res 1b"
                    elif "resolved_2b" in category:  cat_name = r"Res 2b"
                    plt.plot(x, y, label=cat_name, linewidth=3, color=cmap(i))
                LS_file = maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/higgsCombineTest.MultiDimFit.mH120.root'
                x, y = GetDeltaLL(LS_file)
                plt.plot(x, y, label='Combination', linewidth=3, color=cmap(i+1))
                LS_file = maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/higgsCombine.scan.with_syst.statonly_correct.MultiDimFit.mH120.root'
                x_stat, y_stat = GetDeltaLL(LS_file)
                plt.plot(x_stat, y_stat, label='Stat-only', linewidth=3, linestyle='--', color=cmap(i+1))
                plt.legend(loc='upper right', fontsize=18, frameon=True)
                year = version.split("ul_")[1].split("_Z")[0]
                SetStyle(fig, x, year, "", 8)
                WriteResults(fig, x, y, x_stat, y_stat, maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/higgsCombineTest.Significance.mH120.root')
                ver_short = version.split("ul_")[1].split("_Z")[0]
                plt.savefig(maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/DeltaNLL_{ver_short}.png')
                plt.savefig(maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/DeltaNLL_{ver_short}.pdf')

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    def run_comb_years(feature):

        combdir = maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}'
        print(" ### INFO: Saving combination in ", combdir)
        if run: os.system('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for version in versions:
            year = version.split("ul_")[1].split("_Z")[0]
            ver_file = maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/{version}_{feature}_os_iso.txt'
            if os.path.exists(ver_file):
                cmd += f' Y{year}={ver_file}'
        cmd += f' > FullRun2_{o_name}_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        print(cmd)
        if run: os.system(cmd)

        if "KinFit" in feature: r_range = r_range_comb_KinFit ; r_range_setPR = r_range_setPR_KinFit
        else:                   r_range = r_range_comb ; r_range_setPR = r_range_setPR_comb
    
        cmd = f'text2workspace.py FullRun2_{o_name}_{feature}_os_iso.txt -o model.root'
        print(cmd)
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit model.root --algo=singles {r_range} --preFitValue 1 --expectSignal 1 -t -1'
        print(cmd)
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit model.root --algo=grid --points 100 {r_range} --preFitValue 1 --expectSignal 1 -t -1'
        print(cmd)
        if run: os.system(cmd)
        cmd = f'combine -M Significance FullRun2_{o_name}_{feature}_os_iso.txt -t -1 --expectSignal=1 &> Significance_{feature}.log'
        print(cmd)
        if run: os.system(cmd)

        cmd = f'combine -M MultiDimFit model.root -m 125 -n .bestfit.with_syst {r_range_setPR} --saveWorkspace'\
            ' --preFitValue 1 --expectSignal 1 -t -1 &>MultiDimFit.log'
        print(cmd)
        if run: os.system(cmd)
        cmd = f'combine -M MultiDimFit higgsCombine.bestfit.with_syst.MultiDimFit.mH125.root {r_range_setPR} '\
            '--saveWorkspace --preFitValue 1 --expectSignal 1 -t -1 -n .scan.with_syst.statonly_correct --algo grid '\
            '--points 100 --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances &>MultiDimFit_statOnly.log'
        print(cmd)
        if run: os.system(cmd)

        print(" ### INFO: Produce Full Run 2 Impact Plots")

        cmd = f'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 {r_range_setPR} --doInitialFit --robustFit 1 --parallel 50 ' 
        if run: os.system(cmd)
        cmd = f'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 {r_range_setPR} --doFits --robustFit 1 --parallel 50'
        if run: os.system(cmd)
        cmd = 'combineTool.py -M Impacts -d model.root -m 125 -o impacts.json --parallel 50'
        if run: os.system(cmd)
        cmd = f'plotImpacts.py -i impacts.json -o Impacts_{o_name}_{feature}'
        if run: os.system(cmd)
        if run: os.system('mkdir -p impacts')
        if run: os.system('mv higgsCombine_paramFit* higgsCombine_initialFit* impacts')

        cmd = f'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 {r_range_setPR} --doInitialFit --robustFit 1 --parallel 50 ' +  r" --exclude 'rgx{prop_bin.+}'"
        if run: os.system(cmd)
        cmd = f'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 {r_range_setPR} --doFits --robustFit 1 --parallel 50'+  r" --exclude 'rgx{prop_bin.+}'"
        if run: os.system(cmd)
        cmd = 'combineTool.py -M Impacts -d model.root -m 125 -o impacts_noMCstats.json --parallel 50'+  r" --exclude 'rgx{prop_bin.+}'"
        if run: os.system(cmd)
        cmd = f'plotImpacts.py -i impacts_noMCstats.json -o Impacts_{o_name}_{feature}_NoMCstats'
        if run: os.system(cmd)
        if run: os.system('mkdir -p impacts_noMCstats')
        if run: os.system('mv higgsCombine_paramFit* higgsCombine_initialFit* impacts_noMCstats')

    if run_year:

        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Combination of years \n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                    run_comb_years(feature)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        futures.append(exc.submit(run_comb_years, feature))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

        ############################################################################
        print("\n ### INFO: Plot Combination of years \n")
        ############################################################################

        for feature in features:
            fig = plt.figure(figsize=(10, 10))
            cmap = plt.get_cmap('tab10')
            for i, version in enumerate(versions):
                LS_file = maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/higgsCombineTest.MultiDimFit.mH120.root'
                x, y = GetDeltaLL(LS_file)
                plt.plot(x, y, label=version.split("ul_")[1].split("_Z")[0], linewidth=3, color=cmap(i))
            LS_file = maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/higgsCombineTest.MultiDimFit.mH120.root'
            x, y = GetDeltaLL(LS_file)
            plt.plot(x, y, label='Combination', linewidth=3, color=cmap(i+1))
            LS_file = maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/higgsCombine.scan.with_syst.statonly_correct.MultiDimFit.mH120.root'
            x_stat, y_stat = GetDeltaLL(LS_file)
            plt.plot(x_stat, y_stat, label='Stat-only', linewidth=3, linestyle='--', color=cmap(i+1))
            plt.legend(loc='upper right', fontsize=18, frameon=True)
            SetStyle(fig, x, fancy_name, "", 8)
            WriteResults(fig, x, y, x_stat, y_stat, maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/higgsCombineTest.Significance.mH120.root')
            plt.savefig(maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/DeltaNLL_FullRun2_{o_name}.png')
            plt.savefig(maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/DeltaNLL_FullRun2_{o_name}.pdf')

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    if options.move_eos:

        eos_dir = f'/eos/user/e/evernazz/www/ZZbbtautau/B2GPlots/2024_06_14/{o_name}/Limits/NonRes'
        user = options.user_eos
        print(f" ### INFO: Copy results to {user}@lxplus.cern.ch")
        print(f"           Inside directory {eos_dir}\n")

        # [FIXME] Work-around for mkdir on eos
        os.system(f'mkdir -p TMP_RESULTS_NONRES && cp index.php TMP_RESULTS_NONRES')
        for feature in features:
            os.system(f'cp ' + maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/DeltaNLL_*.p* TMP_RESULTS_NONRES')
            os.system(f'cp ' + maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/*_os_iso.txt TMP_RESULTS_NONRES')
            os.system(f'cp ' + maindir + f'/NonRes/FullRun2_{o_name}/{prd}/{feature}/Impacts* TMP_RESULTS_NONRES')
            for version in versions:
                ver_short = version.split("ul_")[1].split("_Z")[0]
                os.system(f'mkdir -p TMP_RESULTS_NONRES/{ver_short} && cp index.php TMP_RESULTS_NONRES/{ver_short}')
                os.system(f'cp ' + maindir + f'/NonRes/{version}/{prd}/{feature}/Combination_Cat/DeltaNLL*.p* TMP_RESULTS_NONRES/{ver_short}')
                for category in categories:
                    os.system(f'cp ' + maindir + f'/NonRes/{version}/{prd}/{feature}/{category}/Combination_Ch/DeltaNLL*.p* TMP_RESULTS_NONRES/{ver_short}')
        os.system(f'rsync -rltv TMP_RESULTS_NONRES/* {user}@lxplus.cern.ch:{eos_dir}')
        os.system(f'rm -r TMP_RESULTS_NONRES')

