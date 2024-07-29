import os,json,concurrent.futures,itertools,sys
import pdb, csv, subprocess
import numpy as np

import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.CMS)

import warnings, logging
warnings.filterwarnings("ignore", message=".*Type 3 font.*")
logging.getLogger('matplotlib').setLevel(logging.ERROR)

# IF YOU ONLY NEED TO CHANGE THE PLOTTING STYLE, RUN WITH THE OPTION: --plot_only

'''
ulimit -s unlimited

python3 RunAsymptoticLimits.py --ver ul_2016_ZZ_v12,ul_2016_HIPM_ZZ_v12,ul_2017_ZZ_v12,ul_2018_ZZ_v12 \
    --cat cat_ZZ_elliptical_cut_90_resolved_1b,cat_ZZ_elliptical_cut_90_resolved_2b,cat_ZZ_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZZbbtt_kl_1 --featureDependsOnMass --prd prod_240528 --grp datacard_zz_res --move_eos --user_eos evernazz \
    --mass 200,210,220,230,240,250,260,280,300,320,350,360,400,450,500,550,600,650,700,750,800,850,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2500,2600,2800,3000,3500,4000,4500,5000 \
    --num 0

python3 RunAsymptoticLimits.py --ver ul_2016_HIPM_ZbbHtt_v12,ul_2016_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_resolved_1b,cat_ZbbHtt_elliptical_cut_90_resolved_2b,cat_ZbbHtt_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZbbHtt_kl_1 --featureDependsOnMass --prd prod_... --grp datacard_zbbhtt_res \
    --mass ... \
    --move_eos --user_eos cuisset

python3 RunAsymptoticLimits.py --ver ul_2016_HIPM_ZttHbb_v12,ul_2016_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_resolved_1b,cat_ZttHbb_elliptical_cut_90_resolved_2b,cat_ZttHbb_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZttHbb_kl_1 --featureDependsOnMass --prd prod_... --grp datacard_ztthbb_res \
    --mass ... \
    --move_eos --user_eos cuisset
'''

# comb_options = '--minimizerAlgo Minuit2'
# comb_options = '--cminDefaultMinimizerType Minuit2'
comb_options = ''

def run_cmd(cmd, run=True, check=True):
    if run:
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            if check:
                raise RuntimeError(f"Command {cmd} failed with exit code {e.returncode}. Working directory : {os.getcwd()}") from None
            else:
                print(f"## ERROR : Command {cmd} failed with exit code {e.returncode}. Working directory : {os.getcwd()}")

#######################################################################
######################### SCRIPT BODY #################################
#######################################################################

if __name__ == "__main__" :

    import matplotlib,argparse
    matplotlib.use('Agg')
    parser = argparse.ArgumentParser("RunAsymptoticLimits")
    
    def makeFlag(arg_name, **kwargs):
        if arg_name.startswith("--"):
            arg_name = arg_name[2:]
        parser.add_argument(f"--{arg_name}", action='store_true', **kwargs)
        parser.add_argument(f"--no_{arg_name}", action='store_false', **kwargs)

    parser.add_argument("--ver",          dest="ver",                   default='')
    parser.add_argument("--cat",          dest="cat",                   default='')
    parser.add_argument("--prd",          dest="prd",                   default='')
    parser.add_argument("--feat",         dest="feat",                  default='dnn_ZZbbtt_kl_1')
    parser.add_argument("--grp",          dest="grp",                   default='datacard_zz')
    parser.add_argument("--channels",     dest="channels",              default="etau,mutau,tautau")
    parser.add_argument("--num",          dest="num",                   default='',               help='Assign number to output directory for versioning')
    parser.add_argument("--user_eos",     dest="user_eos",              default='evernazz',       help='User Name for lxplus account')
    parser.add_argument("--user_cmt",     dest="user_cmt",              default='vernazza',       help='User Name for cmt folder')
    parser.add_argument("--mass",         dest="mass",                  default='200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,3000')
    makeFlag("--run",                     dest="run",                   default=True,             help='Run commands or do a dry-run')
    makeFlag("--comb_2016",               dest="comb_2016",             default=True,             help='Combine 2016 and 2016_HIPM')
    makeFlag("--run_copy",                dest="run_copy",              default=True,             help='Run copy of datacards')
    makeFlag("--run_one",                 dest="run_one",               default=True,             help='Run each channel or not')
    makeFlag("--run_ch",                  dest="run_ch",                default=True,             help='Combine channels or not')
    makeFlag("--run_cat",                 dest="run_cat",               default=True,             help='Combine categories or not')
    makeFlag("--run_zh_comb_cat",         dest="run_zh_comb_cat",       default=False,            help='Run ZbbHtt & ZttHbb combination for each year separately')
    makeFlag("--run_year",                dest="run_year",              default=True,             help='Combine years or not')
    makeFlag("--run_zh_comb_year",        dest="run_zh_comb_year",      default=False,            help='Run ZbbHtt & ZttHbb combination for Full Run2')
    makeFlag("--run_impacts",             dest="run_impacts",           default=True,             help='Make impact plots')
    makeFlag("--run_impacts_noMCStat",    dest="run_impacts_noMCStat",  default=False,            help='Make impact plots, but only without MC stat uncertainties (faster)')
    makeFlag("--plot_only",               dest="plot_only",             default=False,            help='Skip all combine commands and plot only')
    makeFlag("--move_eos",                dest="move_eos",              default=False,            help='Move results to eos')
    makeFlag("--singleThread",            dest="singleThread",          default=False,            help="Don't run in parallel, disable for debugging")
    makeFlag("--featureDependsOnMass",    dest='featureDependsOnMass',  default=False,            help="Add _$MASS to name of feature for each mass for parametrized DNN")
    options = parser.parse_args()

    if ',' in options.ver:  versions = options.ver.split(',')
    else:                   versions = [options.ver]
    
    if ',' in options.cat:  categories = options.cat.split(',')
    else:                   categories = [options.cat]

    if ',' in options.feat: features = options.feat.split(',')
    else:                   features = [options.feat]

    if ',' in options.channels: channels = options.channels.split(',')
    else:                       channels = [options.channels]

    if ',' in options.mass:     mass_points = options.mass.split(',')
    else:                       mass_points = [options.mass]

    prd = options.prd
    grp = options.grp
    run = options.run
    run_copy = options.run_copy
    run_one = options.run_one
    run_ch = options.run_ch
    run_cat = options.run_cat
    run_year = options.run_year
    featureDependsOnMass = options.featureDependsOnMass
    comb_2016 = options.comb_2016

    cmtdir = '/data_CMS/cms/' + options.user_cmt + '/cmt/CreateDatacards/'
    maindir = os.getcwd() + f'/Res{options.num}/'

    if "ZZ" in options.ver:       o_name = 'ZZbbtt'; process_tex = r"$X\rightarrow ZZ\rightarrow bb\tau\tau$";   x_axis = r"$m_{X}$ [GeV]"
    elif "ZbbHtt" in options.ver: o_name = 'ZbbHtt'; process_tex = r"$Z'\rightarrow ZH$";  x_axis = r"$m_{Z'}$ [GeV]" # \rightarrow bb\tau\tau
    elif "ZttHbb" in options.ver: o_name = 'ZttHbb'; process_tex = r"$Z'\rightarrow ZH$"; x_axis = r"$m_{Z'}$ [GeV]" # \rightarrow \tau\tau bb

    dict_ch_name = {"etau": "$\\tau_{e}\\tau_{h}$", "mutau": "$\\tau_{\\mu}\\tau_{h}$", "tautau": "$\\tau_{h}\\tau_{h}$"}

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    def parseFile(filename, CL='50.0', exp=True):
        f = open(filename)
        matches = []
        for line in f:
            search = 'Expected {}%: r <'.format(CL)
            if not exp:                search = 'Observed Limit: r <'
            if not search in line:     continue
            val = line.replace(search, '')
            matches.append(float(val))
        if len(matches) == 0:
            mes = 'Did not find any expected in file: {}, CL={}, exp?={}'
            print(mes.format(filename, CL, exp))
            return -1.0
        else:
            return matches[-1]

    def SaveResults(odir, mass):
        log_file_name = odir + '/combine.log'
        limis_dict = {
            '{}'.format(mass): {
                'exp'  : parseFile(log_file_name),
                'm1s_t': parseFile(log_file_name, CL='16.0'),
                'p1s_t': parseFile(log_file_name, CL='84.0'),
                'm2s_t': parseFile(log_file_name, CL=' 2.5') ,
                'p2s_t': parseFile(log_file_name, CL='97.5'),
            }
        }
        json_file_name = odir + '/limits.json'
        with open(json_file_name, 'w') as json_file:
            json.dump(limis_dict, json_file, indent=2)

    def GetLimits(limit_file_list, remove_invalid=False):
        """ Read from a list of json files the limit results 
        If remove_invalid=True, will exclude the points with limit=-1 from the result
        """
        mass = []; exp = []; m1s_t = []; p1s_t = []; m2s_t = []; p2s_t = []
        for limit_file in limit_file_list:
            try:
                with open(limit_file, 'r') as json_file:
                    mass_dict = json.load(json_file)
                first_key = list(mass_dict.keys())[0]
                mass.append(float(first_key))
                exp.append(mass_dict[first_key]['exp'])
                m1s_t.append(mass_dict[first_key]['m1s_t'])
                p1s_t.append(mass_dict[first_key]['p1s_t'])
                m2s_t.append(mass_dict[first_key]['m2s_t'])
                p2s_t.append(mass_dict[first_key]['p2s_t'])
            except FileNotFoundError:
                print(f"## INFO : plot skipping non-existent limit file {limit_file}")
        if remove_invalid:
            good_values = np.array(exp) > 0.
            mass = np.array(mass)[good_values].tolist()
            exp = np.array(exp)[good_values].tolist()
            m1s_t = np.array(m1s_t)[good_values].tolist()
            p1s_t = np.array(p1s_t)[good_values].tolist()
            m2s_t = np.array(m2s_t)[good_values].tolist()
            p2s_t = np.array(p2s_t)[good_values].tolist()
        return mass, exp, m1s_t, p1s_t, m2s_t, p2s_t
    
    def CheckLimits(limit_file):
        """ Ensure limits file exists and has a valid limit """
        try:
            with open(limit_file, 'r') as json_file:
                mass_dict = json.load(json_file)
            first_key = list(mass_dict.keys())[0]
            return mass_dict[first_key]["exp"] >= 0
        except FileNotFoundError:
            return False  

    def SetStyle(p2s_t, x_axis, y_axis, version, line1=""):
        if '2018' in version:          text_year = r'2018 - 59.7 fb$^{-1}$ (13 TeV)'
        elif '2017' in version:        text_year = r'2017 - 41.5 fb$^{-1}$ (13 TeV)'
        elif '2016_HIPM' in version:   text_year = r'2016 - 19.5 fb$^{-1}$ (13 TeV)'
        elif '2016' in version:        text_year = r'2016 - 16.8 fb$^{-1}$ (13 TeV)'
        elif 'FullRun2' in version:    text_year = r'137.1 fb$^{-1}$ (13 TeV)'
        plt.text(0.03, 0.97, line1, ha="left", va="top", transform=plt.gca().transAxes, color="black", bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
        mplhep.cms.label(data=False, rlabel=text_year)
        plt.xlabel(x_axis)
        plt.ylabel(r"95% CL on $\sigma \times \mathbf{\it{B}}$(" + y_axis + r") [pb]")
        plt.title("")
        plt.ylim(0.003,2*max(p2s_t))
        plt.grid(True, zorder = 4)
        plt.legend(loc='upper right', fontsize=18, frameon=True)
        plt.yscale('log')
        ax = plt.gca()
        ax.set_axisbelow(False)
    
    ############################################################################
    print("\n ### INFO: Copy all datacards \n")
    ############################################################################

    if run_copy:
        for feature in features:
            for version in versions:
                for category in categories:
                    for channel in channels:
                        for mass in mass_points:
                            odir = maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/M{mass}'
                            if run: run_cmd('mkdir -p ' + odir)
                            if featureDependsOnMass: feat_name = f'{feature}_{mass}'
                            else:                    feat_name = f'{feature}'
                            ch_file = cmtdir + f'/{version}/{category}/{prd}_M{mass}/{feat_name}_{grp}_{channel}_os_iso.txt'
                            ch_root = cmtdir + f'/{version}/{category}/{prd}_M{mass}/{feat_name}_{grp}_{channel}_os_iso.root'
                            run_cmd(f'cp {ch_file} {odir}/{version}_{category}_{feat_name}_{grp}_{channel}_os_iso.txt')
                            run_cmd(f'cp {ch_root} {odir}')

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    # run combination of 2016 and 2016_HIPM

    if comb_2016:

        if any("2016_Z" in s for s in versions) and any("2016_HIPM_Z" in s for s in versions):
        
            prefix = "ul_"
            suffix = "_Z" + versions[0].split("_Z")[1]
            v_2016 = prefix + "2016" + suffix
            v_2016_HIPM = prefix + "2016_HIPM" + suffix
            v_combined = prefix + "2016_ALL" + suffix

            def run_comb_2016(feature, category, channel, mass):
            
                combdir = maindir + f'/{v_combined}/{prd}/{feature}/{category}/{channel}/M{mass}'
                print(" ### INFO: Saving combination in ", combdir)
                if run: run_cmd('mkdir -p ' + combdir)

                if featureDependsOnMass: feat_name = f'{feature}_{mass}'
                else:                    feat_name = f'{feature}'
                cmd = f'combineCards.py'
                for version in [v_2016, v_2016_HIPM]:
                    year_file = maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/M{mass}/{version}_{category}_{feat_name}_{grp}_{channel}_os_iso.txt'
                    year = version.split("ul_")[1].split("_Z")[0]
                    cmd += f' Y{year}={year_file}'
                    cmd += f' > {v_combined}_{category}_{feat_name}_{grp}_{channel}_os_iso.txt'
                if run: os.chdir(combdir)
                if run: run_cmd(cmd)

            if run_one:
                if not options.plot_only:

                    ############################################################################
                    print("\n ### INFO: Run Combination of 2016 and 2016_HIPM \n")
                    ############################################################################

                    if options.singleThread:
                        for feature in features:
                            for category in categories:
                                for channel in channels:
                                    for mass in mass_points:
                                        run_comb_2016(feature, category, channel, mass)
                    else:
                        with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                            futures = []
                            for feature in features:
                                for category in categories:
                                    for channel in channels:
                                        for mass in mass_points:
                                            run_comb_2016(feature, category, channel, mass)
                            for res in concurrent.futures.as_completed(futures):
                                res.result()

            versions.remove(v_2016) 
            versions.remove(v_2016_HIPM)
            versions.insert(0, v_combined)

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    ##########################################################            
    # RUN SINGLE LIMIT
    ########################################################## 

    def run_single_limit(maindir, feature, version, prd, category, mass, channel, featureDependsOnMass):

        if featureDependsOnMass: feat_name = f'{feature}_{mass}'
        else:                    feat_name = f'{feature}'

        odir = maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/M{mass}'
        run_cmd('mkdir -p ' + odir)
        datafile = odir + f'/{version}_{category}_{feat_name}_{grp}_{channel}_os_iso.txt'

        if not os.path.isfile(datafile):
            print("#### WARNING : could not find datacard " + datafile + ", ignoring")
            return

        print(" ### INFO: Run Asymptotic limit")
        cmd = f'cd {odir} && combine -M AsymptoticLimits {datafile} --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: run_cmd(cmd)

        SaveResults(odir, mass)
        
    if run_one:

        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Single Channel\n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                    for version in versions:
                        for category in categories:
                            for mass in mass_points:
                                for channel in channels:
                                    run_single_limit(maindir, feature, version, prd, category, mass, channel, featureDependsOnMass)

            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=50) as exc:
                    futures = []
                    for feature in features:
                        for version in versions:
                            for category in categories:
                                for mass in mass_points:
                                    for channel in channels:
                                        futures.append(exc.submit(run_single_limit, \
                                            maindir, feature, version, prd, category, mass, channel, featureDependsOnMass))
                    for res in concurrent.futures.as_completed(futures):
                        try:
                            res.result()
                        except RuntimeError as e:
                            print("#### ERROR : " + str(e))

        ############################################################################
        print("\n ### INFO: Plot Single Channel \n")
        ############################################################################

        for feature in features:
            for version in versions:
                for category in categories:

                    if "boosted" in category:        cat_name = r"Boosted"
                    elif "resolved_1b" in category:  cat_name = r"Res 1b"
                    elif "resolved_2b" in category:  cat_name = r"Res 2b"
                    else:                            cat_name = category

                    for channel in channels:
                        limit_file_list = [maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/M{mass}/limits.json'
                            for mass in mass_points]
                        mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                        if len(mass) == 0:
                            print(f"## INFO : skipping plot for {maindir}/{version}/{prd}/{feature}/{category}/{channel}/Limits_{ver_short}_{cat_short}_{channel} due to no limits available")
                            continue
                        
                        fig, ax = plt.subplots(figsize=(12,10))
                        plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
                        plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                            color = '#FFDF7Fff', label = "68% expected", zorder=2)
                        plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                            color = '#85D1FBff', label = "95% expected", zorder=1)
                        SetStyle(p2s_t, x_axis, process_tex, version, line1=cat_name+"\n"+dict_ch_name[channel])
                        ver_short = version.split("ul_")[1].split("_Z")[0] ; cat_short = category.split("_cut_90_")[1]
                        plt.savefig(maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/Limits_{ver_short}_{cat_short}_{channel}.pdf')
                        plt.savefig(maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/Limits_{ver_short}_{cat_short}_{channel}.png')
                        # print(maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/Limits_{ver_short}_{cat_short}_{channel}.png')
                        plt.close()

    ##########################################################
    # RUN COMBINATION OF CHANNELS
    ##########################################################

    def run_comb_channels(maindir, cmtdir, feature, version, prd, category, mass, featureDependsOnMass):

        if featureDependsOnMass: feat_name = f'{feature}_{mass}'
        else:                    feat_name = f'{feature}'

        combdir = maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}'
        print(" ### INFO: Saving combination in ", combdir)
        if run: run_cmd('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for ch in channels:
            # FIrst check if we have limits for individual channel
            ch_file = cmtdir + f'/{version}/{category}/{prd}_M{mass}/{feat_name}_{grp}_{ch}_os_iso.txt'
            if CheckLimits(maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/M{mass}/limits.json') and os.path.isfile(ch_file): # check expected
                cmd += f' {ch}={ch_file}'
            else:
                print(f"## WARNING : comb_channels : skipping {version}/{prd}/{category}/{channel}/M{mass}")
        if not "=" in cmd:
            print(f"## WARNING Skipping combination {version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}")
            return
        cmd += f' > {version}_{feature}_{category}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd)

        cmd = f'combine -M AsymptoticLimits {version}_{feature}_{category}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        
        if run: os.chdir(combdir)
        if run: run_cmd(cmd, check=False)

        SaveResults(combdir, mass)
        
    if run_ch:

        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Combination of channels \n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                    for version in versions:
                        for category in categories:
                            for mass in mass_points:
                                run_comb_channels(maindir, cmtdir, feature, version, prd, category, mass, featureDependsOnMass)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for version in versions:
                            for category in categories:
                                for mass in mass_points:
                                    futures.append(exc.submit(run_comb_channels, maindir, cmtdir, feature, version, prd, category, mass, featureDependsOnMass))
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

                    limit_file_list = [maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}/limits.json'
                        for mass in mass_points]
                    mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)

                    if len(mass) == 0 or np.max(p1s_t) <= 0:
                        print(f"## INFO : skipping plot for {maindir}/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}_split due to no limits available")
                        continue
                    
                    fig, ax = plt.subplots(figsize=(12,10))
                    plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
                    plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                        color = '#FFDF7Fff', label = "68% expected", zorder=2)
                    plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                        color = '#85D1FBff', label = "95% expected", zorder=1)
                    SetStyle(p2s_t, x_axis, process_tex, version, line1=cat_name)
                    ver_short = version.split("ul_")[1].split("_Z")[0] ; cat_short = category.split("_cut_90_")[1]
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}.pdf')
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}.png')
                    # print(maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}.png')

                    cmap = plt.get_cmap('tab10')
                    for i, channel in enumerate(channels):
                        limit_file_list = [maindir + f'/{version}/{prd}/{feature}/{category}/{channel}/M{mass}/limits.json'
                            for mass in mass_points]
                        mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                        plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {dict_ch_name[channel]}", zorder=3, color=cmap(i))
                    plt.legend(loc='upper right', fontsize=18, frameon=True)
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}_split.pdf')
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}_split.png')
                    # print(maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}_split.png')
                    plt.close()

    ##########################################################
    # RUN COMBINATION OF CATEGORIES
    ##########################################################

    def run_comb_categories(maindir, cmtdir, feature, version, prd, mass, featureDependsOnMass):

        if featureDependsOnMass:    feat_name = f'{feature}_{mass}'
        else:                       feat_name = f'{feature}'

        combdir = maindir + f'/{version}/{prd}/{feature}/Combination_Cat/M{mass}'
        print(" ### INFO: Saving combination in ", combdir)
        if run: run_cmd('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for category in categories:
            cat_file = maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}/{version}_{feature}_{category}_os_iso.txt'
            if CheckLimits(maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}/limits.json') and os.path.isfile(cat_file): # check expected
                cat_short = category.split("_cut_90_")[1]
                cmd += f' {cat_short}={cat_file}'
            else:
                print(f"## WARNING : comb_categories: skipping {version}/{prd}/{category}/Combination_Ch/M{mass}")
        if not "=" in cmd:
            print(f"## WARNING Skipping combination {version}/{prd}/{feature}/Combination_Cat/M{mass}")
            return
        cmd += f' > {version}_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd)

        cmd = f'combine -M AsymptoticLimits {version}_{feature}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd, check=False)

        SaveResults(combdir, mass)

    if run_cat:

        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Combination of categories \n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                    for version in versions:
                        for mass in mass_points:
                            run_comb_categories(maindir, cmtdir, feature, version, prd, mass, featureDependsOnMass)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for version in versions:
                            for mass in mass_points:
                                futures.append(exc.submit(run_comb_categories, maindir, cmtdir, feature, version, prd, mass, featureDependsOnMass))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

        ############################################################################
        print("\n ### INFO: Plot Combination of categories \n")
        ############################################################################

        for feature in features:
            for version in versions:

                limit_file_list = [maindir + f'/{version}/{prd}/{feature}/Combination_Cat/M{mass}/limits.json'
                    for mass in mass_points]
                mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)

                if len(mass) == 0:
                    print(f"## INFO : skipping plot for {maindir}/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split due to no limits available")
                    continue
                
                fig, ax = plt.subplots(figsize=(12,10))
                plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
                plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                    color = '#FFDF7Fff', label = "68% expected", zorder=2)
                plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                    color = '#85D1FBff', label = "95% expected", zorder=1)
                SetStyle(p2s_t, x_axis, process_tex, version)
                ver_short = version.split("ul_")[1].split("_Z")[0]
                try:
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.pdf')
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.png')
                except ValueError:
                    pass # Math domain error due log scale
                # print(maindir + f'/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.png')

                cmap = plt.get_cmap('tab10')
                for i, category in enumerate(categories):

                    if "boosted" in category:        cat_name = r"Boosted"
                    elif "resolved_1b" in category:  cat_name = r"Res 1b"
                    elif "resolved_2b" in category:  cat_name = r"Res 2b"
                    else:                            cat_name = category

                    limit_file_list = [maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}/limits.json'
                        for mass in mass_points]
                    mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                    plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {cat_name}", zorder=3, color=cmap(i))
                plt.legend(loc='upper right', fontsize=18, frameon=True)
                try:
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.pdf')
                    plt.savefig(maindir + f'/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.png')
                except ValueError:
                    pass # Math domain error due log scale
                # print(maindir + f'/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.png')
                plt.close()

    ##########################################################
    # RUN COMBINATION OF ZbbHtt & ZttHbb for a single year
    ##########################################################

    def run_comb_ZH(maindir, cmtdir, feature, version_ZbbHtt, version_ZttHbb, prd, mass, featureDependsOnMass):
        """ Combination of ZbbHtt & ZttHbb for a single year and a given mass hypothesis """
        version_comb = version_ZbbHtt.replace("ZbbHtt", "ZHComb")
        combdir = maindir + f'/{version_comb}/{prd}/{feature}/M{mass}'
        print(" ### INFO: Saving ZH combination in ", combdir)
        if run: run_cmd('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for version, short_name in [(version_ZbbHtt, "ZbbHtt"), (version_ZttHbb, "ZttHbb")]:
            cat_file = maindir + f'/{version}/{prd}/{feature}/Combination_Cat/M{mass}/{version}_{feature}_{category}_os_iso.txt'
            if not os.path.isfile(cat_file):
                print("## WARNING : ignoring missing datacard " + cat_file + " in combination")
                continue
            cmd += f' {short_name}={cat_file}'
        if not "=" in cmd:
            print(f"## WARNING Skipping combination {version_comb}/{prd}/{feature}/M{mass}")
            return
        cmd += f' > {version_comb}_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd)

        cmd = f'combine -M AsymptoticLimits {version_comb}_{feature}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd, check=False)

        SaveResults(combdir, mass)

    def split_versions_ZH():
        """ Split the input versions into ZbbHtt versions and ZttHbb versions """
        versions_ZbbHtt = []
        versions_ZttHbb = []
        for version in versions:
            if "ZbbHtt" in version:
                ver_ZttHbb = version.replace("ZbbHtt", "ZttHbb")
                if ver_ZttHbb not in versions:
                    raise RuntimeError(f"Whilst running ZH combination : {version} was found but not corresponding {ver_ZttHbb}")
                versions_ZbbHtt.append(version)
                versions_ZttHbb.append(ver_ZttHbb)
        return versions_ZbbHtt, versions_ZttHbb
    
    if options.run_zh_comb_cat:
        versions_ZbbHtt, versions_ZttHbb = split_versions_ZH()

        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Combination ZH of categories \n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                    for ver_ZbbHtt, ver_ZttHbb in zip(versions_ZbbHtt, versions_ZttHbb):
                        for mass in mass_points:
                            run_comb_ZH(maindir, cmtdir, feature, ver_ZbbHtt, ver_ZttHbb, prd, mass, featureDependsOnMass)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for ver_ZbbHtt, ver_ZttHbb in zip(versions_ZbbHtt, versions_ZttHbb):
                            for mass in mass_points:
                                futures.append(exc.submit(run_comb_ZH, maindir, cmtdir, feature, ver_ZbbHtt, ver_ZttHbb, prd, mass, featureDependsOnMass))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

        ############################################################################
        print("\n ### INFO: Plot Combination ZH of categories \n")
        ############################################################################

        for feature in features:
            for version_ZbbHtt, version_ZttHbb in zip(versions_ZbbHtt, versions_ZttHbb):
                version_comb = version_ZbbHtt.replace("ZbbHtt", "ZHComb")
                limit_file_list = [maindir + f'/{version_comb}/{prd}/{feature}/M{mass}/limits.json'
                    for mass in mass_points]
                mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                
                if len(mass) == 0:
                    print(f"## INFO : skipping plot for {maindir}/{version_comb}/{prd}/{feature} due to no limits available")
                    continue

                fig, ax = plt.subplots(figsize=(12,10))
                plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
                plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                    color = '#FFDF7Fff', label = "68% expected", zorder=2)
                plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                    color = '#85D1FBff', label = "95% expected", zorder=1)
                SetStyle(p2s_t, x_axis, process_tex, version_comb)
                ver_short = version_comb.split("ul_")[1].split("_Z")[0]
                try:
                    plt.savefig(maindir + f'/{version_comb}/{prd}/{feature}/Limits_{ver_short}.pdf')
                    plt.savefig(maindir + f'/{version_comb}/{prd}/{feature}/Limits_{ver_short}.png')
                except ValueError:
                    pass # Math domain error due log scale
                # print(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.png')

                cmap = plt.get_cmap('tab10')
                for i, (version, short_name) in enumerate([(version_ZbbHtt, "ZbbHtt"), (version_ZttHbb, "ZttHbb")]):
                    limit_file_list = [maindir + f'/{version}/{prd}/{feature}/Combination_Cat/M{mass}/limits.json'
                        for mass in mass_points]
                    mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                    plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {short_name}", zorder=3, color=cmap(i))
                plt.legend(loc='upper right', fontsize=18, frameon=True)
                try:
                    plt.savefig(maindir + f'/{version_comb}/{prd}/{feature}/Limits_{ver_short}_split.pdf')
                    plt.savefig(maindir + f'/{version_comb}/{prd}/{feature}/Limits_{ver_short}_split.png')
                    plt.close()
                except ValueError:
                    pass # Math domain error due log scale
                # print(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.png')

    ##########################################################
    # RUN COMBINATION OF YEARS
    ##########################################################

    def run_comb_years(maindir, cmtdir, feature, prd, mass, featureDependsOnMass):

        if featureDependsOnMass:    feat_name = f'{feature}_{mass}'
        else:                       feat_name = f'{feature}'

        combdir = maindir + f'/FullRun2_{o_name}/{prd}/{feature}/M{mass}'
        print(" ### INFO: Saving combination in ", combdir)
        if run: run_cmd('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for version in versions:
            ver_file = maindir + f'/{version}/{prd}/{feature}/Combination_Cat/M{mass}/{version}_{feature}_os_iso.txt'
            if CheckLimits(maindir + f'/{version}/{prd}/{feature}/Combination_Cat/M{mass}/limits.json') and os.path.isfile(ver_file): # check expected
                ver_short = version.split("ul_")[1].split("_Z")[0]
                cmd += f' Y{ver_short}={ver_file}'
            else:
                print(f"## WARNING : comb_years: skipping {version}/{prd}/{feature}/Combination_Cat/M{mass}")
            
        if not "=" in cmd:
            print(f"## WARNING Skipping combination FullRun2_{o_name}/{prd}/{feature}/M{mass}")
            return
        cmd += f' > FullRun2_{o_name}_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd)

        cmd = f'combine -M AsymptoticLimits FullRun2_{o_name}_{feature}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd, check=False)

        SaveResults(combdir, mass)

        print(" ### INFO: Produce Full Run 2 Impact Plots")

        cmd = f'text2workspace.py FullRun2_{o_name}_{feature}_os_iso.txt -o model.root &> text2workspace.log'
        if run: os.system(cmd)
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
                    for mass in mass_points:
                        run_comb_years(maindir, cmtdir, feature, prd, mass, featureDependsOnMass)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for mass in mass_points:
                            futures.append(exc.submit(run_comb_years, maindir, cmtdir, feature, prd, mass, featureDependsOnMass))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

        ############################################################################
        print("\n ### INFO: Plot Combination of years \n")
        ############################################################################

        for feature in features:

            limit_file_list = [maindir + f'/FullRun2_{o_name}/{prd}/{feature}/M{mass}/limits.json'
                for mass in mass_points]
            mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)

            if len(mass) == 0:
                print(f"## INFO : skipping plot for {maindir}/FullRun2_{o_name}/{prd}/{feature}/ due to no limits available")
                continue
            
            output_csv_path = maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.csv'
            with open(output_csv_path, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(['mass', 'exp', 'm1s_t', 'p1s_t', 'm2s_t', 'p2s_t'])
                for i in range(len(mass)):
                    csvwriter.writerow([mass[i], exp[i], m1s_t[i], p1s_t[i], m2s_t[i], p2s_t[i]])

            fig, ax = plt.subplots(figsize=(12,10))
            plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
            plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                color = '#FFDF7Fff', label = "68% expected", zorder=2)
            plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                color = '#85D1FBff', label = "95% expected", zorder=1)
            SetStyle(p2s_t, x_axis, process_tex, "FullRun2")
            plt.savefig(maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.pdf')
            plt.savefig(maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.png')
            # print(maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.png')

            cmap = plt.get_cmap('tab10')
            for i, version in enumerate(versions):
                ver_short = version.split("ul_")[1].split("_Z")[0]
                limit_file_list = [maindir + f'/{version}/{prd}/{feature}/Combination_Cat/M{mass}/limits.json'
                    for mass in mass_points]
                mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {ver_short}", zorder=3, color=cmap(i))
            plt.legend(loc='upper right', fontsize=18, frameon=True)
            plt.savefig(maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}_split.pdf')
            plt.savefig(maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}_split.png')
            # print(maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}_split.png')
            plt.close()

    ##########################################################
    # RUN COMBINATION OF ZbbHtt & ZttHbb for full run2
    ##########################################################

    def run_comb_years_ZH(maindir, cmtdir, feature, prd, mass, featureDependsOnMass):
        """ Combination of ZbbHtt & ZttHbb for a single year and a given mass hypothesis """
        versions_ZbbHtt, versions_ZttHbb = split_versions_ZH()
        combdir = maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/M{mass}'
        print(" ### INFO: Saving ZH combination in ", combdir)
        if run: run_cmd('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for version in versions_ZbbHtt:
            year = version.split("ul_")[1].split("_Z")[0]
            cat_file = maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/M{mass}/{version}_{feature}_os_iso.txt'
            if not os.path.isfile(cat_file):
                print("## WARNING : ignoring missing datacard " + cat_file + " in combination")
                continue
            cmd += f' Y{year}_ZbbHtt={cat_file}'
        for version in versions_ZttHbb:
            year = version.split("ul_")[1].split("_Z")[0]
            cat_file = maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/M{mass}/{version}_{feature}_os_iso.txt'
            if not os.path.isfile(cat_file):
                print("## WARNING : ignoring missing datacard " + cat_file + " in combination")
                continue
            cmd += f' Y{year}_ZttHbb={cat_file}'
        if not "=" in cmd:
            print(f"## WARNING Skipping combination FullRun2_ZHComb/{prd}/{feature}/M{mass}")
            return
        cmd += f' > FullRun2_ZHComb_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd)

        cmd = f'combine -M AsymptoticLimits FullRun2_ZHComb_{feature}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.chdir(combdir)
        if run: run_cmd(cmd, check=False)

        SaveResults(combdir, mass)

    if options.run_zh_comb_year:
        if not options.plot_only:

            ############################################################################
            print("\n ### INFO: Run Combination ZH of years \n")
            ############################################################################

            if options.singleThread:
                for feature in features:
                        for mass in mass_points:
                            run_comb_years_ZH(maindir, cmtdir, feature, prd, mass, featureDependsOnMass)
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for mass in mass_points:
                            futures.append(exc.submit(run_comb_years_ZH, maindir, cmtdir, feature, prd, mass, featureDependsOnMass))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

        ############################################################################
        print("\n ### INFO: Plot Combination ZH of years \n")
        ############################################################################

        for feature in features:
            limit_file_list = [maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/M{mass}/limits.json'
                for mass in mass_points]
            mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list, remove_invalid=True)
            
            if len(mass) == 0:
                print(f"## INFO : skipping plot for /Res/FullRun2_ZHComb/{prd}/{feature} due to no limits available")
                continue

            fig, ax = plt.subplots(figsize=(12,10))
            plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
            plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                color = '#FFDF7Fff', label = "68% expected", zorder=2)
            plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                color = '#85D1FBff', label = "95% expected", zorder=1)
            SetStyle(p2s_t, x_axis, process_tex, "FullRun2 ZHComb")
            try:
                plt.savefig(maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/Limits_FullRun2_ZHComb.pdf')
                plt.savefig(maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/Limits_FullRun2_ZHComb.png')
            except ValueError:
                pass # Math domain error due log scale
            # print(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.png')

            # Overlaying with split in years
            cmap = plt.get_cmap('tab10')
            versions_ZbbHtt, versions_ZttHbb = split_versions_ZH()
            for i, (version_ZbbHtt, version_ZttHbb) in enumerate(zip(versions_ZbbHtt, versions_ZttHbb)):
                version_comb = version_ZbbHtt.replace("ZbbHtt", "ZHComb")
                short_name = version_ZbbHtt.split("ul_")[1].split("_Z")[0]
                limit_file_list = [maindir + f'/Res/{version_comb}/{prd}/{feature}/M{mass}/limits.json'
                    for mass in mass_points]
                mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {short_name}", zorder=3, color=cmap(i))
            plt.legend(loc='upper right', fontsize=18, frameon=True)
            try:
                plt.savefig(maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/Limits_FullRun2_ZHComb_split.pdf')
                plt.savefig(maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/Limits_FullRun2_ZHComb_split.png')
            except ValueError:
                pass # Math domain error due log scale
            # print(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.png')

            # Overlaying with ZbbHtt / ZttHbb
            fig, ax = plt.subplots(figsize=(12,10))
            plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
            plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                color = '#FFDF7Fff', label = "68% expected", zorder=2)
            plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                color = '#85D1FBff', label = "95% expected", zorder=1)
            SetStyle(p2s_t, x_axis, process_tex, "FullRun2 ZHComb")
            versions_ZbbHtt, versions_ZttHbb = split_versions_ZH()
            for i, (proc, short_name) in enumerate(zip(["ZbbHtt", "ZttHbb"], [r"$Z_{bb}H_{\tau\tau}$", r"$Z_{\tau\tau}H_{bb}$"])):
                limit_file_list = [maindir + f'/Res/FullRun2_{proc}/{prd}/{feature}/M{mass}/limits.json'
                    for mass in mass_points]
                mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {short_name}", zorder=3, color=cmap(i))
            plt.legend(loc='upper right', fontsize=18, frameon=True)
            try:
                plt.savefig(maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/Limits_FullRun2_ZHComb_split_proc.pdf')
                plt.savefig(maindir + f'/Res/FullRun2_ZHComb/{prd}/{feature}/Limits_FullRun2_ZHComb_split_proc.png')
            except ValueError:
                pass # Math domain error due log scale
    
    if options.move_eos:

        eos_dir = f'/eos/user/e/evernazz/www/ZZbbtautau/B2GPlots/2024_06_14/{o_name}/Limits/Res'
        user = options.user_eos
        print(f" ### INFO: Copy results to {user}@lxplus.cern.ch")
        print(f"           Inside directory {eos_dir}\n")

        # [FIXME] Work-around for mkdir on eos
        run_cmd(f'mkdir -p TMP_RESULTS_RES && cp index.php TMP_RESULTS_RES')
        for feature in features:
            run_cmd(f'cp ' + maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_*.p* TMP_RESULTS_RES')
            run_cmd(f'cp ' + maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Limits_*.csv TMP_RESULTS_RES')
            # run_cmd(f'cp ' + maindir + f'/FullRun2_{o_name}/{prd}/{feature}/*_os_iso.txt TMP_RESULTS_RES')
            # run_cmd(f'cp ' + maindir + f'/FullRun2_{o_name}/{prd}/{feature}/Impacts* TMP_RESULTS_RES')
            for version in versions:
                ver_short = version.split("ul_")[1].split("_Z")[0]
                run_cmd(f'mkdir -p TMP_RESULTS_RES/{ver_short} && cp index.php TMP_RESULTS_RES/{ver_short}')
                run_cmd(f'cp ' + maindir + f'/{version}/{prd}/{feature}/Combination_Cat/Limits_*.p* TMP_RESULTS_RES/{ver_short}')
                for category in categories:
                    run_cmd(f'cp ' + maindir + f'/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_*.p* TMP_RESULTS_RES/{ver_short}')        
        run_cmd(f'rsync -rltv TMP_RESULTS_RES/* {user}@lxplus.cern.ch:{eos_dir}')
        run_cmd(f'rm -r TMP_RESULTS_RES')
