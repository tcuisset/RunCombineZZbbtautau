import os,json,concurrent.futures,itertools
import pdb, csv
import numpy as np

import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.CMS)

import warnings, logging
warnings.filterwarnings("ignore", message=".*Type 3 font.*")
logging.getLogger('matplotlib').setLevel(logging.ERROR)


'''
ulimit -s unlimited

python3 RunAsymptoticLimits.py --ver ul_2016_ZZ_v12,ul_2016_HIPM_ZZ_v12,ul_2017_ZZ_v12,ul_2018_ZZ_v12 \
    --cat cat_ZZ_elliptical_cut_90_resolved_1b,cat_ZZ_elliptical_cut_90_resolved_2b,cat_ZZ_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZZbbtt_kl_1 --featureDependsOnMass --prd prod_240528 --grp datacard_zz_res --move_eos --user_eos evernazz \
    --mass 200,210,220,230,240,250,260,280,300,320,350,360,400,450,500,550,600,650,700,750,800,850,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2500,2600,2800,3000,3500,4000,4500,5000

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
    parser.add_option("--mass",         dest="mass",         default='200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,3000')
    parser.add_option("--singleThread", action="store_true", help="Don't run in parallel, disable for debugging")
    parser.add_option("--run_one",      dest="run_one",      default=True,             help='Run each channel or not')
    parser.add_option("--run_ch",       dest="run_ch",       default=True,             help='Combine channels or not')
    parser.add_option("--run_cat",      dest="run_cat",      default=True,             help='Combine categories or not')
    parser.add_option("--run_year",     dest="run_year",     default=True,             help='Combine years or not')
    parser.add_option("--plot_only",    dest="plot_only",    default=False,            action='store_true')
    parser.add_option("--move_eos",     dest="move_eos",     default=False,            action='store_true')
    parser.add_option("--user_eos",     dest="user_eos",     default='evernazz',       help='User Name for lxplus account')
    parser.add_option("--featureDependsOnMass",              default=False,            help="Add _$MASS to name of feature for each mass for parametrized DNN",     action="store_true")
    (options, args) = parser.parse_args()

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
    run = int(options.run) == 1
    run_one = options.run_one == True
    run_ch = options.run_ch == True
    run_cat = options.run_cat == True
    run_year = options.run_year == True
    featureDependsOnMass = options.featureDependsOnMass

    if os.environ["USER"] == 'evernazza':
        cmtdir = '/data_CMS/cms/' + os.environ["USER"][1:] + '/cmt/CreateDatacards/'
    else:
        cmtdir = '/data_CMS/cms/' + os.environ["USER"] + '/cmt/CreateDatacards/'

    maindir = os.getcwd() 

    if "ZZ" in options.ver:       o_name = 'ZZbbtt'; process_tex = r"$X\rightarrow ZZ\rightarrow bb\tau\tau$";   x_axis = r"$m_{X}$ [GeV]"
    elif "ZbbHtt" in options.ver: o_name = 'ZbbHtt'; process_tex = r"$Z'\rightarrow ZH\rightarrow bb\tau\tau$";  x_axis = r"$m_{Z'}$ [GeV]"
    elif "ZttHbb" in options.ver: o_name = 'ZttHbb'; process_tex = r"$Z'\rightarrow ZH\rightarrow \tau\tau bb$"; x_axis = r"$m_{Z'}$ [GeV]"

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

    def GetLimits(limit_file_list):
        mass = []; exp = []; m1s_t = []; p1s_t = []; m2s_t = []; p2s_t = []
        for limit_file in limit_file_list:
            with open(limit_file, 'r') as json_file:
                mass_dict = json.load(json_file)
            first_key = list(mass_dict.keys())[0]
            mass.append(float(first_key))
            exp.append(mass_dict[first_key]['exp'])
            m1s_t.append(mass_dict[first_key]['m1s_t'])
            p1s_t.append(mass_dict[first_key]['p1s_t'])
            m2s_t.append(mass_dict[first_key]['m2s_t'])
            p2s_t.append(mass_dict[first_key]['p2s_t'])
        return mass, exp, m1s_t, p1s_t, m2s_t, p2s_t
        
    def SetStyle(p2s_t, x_axis, y_axis, version, line1=""):
        if '2018' in version:          text_year = r'2018 - 59.7 fb$^{-1}$ (13 TeV)'
        elif '2017' in version:        text_year = r'2017 - 41.5 fb$^{-1}$ (13 TeV)'
        elif '2016' in version:        text_year = r'2016 - 16.8 fb$^{-1}$ (13 TeV)'
        elif '2016_HIPM' in version:   text_year = r'2016 - 19.5 fb$^{-1}$ (13 TeV)'
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

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    ##########################################################            
    # RUN SINGLE LIMIT
    ########################################################## 

    def run_single_limit(maindir, cmtdir, feature, version, prd, category, mass, channel, featureDependsOnMass):

        if featureDependsOnMass: feat_name = f'{feature}_{mass}'
        else:                    feat_name = f'{feature}'

        odir = maindir + f'/Res/{version}/{prd}/{feature}/{category}/{channel}/M{mass}'        
        os.system('mkdir -p ' + odir)
        datadir = cmtdir + f'/{version}/{category}/{prd}_M{mass}'
        datafile = datadir + f'/{feat_name}_{grp}_{channel}_os_iso.txt'

        print(" ### INFO: Run Asymptotic limit")
        cmd = f'cd {odir} && combine -M AsymptoticLimits {datafile} --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.system(cmd)

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
                                    run_single_limit(maindir, cmtdir, feature, version, prd, category, mass, channel, featureDependsOnMass)

            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=15) as exc:
                    futures = []
                    for feature in features:
                        for version in versions:
                            for category in categories:
                                for mass in mass_points:
                                    for channel in channels:
                                        futures.append(exc.submit(run_single_limit, \
                                            maindir, cmtdir, feature, version, prd, category, mass, channel, featureDependsOnMass))
                    for res in concurrent.futures.as_completed(futures):
                        res.result()

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
                        limit_file_list = [maindir + f'/Res/{version}/{prd}/{feature}/{category}/{channel}/M{mass}/limits.json'
                            for mass in mass_points]
                        mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                        
                        fig, ax = plt.subplots(figsize=(12,10))
                        plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
                        plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                            color = '#FFDF7Fff', label = "68% expected", zorder=2)
                        plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                            color = '#85D1FBff', label = "95% expected", zorder=1)
                        SetStyle(p2s_t, x_axis, process_tex, version, line1=cat_name+"\n"+dict_ch_name[channel])
                        ver_short = version.split("ul_")[1].split("_Z")[0] ; cat_short = category.split("_cut_90_")[1]
                        plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/{category}/{channel}/Limits_{ver_short}_{cat_short}_{channel}.pdf')
                        plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/{category}/{channel}/Limits_{ver_short}_{cat_short}_{channel}.png')
                        # print(maindir + f'/Res/{version}/{prd}/{feature}/{category}/{channel}/Limits_{ver_short}_{cat_short}_{channel}.png')
                        plt.close()

    ##########################################################
    # RUN COMBINATION OF CHANNELS
    ##########################################################

    def run_comb_channels(maindir, cmtdir, feature, version, prd, category, mass, featureDependsOnMass):

        if featureDependsOnMass: feat_name = f'{feature}_{mass}'
        else:                    feat_name = f'{feature}'

        combdir = maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}'
        print(" ### INFO: Saving combination in ", combdir)
        if run: os.system('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for ch in channels:
            ch_file = cmtdir + f'/{version}/{category}/{prd}_M{mass}/{feat_name}_{grp}_{ch}_os_iso.txt'
            cmd += f' {ch}={ch_file}'
        cmd += f' > {version}_{feature}_{category}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

        cmd = f'combine -M AsymptoticLimits {version}_{feature}_{category}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

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

                    limit_file_list = [maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}/limits.json'
                        for mass in mass_points]
                    mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                    
                    fig, ax = plt.subplots(figsize=(12,10))
                    plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
                    plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                        color = '#FFDF7Fff', label = "68% expected", zorder=2)
                    plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                        color = '#85D1FBff', label = "95% expected", zorder=1)
                    SetStyle(p2s_t, x_axis, process_tex, version, line1=cat_name)
                    ver_short = version.split("ul_")[1].split("_Z")[0] ; cat_short = category.split("_cut_90_")[1]
                    plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}.pdf')
                    plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}.png')
                    # print(maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}.png')

                    cmap = plt.get_cmap('tab10')
                    for i, channel in enumerate(channels):
                        limit_file_list = [maindir + f'/Res/{version}/{prd}/{feature}/{category}/{channel}/M{mass}/limits.json'
                            for mass in mass_points]
                        mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                        plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {dict_ch_name[channel]}", zorder=3, color=cmap(i))
                    plt.legend(loc='upper right', fontsize=18, frameon=True)
                    plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}_split.pdf')
                    plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}_split.png')
                    # print(maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_{ver_short}_{cat_short}_split.png')

    ##########################################################
    # RUN COMBINATION OF CATEGORIES
    ##########################################################

    def run_comb_categories(maindir, cmtdir, feature, version, prd, mass, featureDependsOnMass):

        if featureDependsOnMass:    feat_name = f'{feature}_{mass}'
        else:                       feat_name = f'{feature}'

        combdir = maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/M{mass}'
        print(" ### INFO: Saving combination in ", combdir)
        if run: os.system('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for category in categories:
            cat_file = maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}/{version}_{feature}_{category}_os_iso.txt'
            cat_short = category.split("_cut_90_")[1]
            cmd += f' {cat_short}={cat_file}'
        cmd += f' > {version}_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

        cmd = f'combine -M AsymptoticLimits {version}_{feature}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

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

                limit_file_list = [maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/M{mass}/limits.json'
                    for mass in mass_points]
                mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                
                fig, ax = plt.subplots(figsize=(12,10))
                plt.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=3)
                plt.fill_between(np.asarray(mass), np.asarray(p1s_t), np.asarray(m1s_t), 
                    color = '#FFDF7Fff', label = "68% expected", zorder=2)
                plt.fill_between(np.asarray(mass), np.asarray(p2s_t), np.asarray(m2s_t), 
                    color = '#85D1FBff', label = "95% expected", zorder=1)
                SetStyle(p2s_t, x_axis, process_tex, version)
                ver_short = version.split("ul_")[1].split("_Z")[0]
                plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.pdf')
                plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.png')
                # print(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}.png')

                cmap = plt.get_cmap('tab10')
                for i, category in enumerate(categories):

                    if "boosted" in category:        cat_name = r"Boosted"
                    elif "resolved_1b" in category:  cat_name = r"Res 1b"
                    elif "resolved_2b" in category:  cat_name = r"Res 2b"
                    else:                            cat_name = category

                    limit_file_list = [maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/M{mass}/limits.json'
                        for mass in mass_points]
                    mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                    plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {cat_name}", zorder=3, color=cmap(i))
                plt.legend(loc='upper right', fontsize=18, frameon=True)
                plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.pdf')
                plt.savefig(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.png')
                # print(maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_{ver_short}_split.png')

    ##########################################################
    # RUN COMBINATION OF YEARS
    ##########################################################

    def run_comb_years(maindir, cmtdir, feature, prd, mass, featureDependsOnMass):

        if featureDependsOnMass:    feat_name = f'{feature}_{mass}'
        else:                       feat_name = f'{feature}'

        combdir = maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/M{mass}'
        print(" ### INFO: Saving combination in ", combdir)
        if run: os.system('mkdir -p ' + combdir)

        cmd = f'combineCards.py'
        for version in versions:
            ver_file = maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/M{mass}/{version}_{feature}_os_iso.txt'
            ver_short = version.split("ul_")[1].split("_Z")[0]
            cmd += f' Y{ver_short}={ver_file}'
        cmd += f' > FullRun2_{o_name}_{feature}_os_iso.txt'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

        cmd = f'combine -M AsymptoticLimits FullRun2_{o_name}_{feature}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        if run: os.chdir(combdir)
        if run: os.system(cmd)

        SaveResults(combdir, mass)

        # cmd = f'combineTool.py -M Impacts -d FullRun2_{o_name}_{feature}_os_iso.txt -m 125 --run blind --noFitAsimov {comb_options}'
        # [FIXME] Missing impact plots
                
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

            limit_file_list = [maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/M{mass}/limits.json'
                for mass in mass_points]
            mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
            
            output_csv_path = maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.csv'
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
            plt.savefig(maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.pdf')
            plt.savefig(maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.png')
            # print(maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}.png')

            cmap = plt.get_cmap('tab10')
            for i, version in enumerate(versions):
                ver_short = version.split("ul_")[1].split("_Z")[0]
                limit_file_list = [maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/M{mass}/limits.json'
                    for mass in mass_points]
                mass, exp, m1s_t, p1s_t, m2s_t, p2s_t = GetLimits(limit_file_list)
                plt.plot(mass, exp, marker='o', linestyle='--', label = f"Expected {ver_short}", zorder=3, color=cmap(i))
            plt.legend(loc='upper right', fontsize=18, frameon=True)
            plt.savefig(maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}_split.pdf')
            plt.savefig(maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}_split.png')
            # print(maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_FullRun2_{o_name}_split.png')

    if options.move_eos:

        eos_dir = f'/eos/user/e/evernazz/www/ZZbbtautau/B2GPlots/2024_06_14/{o_name}/Limits/Res'
        user = options.user_eos
        print(f" ### INFO: Copy results to {user}@lxplus.cern.ch")
        print(f"           Inside directory {eos_dir}\n")

        # [FIXME] Work-around for mkdir on eos
        os.system(f'mkdir -p TMP_RESULTS_RES && cp index.php TMP_RESULTS_RES')
        for feature in features:
            os.system(f'cp ' + maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_*.p* TMP_RESULTS_RES')
            os.system(f'cp ' + maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Limits_*.csv TMP_RESULTS_RES')
            # os.system(f'cp ' + maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/*_os_iso.txt TMP_RESULTS_RES')
            # os.system(f'cp ' + maindir + f'/Res/FullRun2_{o_name}/{prd}/{feature}/Impacts* TMP_RESULTS_RES')
            for version in versions:
                ver_short = version.split("ul_")[1].split("_Z")[0]
                os.system(f'mkdir -p TMP_RESULTS_RES/{ver_short} && cp index.php TMP_RESULTS_RES/{ver_short}')
                os.system(f'cp ' + maindir + f'/Res/{version}/{prd}/{feature}/Combination_Cat/Limits_*.p* TMP_RESULTS_RES/{ver_short}')
                for category in categories:
                    os.system(f'cp ' + maindir + f'/Res/{version}/{prd}/{feature}/{category}/Combination_Ch/Limits_*.p* TMP_RESULTS_RES/{ver_short}')        
        os.system(f'rsync -rltv TMP_RESULTS_RES/* {user}@lxplus.cern.ch:{eos_dir}')
        os.system(f'rm -r TMP_RESULTS_RES')
