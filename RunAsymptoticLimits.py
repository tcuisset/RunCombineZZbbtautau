import os,json,concurrent.futures,itertools
import pdb
import numpy as np

import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.CMS)

'''
python3 RunAsymptoticLimits.py --cfg ul_2016_HIPM_ZZ_v12,ul_2016_ZZ_v12,ul_2017_ZZ_v12,ul_2018_ZZ_v12 \
    --feat dnn_ZZbbtt_kl_1 --featureDependsOnMass --prd prod_240513 \
    --mass 200,210,220,230,240,250,260,280,300,320,350,360,400,450,500,550,600,650,700,750,800,850,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2500,2600,2800,3000,3500,4000,4500,5000 \
    --cat ZZ_elliptical_cut_90 --usr vernazza --grp datacard_zz_res

python3 RunAsymptoticLimits.py --cfg ul_2016_HIPM_ZbbHtt_v12,ul_2016_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --feat dnn_ZbbHtt_kl_1 --featureDependsOnMass --prd ... \
    --mass ... \
    --cat ZZ_elliptical_cut_90 --grp datacard_zz_res

python3 RunAsymptoticLimits.py --cfg ul_2016_HIPM_ZttHbb_v12,ul_2016_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --feat dnn_ZttHbb_kl_1 --featureDependsOnMass --prd ... \
    --mass ... \
    --cat ZZ_elliptical_cut_90 --grp datacard_zz_res
'''


# comb_options = '--minimizerAlgo Minuit2'
# comb_options = '--cminDefaultMinimizerType Minuit2'
comb_options = ''

def parseFile(filename, CL='50.0', exp=True):
    f = open(filename)
    matches = []
    for line in f:
        search = 'Expected {}%: r <'.format(CL)
        if not exp:
            search = 'Observed Limit: r <'
        if not search in line:
            continue
        val = line.replace(search, '')
        val = float(val)
        matches.append(val)
    if len(matches) == 0:
        mes = 'Did not find any expected in file: {}, CL={}, exp?={}'
        print(mes.format(filename, CL, exp))
        return -1.0
    else:
        return matches[-1]

def run_single_limit(cfg_dir, cfg, cat, grp, prd, ch, feat_ver, mass, run, usr):

    if usr: basedir = '/data_CMS/cms/' + usr + '/cmt/CreateDatacards/'
    else:   basedir = '/data_CMS/cms/' + os.environ["USER"] + '/cmt/CreateDatacards/'
    datadir = basedir + f'/{cfg}/cat_{cat}/{prd}_M{mass}/' # here change depending on whether catgeory has channel inside it or not
    datafile = datadir + f'/{feat_ver}_{grp}_{ch}_os_iso.txt'
    rootfile = datadir + f'/{feat_ver}_{grp}_{ch}_os_iso.root'

    # Define output directory
    odir = cfg_dir + f'/{prd}/{ch}/{feat_ver}'
    print(" ### INFO: Saving output in ", odir)
    os.system('mkdir -p ' + odir)

    # Get datacards
    os.system(f'cp {datafile} {odir}/{feat_ver}_{ch}_os_iso.txt')
    os.system(f'cp {rootfile} {odir}/{feat_ver}_{ch}_os_iso.root')

    if run:
        # Run asymptotic limit
        print(" ### INFO: Run Asymptotic limit")
        cmd = f'combine -M AsymptoticLimits {feat_ver}_{ch}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        os.chdir(odir)
        os.system(cmd)

    # Save results
    log_file_name = odir + '/combine.log'
    exp   = parseFile(log_file_name)
    m1s_t = parseFile(log_file_name, CL='16.0')
    p1s_t = parseFile(log_file_name, CL='84.0')
    m2s_t = parseFile(log_file_name, CL=' 2.5') 
    p2s_t = parseFile(log_file_name, CL='97.5')

    limis_dict = {
        '{}'.format(mass): {
            'exp'  : exp,
            'm1s_t': m1s_t,
            'p1s_t': p1s_t,
            'm2s_t': m2s_t,
            'p2s_t': p2s_t,
        }
    }

    json_file_name = odir + '/limits.json'
    with open(json_file_name, 'w') as json_file:
        json.dump(limis_dict, json_file, indent=2)

def run_comb_limit_singleYear(cfg_dir, feat_ver, prd, mass, run):

    # Define output directory
    odir = cfg_dir + f'/{prd}/combination/{feat_ver}'
    print(" ### INFO: Saving output in ", odir)
    os.system('mkdir -p ' + odir)

    # Get datacards
    for ch in ['etau', 'mutau', 'tautau']:
        datafile = cfg_dir + f'/{prd}/{ch}/{feat_ver}' + f'/{feat_ver}_{ch}_os_iso.txt'
        rootfile = cfg_dir + f'/{prd}/{ch}/{feat_ver}' + f'/{feat_ver}_{ch}_os_iso.root'
        os.system('cp '+datafile+' '+rootfile+' '+odir)
        os.system('cp '+datafile+' '+rootfile+' '+odir)

    # Combine datacards
    cmd = f'combineCards.py ch_etau={feat_ver}_etau_os_iso.txt ch_mutau={feat_ver}_mutau_os_iso.txt ch_tautau={feat_ver}_tautau_os_iso.txt > {feat_ver}_comb_os_iso.txt'
    os.chdir(odir)
    os.system(cmd)

    if run:
        # Run asymptotic limit
        print(" ### INFO: Run Asymptotic limit")
        cmd = f'combine -M AsymptoticLimits {feat_ver}_comb_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        os.chdir(odir)
        os.system(cmd)

    log_file_name = odir + '/combine.log'
    exp   = parseFile(log_file_name)
    m1s_t = parseFile(log_file_name, CL='16.0')
    p1s_t = parseFile(log_file_name, CL='84.0')
    m2s_t = parseFile(log_file_name, CL=' 2.5') 
    p2s_t = parseFile(log_file_name, CL='97.5')

    limis_dict = {
        '{}'.format(mass): {
            'exp'  : exp,
            'm1s_t': m1s_t,
            'p1s_t': p1s_t,
            'm2s_t': m2s_t,
            'p2s_t': p2s_t,
        }
    }

    json_file_name = odir + '/limits.json'
    with open(json_file_name, 'w') as json_file:
        json.dump(limis_dict, json_file, indent=2)  

def run_comb_limit_allYears(base_dir, prefix, feat_ver, prd, mass, config, run ):

    # Define output directory
    odir = base_dir + f'/{prefix}_FullRun2/{prd}/{feat_ver}'
    print(" ### INFO: Saving output in ", odir)
    os.system('mkdir -p ' + odir)

    # Get datacards
    cmd = 'combineCards.py'
    for cfg in config:
        year = cfg.split('_')[1] if "HIPM" not in cfg else "2016_HIPM"
        for ch in ['etau', 'mutau', 'tautau']:
            datafile = base_dir + f'/{cfg}/{prd}/{ch}/{feat_ver}' + f'/{feat_ver}_{ch}_os_iso.txt'
            rootfile = base_dir + f'/{cfg}/{prd}/{ch}/{feat_ver}' + f'/{feat_ver}_{ch}_os_iso.root'

            os.system('cp '+datafile+' '+odir+f"/{feat_ver}_{year}_{ch}_os_iso.txt")
            os.system('cp '+rootfile+' '+odir+f"/{feat_ver}_{year}_{ch}_os_iso.root")
            cmd += f' Year{year}_{ch}={feat_ver}_{year}_{ch}_os_iso.txt'

    # Combine datacards
    cmd += f' > FullRun2_{feat_ver}_os_iso.txt'
    os.chdir(odir)
    os.system(cmd)

    if run:
        # Run asymptotic limit
        print(" ### INFO: Run Asymptotic limit")
        cmd = f'combine -M AsymptoticLimits FullRun2_{feat_ver}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
        os.chdir(odir)
        os.system(cmd)

    log_file_name = odir + '/combine.log'
    exp   = parseFile(log_file_name)
    m1s_t = parseFile(log_file_name, CL='16.0')
    p1s_t = parseFile(log_file_name, CL='84.0')
    m2s_t = parseFile(log_file_name, CL=' 2.5') 
    p2s_t = parseFile(log_file_name, CL='97.5')

    limis_dict = {
        '{}'.format(mass): {
            'exp'  : exp,
            'm1s_t': m1s_t,
            'p1s_t': p1s_t,
            'm2s_t': m2s_t,
            'p2s_t': p2s_t,
        }
    }

    json_file_name = odir + '/limits.json'
    with open(json_file_name, 'w') as json_file:
        json.dump(limis_dict, json_file, indent=2)


def PlotAsymptoticLimits(base_dir, mass_points, prd, cfg, feat, prefix):

    mass_dict = []
    
    mass  = []
    exp   = []
    m1s_t = []
    p1s_t = []
    m2s_t = []
    p2s_t = []

    if cfg == 'FullRun2':
        folder_name = prefix + f'_FullRun2/{prd}/'
        out_folder = folder_name
    else:
        folder_name = f'{cfg}/{prd}/combination'
        out_folder = f'{cfg}/{prd}'

    if '2018' in cfg:          text_year = r'2018 - 59.7 fb$^{-1}$ (13 TeV)'
    elif '2017' in cfg:        text_year = r'2017 - 41.5 fb$^{-1}$ (13 TeV)'
    elif '2016' in cfg:        text_year = r'2016 PostVFP - 16.8 fb$^{-1}$ (13 TeV)'
    elif '2016_HIPM' in cfg:   text_year = r'2016 PreVFP - 19.5 fb$^{-1}$ (13 TeV)'
    elif 'FullRun2' in cfg: text_year = r'137.1 fb$^{-1}$ (13 TeV)'
    
    for mass_value in mass_points:
        limit_file = base_dir + f'/{folder_name}/{feat}_{mass_value}/limits.json'
        with open(limit_file, 'r') as json_file:
            mass_dict = json.load(json_file)
        m = list(mass_dict.keys())[0]
        first_key = list(mass_dict.keys())[0]
        mass.append(float(m))
        exp.append(mass_dict[first_key]['exp'])
        m1s_t.append(mass_dict[first_key]['m1s_t'])
        p1s_t.append(mass_dict[first_key]['p1s_t'])
        m2s_t.append(mass_dict[first_key]['m2s_t'])
        p2s_t.append(mass_dict[first_key]['p2s_t'])

    # Setup matplotlib figure
    fig, ax = plt.subplots(figsize=(12,10))

    ax.plot(mass, exp, color='k', marker='o', label = "Expected", zorder=10)
    ax.fill_between(np.asarray(mass), 
                    np.asarray(p1s_t), 
                    np.asarray(m1s_t), color = '#FFDF7Fff', label = "68% expected", zorder=3)
    ax.fill_between(np.asarray(mass), 
                    np.asarray(p2s_t), 
                    np.asarray(m2s_t), color = '#85D1FBff', label = "95% expected")

    mplhep.cms.label(data=False, rlabel=text_year, fontsize=22)

    if "ZZ" in prefix:
        process_tex = r"$X\rightarrow ZZ\rightarrow bb\tau\tau$"
        ax.set_xlabel(r"$m_{X}$ [GeV]")
    elif "ZbbHtt" in prefix:
        process_tex = r"$Z'\rightarrow ZH\rightarrow bb\tau\tau$"
        ax.set_xlabel(r"$m_{Z'}$ [GeV]")
    elif "ZttHbb" in prefix:
        process_tex = r"$Z'\rightarrow ZH\rightarrow \tau\tau bb$"
        ax.set_xlabel(r"$m_{Z'}$ [GeV]")
    else:
        process_tex = "unknown process"
    ax.set_ylabel(r"95% CL on $\sigma \times \mathbf{\it{B}}$(" + process_tex + r") [pb]")
    
    # Style
    ax.legend()
    ax.set_yscale('log')
    ax.grid()

    fig.savefig(base_dir + f'/{out_folder}/Limits_'+feat+'_'+prd+'_'+cfg+'.png')
    fig.savefig(base_dir + f'/{out_folder}/Limits_'+feat+'_'+prd+'_'+cfg+'.pdf')
    print(base_dir + f'/{out_folder}/Limits_'+feat+'_'+prd+'_'+cfg+'.png')


if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run",     dest="run",      default=True,             help='Run or not')
    parser.add_option("--cfg",     dest="cfg",      default='',               help='Comma separated list of config names (ex: ul_2016_ZZ_v12,ul_2016_HIPM_ZZ_v12,ul_2017_ZZ_v12,ul_2018_ZZ_v12)')
    parser.add_option("--cat",     dest="cat",      default='',               help='Category name to fetch Datacards from (ex : ZZ_elliptical_cut_90)')
    parser.add_option("--prd",     dest="prd",      default='',               help='Prod version')
    parser.add_option("--grp",     dest="grp",      default='datacard_zz',    help='Process group name (ex : datacard_ZZ_res)')
    parser.add_option("--feat",    dest="feat",     default='dnn_ZZbbtt_kl_1',help='Comma separated list of features')
    parser.add_option("--usr",     dest="usr",      default=None,             help='USER in case it is different from default')
    parser.add_option("--mass",    dest="mass",     default='200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,3000')
    parser.add_option("--featureDependsOnMass",     default=False,            help="Add _$MASS to name of feature for each mass for parametrized DNN", action="store_true")
    parser.add_option("--singleThread",             default=False,            help="Don't run in parallel, disable for debugging", action="store_true")
    parser.add_option("--run_ch",  dest="run_ch",   default=True,             help='Run each channel or not')
    parser.add_option("--run_year",dest="run_year", default=True,             help='Run each year or not')
    parser.add_option("--run_comb",dest="run_comb", default=True,             help='Run combination or not')
    (options, args) = parser.parse_args()

    run = options.run == True
    run_ch = options.run_ch == True
    run_year = options.run_year == True
    run_comb = options.run_comb == True
    cat = options.cat
    prd = options.prd
    grp = options.grp
    usr = options.usr
    if ',' in options.feat: features = options.feat.split(',')
    else:                   features = [options.feat]
    if ',' in options.cfg:  config = options.cfg.split(',')
    else:                   config = [options.cfg]
    if ',' in options.mass: mass_points = options.mass.split(',')
    else:                   mass_points = [options.mass]

    maindir = os.getcwd() 
    base_dir = maindir + '/ResLimits/'
    prefix = config[0].split('_')[-2]

    ##########################################################            
    # RUN SINGLE LIMIT
    ########################################################## 

    def processMapPoint(mass, cfg, feat):
        cfg_dir = f'{base_dir}/{cfg}'
        os.system('mkdir -p ' + cfg_dir)
        mass_ver = prd + f'_M{mass}'
        feat_ver = feat + f'_{mass}' if options.featureDependsOnMass else feat
        if run_ch:
            for ch in ['etau', 'mutau', 'tautau']:
                run_single_limit(cfg_dir, cfg, cat, grp, prd, ch, feat_ver, mass, run, usr)

        if run_year:
            run_comb_limit_singleYear(cfg_dir, feat_ver, prd, mass, run)

    if options.singleThread:
        for feat in features:
            for cfg in config:
                for mass in mass_points:
                    processMapPoint(mass, cfg, feat)
    
    else:
        # set max_workers=1 for debugging
        with concurrent.futures.ProcessPoolExecutor(max_workers=30) as exc:
            mass_cfg_feat = list(itertools.product(mass_points, config, features))
            exc.map(processMapPoint, [x[0] for x in mass_cfg_feat], [x[1] for x in mass_cfg_feat], [x[2] for x in mass_cfg_feat])

    for feat in features:
        for cfg in config:
            PlotAsymptoticLimits(base_dir, mass_points, prd, cfg, feat, prefix)

    ##########################################################            
    # RUN COMBINATION
    ##########################################################

    def processMapPointAllYears(mass, feat):
        feat_ver = feat + f'_{mass}' if options.featureDependsOnMass else feat
        if run_comb:
            run_comb_limit_allYears(base_dir, prefix, feat_ver, prd, mass, config, run)

    if options.singleThread:
        for feat in features:
            for mass in mass_points:
                processMapPointAllYears(mass, feat)

    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=30) as exc:
            mass_feat = list(itertools.product(mass_points, features))
            exc.map(processMapPointAllYears, [x[0] for x in mass_feat], [x[1] for x in mass_feat])   

    for feat in features:
        PlotAsymptoticLimits(base_dir, mass_points, prd, 'FullRun2', feat, prefix)
