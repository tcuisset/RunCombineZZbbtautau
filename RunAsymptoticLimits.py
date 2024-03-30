import os,json,concurrent.futures

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

def run_single_limit(maindir, config, base_cat, process_group_name, ch, feat, ver, mass, run):

    basedir = '/data_CMS/cms/' + os.environ["USER"] + '/cmt/CreateDatacards/'
    datadir = basedir + f'/{config}/cat_{base_cat}_{ch}/{ver}'
    datafile = datadir + f'/{feat}_{process_group_name}_{ch}_os_iso.txt'
    rootfile = datadir + f'/{feat}_{process_group_name}_{ch}_os_iso.root'

    # Define output directory
    odir = maindir + f'/{ver}/{ch}/{feat}'
    print(" ### INFO: Saving output in ", odir)
    os.system('mkdir -p ' + odir)

    # Get datacards
    os.system(f'cp {datafile} {odir}/{feat}_{ch}_os_iso.txt')
    os.system(f'cp {rootfile} {odir}/{feat}_{ch}_os_iso.root')

    if run:
        # Run asymptotic limit
        print(" ### INFO: Run Asymptotic limit")
        cmd = f'combine -M AsymptoticLimits {feat}_{ch}_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
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

def run_comb_limit(maindir, feat, ver, mass, run):

    # Define output directory
    odir = maindir + f'/{ver}/combination/{feat}'
    print(" ### INFO: Saving output in ", odir)
    os.system('mkdir -p ' + odir)

    # Get datacards
    for ch in ['etau', 'mutau', 'tautau']:
        datafile = maindir + f'/{ver}/{ch}/{feat}' + f'/{feat}_{ch}_os_iso.txt'
        rootfile = maindir + f'/{ver}/{ch}/{feat}' + f'/{feat}_{ch}_os_iso.root'
        os.system('cp '+datafile+' '+rootfile+' '+odir)
        os.system('cp '+datafile+' '+rootfile+' '+odir)

    # Combine datacards
    cmd = f'combineCards.py ch_etau={feat}_etau_os_iso.txt ch_mutau={feat}_mutau_os_iso.txt ch_tautau={feat}_tautau_os_iso.txt > {feat}_comb_os_iso.txt'
    os.chdir(odir)
    os.system(cmd)

    # Run asymptotic limit
    print(" ### INFO: Run Asymptotic limit")
    cmd = f'combine -M AsymptoticLimits {feat}_comb_os_iso.txt --run blind --noFitAsimov {comb_options} &> combine.log'
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

if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run",     dest="run",      default=True)
    parser.add_option("--config", help="Configuration name (ex : ul_2018_ZZ_v12)")
    parser.add_option("--category", help="Category name (ex : ZbbHtt_elliptical_cut_90) to fetch Datacards from")
    parser.add_option("--process-group-name", help="Value of process_group_name used to make datacards (ex : datacard_ZbbHtt_res)")
    parser.add_option("--feat",    dest="feat",     default='ZZKinFit_mass')
    parser.add_option("--featureDependsOnMass", help="Add _$MASS to name of feature for each mass -> for parametrized DNN", action="store_true", default=False)
    parser.add_option("--ver",     dest="ver",      default='prod_231129')
    parser.add_option("--mass",    dest="mass",     default='200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,3000')
    (options, args) = parser.parse_args()

    run = options.run
    feat = options.feat
    ver = options.ver
    if ',' in options.mass:
        mass_points = options.mass.split(',')
    else:
        mass_points = [options.mass]

    maindir = os.getcwd() + f'/ResLimits/{options.config}'
    os.system('mkdir -p ' + maindir)

    def processMapPoint(mass):
        mass_ver = ver + f'_M{mass}'
        feat_ver = feat + f'_{mass}' if options.featureDependsOnMass else feat
        for ch in ['etau', 'mutau', 'tautau']:
            run_single_limit(maindir, options.config, options.category, options.process_group_name, ch, feat_ver, mass_ver, mass, run)

        run_comb_limit(maindir, feat_ver, mass_ver, mass, run)
    
    # set max_workers=1 for debugging
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as exc:
        exc.map(processMapPoint, mass_points)