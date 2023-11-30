import os,json

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

def run_single_limit(maindir, ch, feat, ver, mass, run):

    basedir = '/data_CMS/cms/vernazza/cmt/CreateDatacards/'
    datadir = basedir + f'/ul_2018_ZZ_v10/cat_ZZ_elliptical_cut_80_{ch}/{ver}'
    datafile = datadir + f'/{feat}_{ch}_os_iso.txt'
    rootfile = datadir + f'/{feat}_{ch}_os_iso.root'

    # Define output directory
    odir = maindir + f'/{ver}/{ch}/{feat}'
    print(" ### INFO: Saving output in ", odir)
    os.system('mkdir -p ' + odir)

    # Get datacards
    os.system('cp '+datafile+' '+rootfile+' '+odir)

    if run:
        # Run asymptotic limit
        print(" ### INFO: Run Asymptotic limit")
        cmd = f'combine -M AsymptoticLimits {feat}_{ch}_os_iso.txt --run blind --noFitAsimov &> combine.log'
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
    cmd = f'combineCards.py ch1={feat}_etau_os_iso.txt ch2={feat}_mutau_os_iso.txt ch3={feat}_tautau_os_iso.txt > {feat}_comb_os_iso.txt'
    os.chdir(odir)
    os.system(cmd)

    # Run asymptotic limit
    print(" ### INFO: Run Asymptotic limit")
    cmd = f'combine -M AsymptoticLimits {feat}_comb_os_iso.txt --run blind --noFitAsimov &> combine.log'
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
    parser.add_option("--feat",    dest="feat",     default='ZZKinFit_mass')
    parser.add_option("--ver",     dest="ver",      default='prod_231129_M300')
    parser.add_option("--mass",    dest="mass",     default='300')
    (options, args) = parser.parse_args()

    run = options.run
    feat = options.feat
    ver = options.ver
    mass = options.mass

    maindir = os.getcwd() + '/ResLimits'
    os.system('mkdir -p ' + maindir)

    for ch in ['etau', 'mutau', 'tautau']:
        run_single_limit(maindir, ch, feat, ver, mass, run)
    
    run_comb_limit(maindir, feat, ver, mass, run)