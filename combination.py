import os,sys
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.CMS)
ROOT.gROOT.SetBatch(True)

#######################################################################
######################### SCRIPT BODY #################################
#######################################################################

if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--run",     dest="run",      default=True)
    parser.add_option("--feat",    dest="feat",     default='dnn_ZZbbtt_kl_1')
    parser.add_option("--ver",     dest="ver",      default='prod_231019')
    (options, args) = parser.parse_args()

    feat = options.feat
    ver = options.ver
    run = int(options.run) == 1

    basedir = '/data_CMS/cms/vernazza/FrameworkNanoAOD/CombineTests/CMSSW_11_3_4/src/SignalStrengthPlotting/'

    #######################################################################
    #######################################################################
    #######################################################################

# combineCards.py Name1=dnn_new_zzbbtt_kl_1_etau_os_iso.txt Name2=dnn_new_zzbbtt_kl_1_mutau_os_iso.txt Name3=dnn_new_zzbbtt_kl_1_tautau_os_iso.txt  > dnn_comb.txt
# text2workspace.py dnn_comb.txt -o model.root
# combine -M MultiDimFit model.root --algo=singles --rMin 0 --rMax 2 --preFitValue 1 --expectSignal 1 -t -1
# combine -M MultiDimFit model.root --algo=grid --points 100 --rMin 0 --rMax 2 --preFitValue 1 --expectSignal 1 -t -1
# combine -M Significance dnn_comb.txt -t -1 --expectSignal=1
# combine -M MultiDimFit model.root -m 125 -n .bestfit.with_syst --setParameterRanges r=0,2 --saveWorkspace --preFitValue 1 --expectSignal 1 -t -1
# combine -M MultiDimFit higgsCombine.bestfit.with_syst.MultiDimFit.mH125.root --setParameterRanges r=0,2 --saveWorkspace --preFitValue 1 --expectSignal 1 -t -1 -n .scan.with_syst.statonly_correct --algo grid --points 100 --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances

    combdir = basedir + f'/{ver}/combination/{feat}'
    print(" ### INFO: Saving combination in ", combdir)
    if run: os.system('mkdir -p ' + combdir)

    cmtdir = '/data_CMS/cms/vernazza/cmt/CreateDatacards/'
    etau_file = cmtdir + f'/ul_2018_ZZ_v10/cat_ZZ_elliptical_cut_80_etau/{ver}/{feat}_etau_os_iso.txt'
    mutau_file = cmtdir + f'/ul_2018_ZZ_v10/cat_ZZ_elliptical_cut_80_mutau/{ver}/{feat}_mutau_os_iso.txt'
    tautau_file = cmtdir + f'/ul_2018_ZZ_v10/cat_ZZ_elliptical_cut_80_tautau/{ver}/{feat}_tautau_os_iso.txt'
    etau_root = cmtdir + f'/ul_2018_ZZ_v10/cat_ZZ_elliptical_cut_80_etau/{ver}/{feat}_etau_os_iso.root'
    mutau_root = cmtdir + f'/ul_2018_ZZ_v10/cat_ZZ_elliptical_cut_80_mutau/{ver}/{feat}_mutau_os_iso.root'
    tautau_root = cmtdir + f'/ul_2018_ZZ_v10/cat_ZZ_elliptical_cut_80_tautau/{ver}/{feat}_tautau_os_iso.root'

    cmd = 'combineCards.py'
    n_comb = 1
    if os.path.exists(etau_file):
        os.system('cp {} {} {}'.format(etau_file, etau_root, combdir))
        cmd += f' Name{n_comb}={feat}_etau_os_iso.txt'
        n_comb += 1
    if os.path.exists(mutau_file):
        os.system('cp {} {} {}'.format(mutau_file, mutau_root, combdir))
        cmd += f' Name{n_comb}={feat}_mutau_os_iso.txt'
        n_comb += 1
    if os.path.exists(tautau_file):
        os.system('cp {} {} {}'.format(tautau_file, tautau_root, combdir))
        cmd += f' Name{n_comb}={feat}_tautau_os_iso.txt'
        n_comb += 1
    cmd += ' > dnn_comb.txt'
    print(cmd)
    if run: os.chdir(combdir)
    if run: os.system(cmd)

    cmd = 'text2workspace.py dnn_comb.txt -o model.root'
    if run: os.system(cmd)
    cmd = 'combine -M MultiDimFit model.root --algo=singles --rMin 0 --rMax 2 --preFitValue 1 --expectSignal 1 -t -1'
    if run: os.system(cmd)
    cmd = 'combine -M MultiDimFit model.root --algo=grid --points 100 --rMin 0 --rMax 2 --preFitValue 1 --expectSignal 1 -t -1'
    if run: os.system(cmd)
    cmd = 'combine -M Significance dnn_comb.txt -t -1 --expectSignal=1'
    if run: os.system(cmd)

    cmd = 'combine -M MultiDimFit model.root -m 125 -n .bestfit.with_syst --setParameterRanges r=0,2 --saveWorkspace'\
        ' --preFitValue 1 --expectSignal 1 -t -1'
    if run: os.system(cmd)
    cmd = 'combine -M MultiDimFit higgsCombine.bestfit.with_syst.MultiDimFit.mH125.root --setParameterRanges r=0,2 '\
        '--saveWorkspace --preFitValue 1 --expectSignal 1 -t -1 -n .scan.with_syst.statonly_correct --algo grid '\
        '--points 100 --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances'
    if run: os.system(cmd)

    #######################################################################
    #######################################################################
    #######################################################################

    LS_file_etau = basedir + f'/{ver}/etau/{feat}/higgsCombineTest.MultiDimFit.mH120.root'
    if os.path.exists(LS_file_etau):
        f_etau = ROOT.TFile(LS_file_etau)
        limit_etau = f_etau.Get("limit")
        to_draw_etau = ROOT.TString("2*deltaNLL:r")
        n_etau = limit_etau.Draw( to_draw_etau.Data(), "", "l")

        x_etau = np.array(np.ndarray((n_etau), 'd', limit_etau.GetV2())[30:70])
        y_etau = np.array(np.ndarray((n_etau), 'd', limit_etau.GetV1())[30:70])

        LS_file = basedir + f'/{ver}/etau/{feat}/higgsCombineTest.Significance.mH120.significance.root'
        f = ROOT.TFile(LS_file)
        limit = f.Get("limit")
        limit.GetEntry(0)
        sig_etau = limit.limit

    else:
        x_etau, y_etau = np.array([]), np.array([])
        sig_etau = None
        print("File not found:", LS_file_etau)

    LS_file_mutau = basedir + f'/{ver}/mutau/{feat}/higgsCombineTest.MultiDimFit.mH120.root'
    if os.path.exists(LS_file_mutau):
        f_mutau = ROOT.TFile(LS_file_mutau)
        limit_mutau = f_mutau.Get("limit")
        to_draw_mutau = ROOT.TString("2*deltaNLL:r")
        n_mutau = limit_mutau.Draw( to_draw_mutau.Data(), "", "l")

        x_mutau = np.array(np.ndarray((n_mutau), 'd', limit_mutau.GetV2())[30:70])
        y_mutau = np.array(np.ndarray((n_mutau), 'd', limit_mutau.GetV1())[30:70])

        LS_file = basedir + f'/{ver}/mutau/{feat}/higgsCombineTest.Significance.mH120.significance.root'
        f = ROOT.TFile(LS_file)
        limit = f.Get("limit")
        limit.GetEntry(0)
        sig_mutau = limit.limit

    else:
        x_mutau, y_mutau = np.array([]), np.array([])
        sig_mutau = None
        print("File not found:", LS_file_mutau)

    LS_file_tautau = basedir + f'/{ver}/tautau/{feat}/higgsCombineTest.MultiDimFit.mH120.root'
    if os.path.exists(LS_file_tautau):
        f_tautau = ROOT.TFile(LS_file_tautau)
        limit_tautau = f_tautau.Get("limit")
        to_draw_tautau = ROOT.TString("2*deltaNLL:r")
        n_tautau = limit_tautau.Draw( to_draw_tautau.Data(), "", "l")

        x_tautau = np.array(np.ndarray((n_tautau), 'd', limit_tautau.GetV2())[30:70])
        y_tautau = np.array(np.ndarray((n_tautau), 'd', limit_tautau.GetV1())[30:70])

        LS_file = basedir + f'/{ver}/tautau/{feat}/higgsCombineTest.Significance.mH120.significance.root'
        f = ROOT.TFile(LS_file)
        limit = f.Get("limit")
        limit.GetEntry(0)
        sig_tautau = limit.limit

    else:
        x_tautau, y_tautau = np.array([]), np.array([])
        sig_tautau = None
        print("File not found:", LS_file_tautau)

    LS_file_comb = combdir + f'/higgsCombineTest.MultiDimFit.mH120.root'
    f_comb= ROOT.TFile(LS_file_comb)
    limit_comb = f_comb.Get("limit")
    to_draw_comb = ROOT.TString("2*deltaNLL:r")
    n_comb = limit_comb.Draw( to_draw_comb.Data(), "", "l")

    x_comb = np.array(np.ndarray((n_comb), 'd', limit_comb.GetV2())[30:70])
    y_comb = np.array(np.ndarray((n_comb), 'd', limit_comb.GetV1())[30:70])

    LS_file_comb_stat = combdir + f'/higgsCombine.scan.with_syst.statonly_correct.MultiDimFit.mH120.root'
    f_comb_stat= ROOT.TFile(LS_file_comb_stat)
    limit_comb_stat = f_comb_stat.Get("limit")
    to_draw_comb_stat = ROOT.TString("2*deltaNLL:r")
    n_comb_stat = limit_comb_stat.Draw( to_draw_comb_stat.Data(), "", "l")

    x_comb_stat = np.array(np.ndarray((n_comb_stat), 'd', limit_comb_stat.GetV2())[30:70])
    y_comb_stat = np.array(np.ndarray((n_comb_stat), 'd', limit_comb_stat.GetV1())[30:70])

    def GetLegend(x, y, x_stat, y_stat, round=3):
        central = x[np.argmin(y)]
        interval_1sigma = x[np.where(y < 1)]
        min_1sigma = np.abs(min(interval_1sigma)-central)
        max_1sigma = np.abs(max(interval_1sigma)-central)
        interval_1sigma_stat = x_stat[np.where(y_stat < 1)]
        min_1sigma_stat = np.abs(min(interval_1sigma_stat)-central)
        max_1sigma_stat = np.abs(max(interval_1sigma_stat)-central)
        r = np.round(central, round)
        up = np.round(max_1sigma, round)
        down = np.round(min_1sigma, round)
        up_stat = np.round(max_1sigma_stat, round)
        down_stat = np.round(min_1sigma_stat, round)
        up_syst = np.round(np.sqrt(max_1sigma**2 - max_1sigma_stat**2), round)
        down_syst = np.round(np.sqrt(min_1sigma**2 - min_1sigma_stat**2), round)
        return r, up, down, up_stat, down_stat, up_syst, down_syst

    f, ax = plt.subplots(figsize = [10,10])
    plt.plot(x_etau, y_etau, linewidth=2, label=r'$\tau_{e}\tau_{h}$')
    plt.plot(x_mutau, y_mutau, linewidth=2, label=r'$\tau_{\mu}\tau_{h}$')
    plt.plot(x_tautau, y_tautau, linewidth=2, label=r'$\tau_{h}\tau_{h}$')
    plt.plot(x_comb, y_comb, linewidth=2, color='firebrick', label='Combination')
    plt.plot(x_comb_stat, y_comb_stat, linewidth=2, color='firebrick', linestyle='--', label='Stat-only')
    # central = round(x_comb[np.argmin(y_comb)], 3)
    r, up, down, up_stat, down_stat, up_syst, down_syst = GetLegend(x_comb, y_comb, x_comb_stat, y_comb_stat)
    text1 = fr"$\mu = 1.00^{{+{up}}}_{{-{down}}}$"
    plt.text(.02, .98, text1, ha='left', va='top', transform=ax.transAxes, fontsize='small', 
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
    text2 = fr"$\mu = 1.00^{{+{up_syst}}}_{{-{down_syst}}}(syst)^{{+{up_stat}}}_{{-{down_stat}}}(stat)$"
    plt.text(.02, .92, text2, ha='left', va='top', transform=ax.transAxes, fontsize='small', 
             bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'-2$\Delta$ LL')
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    plt.axhline(y=3.84, color='black', linestyle='--', alpha=0.5)
    plt.text(x=x_etau[0], y=1.1, s="68% C.L.", fontsize="x-small")
    plt.text(x=x_etau[0], y=3.94, s="95% C.L.", fontsize="x-small")
    plt.legend(fontsize=20, loc='upper right', frameon=True)
    plt.ylim(-0.05, 1.1*np.max(y_comb))
    mplhep.cms.label(data=False, rlabel='2018, (13.6 TeV) 59.7 $fb^{-1}$', fontsize=20)
    plt.grid()
    savefile = combdir + '/Combination_Mu'
    plt.savefig(savefile+'.png')
    plt.savefig(savefile+'.pdf')
    print(savefile+'.png')
    plt.close() 

    f, ax = plt.subplots(figsize = [10,10])
    xs=5.52*0.033658*0.1512
    plt.plot(xs*x_etau, y_etau, linewidth=2, label=r'$\tau_{e}\tau_{h}$')
    plt.plot(xs*x_mutau, y_mutau, linewidth=2, label=r'$\tau_{\mu}\tau_{h}$')
    plt.plot(xs*x_tautau, y_tautau, linewidth=2, label=r'$\tau_{h}\tau_{h}$')
    plt.plot(xs*x_comb, y_comb, linewidth=2, color='firebrick', label='Combination')
    plt.plot(xs*x_comb_stat, y_comb_stat, linewidth=2, color='firebrick', linestyle='--', label='Stat-only')
    r, up, down, up_stat, down_stat, up_syst, down_syst = GetLegend(xs*x_comb, y_comb, xs*x_comb_stat, y_comb_stat)
    text1 = fr"$\sigma = {{{r}}}^{{+{up}}}_{{-{down}}}$ pb"
    plt.text(.02, .98, text1, ha='left', va='top', transform=ax.transAxes, fontsize='small', 
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
    text2 = fr"$\sigma = {{{r}}}^{{+{up_syst}}}_{{-{down_syst}}}(syst)^{{+{up_stat}}}_{{-{down_stat}}}(stat)$ pb"
    plt.text(.02, .92, text2, ha='left', va='top', transform=ax.transAxes, fontsize='small', 
             bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.xlabel(r'$\sigma\;(ZZ \rightarrow bb\tau\tau)$ [pb]')
    plt.ylabel(r'-2$\Delta$ LL')
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    plt.axhline(y=3.84, color='black', linestyle='--', alpha=0.5)
    plt.text(x=xs*x_etau[0], y=1.1, s="68% C.L.", fontsize="x-small")
    plt.text(x=xs*x_etau[0], y=3.94, s="95% C.L.", fontsize="x-small")
    plt.legend(fontsize=20, loc='upper right', frameon=True)
    plt.ylim(-0.05, 1.1*np.max(y_comb))
    mplhep.cms.label(data=False, rlabel='2018, (13.6 TeV) 59.7 $fb^{-1}$', fontsize=20)
    plt.grid()
    savefile = combdir + '/Combination_Sigma'
    plt.savefig(savefile+'.png')
    plt.savefig(savefile+'.pdf')
    print(savefile+'.png')
    plt.close() 

    LS_file = f'{combdir}/higgsCombineTest.Significance.mH120.root'
    f = ROOT.TFile(LS_file)
    limit = f.Get("limit")
    limit.GetEntry(0)
    b = limit.limit

    print(" ### INFO: Produce impact plots")
    cmd = 'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 --setParameterRanges r=0,2 --doInitialFit --robustFit 1'
    if run: os.system(cmd)
    cmd = 'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 --setParameterRanges r=0,2 --doFits --robustFit 1'
    if run: os.system(cmd)
    cmd = 'combineTool.py -M Impacts -d model.root -m 125 -o impacts.json'
    if run: os.system(cmd)
    cmd = 'plotImpacts.py -i impacts.json -o impacts'
    if run: os.system(cmd)
    if run: os.system('mkdir -p impacts')
    if run: os.system('mv higgsCombine_paramFit* higgsCombine_initialFit* impacts')

    print(" ### INFO: Results for combination")
    print(" ### significance etau = ", sig_etau)
    print(" ### significance mutau = ", sig_mutau)
    print(" ### significance tautau = ", sig_tautau)
    print(" ### significance combination = ", b)

    print(GetLegend(x_etau, y_etau, x_etau, y_etau, round=5))
    print(GetLegend(x_mutau, y_mutau, x_mutau, y_mutau, round=5))
    print(GetLegend(x_tautau, y_tautau, x_tautau, y_tautau, round=5))
    
    # to plot sigma, we can just multiply mu times sigma theoretical

    # combine -M MultiDimFit model.root -m 125 -n .bestfit.with_syst --setParameterRanges r=0,2 --saveWorkspace --preFitValue 1 --expectSignal 1 -t -1

# combine -M MultiDimFit higgsCombine.bestfit.with_syst.MultiDimFit.mH125.root --setParameterRanges r=0,2 --saveWorkspace --preFitValue 1 --expectSignal 1 -t -1 \
#  -n .scan.with_syst.statonly_correct --algo grid --points 100 --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances

# combine -M FitDiagnostics -d combined_datacard.txt --setParameterRanges signal=0,2 --saveShapes --saveWithUncertainties --expectSignal 1 -n FitResult
# model.root

# 'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 --setParameterRanges r=0,2 --doInitialFit --robustFit 1'
# 'combineTool.py -M Impacts -d model.root -m 125 --expectSignal 1 -t -1 --preFitValue 1 --setParameterRanges r=0,2 --doFits --robustFit 1'
# 'combineTool.py -M Impacts -d model.root -m 125 -o impacts.json'
# 'plotImpacts.py -i impacts.json -o impacts'
