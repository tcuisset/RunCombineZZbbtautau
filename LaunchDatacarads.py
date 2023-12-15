import os, sys

if __name__ == "__main__" :

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--mass",    dest="mass",     default='200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,3000')
    parser.add_option("--feat",    dest="feat",     default='ZZKinFit_mass,ZZKinFit_highmass,dnn_ZZbbtt_kl_1')
    parser.add_option("--ver",     dest="ver",      default='prod_231129')
    parser.add_option("--ch",      dest="ch",       default='combination')
    (options, args) = parser.parse_args()

    mass_points = options.mass.split(',')

    if options.ch == 'combination': ch_list = ['etau', 'mutau', 'tautau']
    else: 
        if ',' in options.ch:
            ch_list = options.ch.split(',')
        else:
            ch_list = [options.ch]

    print(" ### INFO: Running over mass points ", mass_points)
    print(" ### INFO: Running over channels ", ch_list)

    for ch in ch_list:

        for mass in mass_points:

            sig = f'ggXZZbbtt_M{mass}'
            bkg = 'zz_bbtt,dy,tt_dl,tt_sl,tt_fh,wjets,ttw_lnu,ttw_qq,ttww,ttwz,ttwh,ttzh,ttz_llnunu,ttz_qq,ttzz,tth_bb,tth_nonbb,tth_tautau,' \
                'zh_hbb_zll,wplush_htt,wminush_htt,ggH_ZZ,ggf_sm,zz_dl,zz_sl_background,zz_lnu,zz_qnu,wz_lllnu,wz_lnuqq,wz_llqq,wz_lnununu,ww_llnunu,ww_lnuqq,ww_qqqq,' \
                'zzz,wzz,www,wwz,ewk_z,ewk_wplus,ewk_wminus,st_tw_antitop,st_tw_top,st_antitop,st_top'
            
            if ch == 'tautau': data_ch = 'tau'
            else: data_ch = ch
            data = f'data_{data_ch}_a,data_{data_ch}_b,data_{data_ch}_c,data_{data_ch}_d'

            cmd = 'law run CreateDatacards'
            cmd += ' --version ' + options.ver + f'_M{mass}'
            cmd += ' --FeaturePlot-version ' + options.ver
            cmd += ' --PrePlot-version ' + options.ver + ' --PrePlot-skip-merging'
            cmd += ' --config-name ul_2018_ZZ_v10'
            cmd += ' --feature-names ' + options.feat
            cmd += ' --dataset-names ' + sig + ',' + bkg + ',' + data
            cmd += ' --workers 50'
            cmd += ' --MergeCategorizationStats-version prod_231005'
            cmd += ' --Categorization-version prod_231005'
            cmd += ' --process-group-name datacard_zz_res'
            cmd += ' --category-name ZZ_elliptical_cut_80_' + ch
            cmd += ' --region-name ' + ch + '_os_iso'
            cmd += ' --do-qcd'

            print('\n',cmd)
            os.system(cmd)

