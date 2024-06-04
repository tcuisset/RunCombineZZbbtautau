# Non resonant 2018 split in 3 categories (res1B, etc)

cd ../CMSSW_11_3_4/; cmsenv; cd -

#  ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,
#ZHKinFit_mass

############ ZbbHtt # 
python3 ../RunNonResLimits.py --ver ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_resolved_1b,cat_ZbbHtt_elliptical_cut_90_resolved_2b,cat_ZbbHtt_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240517_nonres --grp datacard_ZbbHtt


python3 ../PlotNonResLimits.py --ver ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_resolved_1b,cat_ZbbHtt_elliptical_cut_90_resolved_2b,cat_ZbbHtt_elliptical_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240517_nonres --grp datacard_ZbbHtt

# testing
python3 ../RunNonResLimits.py --ver ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_boosted_noPNet --channels etau \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240517_nonres_testQCD --grp datacard_ZbbHtt

############## ZttHbb
python3 ../RunNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_etau,cat_ZttHbb_elliptical_cut_90_mutau,cat_ZttHbb_elliptical_cut_90_tautau \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240418_nonres --grp datacard_ZttHbb


python3 ../PlotNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_etau,cat_ZttHbb_elliptical_cut_90_mutau,cat_ZttHbb_elliptical_cut_90_tautau \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240418_nonres --grp datacard_ZttHbb




########## ZHKinFit_mass

python3 RunNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_etau,cat_ZbbHtt_elliptical_cut_90_mutau,cat_ZbbHtt_elliptical_cut_90_tautau \
    --feat ZHKinFit_mass --prd prod_240403_nonres --grp datacard_ZbbHtt


python3 PlotNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_etau,cat_ZbbHtt_elliptical_cut_90_mutau,cat_ZbbHtt_elliptical_cut_90_tautau \
    --feat ZHKinFit_mass --prd prod_240403_nonres --grp datacard_ZbbHtt




############## ZttHbb
python3 RunNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_etau,cat_ZttHbb_elliptical_cut_90_mutau,cat_ZttHbb_elliptical_cut_90_tautau \
    --feat ZHKinFit_mass --prd prod_240403_nonres --grp datacard_ZttHbb


python3 PlotNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_etau,cat_ZttHbb_elliptical_cut_90_mutau,cat_ZttHbb_elliptical_cut_90_tautau \
    --feat ZHKinFit_mass --prd prod_240403_nonres --grp datacard_ZttHbb

#########################
####################### Move to EOS
rsync -rltv /grid_mnt/data__data.polcms/cms/cuisset/ZHbbtautau/combine/RunCombineZZbbtautau/NonRes tcuisset@lxplus.cern.ch:/eos/user/t/tcuisset/www/zzbbtt/CombineLimits --delete --delete-excluded \
 --filter '- impacts/' --filter '+ */' \
 --filter '+ **.png' --filter '+ **.pdf' --filter '+ **/Significance.log' --filter '+ **/*_datacard_*.txt' --filter '- *' 


