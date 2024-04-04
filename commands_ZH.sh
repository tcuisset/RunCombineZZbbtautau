cd ../CMSSW_11_3_4/; cmsenv; cd -

#  ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,
#ZHKinFit_mass

############ ZbbHtt
python3 RunNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_etau,cat_ZbbHtt_elliptical_cut_90_mutau,cat_ZbbHtt_elliptical_cut_90_tautau \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240403_nonres --grp datacard_ZbbHtt


python3 PlotNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_elliptical_cut_90_etau,cat_ZbbHtt_elliptical_cut_90_mutau,cat_ZbbHtt_elliptical_cut_90_tautau \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240403_nonres --grp datacard_ZbbHtt




############## ZttHbb
python3 RunNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_etau,cat_ZttHbb_elliptical_cut_90_mutau,cat_ZttHbb_elliptical_cut_90_tautau \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240403_nonres --grp datacard_ZttHbb


python3 PlotNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_elliptical_cut_90_etau,cat_ZttHbb_elliptical_cut_90_mutau,cat_ZttHbb_elliptical_cut_90_tautau \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240403_nonres --grp datacard_ZttHbb




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
rsync -rltv /grid_mnt/data__data.polcms/cms/cuisset/ZHbbtautau/combine/RunCombineZZbbtautau/NonRes tcuisset@lxplus.cern.ch:/eos/user/t/tcuisset/www/zzbbtt/CombineLimits --delete \
 --filter '- impacts/' --filter '+ */' \
 --filter '+ **.png' --filter '+ **.pdf' --filter '+ **/Significance.log' --filter '+ **/*_datacard_*.txt' --filter '- *' 


###############################################################################################################################################
############################################################  RESONANT  #################################################################
###############################################################################################################################################
# MASSES=500,1000,2000,3000
MASSES=500,1000,2000


#################### parametrized DNN
python3 RunAsymptoticLimits.py --config ul_2018_ZbbHtt_v12 --category ZbbHtt_elliptical_cut_90 --process-group-name datacard_ZbbHtt_res \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --mass $MASSES

for CHANNEL in etau mutau tautau combination; do
python3 PlotAsymptoticLimits.py --mass $MASSES --config ul_2018_ZbbHtt_v12  \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --ch $CHANNEL &
done

### ZttHbb
python3 RunAsymptoticLimits.py --config ul_2018_ZttHbb_v12 --category ZttHbb_elliptical_cut_90 --process-group-name datacard_ZttHbb_res \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --mass $MASSES

for CHANNEL in etau mutau tautau combination; do
python3 PlotAsymptoticLimits.py --mass $MASSES --config ul_2018_ZttHbb_v12  \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --ch $CHANNEL &
done

#################### non-parametrized DNN
python3 RunAsymptoticLimits.py --config ul_2018_ZbbHtt_v12 --category ZbbHtt_elliptical_cut_90 --process-group-name datacard_ZbbHtt_res \
    --feat dnn_ZHbbtt_kl_1 --ver prod_240329 --mass $MASSES

for CHANNEL in etau mutau tautau combination; do
python3 PlotAsymptoticLimits.py --mass $MASSES --config ul_2018_ZbbHtt_v12  \
    --feat dnn_ZHbbtt_kl_1ZttHbb --ver prod_240329 --ch $CHANNEL &
done

### ZttHbb
python3 RunAsymptoticLimits.py --config ul_2018_ZttHbb_v12 --category ZttHbb_elliptical_cut_90 --process-group-name datacard_ZttHbb_res \
    --feat dnn_ZHbbtt_kl_1ZttHbb --ver prod_240329 --mass $MASSES

for CHANNEL in etau mutau tautau combination; do
python3 PlotAsymptoticLimits.py --mass $MASSES --config ul_2018_ZttHbb_v12  \
    --feat dnn_ZHbbtt_kl_1ZttHbb --ver prod_240329 --ch $CHANNEL &
done


#########################
####################### Move to EOS
rsync -rltv /grid_mnt/data__data.polcms/cms/cuisset/ZHbbtautau/combine/RunCombineZZbbtautau/ResLimits/ tcuisset@lxplus.cern.ch:/eos/user/t/tcuisset/www/zzbbtt/CombineLimits/ResLimits --delete --delete-excluded \
  --filter '+ */' \
 --filter '+ **.png' --filter '+ **.pdf' --filter '+ **/Significance.log' --filter '+ **/*_datacard_*.txt' --filter '- *' 