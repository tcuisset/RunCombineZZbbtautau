# Non resonant 2018 split in 3 categories (res1B, etc) prod_240528 *with ZH combination*

cd ../CMSSW_11_3_4/ && cmsenv && cd - && ulimit -s unlimited && cd NonRes_v10/

#  ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,
#ZHKinFit_mass

############ ZbbHtt # 
python3 ../RunNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard --run_year True # --move_eos --user_eos tcuisset

python3 ../RunNonResLimits.py --ver ul_2017_ZbbHtt_v12 \
    --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard --run_year False

# Impacts
python3 ../RunNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12  \
  --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet   \
  --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard --no_run_one --no_run_ch --no_run_cat --run_year --no_run_impacts --run_impacts_noMCStat

############## ZttHbb
python3 ../RunNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_orthogonal_cut_90_resolved_1b,cat_ZttHbb_orthogonal_cut_90_resolved_2b,cat_ZttHbb_orthogonal_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard 

# Impacts
python3 ../RunNonResLimits.py --ver ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12  \
  --cat cat_ZttHbb_orthogonal_cut_90_resolved_1b,cat_ZttHbb_orthogonal_cut_90_resolved_2b,cat_ZttHbb_orthogonal_cut_90_boosted_noPNet   \
  --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard --no_run_one --no_run_ch --no_run_cat --run_zh_comb_year --run_year --no_run_impacts --run_impacts_noMCStat

############## Combination
python3 ../RunNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12,ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12  \
  --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet   \
  --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard --no_run_one --no_run_ch --no_run_cat --run_zh_comb_cat --no_run_year --run_zh_comb_year --no_run_impacts --no_run_impacts_noMCStat

# Impacts
python3 ../RunNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12,ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12  \
  --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet   \
  --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard --no_run_one --no_run_ch --no_run_cat --no_run_zh_comb_cat --no_run_year --run_zh_comb_year --no_run_impacts --run_impacts_noMCStat

#########################
####################### Move to EOS
# we use all categories (Zbb&Ztt), this will fail for some but it is expected
python3 ../RunNonResLimits.py --ver ul_2016_ZbbHtt_v12,ul_2016_HIPM_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12,ul_2016_ZttHbb_v12,ul_2016_HIPM_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12  \
    --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet,cat_ZttHbb_orthogonal_cut_90_resolved_1b,cat_ZttHbb_orthogonal_cut_90_resolved_2b,cat_ZttHbb_orthogonal_cut_90_boosted_noPNet  \
    --feat dnn_ZHbbtt_kl_1 --prd prod_240528_nonres --grp datacard --no_run --no_run_one --no_run_ch --no_run_cat --no_run_year --run_zh_comb_year --plot_only --move_eos --user_eos tcuisset


rsync -rltv /grid_mnt/data__data.polcms/cms/cuisset/ZHbbtautau/combine/RunCombineZZbbtautau/NonRes tcuisset@lxplus.cern.ch:/eos/user/t/tcuisset/www/zzbbtt/CombineLimits --delete --delete-excluded \
 --filter '- impacts/' --filter '+ */' \
 --filter '+ **.png' --filter '+ **.pdf' --filter '+ **/Significance.log' --filter '+ **/*_datacard_*.txt' --filter '- *' 



#########################
####################### Resonant limits
cd ../CMSSW_11_3_4/ && cmsenv && cd - && ulimit -s unlimited && cd ResLimits_v10/


python3 ../RunAsymptoticLimits.py --ver ul_2016_HIPM_ZbbHtt_v12,ul_2016_ZbbHtt_v12,ul_2017_ZbbHtt_v12,ul_2018_ZbbHtt_v12 \
    --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --prd prod_240528 --grp datacard_res \
    --mass 600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000 

python3 ../RunAsymptoticLimits.py --ver ul_2016_HIPM_ZttHbb_v12,ul_2016_ZttHbb_v12,ul_2017_ZttHbb_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZttHbb_orthogonal_cut_90_resolved_1b,cat_ZttHbb_orthogonal_cut_90_resolved_2b,cat_ZttHbb_orthogonal_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --prd prod_240528 --grp datacard_res \
    --mass 600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000 

python3 ../RunAsymptoticLimits.py --ver ul_2018_ZbbHtt_v12,ul_2018_ZttHbb_v12 \
    --cat cat_ZbbHtt_orthogonal_cut_90_resolved_1b,cat_ZbbHtt_orthogonal_cut_90_resolved_2b,cat_ZbbHtt_orthogonal_cut_90_boosted_noPNet \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --prd prod_240528 --grp datacard_res \
    --mass 600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000 --no_run_one --no_run_ch --no_run_cat --run_zh_comb_cat

#     --move_eos --user_eos cuisset
