cd /grid_mnt/data__data.polcms/cms/cuisset/ZHbbtautau/combine/RunCombineZZbbtautau/
cd ../CMSSW_11_3_4/ && cmsenv && cd - && ulimit -s unlimited
RESONANT_MASSES="200,210,220,230,240,250,260,280,300,320,350,360,400,450,500,550,600,650,700,750,800,850,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2500,2600,2800,3000,3500,4000,4500,5000"


cd ResLimits_BT_v1

python3 ../RunAsymptoticLimits.py --ver bul_2018_ZZ_v12 --user_cmt cuisset --num 3 \
    --cat cat_ZZ_EC90_boosted_bb_boostedTau \
    --feat dnn_ZZbbtt_M --featureDependsOnMass --prd prod_241119e --grp datacard_zz_res_reduced \
    --mass $RESONANT_MASSES --no_run_cat --no_run_year --no_run_impacts

python3 ../RunAsymptoticLimits.py --ver bul_2018_ZZ_v12 --user_cmt cuisset \
    --cat cat_ZZ_EC90_boosted_bb_boostedTau \
    --feat dnn_ZZbbtt_M --featureDependsOnMass --prd prod_241119e --grp datacard_zz_res_reduced \
    --mass $RESONANT_MASSES  --no_run_copy --no_run_cat --no_run_year --no_run_impacts --plot_only

## debugging

combine -M AsymptoticLimits bul_2018_ZZ_v12_cat_ZZ_EC90_boosted_bb_boostedTau_dnn_ZZbbtt_M_3000_datacard_zz_res_reduced_mutau_os_iso.txt --run blind --noFitAsimov


python3 ../RunAsymptoticLimits.py --ver bul_2018_ZZ_v12 --user_cmt cuisset --num 4 \
    --cat cat_ZZ_EC90_boosted_bb_boostedTau,cat_ZZ_EC90_resolved_1b_boostedTau,cat_ZZ_EC90_resolved_2b_boostedTau \
    --feat dnn_ZZbbtt_M --featureDependsOnMass --prd prod_241119e --grp datacard_zz_res_reduced \
    --mass $RESONANT_MASSES --no_run_year --no_run_impacts
