cd ../CMSSW_11_3_4/; cmsenv; cd -

###############################################################################################################################################
############################################################  RESONANT  #################################################################
###############################################################################################################################################

python3 RunAsymptoticLimits.py --config ul_2018_ZbbHtt_v12 --category ZbbHtt_elliptical_cut_90 --process-group-name datacard_ZbbHtt_res \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --mass 500,1000,2000,3000

for CHANNEL in etau mutau tautau combination; do
python3 PlotAsymptoticLimits.py --mass 500,1000,2000,3000 --config ul_2018_ZbbHtt_v12  \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --ch $CHANNEL &
done

# python3 PlotAsymptoticLimits.py --mass 500,1000,2000,3000 --config ul_2018_ZbbHtt_v12  \
#     --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --ch combination


### ZttHbb
python3 RunAsymptoticLimits.py --config ul_2018_ZttHbb_v12 --category ZttHbb_elliptical_cut_90 --process-group-name datacard_ZttHbb_res \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --mass 500,1000,2000,3000

for CHANNEL in etau mutau tautau combination; do
python3 PlotAsymptoticLimits.py --mass 500,1000,2000,3000 --config ul_2018_ZttHbb_v12  \
    --feat dnn_ZHbbtt_kl_1 --featureDependsOnMass --ver prod_240329 --ch $CHANNEL &
done
