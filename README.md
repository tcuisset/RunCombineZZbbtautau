# RunCombineZZbbtautau

## Install

Install combine:
```
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit

cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.1.0
scramv1 b clean; scramv1 b # always make a clean build
```

Install combine harvester:
```
cd ../../
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
git checkout v2.0.0
scram b
```

Clone this repository:
```
git clone git@github.com:elenavernazza/RunCombineZZbbtautau.git
```

## Description

## Run

```
cmsenv 

python3 run_combine.py --ver prod_231128 --ch etau --feat dnn_ZZbbtt_kl_1
python3 run_combine.py --ver prod_231128 --ch mutau --feat dnn_ZZbbtt_kl_1
python3 run_combine.py --ver prod_231128 --ch tautau --feat dnn_ZZbbtt_kl_1
python3 combination.py --ver prod_231128 --feat dnn_ZZbbtt_kl_1

python3 run_combine.py --ver prod_231128 --ch etau --feat ZZKinFit_mass
python3 run_combine.py --ver prod_231128 --ch mutau --feat ZZKinFit_mass
python3 run_combine.py --ver prod_231128 --ch tautau --feat ZZKinFit_mass
python3 combination.py --ver prod_231128 --feat ZZKinFit_mass
```