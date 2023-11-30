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
git@github.com:elenavernazza/RunCombineZZbbtautau.git
```

## Description

