# Set up

(1) Navigate to your work folder, let's say myWJPsiWorkFolder/, do:
```
setupATLAS -q
lsetup git
mkdir -p build run
```

(2) Clone the framework, here we use the branch called main-WJPsi-dev, do:
```
git clone -b main-WJPsi-dev https://github.com/lianglianghan12/TruthAnalysis.git source
```

(3) Compile the framework, do:
```
cd build
asetup AnalysisBase,21.2.104
cmake ../source
make
source x86*/setup.sh
cd ../
```

(4) Run the framework, do:
```
cd run
RunTruthAnalysis.py -s submitDir -i /afs/cern.ch/work/l/lihan/public/W_JPsi_study/Pythia8_Jpsipi_2mu_truth_derivation -p *W2JPsiPi2mumuPi.pool.root -n 10

Explanation of the arguments:
-s: specify the folder name where your output file will be stored
-i: specify the folder name where your input file is located
-p: specify the input file name. In principle, all the files with their names containing "W2JPsiPi2mumuPi.pool.root" will be considered. In our case, we only need to use one file.
-n: specify the number of events that you want the framework to process. This argument is only used for test purpose. For processing all the events, you just need to remove this argument.
```

(5) Check your output:
```
  1. We got those debugging messages printed out when executing the command RunTruthAnalysis.py. Those messages can give you more details.
  2. Those variables that we are interested in are filled in the TTree, which is further stored in the root file: submitDir-2023-*/data-TruthAna/Pythia8_Jpsipi_2mu_truth_derivation.root
```

# Next time you login, you don't need to repeat all the steps above. What you need to do is:
```
Inside myWJPsiWorkFolder/, execute: 
setupATLAS -q && pushd build && asetup && source x86*/setup.sh && popd
```

# How to understand the framework:
```
Basically, all the operations are done in the 2 files:
source/MyTruthAnalysis/MyTruthAnalysis/TruthAnaHHbbtautau.h
source/MyTruthAnalysis/Root/TruthAnaHHbbtautau.cxx

I simplified the two files, so you can easily understand how they work. As you can see in the code, 4 variables are added already. You can think about how to add more variables, which can be done in a similar way.
```

# Info
What/Why truth level analysis and how to do it in ATLAS software: [slides](https://indico.cern.ch/event/472469/contributions/1982685/attachments/1222751/1789718/truth_tutorial.pdf)

# To do
Add more variables to the code
