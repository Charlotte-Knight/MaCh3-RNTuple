```
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc15-opt/setup.sh
```

```
g++ -O3 -o plot_rntuple_pointer_nojit plot_rntuple_pointer_nojit.cpp `root-config --cflags --libs`
```

```
./plot_rntuple_pointer_nojit RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root BinnedSplinesTutorialInputs2D.root
```

