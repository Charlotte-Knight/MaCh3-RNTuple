# MaCh3-RNTuple

I wrote this repository to experiment with inserting RDataFrame into [MaCh3](https://github.com/mach3-software/MaCh3). There are several scripts here that represent different stages of my testing. The majority are still here for prosperity. For my final conclusions, the important files are:
- [complexity_test.cpp](complexity_test.cpp)
- [FastTSpline3Eval.h](FastTSpline3Eval.h)
- [optimised_splines.cpp](optimised_splines.cpp)

The following tests were run using the LCG 108 (x86_64-el9-gcc15-opt) release which comes with ROOT 5.36.02.  

## Complexity Test

The aim of [complexity_test.cpp](complexity_test.cpp) was to ascertain how much of an overhead using RDataFrame introduces compared to an approach running over C++ std vectors. I was curious to know how complex the operation on an event must be before the overhead becomes negligible. 

For this test, the script reads [RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root](RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root) which is a file of ~50K events, and creates a histogram based on ELep. The 'operations' performed on the event are to multiply the event weight by an arbitary number of sin functions (see [here](complexity_test.cpp#L21)). By varying the number of sin functions, I could change the complexity of the operation, and then time both approaches (RDataFrame vs std vectors).

To run this test:
```
g++ -O3 $(root-config --cflags --libs) -o complexity_test.out complexity_test.cpp
./complexity_test.out RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root
```

On my machine (Intel(R) Xeon(R) Gold 6240 CPU), I got the following results:
```
Number of events: 49996
Running with n_sin = 1
  Running vectors
    Average time per trial: 1 milliseconds
    Time per trial per event: 35 nanoseconds
  Running RDataFrame
    Average time per trial: 3 milliseconds
    Time per trial per event: 74 nanoseconds
Running with n_sin = 5
  Running vectors
    Average time per trial: 4 milliseconds
    Time per trial per event: 87 nanoseconds
  Running RDataFrame
    Average time per trial: 6 milliseconds
    Time per trial per event: 124 nanoseconds
Running with n_sin = 10
  Running vectors
    Average time per trial: 7 milliseconds
    Time per trial per event: 152 nanoseconds
  Running RDataFrame
    Average time per trial: 9 milliseconds
    Time per trial per event: 192 nanoseconds
Running with n_sin = 50
  Running vectors
    Average time per trial: 35 milliseconds
    Time per trial per event: 702 nanoseconds
  Running RDataFrame
    Average time per trial: 38 milliseconds
    Time per trial per event: 779 nanoseconds
Running with n_sin = 100
  Running vectors
    Average time per trial: 70 milliseconds
    Time per trial per event: 1411 nanoseconds
  Running RDataFrame
    Average time per trial: 73 milliseconds
    Time per trial per event: 1476 nanoseconds
```

which teaches us that by the time the operations equate to ~1 microsecond, the two approaches become similar in speed. 

To test whether the total number of events being ran over made any difference, I hadded together one of the root files to create a bigger dataset and ran on that 
```
hadd big_file.root RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root{,,,,,,,,,}
./complexity_test.out RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root
```
but the time per events were were the same.

Testing a DUNE MaCh3 fit, I found a step time of about 1s running over 1.24M events, which corresponds to 0.77 microseconds per event. Therefore, we are in the territory where the two approaches have similar performance.

However, in [MaCh3Tutorial](https://github.com/mach3-software/MaCh3Tutorial/tree/main), which is intentionally simple, the time per step of 0.04 microseconds, which is equivalent to the `n_sin = 1` scenario, where RDataFrame is ~2x slower than std vectors. Given that MaCh3Tutorial is often used to benchmark performance, this is undesirable. Furthermore, MaCh3 developments may bring the time per events down further than it is currently, in which case, these results suggest we would begin to be hampered by RDataFrame.

There are some advantages still to using RDataFrame. They will be discussed later.

## Replicating MaCh3Tutorial with RDataFrame

To get a feeling for using RDataFrame, I set out to implement a 'fit' similar to MaCh3Tutorial. Here, I skip the likelihood calculating and MCMC parts, and just focus on the histogram reweighting. The steps that are implemented and benchmarked are:
- evaluating the binned splines
- shifting ELep
- a norm weight calculation based on q2
- applying a selection on Enu_true
- histogramming based on ELep

I followed MaCh3Tutorial fairly closely, with the exception of the splines which I did not load completely. However, when needed, I created copies of the splines I did load to introduce complexity to the 'fit'. 

The bulk of the implementation can be found in [optimised_splines.cpp](optimised_splines.cpp). The rest can be found in [FastTSpline3Eval.h](FastTSpline3Eval.h) which implements a lot of the spline optimisations found in MaCh3. I did not implement the fast rebinning seen in MaCh3 which caches the bin indices from the last MCMC step.

The fit runs over just one sample from MaCh3Tutorial [RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root](RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root) which is about 1/4 of the events, and since I do not load in all splines, the complexity is less per event. 

As before, I implemented this in RDataFrame and in C++ std vectors to compare. To run the fits:
```
g++ -O3 $(root-config --cflags --libs) -o optimised_splines.out optimised_splines.cpp
./optimised_splines.out
```

I get the following results:
```
Running vectors
Total time (RNTuple - Fast): 137 ms
Average time per trial (RNTuple - Fast): 1.37 ms
Running dataframe
Total time (RDF - Fast): 663 ms
Average time per trial (RDF - Fast): 6.63 ms
```
for comparison, MaCh3Tutorial runs in 8ms (all of these tests are on single threads). Therefore, things are in the right ballpark. However, the RDataFrame implementation is against noticeably slower. 

By editing this [line](optimised_splines.cpp#L302), I can create copies of the splines to increase complexity. When creating 100 copies, I get:
```
Running vectors
Total time (RNTuple - Fast): 557 ms
Average time per trial (RNTuple - Fast): 5.57 ms
Running dataframe
Total time (RDF - Fast): 1153 ms
Average time per trial (RDF - Fast): 11.53 ms
```
and for 1000 copies:
```
Running vectors
Total time (RNTuple - Fast): 5063 ms
Average time per trial (RNTuple - Fast): 50.63 ms
Running dataframe
Total time (RDF - Fast): 5879 ms
Average time per trial (RDF - Fast): 58.79 ms
```
When turning on RDataFrame multithreading by uncommenting this [line]([line](optimised_splines.cpp#L296)), I then get:
```
Running vectors
Total time (RNTuple - Fast): 5276 ms
Average time per trial (RNTuple - Fast): 52.76 ms
Running dataframe
Total time (RDF - Fast): 1010 ms
Average time per trial (RDF - Fast): 10.1 ms
```
which shows that RDataFrame does use multiple threads efficiently. However, when considering whether to use RDataFrame in MaCh3, doesn't help much because multithreading is already implemented in MaCh3. 

## RDataFrame vs C++ std vectors

When it comes to speed, in the regime that MaCh3 is operating in, RDataFrame currently doesn't make a lot of sense. It is possible that it scales better with multithreading, because it splits by events, rather than by operations which MaCh3 does currently, but I would not say that's a good enough argument to switch.

There are plenty of nice features that could come with RDataFrame, like natural/nice interface for reading in files which could lead to a cleaner code base. However, in practice I found that a lot of the codebase was spent on implementing e.g. the spline optimisations, which were outside of the RDataFrame world anyway. There is the snapshotting dataset feature, which could be useful for debugging, but this alone does not feel like a suitable enough reason to switch. Furthermore, given that a lot of RDataFrame is abstracted from the user (for good reason), I found it difficult to profile and diagnose speed issues. For the MaCh3 use case, I think it is preferable that we have full control and visibility.

## Conclusion

Sticking with our current reweighting and histogramming implementations in MaCh3 makes sense. 

RNTuples promise speed improvements for read in which we should take advantage of, plus some of the ROOT tools surronding RNTuple read-in may come in handy. In principle, we could use RDataFrame for initial read-in and preprocessing, before 'taking' std vectors from the RDataFrame for the speedy MCMC later on.
