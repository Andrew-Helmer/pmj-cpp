# Progressive Multi-Jittered Sample Sequences

This is a C++ implementation of Progressive Multi-Jittered Sample Sequences (in 2D), based off two papers:

* [Progressive Multi-Jittered Sample Sequences (2018) by Christensen, Kensler, and Kilpatrick](https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/paper.pdf)
* [Efficient Generation of Points that Satisfy
Two-Dimensional Elementary Intervals (2019) by Matt Pharr](http://jcgt.org/published/0008/01/04/)

Much thanks to Per Christensen and Matt Pharr for help and encouragement with this implementation. For a Rust implementation, see [Simon Brown's repository](https://github.com/sjb3d/pmj). Simon had a key insight to get optimal integration performance, something that I never would've figured out myself.

These sample sequences are really great for certain types of Monte Carlo integration problems, especially for use in computer graphics rendering. Here's one such sequence, with a subdividing grid to illustrate the most basic progressive stratification.

<p align="center">
<img src='https://github.com/Andrew-Helmer/pmj-cpp/blob/master/docs/pmj02bn.gif'>
</p>

## Using the Sample Sequences for Rendering

See the [Sample Usage Documentation](docs/sample_usage.md). Also provided are some precomputed sample sequences in the [sample_sequences directory](/sample_sequences), so you don't even need to run any of this code.

## Using the Code (Building and Generating Your Own Samples)

See the [Code Usage Documentation](docs/code_usage.md).

## Sequence Properties and Convergence

The PMJ(0,2) or pmj02 sample sequence has the property that any prefix of the samples are stratified on every elementary (0,2) interval.

<p align="center">
<img src='https://github.com/Andrew-Helmer/pmj-cpp/blob/master/docs/pmj02_intervals.svg'>
</p>

This yields extremely good convergence on test integrals, especially at powers of two for the PMJ(0,2) sequence. Note that both uniform random sampling, also known as "jittered" sampling, and [best-candidate sampling](http://www.pbr-book.org/3ed-2018/Sampling_and_Reconstruction/Maximized_Minimal_Distance_Sampler.html) converge at a rate of approximately N<sup>-.5</sup>. Refer to Christensen et al. for comparison against more sample sequences.

<p align="center">
<img src='https://github.com/Andrew-Helmer/pmj-cpp/blob/master/docs/error_analysis.svg'>
</p>

It's also very easy to shuffle the ordering of points in a PMJ(0,2) sequence to generate another PMJ(0,2) sequence. Check out [the section on shuffling](docs/sample_usage.md#Shuffling).

The Progressive Multi-Jittered Sequence with Blue Noise, or pmjbn, doesn't have as good convergence as pmj02, but it does have very nice blue noise properties, while still having better convergence than best-candidate sampling. This animation shows how the samples are (mostly) distributed far away from each other at any given number of samples.

<p align="center">
<img src='https://github.com/Andrew-Helmer/pmj-cpp/blob/master/docs/pmjbn.gif'>
</p>

## Licensing

See the [LICENSE](/LICENSE).
