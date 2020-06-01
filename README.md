# Progressive Multi-Jittered Sample Sequences in C++

This is a C++ implementation of Progressive Multi-Jittered Sample Sequences (in 2D), based off two papers:

* [Progressive Multi-Jittered Sample Sequences (2018) by Christensen, Kensler, and Kilpatrick](https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/paper.pdf)
* [Efficient Generation of Points that Satisfy
Two-Dimensional Elementary Intervals (2019) by Matt Pharr](http://jcgt.org/published/0008/01/04/)

Much thanks to Per and Matt for their help and encouragement with this implementation.

These sample sequences are really great for certain types of Monte Carlo integration problems, especially for use in computer graphics rendering. They're currently used in Pixar's Renderman. Here's one such sequence.

<p align="center">
<img src='https://github.com/Andrew-Helmer/pmj-cpp/blob/master/docs/pmj02bn.gif'>
</p>

## Sequence Properties and Error Convergence

The PMJ(0,2) or pmj02 sample sequence has the property that any prefix of the samples are stratified on every elementary (0,2) interval.

<p align="center">
<img src='https://github.com/Andrew-Helmer/pmj-cpp/blob/master/docs/pmj02_intervals.svg'>
</p>

This yields extremely good convergence on test integrals, especially at powers of two for the PMJ(0,2) sequence. Note that uniform random sampling and typical best-candidate sampling converge at a rate of approximately N<sup>-.5</sup>.

<p align="center">
<img src='https://github.com/Andrew-Helmer/pmj-cpp/blob/master/docs/error_analysis.svg'>
</p>

It's also *very* easy to shuffle the points in a PMJ(0,2) sequence to generate another PMJ(0,2) sequence. Check out the section on Shuffling within "Using the Samples in Rendering".

## License and Attribution

See the LICENSE.

# Using the Samples in Rendering

There are three sample sequences that are useful for rendering.

* pmj02, also known as Progressive Multi-Jittered Sequences (0,2). This sequence has extremely good convergence properties - for a given number of samples, it will typically have the lowest error. 
* pmjbn, also known as Progressive Multi-Jittered Sequences with Blue Noise (or best candidate sampling). This has the best "blue noise" characteristics, which creates noise that looks very pleasing to the human eye. This might be useful for camera rays in a path-tracer, for example, especially if you're not using it in tandem with any denoising. The error convergence is not as good as pmj02, but it's still much better than uniform random sampling or best-candidate sampling.
* pmj02bn, also known as Progressive Multi-Jittered Sequences (0,2) with Blue Noise. These sequences don't have as good blue-noise properties as pmjbn, but effectively maintain the error convergence of pmj02 with better blue noise.

I've also implemented the basic Progressive Jittered samples, but that's mostly for educational purposes.

## Precomputed Sample Sequences

There are already sample sets that can be used in the sample_sets directly. You can go ahead and download them! They're formatted as (X, Y) per line.

## Generating Your Sequences

You can also generate your own sequences, using the generate_samples utility. See the "Generating Your Own Sample Sequences" section under Using the Code.

## Shuffling

It's extremely easy to shuffle a balanced\* PMJ(0,2) sequence and still have it be a balanced PMJ(0,2) sequence, maintaining the good convergence. For a PMJ(0,2) sequence of N samples, where N is a power of two:
1. Generate a random integer **r** in the range [0, N)
2. To get the **i**'th sample in a shuffled sequence, get the *i^r*th sample from the original sequence, where *^* is the bit-wise xor operator.

This property is fantastic for a renderer, because you can compute a hash from the pixel coordinates and ray-bounce, and use that to index into your table, for decorrelated samples. You *might* even be able to get away with a single very large table for your entire renderer! I don't know if anyone has tried this though. In Renderman, they store hundreds of 4096-sample tables and index into them.

If you shuffle a sequence with blue-noise characteristics, you'll likely lose those blue noise characteristics.

<sub>\* "Balanced" here refers to the property of sub-sequence stratification. When generating the PMJ(0,2) samples, after generating N samples, where N is an odd power of two, the next N/2 samples and the N/2 samples following that should each be stratified in the (0,2) intervals. Refer to Christensen et al. for more information.*</sub>

## (Do not use) Cranley-Patterson Rotations

Cranley-Patterson Rotations (i.e. jittering the entire sample sequence and wrapping it around the edges) significantly harms the convergence of the PMJ(0,2) sequence.

## Stratification in More Dimensions

Christensen et al. describe an algorithm to take 2D PMJ(0,2) sequences and shuffle them to have good stratification in higher dimensions. [See Section 11 (Discussion) of their paper](https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/paper.pdf).

# Using the Code

## Generating Your Own Sample Sequences

## Building

If you want to use this code simply to generate your own samples, see the section "Generating your own samples". You can build the sample generation using a simple "make release" command.

If you want to run the performance testing or error analysis code, you will need to install the [Bazel build tool](https://bazel.build/).

## Code Performance

Most likely, you'll want to use precomputed tables of these samples for optimal performance. See "Sample Usage". But a couple of notes on performance here:

* Because of Matt Pharr's optimization, the generation of pmj02 sequences is quite fast. On a single-thread on my 2017 Macbook Pro, I could generate 65536 samples in ~40ms, for approximately 1.65 million samples/sec.
* The best candidate sampling could most likely be a lot faster. The sample_generation/pmj.cc and sample_generation/pmj02.cc files use a constant 100 candidates per sample to generate high-quality samples, but this slows them down considerably. If you wanted to generate these samples during rendering, you'd want to use fewer candidates.
* The PMJ sequence is not well optimized. It scans a linear array of strata, but the strata should really be kept in a binary tree to be faster. Then again, there isn't much use for the PMJ sequence on its own, most likely you'd only want to use the PMJBN sequence, which you'd want to precompute to have the best blue-noise characteristics.
