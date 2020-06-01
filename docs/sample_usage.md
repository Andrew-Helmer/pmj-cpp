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

<sub>\* "Balanced" here refers to the property of sub-sequence stratification. When generating the PMJ(0,2) samples, after generating N samples, where N is an odd power of two, the next N/2 samples, and the N/2 samples after that, should each be stratified in the (0,2) intervals. Refer to Christensen et al. for more information.</sub>

## (Do not use) Cranley-Patterson Rotations

Cranley-Patterson Rotations (i.e. jittering the entire sample sequence and wrapping it around the edges) significantly harms the convergence of the PMJ(0,2) sequence.

## Stratification in More Dimensions

Christensen et al. describe an algorithm to take 2D PMJ(0,2) sequences and shuffle them to have good stratification in higher dimensions. [See Section 11 (Discussion) of their paper](https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/paper.pdf).
