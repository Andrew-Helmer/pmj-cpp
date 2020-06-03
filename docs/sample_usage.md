# Using the Samples in Rendering

There are three sample sequences that could be useful for rendering.

* ***pmj02***, also known as Progressive Multi-Jittered Sequences (0,2). This sequence has extremely good convergence properties - for a given number of samples, it will typically have the lowest error. 
* ***pmjbn***, also known as Progressive Multi-Jittered Sequences with Blue Noise (or best candidate sampling). This has the best "blue noise" characteristics, which creates noise that looks very pleasing to the human eye. So it might be useful for camera rays in a path-tracer, for example, especially if you're not using it in tandem with any denoising. The error convergence is not as good as pmj02, but it's still much better than uniform random sampling or best-candidate sampling.
* ***pmj02bn***, also known as Progressive Multi-Jittered Sequences (0,2) with Blue Noise. These sequences don't have as good blue-noise properties as pmjbn, but better blue noise properties than pmj02, while still maintaining the error convergence of pmj02.

I've also implemented the basic Progressive Jittered samples, but that's mostly for educational purposes.

## Precomputed Sample Sequences

This repository contains a set of precomputed sample sequences you can download and use directly, in the [sample_sequences directory](/sample_sequences). They're formatted as one (X, Y) point per line in plain text.

## Generating Your Sequences

You can also generate your own sequences, using the generate_samples utility. See the ["Generating Your Own Sample Sequences" section](/docs/code_usage.md#generating-your-own-sample-sequences) of Using the Code.

## References

You may want to consult the Discussion section (section 11) of [Christensen et al](https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/paper.pdf) for some information of how these samples are used in RenderMan.

Another good reference might be the [(0,2) Sequence](http://www.pbr-book.org/3ed-2018/Sampling_and_Reconstruction/(0,_2)-Sequence_Sampler.html) section of PBRT.

## Shuffling

It's easy to shuffle a balanced\* PMJ(0,2) sequence and still have it be a balanced (0,2) sequence, maintaining the good convergence, with at most 1 additional integer of memory. For a PMJ(0,2) sequence of N samples, where N is a power of two:
1. Generate a random integer <code>r</code> in the range [0, N)
2. To get the <code>i</code>'th sample in a shuffled sequence, get the <code>(i^r)</code>'th sample from the original sequence, where <code>^</code> is the bit-wise xor operator.

This property may be useful to a renderer, because you can compute a hash from the pixel coordinates and ray-bounce, and use that to index into your sequence, which may help to decorrelate precomputed sequences. However, if you shuffle a sequence with progressive blue-noise characteristics, you'll likely lose the progressive blue-noise characteristics.

According to Christensen et al., in Renderman they store hundreds of 4096-sample tables and index into one of the tables based on a hash of pixel coordinates and bounce.

It's possible to do a somewhat better shuffle, and it's fast, but it's not so memory-efficient. I think you would need to store the whole list of indices. The procedure is:
1. Iterate over sequential pairs in the sequence. For every pair, randomly decide (with 50% probability) whether to swap it. E.g. randomly swap index 0 with index 1, randomly swap 2 with 3, randomly swap 4 with 5, etc.
2. Iterate over sequential pairs of two values. For every pair, randomly decide whether to swap *both* values. E.g. randomly swap indices 0,1 with 2,3. Randomly swap 4,5 with 6,7, etc.
3. Iterate over sequential pairs of 4 values, and randomly swap the 4 values at a time. For instance, randomly swap 0,1,2,3 with 4,5,6,7. 
4. Continue multiplying each length by two, until you've randomly swapped the first half of the entire sequence with the second half.

Both of these methods work because in a balanced PMJ(0,2) sequences, any sub-sequence with a power of two length, and starting at an integer multiple of its length, is itself a balanced progressive (0,2) sequence. So samples 1-4 (indices 0-3) are an (0,2) sequence, as are samples 5-8, 9-13, etc. Samples 1-16 are an (0,2) sequence, as are 17-32, 33-48, etc.

<sub>\* "Balanced" here refers to the property of sub-sequence stratification. When generating the PMJ(0,2) samples, after generating N samples, where N is an odd power of two, the next N/2 samples, and the N/2 samples after that, should each be (0,2) sequences themselves. Refer to Christensen et al. for more information.</sub>

## (Do not use) Cranley-Patterson Rotations

Cranley-Patterson Rotations (i.e. jittering the entire sample sequence and wrapping it around the edges) harm the convergence of the PMJ(0,2) sequence on test integrals, so probably better to not use them in conjunction.

## Stratification in More Dimensions

In Christensen et al. they describe an offline algorithm to take precomputed 2D PMJ(0,2) sequences and shuffle them to have good stratification in higher dimensions. [See Section 11 (Discussion) of their paper](https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/paper.pdf).
