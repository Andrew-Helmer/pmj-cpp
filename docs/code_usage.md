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