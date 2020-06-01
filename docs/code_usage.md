# Using the Code

## Generating Your Own Sample Sequences

To generate your own sample sequence, you need to build the generate_samples too. This can be done simply with:
<pre><code>make release</code></pre>

And then running the generate_samples command-line utility, redirecting the standard output to a file:
<pre><code>./generate_samples --n=4096 --algorithm=pmj02 > pmj02_samples.txt</code></pre>

This outputs the sample sequence into the pmj02_samples.txt file.

You can also use [Bazel](https://bazel.build/).

<pre><code>bazel run -c opt generate_samples -- --n=4096 --algorithm=pmj02 > pmj02_samples.txt</code></pre>

Useful options for the algorithm flag are "pmjbn", "pmj02", and "pmj02bn". See the [sample usage documentation](docs/sample_usage.md) for the differences.

## Building

If you want to use this code simply to generate your own samples, see the section "Generating your own samples". You can build the sample generation using a simple "make release" command. The "release" isn't necessary, but it does significantly speed up sample generation.

If you want to run the performance testing or error analysis code, you will need to install the [Bazel build tool](https://bazel.build/).

Note that this C++ code uses one C++14 feature as of now: std::make_unique. If you want to build this with C++11, you could change the few calls to that.

## Code Performance

Most likely, you'll want to use precomputed tables of these samples for optimal performance. See the [sample usage documentation](docs/sample_usage.md). But a couple of notes on performance here:

* Because of Matt Pharr's optimization, the generation of pmj02 sequences is quite fast. On a single-thread on my 2017 Macbook Pro, I could generate 65536 samples in ~40ms, for approximately 1.65 million samples/sec.
* The best candidate sampling could most likely be a lot faster. The sample_generation/pmj.cc and sample_generation/pmj02.cc files use a constant 100 candidates per sample to generate high-quality samples, but this slows them down considerably. If you wanted to generate these samples during rendering, you'd want to use fewer candidates.
* The PMJ sequence is not well optimized. It scans a linear array of strata, but the strata should really be kept in a binary tree to be faster. Then again, there isn't much use for the PMJ sequence on its own, most likely you'd only want to use the PMJBN sequence, which you'd want to precompute to have the best blue-noise characteristics.

## Performance Testing

If you want to evaluate the performance of the different algorithms, you can do so with the [test_performance](/test_performance.cc) tool. You'll need to use Bazel for this.

<pre><code>bazel run -c opt test_performance -- --n=65536 --runs=64 --algorithms=pj,pmj,pmj02</code></pre>
