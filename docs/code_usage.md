# Using the Code

## Generating Your Own Sample Sequences

To generate your own sample sequence, you need to build the generate_samples too. This can be done simply with:
<pre><code>make release</code></pre>

And then running the generate_samples command-line tool:
<pre><code>./generate_samples --n=4096 --algorithm=pmj02 --out$PWD/pmj02_samples.txt</code></pre>

This outputs the sample sequence into a pmj02_samples.txt file in the current working directory.

You can also use [Bazel](https://bazel.build/).

<pre><code>bazel run -c opt generate_samples -- --n=4096 --algorithm=pmj02 --out=$PWD/pmj02_samples.txt</code></pre>

Useful options for the algorithm flag are "pmjbn", "pmj02", and "pmj02bn". See the [sample usage documentation](/docs/sample_usage.md) for the differences.

## C++14

This code uses one **C++14** feature as of now: std::make_unique. If you want to build this with C++11, you could change the few calls to that.

## Calling the code directly

The [pmj.h](/sample_generation/pmj.h) and [pmj02.h](/sample_generation/pmj02.h) header files contain the functions which generate the samples. They return an array (wrapped in a unique_ptr) of "Point" structs, which are simply two doubles: x and y.

## Error Analysis

If you want to evaluate the error of different sampling algorithms, you can do so with the [analyze_error](/analyze_error.cc) tool. You'll need to use Bazel for this.

<pre><code>bazel run -c opt analyze_error -- --algorithms=uniform,pj,pmj,pmj02 --pyfile=$PWD/analyses/error_analysis.py</code></pre>

This will output a python file that can be read, used to generate graphs of error convergence. [Here's an example Colab notebook](/analyses/PMJ(0%2C2)_Error_Analysis.ipynb) using the error_analysis.py file. If you don't supply the --pyfile flag, you can just get the final error for a given number of samples (averaged over many runs):

<pre><code>bazel run -c opt analyze_error -- --max_n=256 --runs=1024</code></pre>

## Performance Testing

If you want to evaluate the performance of the different algorithms, you can do so with the [test_performance](/test_performance.cc) tool. You'll need to use Bazel for this.

<pre><code>bazel run -c opt test_performance -- --n=65536 --runs=64 --algorithms=pj,pmj,pmj02</code></pre>

## Code Performance

Most likely, you'll want to use precomputed tables of these samples for optimal performance. See the [sample usage documentation](docs/sample_usage.md). But a couple of notes on performance here:

* Because of Matt Pharr's optimization, the generation of pmj02 sequences is quite fast. On a single-thread on my 2017 Macbook Pro, I could generate 65536 samples in ~40ms, for approximately 1.65 million samples/sec.
* The best candidate sampling could most likely be a lot faster. The sample_generation/pmj.cc and sample_generation/pmj02.cc files use a constant 100 candidates per sample to generate high-quality samples, but this slows them down considerably. If you wanted to generate these samples during rendering, you'd want to use fewer candidates.
* The PMJ sequence is not well optimized. It scans a linear array of strata, but the strata should really be kept in a binary tree to be faster. Then again, there isn't much use for the PMJ sequence on its own, most likely you'd only want to use the PMJBN sequence, which you'd want to precompute to have the best blue-noise characteristics.
