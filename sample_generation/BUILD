package(default_visibility = ["//visibility:public"])

# If you want symbols on a Mac, build with -c dbg and --spawn_strategy=local. 
# This is annoyingly necessary.
genrule(
    name = "sample_generation_dsym",
    srcs = [":sample_generation"],
    outs = ["sample_generation.dSYM"],
    output_to_bindir = True,
    cmd = "dsymutil $(location :sample_generation) -o $@"
)

cc_library(
    name = "sample_generation",
    srcs = glob(["*.cc"]),
    hdrs = glob(["*.h"]),
)

cc_library(
    name = "pj",
    srcs = ["pj.cc"],
    hdrs = ["pj.h"],
    deps = [":util"]
)

cc_library(
    name = "pmj",
    srcs = ["pmj.cc"],
    hdrs = ["pmj.h"],
    deps = [
        ":pmj_util",
        ":util",
    ]
)

cc_library(
    name = "pmj_util",
    srcs = ["pmj_util.cc"],
    hdrs = ["pmj_util.h"],
    deps = [":util"]
)

cc_library(
    name = "pmj02",
    srcs = ["pmj02.cc"],
    hdrs = ["pmj02.h"],
    deps = [":util"]
)

cc_library(
    name = "util",
    srcs = ["util.cc"],
    hdrs = ["util.h"],
)