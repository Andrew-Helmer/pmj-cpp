package(default_visibility = ["//visibility:public"])

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