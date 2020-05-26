cc_binary(
    name = "generate_samples",
    srcs = ["generate_samples.cc"],
    deps = [
    	"//sample_generation:sample_generation",
    ],
)

cc_binary(
    name = "test_performance",
    srcs = ["test_performance.cc"],
    deps = [
    	"//sample_generation:sample_generation",
    	"@com_google_absl//absl/strings",
    	"@com_github_gflags_gflags//:gflags",
    ],
)