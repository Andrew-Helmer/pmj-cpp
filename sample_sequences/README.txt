The license also applies to use of these sample sets. If you use these precomputed sample sets in your renderer, I'd appreciate a little attribution in your source code.

There are three algorithms of precomputed sample sets: pmjbn, pmj02, and pmj02bn.

pmjbn has the best blue-noise characteristics. Blue-noise causes errors that are more pleasing and subtle to the eye. It would especially be useful, in ray tracing, for primary samples (i.e. camera rays), as well as first bounce samples, e.g. direct light shadow rays cast from a point of primary sample visiblity.

pmj02 has the best error convergence overall. It would probably be good for something like multi-bounce illumination. It would probably also be best for primary rays if you knew, for example, that you were going to use a denoiser.

pmj02bn is sort of in between. It might be good for the first bounce of indirect illumination.