# Note 9: Geographic inference

Modern catalog coordinates enter as present-day reference points. Historical
varieties do not inherit those points as ancient observations. EDH records
retain reported origin or findspot, source dates, locality precision and
unresolved language assignment. Their current role is evidence inspection;
they do not enter the IE-CoR geographic likelihood.

The planar model integrates Gaussian lineage locations by pruning. Its units
are degrees squared per year, so it is a local approximation with latitude
and dateline limitations. A fixed tip observation error can be added in
degrees. The regional alternative uses a distance-decayed reversible CTMC
over three or six fixed anchors. Anchor choice and reference error are model
assumptions that need sensitivity checks.

Offline workbench export samples conditional Gaussian locations at the root
and frequently sampled clades from retained tree and diffusion states. These
are model draws, rather than decorative jitter. Independent precision-matrix
calculations check the ancestral Gaussian marginals. The workbench separates
these clouds from modern reference marks and sourced historical localities.

A broad cloud is not proof of historical spread. Several clade or root modes
can also reflect unconverged chains. Current observed fits fail convergence,
so the exported clouds are diagnostic model distributions. A historical
multimodality claim requires chain agreement, a specified mode analysis and
sensitivity to the location model and influential observations. That claim
remains open.
