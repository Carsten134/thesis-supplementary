# Supplementary Material, Master's Thesis
This repository contains the code developed for this thesis. It includes all explored research avenues, even those not featured in the final submission, to maintain full transparency regarding the project's evolution.

## Code for the Thesis

### `simulation-study`
The file `simulation_runs.Rmd` contains all simulation runs, that were conducted for the thesis. The results from multiple runs were also stored in the `sim_results` folder of the sub directory `plotting-simulation-results`.

### `case-study`
This code contains the final versions of the "Application" Section of the thesis, from data fetching through the API, over data processing to the final testing file in `.Rmd`.

### `educational-plots`
Some plots in the thesis (such as the kernel on torus presentation) have been soley created for illustrative purposes. This folder contains code for these plots.

## About the Process
This section will talk about some of the ideas behind the thesis, how they have been explored in early stages, and how they ended up contributing to the overall work.

### Early (failed) Implementations
In the `early-implentations` section, we can see first attempts of implementing the testing procedure $\varphi_n^*$ for $q=2$ and $p=1$. These early implementations did not scale well, were constrained in the data they could handle, and, most importantly, did not exhibit good simulation results.

The sub directories `debugging-functional-shape`, and `periodogram` show attempts of understanding, how the functional shape changes under permutation and checking whether the estimators actually worked. This is also where I developed an intuition for the isotropy heuristic.

### Python Implementations: Chasing Runtime Efficiency
After discovering the computational efficiency gains through convolution, I wanted to push this concept to the limit. In `python-implementation` first attempts can be found at transitioning to [`pytorch`](https://pytorch.org/) to run all computations vectorized on the GPU. Doing so, lead to runtime decreases when involving the GPU (because the CPU needs to write each tensor into VRAM, which takes too much time). But when running on CPU only and vectorizing sampling and randomization, we got runtime increases of a factor of $\approx 2000$.

Convinced by the imense increases in performance achieved through convolution, I rewrote the AR(1) sampler into a convolution operation on a gaussian lattice process utilizing the wold decomposition (or impulse response function). This is why we see the recursive solution for the wold weights in the thesis.

Because of these python enabled efficiencies, we ran the first simulations in in this environment. For the thesis this felt as if we would provide the reader with one implementation, but provide simulation results for an entirely different language, and computational framework. Subsequent refinements of simulations were thus carried out for `rspsp` functions and the python route is now just a case study of how far one can push algorithmic innovation for this method.

Beyond chasing runtime efficiency, first draft ideas of the periodogram based test as well as simulation studies are also contained in the repository. 


### Phoenix Test ($\varphi_n^p$): How it Challenged the Structure, Forced Innovation, and Found its Place...
The development of the periodogram-based test, $\varphi_n^p$, was less a targeted hunt and more a fortunate accident. While "tinkering" with established research and deconstructing existing methodologies, we stumbled upon a realization that felt like a breakthrough. At the time, the field was dominated by expensive, complex smoothing processes and the constant frustration of tuning hyperparameters. By simply exploring the boundaries of other people’s work, we inadvertently found a way to bypass these hurdles entirely, allowing a clean, stochastically elegant method to arise like a Phoenix from the ashes of inconsistent asymptotic behavior.

However, this serendipitous discovery created an immediate structural problem for the thesis. While $\varphi_n^p$ was a powerful and "parameter-free" specialist, the more general framework of the work, $\varphi_n^*$, was still in its infancy, implemented only for simple univariate two sample cases. I found myself at a crossroads: I had stumbled upon a highly polished, specific tool that threatened to overshadow the broader mathematical proposal. On one side was this elegant "accident," and on the other was a more general, yet less refined, method that held the actual theoretical weight of the project.

To bridge this gap and create a compelling narrative, I realized I had to put the same level of care into the general framework as I had found in the periodogram discovery. I pivoted the implementation focus back to randomization-based inference, scaling it to handle $q$-samples and developing a new framework for testing weak stationarity. This was the moment the "story" of the thesis truly took shape. Instead of letting the accidental discovery stand alone, I leveraged the generality of the randomized approach to serve as the heart of the work, ensuring the implementation finally lived up to the mathematical promise and the title of the thesis.

In the final arrangement of the material, I chose to position $\varphi_n^p$ as the "appetizer." It serves as a proof of concept—a demonstration of what is possible when we solve the issues of asymptotic behavior—before we move into the more complex, "main course" of the generalized framework. This structure ensures the reader doesn't feel the general method is a secondary alternative, but rather the ultimate destination that the periodogram discovery helped us reach.

Early implementations of $\varphi_n^p$ can be found in the `python-implementation` folder. Further I have written a small article on preliminary results of $\varphi_n^p$,  which can be found in `phoenix-test-tex` folder. This article ended up forming the basis of the *A Periodogram Based Approach* section in the thesis.

### $\varphi_n^{*iso}$: Bouncing Back From a Failed Conjecture
One of the first drafts for a heuristic for testing isotropy has been illustrated in the second presentation of the `persentations` folder. This early version relied on the conjecture, that full isotropy of the ACF is preserved in the spectral representation, meaning that:

$$
C(\|h\|_2) = c\quad \forall h\in \mathbb Z^2 \Longleftrightarrow  f(\|\omega\|_2) = \tilde c \quad \forall \omega\in [-\pi,\pi]^2
$$

Which would justify permuting periodogram ordiantes at Fourier frequency pairs with equal $L_2$ norms. While being an intuitive and to be expected result, this conjecture turned out to be false, being easily disproven by a counterexample. 

In light of this, I adapted the isotropy hypothesis to a weaker varaint that was preserved in the spectrum and simplified the evaluation of the integral to be more in line with the original approach. 