# VP-WBT
A kernel-independent SOE/SOG approximation software with VP-sum and WBT, whose high precision computing is powered by [the Multiprecision Computng Toolbox](https://www.advanpix.com/). 

You could check the [reference](https://arxiv.org/abs/2503.03183) here and the former work [VPMR](https://github.com/ZXGao97/VPMR) by Z.X.Gao. 

For ease to use, the Julia version is coming soon. The Julia version of [VPMR](https://github.com/HPMolSim/SumOfExpVPMR.jl) by X.Z.Gao is available here. 

# BSA-WBT
In our paper, BSA is used for high precision SOE approximation of various inverse power kernels, including Coulomb kernel. 
It is a quadrature method to obtain SOE/SOG as following:

$$
r^{-\alpha} = \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\int_{-\infty}^{\infty}b^{\alpha x}e^{-b^x r} dx \approx \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\sum_{\ell=-\infty}^{+\infty}b^{\alpha \ell}e^{-b^\ell r}
$$