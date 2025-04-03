# VP-WBT
A kernel-independent SOE/SOG approximation software with VP-sum and WBT, whose high precision computing is powered by [the Multiprecision Computng Toolbox](https://www.advanpix.com/). 

You could check the [reference](https://arxiv.org/abs/2503.03183) here and the former work [VPMR](https://github.com/ZXGao97/VPMR) by Z.X.Gao. 

For ease to use, the Julia version is coming soon. The Julia version of [VPMR](https://github.com/HPMolSim/SumOfExpVPMR.jl) by X.Z.Gao is available here. 

## Installation

Follow the steps below to install the `VP-WBT.mlappinstall` app in MATLAB.

### 1. Download from GitHub

- Go to the [Releases](https://github.com/linyuanshen114/VP-WBT/releases) page of this repository.
- Download the file named `VP-WBT.mlappinstall`.

Alternatively, you can click the green **Code** button and download the entire repository as a ZIP file, then extract it and find the `.mlappinstall` file.

### 2. Install in MATLAB

There are two common ways to install a `.mlappinstall` file in MATLAB:

#### Option A: Using MATLAB Command Window

1. Open MATLAB.
2. Run the following command in the Command Window:
   ```matlab
   matlab.apputil.install('path_to/VP-WBT.mlappinstall')
   ```
   Replace `'path_to/VP-WBT.mlappinstall'` with the actual file path.

#### Option B: Using MATLAB GUI

1. Open MATLAB.
2. Go to the **APPS** tab.
3. Click on **Install App**.
4. Select the `VP-WBT.mlappinstall` file and click **Open**.

### 3. Launch the App

Once installed, you can find **VP-WBT** in the **APPS** tab. Click to launch and start using it.

## Requirements

- MATLAB R2020b or newer is recommended.
- [The Multiprecision Computing Toolbox](https://www.advanpix.com/) must be installed and added to your MATLAB path before running the app.
- App Designer must be available (included in recent versions of MATLAB).


# BSA-WBT
In our paper, BSA is used for high precision SOE approximation of various inverse power kernels, including Coulomb kernel. 
It is a quadrature method to obtain SOE/SOG as following:

$$
r^{-\alpha} = \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\int_{-\infty}^{\infty}b^{\alpha x}e^{-b^x r} dx \approx \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\sum_{\ell=-\infty}^{+\infty}b^{\alpha \ell}e^{-b^\ell r} \approx \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\sum_{\ell=N}^{M}b^{\alpha \ell}e^{-b^\ell r}
$$