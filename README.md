# 📚 Table of Contents

- [VP-WBT](#vp-wbt)  
  *A MATLAB App for kernel-independent SOE/SOG approximation.*
  
- [BSA-WBT](#bsa-wbt)  
  *A script-based tool for high-precision SOE approximation of inverse power kernels using the BSA quadrature method.*


# VP-WBT
A kernel-independent SOE/SOG approximation software with VP-sum and WBT, whose high precision computing is powered by [the Multiprecision Computng Toolbox](https://www.advanpix.com/). 

You could check the [reference](https://arxiv.org/abs/2503.03183) here and the former work [VPMR](https://github.com/ZXGao97/VPMR) by Z.X.Gao. 

For ease to use, the Julia version is coming soon. The Julia version of [VPMR](https://github.com/HPMolSim/SumOfExpVPMR.jl) by X.Z.Gao is available here. 


## Requirements

- MATLAB R2020b or newer is recommended.
- ⚠️[The Multiprecision Computing Toolbox](https://www.advanpix.com/) must be installed and added to your MATLAB path before running the app.
- App Designer must be available (included in recent versions of MATLAB).

## Installation

Follow the steps below to install the `VP_WBT.mlappinstall` app in MATLAB.

### 1. Download from GitHub

- Go to the [Releases](https://github.com/linyuanshen114/VP-WBT/releases) page of this repository.
- Download the file named `VP_WBT.mlappinstall`.

Alternatively, you can click the green **Code** button and download the entire repository as a ZIP file, then extract it and find the `.mlappinstall` file.

### 2. Install in MATLAB

There are two common ways to install a `.mlappinstall` file in MATLAB:

#### Option A: Using MATLAB Command Window

1. Open MATLAB.
2. Run the following command in the Command Window:
   ```matlab
   matlab.apputil.install('path_to/VP_WBT.mlappinstall')
   ```
   Replace `'path_to/VP_WBT.mlappinstall'` with the actual file path.

#### Option B: Using MATLAB GUI

1. Open MATLAB.
2. Go to the **APPS** tab.
3. Click on **Install App**.
4. Select the `VP_WBT.mlappinstall` file and click **Open**.

Alternatively, if you prefer to run the code manually, the source file `VP_WBT.m` is available in the `src` folder and can be executed directly from the MATLAB editor or command window.

### 3. Launch the App

Once installed, you can find **VP_WBT** in the **APPS** tab. Click to launch and start using it.

## Usage

For basic software operations and usage examples, you may refer to the following guide:

📄 [VPMR Readme.pdf](https://github.com/ZXGao97/VPMR/blob/main/Readme.pdf)

The following example demonstrates how to use the newly added module: 

![](figure/VP-WBT.png)

- **`T`**: The truncation interval $[0, T]$.
- **`Weight`**: Choose the desired weight function.
- **`Weight Setting`**: Configure parameters for the selected weight function:
  - **`1`**: A trivial weight function $w(t) = 1$, corresponding to Time-limited Balanced Truncation.
  - **`(x + a).^(-1/2)`**: A weight function $w(t) = \frac{1}{\sqrt{t + a}}$. Specify the value of $a$ in `Weight Setting`.
  - **`exp(-s.*x)`**: An exponential weight function $w(t) = e^{-st}$. Specify the value of $s$ in `Weight Setting`.
  - **`Custom`**: Define your own weight function. Enter the expression for $w(t)$ in `Weight Setting` using MATLAB's standard function definition syntax.
- **`MR`**: Output the MR terms. You may enter a single integer to obtain SOE of a specific number of terms, or a vector of integers to compute SOE of multiple numbers of terms in one run.


# BSA-WBT
In our paper, BSA is used for high precision SOE approximation of various inverse power kernels, including Coulomb kernel. 
It is a quadrature method to obtain SOE/SOG as following:

$$
r^{-\alpha} = \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\int_{-\infty}^{\infty}b^{\alpha x}e^{-b^x r} dx \approx \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\sum_{\ell=-\infty}^{+\infty}b^{\alpha \ell}e^{-b^\ell r} \approx \frac{\log(b)}{\Gamma ( \alpha)}\displaystyle\sum_{\ell=N}^{M}b^{\alpha \ell}e^{-b^\ell r}
$$

The main executable script `BSA_WBT.m`, located in the `src` folder, implements the BSA method and can be used to reproduce the results presented in the paper.
