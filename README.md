# Note: If your work uses this algorithm or makes improvements based on it, please be sure to cite this paper. Thank you for your cooperation.

------

# 注意：如果您的工作用到了本算法，或者基于本算法进行了改进，请您务必引用本论文，谢谢配合。

The core algorithm of paper " SODL-IR-FISTA: sparse online dictionary learning with iterative reduction FISTA for cone-beam X-ray luminescence computed tomography"

Xin Cao, Wenlong Tang, Huimin Gao, Yifan Wang, Yi Chen, Chengyi Gao, Fengjun Zhao, and Linzhi Su

Biomedical Optics Express. (2024).

Environments: Matlab R2023a

Reference: 

- IR-FISTA: [HomSolver/irfista at main · wangguojim/HomSolver · GitHub](https://github.com/wangguojim/HomSolver/tree/main/irfista)

- ODL: [GitHub - tiepvupsu/DICTOL: DICTOL - A Dictionary Learning Toolbox in Matlab and Python](https://github.com/tiepvupsu/DICTOL/tree/master)

Usage:

Consider having a system matrix A(m\*n) and surface optical data b(m\*1), we can obtain the result vector x(n\*1) in the following way:

```matlab
lambda = 0.001
[~, ~, x] = ODL(b, A, lambda, [], []); 
```
