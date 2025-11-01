# Baron-Myerson-Using-Rocket-Science
This is set of self-contained Matlab implementation for numerical simulation of models on monopoly regulation, based on [OptimTraj](https://github.com/MatthewPeterKelly/OptimTraj). This is the only dependency of the project.
Without the support and commitment to open source from engineers like @MatthewPeterKelly, it would be intimidating for economic theorist to simulate a policy under mechanism design. But now, it only takes the effort of a single click. 

The basic function of this project is to simulate the policy of a simplified [Baron-Myerson](https://www.jstor.org/stable/1912769?seq=1) model *(consumer value and seller's cost have the same support of [0,1], with or without a fixed cost)* with welfare weight $\alpha \in [0,1]$.

The approach of this project is to transcribe a nonlinear pricing/tax problem into an optimal control. Basic understanding of the optimal control approach to mechanism design is essential. The main file is $\texttt{Baron}\underline{ }\texttt{Myerson.m}$. The dynamics is contained in $\texttt{IC.m}$. The objective function (integrand) is contained in $\texttt{Objective}\underline{ }\texttt{BM.m}$. Several kinds of path constraints are presented in $\texttt{pathConstraint.m}$, including a simple exercise of transfer cap, [Amador and Bagwell (2022)](https://econtheory.org/ojs/index.php/te/article/view/20221719), and [Wei and Zou (2025)](https://sites.google.com/view/dihanzou/research). I suggest that first-time user should focus on the no path-constraint case (by default) and study the behavior of **Baron-Myerson model**, especially for different welfare weights and distributions. 

As for the method part of $\texttt{Baron}\underline{ }\texttt{Myerson.m}$, I strongly suggest users take the $\texttt{trapezoid}$ method or $\texttt{hermiteSimpson}$ method, as they are more robust to discontinuities in $q'$ or even $q$ (demand policy). Higher-order methods such as $\texttt{chebyshev}$ tend to over-smooth the policy, leading to non-convergence. See [here](https://epubs.siam.org/doi/10.1137/16M1062569) for technical details.

The attached figures showcase the optimal policy in Baron-Myerson setting with linear demand and various cost settings for $\alpha = 0$ or $1$. See [here](https://drive.google.com/file/d/110uUldyFOqUwwcxvAnlAUYcNniI8MHJ_/view?usp=sharing) for a quick tutorial on Baron-Myerson model. Additional figures presenting other policy situations are contained in the folder named ``Figures.'' 

I am happy to chat if you have any interesting findings playing with this project. Feel free to shoot an email at $\texttt{dihanzou AT gmail DOT com}$.


**References**

Amador, Manuel and Kyle Bagwell (2022): ``Regulating a Monopolist with Uncertain Costs without Tranfers,'' *Theoretical Economics*, 17 (4): 1719-1760.

Baron, David P. and Roger B. Myerson (1982): ``Regulating a Monopolist with Unknown Costs,'' *Econometrica*, 50 (4), 911 - 930.

Wei, Jiaming and Dihan Zou (2025): ``Monopoly Regulation without Subsidy,'' Working Paper. 
