# Statistics-Toolbox

This reporsitory sums some of the most frequently used statistics tools for every to be used.

ONE WAY ANOVA FUNCTIONS
-----------------------
* def ANOVA1_partition_TSS (df)
	- This function takes a dataframe as an input as returns SStotal, SSw, and SSb. To implement this equation 2.2 is used.

* def ANOVA1_test_equality(df, alpha)
	- This function takes the dataframe and significance level as an input and returns an ANOVA table with addion of p_value and desicion. This function does not require same number of data in each column. So that, it can even work on collections with different sample sizes. I used

* def ANOVA1_is_contrast(c)
	- This function takes a vector input and returns a boolean value whether given vector is contrast or not. Returns true if it is contrast, false otherwise. I used the definition 2.1 to implement this.

* def ANOVA1_is_orthogonal(ns, u, v)
	- This function takes two vectors u, v and data list (n) as an input and return whether given 2 vector is orthogonal or not. To implement this, definition 2.2 and 2.3 is used. This function also returns an error if one of those given vectors is not contrast.

* def Bonferroni_correction (alph, m_)
	- This function takes alpha and number of test as an input and returns the corrected alpha version by using Bonferroni correction. To implement this, theorem 2.3 is used.

* def Sidak_correction (alpha, m_)
	- This function takes alpha and number of test as an input and returns the corrected alpha version by using Sidak correction. To implement this, subtitle B.2.2. Sidak Correction is used

* def is_pairwise (C):
	- This is an additional function that takes the contrast matrix and check whether given matrix can be used as as pairwise matrix or not. It returns true if it is useable, otherwise false.

* def ANOVA1_CI_linear_combs (df, alpha, C, method)
	- This function take a dataframe, alpha, matrix C, and type of method as an input and returns Simultaneous Confidence Intervals. You can either pick Tukey or Scheffe as a method. FWER is controlled while creating simultaneous confidence intervals by using Scheffe or Bonferroni Correction depending on the matrix C:
	* If Scheffe is chosen
	  * If given matrix C is both contrast and orthogonal, apply Scheffe by using Sidak Correction 2.8
	  * If given matrix C is contrast but not orthogonal, apply Scheffe by using Bonferroni Correction 2.8
	  * If given matrix C is not contrast as well as orthogonal, apply Scheffe by using Bonferroni Correction and Theorem 2.7
	* If Tukey is chosen
	  * Check whether given matrix C contains pairwise tests or not.
	    * If not, return an error
	  * If given matrix C is orthogonal, apply Tukey with Bonferroni Correction and Theorem 2.11 (if given dataset does not contain same amount of sample in each feature, then $n_*$ is taken as the minimum sample size)
	  * If given matrix C is not orthogonal but contrast, apply Tukey with Sidak Correction and Theorem 2.11

* def ANOVA1_test_linear_combs (df, alpha, C, method):
	- This function take a dataframe, alpha, matrix C, and type of method as an input and returns Simultaneous Hypothesis Test results. You can either pick Tukey or Scheffe as a method. FWER is controlled while testing Simultaneous Hypothesis Tests by using Scheffe or Bonferroni Correction depending on the matrix C:
	* If Scheffe is chosen
	  * If given matrix C is both contrast and orthogonal, apply Scheffe by using Sidak Correction and 2.14 with Exercise 2.20-c
	  * If given matrix C is contrast but not orthogonal, apply Scheffe by using Bonferroni Correction and 2.14 with Exercise 2.20-c
	  * If given matrix C is not contrast as well as orthogonal, apply Scheffe by using Bonferroni Correction and 2.14 with Exercise 2.20-c
	* If Tukey is chosen
	  * Check whether given matrix C contains pairwise tests or not.
	    * If not, return an error
	  * If given matrix C is orthogonal, apply Tukey with Bonferroni Correction and Theorem 2.11 (if given dataset does not contain same amount of sample in each feature, then n* is taken as the minimum sample size)
	  * If given matrix C is not orthogonal but contrast, apply Tukey with Sidak Correction and Theorem 2.11

TWO WAY ANOVA FUNCTIONS
-----------------------
* def ANOVA2_partition_TSS(Xijk)
	- This function takes samples as an numpy array and returns SStotal, SSA, SSB, SSAB, SSE as an output. The Outputs were calculated according to C.2 Inference / Exercise 2.33-a

* def ANOVA2_MLE (Xijk)
	- This function takes samples as an numpy array and returns  mu_hat,ai_hat,bj_hat,deltaij_hat  as an output. The Outputs were calculated according to C.2 Inference / Exercise 2.31

* def ANOVA2_test_equality (Xijk, alpha, desicion)
	- This function takes numpy array (samples), alpha, and desicion to return Two Way ANOVA Table. According to the desicion
	* If A is chosen:
	  * Returns you Two Way ANOVA Table which only contains dof, SS, MS scores of SSA and SSE with F score and desicion. The table was created by using Exercise 2.35-b
	* If B is chosen:
	  * Returns you Two Way ANOVA Table which only contains dof, SS, MS scores of SSB and SSE with F score and desicion. The table was created by using Exercise 2.35-b
	* If AB is chosen:
	  * Returns you Two Way ANOVA Table which only contains dof, SS, MS scores of SSAB and SSE with F score and desicion. The table was created by using Exercise 2.35-b

LINEAR REGRESSION FUNCTIONS
---------------------------
* def Mult_LR_Least_squares(X, Y)
	- This function takes input and output samples as a numpy array and returns beta_hat, sigma2_hat, Se2_hat values. This function is implemented by using Exercise 3.29

* def Mult_LR_partition_TSS (X, Y)
	- This function takes input and output samples as a numpy array and returns TSS, ResidSs, RegSS by using Corollary 3.5

* def Mult_norm_LR_simul_CI (X, Y, alpha)
	- This function takes input and output samples as a numpy array and alpha as a float variable and returns a vector of simultaneous confidence interval. FWER is preserved by using Bonferroni Correction. Equation 3.47 is used to create confidence intervals

* def Mult_norm_LR_CR (X,Y,C,alpha)
	- This function takes input and output samples as a numpy array and C as a full rank matrix and alpha as a float value and returns parameters to create a Confidence Region. Equation 3.61 is used to implement this function.

* def Multnorm_LR_is_in_CR (X,Y,C,c0,alpha)
	- This function takes input and output values as a numpy array and C as a full rank matrix and c0 vector and alpha as a float variable and returns whether given c0 is inside of Confidence Region or not. Returns True if it is in Critical Region, otherwise False. Equation 3.61 is used to implement this function.

* def Mult_norm_LR_test_general (X, Y, C, c0, alpha)
	- This function takes input and output values as a numpy array and C as a full rank matrix and c0 vector and alpha as a float variable and returns the test result. Equation 3.62 is used to implement this function.

* def Mult_norm_LR_test_comp(X, Y, alpha, j)
	- This function takes input and output as a numpy array, alpha and j as a float variable and returns hypothesis test result. j represents the part of $\beta$ where that part (beta_1, beta_2,...beta_j) is the one to be considered (beta_0 is not considered). To implement this function, Equation 3.63 is used.

* def Mult_norm_LR_test_linear_reg (X,Y,alpha)
	- This function takes input and output as a numpy array, alpha as a float variable and returns hypothesis test that is considerinf all $\beta s$. To implement this function, Equation 3.63 is used.

* def Mult_norm_LR_pred_CI(X,Y,D,alpha,method)
	- This function takes input and output as a numpy array, matrix D, alpha with a float variable and a metod. It returns a list of simultaneos confidence intervals. To implement this task, Theorem 3.19 is used. You need to choose either Bonferroni or Scheffe
  
MAIN SCRIPT
-----------
* This Seciton is separeted into 3. Each section has its own sample dataset. You may try any of them any time you want.
None of them will effect others. 

