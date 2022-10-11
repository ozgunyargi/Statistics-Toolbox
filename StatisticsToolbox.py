# -*- coding: utf-8 -*-
"""
Özgün Yargı*

##ANOVA
---

###ANOVA1 partinon TSS

This function takes a dataframe as an input as returns $SS_{total}, SS_w$, and $SS_b$. To implement this equation 2.2 is used
"""

def ANOVA1_partition_TSS (df):

  SStotal_ = 0
  SSw_ = 0
  SSb_ = 0
  X_bar = sum(df.sum())/(len(df.index)*len(df.columns)-df.isna().sum().sum()) # Calculate the global mean
  Xi_bar = df.mean() # Calculate the column's means
  
  for indx, i in enumerate(df.columns):

    ni = len(df[i])-sum(df[i].isna())
    Xi_b = Xi_bar[indx]
    SSb_ += ni*(Xi_b-X_bar)**2
    
    for Xij in df[i]:

      if np.isnan(Xij) == False:

        SStotal_ += (Xij-X_bar)**2
        SSw_ += (Xij-Xi_b)**2

  return SStotal_, SSw_, SSb_

"""### ANOVA1 test equality

This function takes the dataframe and significance level as an input and returns an ANOVA table with addion of p_value and desicion. This function does not require same number of data in each column. So that, it can even work on collections with different sample sizes. I used
"""

def ANOVA1_test_equality(df, alpha):

  SSt_, SSw_, SSb_ = ANOVA1_partition_TSS (df)
  I = len(df.columns) #Number of columns
  n = (len(df.index)*len(df.columns)-df.isna().sum().sum()) #find the total number of data (subtracted the nan values)
  Ms_b = SSb_/(I-1)
  Ms_w = SSw_/(n-I)
  F = Ms_b/Ms_w
  F_crit = stats.f.ppf(1-alpha, (I-1), (n-I))
  p_val = stats.f.sf(F, (I-1), (n-I))
  test = ""

  if p_val < alpha:

    test = "Reject"

  else:

    test = "Don't Reject"

  table = {"Source": ["Between Groups", "Within Groups", "Total"],
           "df": [I-1, n-I, n-1],
           "SS": [Ms_b, Ms_w, np.nan],
           "F": [F, np.nan, np.nan],
           "F_crit": [F_crit, np.nan, np.nan],
           "p_value": [p_val, np.nan, np.nan],
           "Desicion": [test, np.nan, np.nan]
           }

  df_table = pd.DataFrame(table)
  df_table = df_table.set_index("Source")

  return df_table

"""###ANOVA1 is contrast

This function takes a vector input and returns a boolean value whether given vector is contrast or not. Returns true if it is contrast, false otherwise. I used the definition 2.1 to implement this.
"""

def ANOVA1_is_contrast(c):

  total = sum(c)

  if total == 0:

    return True

  else:

    return False

"""###ANOVA1 is orthogonal

This function takes two vectors u, v and data list (n) as an input and return whether given 2 vector is orthogonal or not. To implement this, definition 2.2 and 2.3 is used. This function also returns an error if one of those given vectors is not contrast.
"""

def ANOVA1_is_orthogonal(ns, u, v):

  flag = True

  if ANOVA1_is_contrast(c1) == False:

    print("Warning! First linear combination is not a contrast")
    flag = False
    return False

  if ANOVA1_is_contrast(c2) == False: 

    print("Warning! First linear combination is not a contrast")
    flag = False
    return False

  uv = 0

  for i in range(len(c1)):

    uv +=  u[i]*v[i]/ns[i] # From definition 2.2

  if (uv == 0):

    return True

  else: 

    return False

"""###Bonferroni correction

This function takes alpha and number of test as an input and returns the corrected alpha version by using Bonferroni correction. To implement this, theorem 2.3 is used.
"""

def Bonferroni_correction (alph, m_):

  corrected_alpha = alph/m_

  return corrected_alpha

"""###Sidak Correction

This function takes alpha and number of test as an input and returns the corrected alpha version by using Sidak correction. To implement this, subtitle B.2.2. Sidak Correction is used
"""

def Sidak_correction (alpha, m_):

  corrected_alpha = 1-(1-alpha)**(1/m)

  return corrected_alpha

"""### is_pairwise

This is an additional function that takes the contrast matrix and check whether given matrix can be used as as pairwise matrix or not. It returns true if it is useable, otherwise false.
"""

def is_pairwise (C):

  for test in C:

    size = 0

    for element in test:

      if element != 0:

        size += 1

    if size != 2:

      return False

    else:

      return True

"""###ANOVA1 CI linear combs

This function take a dataframe, alpha, matrix C, and type of method as an input and returns Simultaneous Confidence Intervals. You can either pick Tukey or Scheffe as a method. FWER is controlled while creating simultaneous confidence intervals by using Scheffe or Bonferroni Correction depending on the matrix C:



* If Scheffe is chosen 
  * If given matrix C is both contrast and orthogonal, apply Scheffe by using Sidak Correction 2.8
  * If given matrix C is contrast but not orthogonal, apply Scheffe by using Bonferroni Correction 2.8
  * If given matrix C is not contrast as well as orthogonal, apply Scheffe by using Bonferroni Correction and Theorem 2.7

* If Tukey is chosen
  * Check whether given matrix C contains pairwise tests or not.
    * If not, return an error
  * If given matrix C is orthogonal, apply Tukey with Bonferroni Correction and Theorem 2.11 (if given dataset does not contain same amount of sample in each feature, then $n_*$ is taken as the minimum sample size)
  * If given matrix C is not orthogonal but contrast, apply Tukey with Sidak Correction and Theorem 2.11
"""

def ANOVA1_CI_linear_combs (df, alpha, C, method):

  CI = []
  SSt_, SSw_, SSb_ = ANOVA1_partition_TSS(df)
  I = len(df.columns)
  ni = [len(df[i])-sum(df[i].isna()) for i in df.columns]
  n_ = len(df.index)*len(df.columns)-df.isna().sum().sum()
  m_ = sum([a for a in range(len(df.columns))])
  means = list(df.sum().values/n_)

#Scheffe
#------------------------------------------------------
  if method == "Scheffe":

    flag = True
    
# Check Contrast
#------------------------------------------------------
    for row in C:

      if ANOVA1_is_contrast(row) == False:
        flag = False
        break

#------------------------------------------------------
# Check orthogonality
#------------------------------------------------------
    if flag:

      flag1 = True

      for indx, row in enumerate(C):     
        a = [i for i in range(indx+1, C.shape[0])]      
        for j in a:         
          if ANOVA1_is_orthogonal(ni, row.tolist(), C[j].tolist()) == False:           
            flag1 = False
            break
        if flag1 == False:         
          break

      Constant = sqrt(SSw_/(n_-I)*sum([row[i]**2/ni[i] for i in range(len(ni))]))
  
#------------------------------------------------------
# Scheffe with Sidak Correction
#------------------------------------------------------
      if flag1:
        alph_s = Sidak_correction(alpha, m_)
        f = stats.f.ppf(1-alph_s, (I-1), (n_-I))
        M = sqrt((I-1)*f)

        for row in C:
          tot_mean = sum([row[i]*means[i] for i in range(len(row))])
          CI.append([tot_mean-M*Constant, tot_mean+M*Constant])

        return CI

#------------------------------------------------------
# Scheffe with Bonferroni Correction
#------------------------------------------------------
      else:
        alph_b = Bonferroni_correction(alpha, m_)
        f = stats.f.ppf(1-alph_b, (I-1), (n_-I))
        M = sqrt((I-1)*f)

        for row in C:
          tot_mean = sum([row[i]*means[i] for i in range(len(row))])
          CI.append([tot_mean-M*Constant, tot_mean+M*Constant])

        return CI

#------------------------------------------------------
# Scheffe with not contrast
#------------------------------------------------------
    else:
      alph_b = Bonferroni_correction(alpha, m_)
      f = stats.f.ppf(1-alph_b, (I), (n_-I))
      M = sqrt((I-1)*f)

      for row in C:
        tot_mean = sum([row[i]*means[i] for i in range(len(row))])
        CI.append([tot_mean-M*Constant, tot_mean+M*Constant])

      return CI

#------------------------------------------------------
#------------------------------------------------------

# Tukey
#------------------------------------------------------
  elif method == "Tukey":

# Check Pairwise
#------------------------------------------------------
      if is_pairwise(C):

#------------------------------------------------------
# Check Orthogonality
#------------------------------------------------------
        flag = True

        for indx, row in enumerate(C):
          a = [i for i in range(indx+1, C.shape[0])]    
          for j in a:          
            if ANOVA1_is_orthogonal(ni, row.tolist(), C[j].tolist()) == False:           
              flag = False
              break
          if flag == False:         
            break

#------------------------------------------------------
# Tukey with Bonferroni Correction
#------------------------------------------------------
      if flag:
           alph_b = Bonferroni_correction(alpha, m_)
           q = qsturng(1-alph_b, I, n_-I)
           Constant = 1/sqrt(2)*sqrt(SSw*2/((n_-I)*min(ni)))
           weight_locs = []

           for test in C:
             weights = []
             for weight in test:
               if weight != 0:
                 weights.append(weight)
             weight_locs.append(weights)

           for test in weight_locs:
             mean_diff = abs(means[test[0]]-means[test[1]])
             CI.append((mean_diff-q*Constant, mean_diff+q*Constant))

           return CI

#------------------------------------------------------
# Tukey with Sidak Correction
#------------------------------------------------------
      else:
        alph_s = Sidak_correction(alpha, m_)
        q = qsturng(1-alph_s, I, n_-I)
        Constant = 1/sqrt(2)*sqrt(SSw*2/((n_-I)*min(ni)))

        weight_locs = []
        for test in C:
          weights = []
          counter = 0 
          for weight in test:
            if weight != 0:
              weights.append(counter)
            counter += 1
          weight_locs.append(weights)

        for test in weight_locs:

          mean_diff = abs(means[test[0]]-means[test[1]])
          CI.append((mean_diff-q*Constant, mean_diff+q*Constant))

        return CI

#------------------------------------------------------
#------------------------------------------------------



#  elif method == "best":

"""###ANOVA1 test linear combs

This function take a dataframe, alpha, matrix C, and type of method as an input and returns Simultaneous Hypothesis Test results. You can either pick Tukey or Scheffe as a method. FWER is controlled while testing Simultaneous Hypothesis Tests by using Scheffe or Bonferroni Correction depending on the matrix C:



* If Scheffe is chosen 
  * If given matrix C is both contrast and orthogonal, apply Scheffe by using Sidak Correction and 2.14 with Exercise 2.20-c
  * If given matrix C is contrast but not orthogonal, apply Scheffe by using Bonferroni Correction and 2.14 with Exercise 2.20-c
  * If given matrix C is not contrast as well as orthogonal, apply Scheffe by using Bonferroni Correction and 2.14 with Exercise 2.20-c

* If Tukey is chosen
  * Check whether given matrix C contains pairwise tests or not.
    * If not, return an error
  * If given matrix C is orthogonal, apply Tukey with Bonferroni Correction and Theorem 2.11 (if given dataset does not contain same amount of sample in each feature, then $n_*$ is taken as the minimum sample size)
  * If given matrix C is not orthogonal but contrast, apply Tukey with Sidak Correction and Theorem 2.11
"""

def ANOVA1_test_linear_combs (df, alpha, C, method):

  result = []
  SSt_, SSw_, SSb_ = ANOVA1_partition_TSS(df)
  I = len(df.columns)
  ni = [len(df[i])-sum(df[i].isna()) for i in df.columns]
  n_ = len(df.index)*len(df.columns)-df.isna().sum().sum()
  m_ = sum([a for a in range(len(df.columns))])
  means = list(df.sum().values/n_)

#Scheffe
#------------------------------------------------------
  if method == "Scheffe":

    flag = True
    
# Check Contrast
#------------------------------------------------------
    for row in C:

      if ANOVA1_is_contrast(row) == False:
        flag = False
        break

#------------------------------------------------------
# Check orthogonality
#------------------------------------------------------
    if flag:

      flag1 = True

      for indx, row in enumerate(C):     
        a = [i for i in range(indx+1, C.shape[0])]      
        for j in a:         
          if ANOVA1_is_orthogonal(ni, row.tolist(), C[j].tolist()) == False:           
            flag1 = False
            break
        if flag1 == False:         
          break

      Constant = sqrt(SSw_/(n_-I)*sum([row[i]**2/ni[i] for i in range(len(ni))]))
  
#------------------------------------------------------
# Scheffe with Sidak Correction
#------------------------------------------------------
      if flag1:
        alph_s = Sidak_correction(alpha, m_)
        f = stats.f.ppf(1-alph_s, (I-1), (n_-I))
        M = sqrt((I-1)*f)

        for row in C:
          ghost_res = []
          tot_mean = sum([row[i]*means[i] for i in range(len(row))])
          CI = [tot_mean-M*Constant, tot_mean+M*Constant]

          if CI[0] < tot_mean and tot_mean < CI[1]:
            ghost_res.append("Reject")

          else:
            ghost_res.append("Don't Reject")

          result.append(ghost_res)

        return result

#------------------------------------------------------
# Scheffe with Bonferroni Correction
#------------------------------------------------------
      else:
        alph_b = Bonferroni_correction(alpha, m_)
        f = stats.f.ppf(1-alph_b, (I-1), (n_-I))
        M = sqrt((I-1)*f)

        for row in C:
          ghost_res = []
          tot_mean = sum([row[i]*means[i] for i in range(len(row))])
          CI = [tot_mean-M*Constant, tot_mean+M*Constant]
        
          if CI[0] < tot_mean and tot_mean < CI[1]:
            ghost_res.append("Reject")

          else:
            ghost_res.append("Don't Reject")

          result.append(ghost_res)

        return result

#------------------------------------------------------
# Scheffe with not contrast
#------------------------------------------------------
    else:
      alph_b = Bonferroni_correction(alpha, m_)
      f = stats.f.ppf(1-alph_b, (I), (n_-I))
      M = sqrt((I-1)*f)

      for row in C:
        tot_mean = sum([row[i]*means[i] for i in range(len(row))])
        CI = [tot_mean-M*Constant, tot_mean+M*Constant]
        
        if CI[0] < tot_mean and tot_mean < CI[1]:
            ghost_res.append("Reject")

        else:
            ghost_res.append("Don't Reject")

        result.append(ghost_res)

      return result

#------------------------------------------------------
#------------------------------------------------------

# Tukey
#------------------------------------------------------
  elif method == "Tukey":

# Check Pairwise
#------------------------------------------------------
      if is_pairwise(C):

#------------------------------------------------------

# Check Orthogonality
#------------------------------------------------------
        flag = True

        for indx, row in enumerate(C):
          a = [i for i in range(indx+1, C.shape[0])]    
          for j in a:          
            if ANOVA1_is_orthogonal(ni, row.tolist(), C[j].tolist()) == False:           
              flag = False
              break
          if flag == False:         
            break

#------------------------------------------------------
# Tukey with Bonferroni Correction
#------------------------------------------------------
      if flag:
           alph_b = Bonferroni_correction(alpha, m_)
           q = qsturng(1-alph_b, I, n_-I)
           Constant = 1/sqrt(2)*sqrt(SSw*2/((n_-I)*min(ni)))
           weight_locs = []

           for test in C:
             weights = []
             for weight in test:
               if weight != 0:
                 weights.append(weight)
             weight_locs.append(weights)

           for test in weight_locs:
             mean_diff = abs(means[test[0]]-means[test[1]])
             Critical_Region = q*Constant

             if mean_diff<Critical_Region:
               result.append("Reject")

             else:
               result.append("Don't Reject")

           return result

#------------------------------------------------------
# Tukey with Sidak Correction
#------------------------------------------------------
      else:
        alph_s = Sidak_correction(alpha, m_)
        q = qsturng(1-alph_s, I, n_-I)
        Constant = 1/sqrt(2)*sqrt(SSw*2/((n_-I)*min(ni)))

        weight_locs = []
        for test in C:
          weights = []
          counter = 0 
          for weight in test:
            if weight != 0:
              weights.append(counter)
            counter += 1
          weight_locs.append(weights)

        for test in weight_locs:
             mean_diff = abs(means[test[0]]-means[test[1]])
             Critical_Region = q*Constant

             if mean_diff<Critical_Region:
               result.append("Reject")

             else:
               result.append("Don't Reject")

        return result

#------------------------------------------------------
#------------------------------------------------------

"""###ANOVA2 partition TSS

This function takes samples as an numpy array and returns $SS_{total}, SS_A, SS_B, SS_{AB}, SS_E$ as an output. The Outputs were calculated according to C.2 Inference / Exercise 2.33-a
"""

def ANOVA2_partition_TSS(Xijk):
  K = len(Xijk[0][0])
  I = len(Xijk)
  J = len(Xijk[0])

  X_bar = np.array(Xijk.sum()).mean()
  Xij_bar = np.array([sum(b)/len(b) for a in Xijk for b in a]).reshape(I,J)
  Xi_bar = [np.array(Xijk[i].sum()).mean() for i in range(I)]
  Xj_bar = [np.array(Xijk[:, i].sum()).mean() for i in range(J)]

  SSA_ = I*K*sum([(Xi_b-X_bar)**2 for Xi_b in Xi_bar])
  SSB_ = I*K*sum([(Xj_b-X_bar)**2 for Xj_b in Xj_bar])
  SSAB_ = K*sum([(Xij_bar[i][j]-Xi_bar[i]-Xj_bar[j]+X_bar)**2 for i in range(I) for j in range(J)])
  SSE_ = sum([(Xijk[i][j][k]-Xij_bar[i][j])**2 for i in range(I) for j in range(J) for k in range(K)])
  SStotal_ = SSA_+SSB_+SSAB_+SSE_

  return SStotal_, SSA_, SSB_, SSAB_, SSE_

"""###ANOVA2 MLE

This function takes samples as an numpy array and returns  $\hat{\mu},\hat{a_i},\hat{b_j},\hat{\delta_{ij}}$  as an output. The Outputs were calculated according to C.2 Inference / Exercise 2.31
"""

def ANOVA2_MLE (Xijk):

  K = len(Xijk[0][0])
  I = len(Xijk)
  J = len(Xijk[0])

  X_bar = np.array(Xijk.sum()).mean()
  Xij_bar = np.array([sum(b)/len(b) for a in Xijk for b in a]).reshape(I,J)
  Xi_bar = [np.array(Xijk[i].sum()).mean() for i in range(I)]
  Xj_bar = [np.array(Xijk[:, i].sum()).mean() for i in range(J)]

  mu_hat_ = X_bar
  ai_hat_ = [Xi_b - X_bar for Xi_b in Xi_bar]
  bj_hat_ = [Xj_b - X_bar for Xj_b in Xj_bar]
  gammaij_hat_ = [Xij_bar[i][j] - Xj_bar[j] - Xi_bar[i] + X_bar for i in range(I) for j in range(J)]

  return mu_hat_, ai_hat_, bj_hat_, gammaij_hat_

"""###ANOVA2 test equality

This function takes numpy array (samples), alpha, and desicion to return Two Way ANOVA Table. According to the desicion 

* If A is chosen:
  * Returns you Two Way ANOVA Table which only contains dof, SS, MS scores of $SS_A$ and $SS_E$ with F score and desicion. The table was created by using Exercise 2.35-b
* If B is chosen:
  * Returns you Two Way ANOVA Table which only contains dof, SS, MS scores of $SS_B$ and $SS_E$ with F score and desicion. The table was created by using Exercise 2.35-b
* If AB is chosen:
  * Returns you Two Way ANOVA Table which only contains dof, SS, MS scores of $SS_{AB}$ and $SS_E$ with F score and desicion. The table was created by using Exercise 2.35-b
"""

def ANOVA2_test_equality (Xijk, alpha, desicion):

  K = len(Xijk[0][0])
  I = len(Xijk)
  J = len(Xijk[0])

  SStotal_, SSA_, SSB_, SSAB_, SSE_ = ANOVA2_partition_TSS(Xijk)

  dof_E = I*J*(K-1)
  MSE = SSE_/dof_E

  if desicion == "A":

    dof = I-1
    SS = SSA_
    MS = SS/dof
    F = MS/MSE
    label = "A"
    f =  stats.f.ppf(1-alpha, dof, dof_E)

    if F > f:
      test = "Reject"
    else:
      test = "Don't Reject"

  elif desicion == "B":

    dof = J-1
    SS = SSB_ 
    MS = SS/dof
    F = MS/MSE
    label = "B"
    f =  stats.f.ppf(1-alpha, dof, dof_E)

    if F > f:
      test = "Reject"
    else:
      test = "Don't Reject"

  elif desicion == "AB":

    dof = (I-1)*(J-1)
    SS = SSAB_   
    MS = SS/dof
    F = MS/MSE
    label = "AB"
    f =  stats.f.ppf(1-alpha, dof, dof_E)

    if F > f:
      test = "Reject"
    else:
      test = "Don't Reject"

  table = {"Source": [label, "within"],
           "degrees of freedom": [dof, dof_E],
           "SS": [SS, SSE_],
           "MS": [MS, MSE],
           "F": [F, np.nan],
           "Test":[test, np.nan] 
           }

  df_table = pd.DataFrame(table)
  df_table = df_table.set_index("Source")

  return df_table

"""##Linear Regression

###Mult LR Least squares

This function takes input and output samples as a numpy array and returns $\hat{\beta}, \hat{\sigma^2}, \hat{Se^2}$ values. This function is implemented by using Exercise 3.29
"""

def Mult_LR_Least_squares(X, Y):

  k = X.shape[1]-1
  n = X.shape[0]

  XTX = X.T.dot(X)
  B_hat_ = np.linalg.inv(XTX).dot(X.T).dot(Y)
  var_hat_ = 1/n*(Y-X.dot(B_hat_)).T.dot(Y-X.dot(B_hat_))
  Se_2_ = 1/(n-k)*(Y-X.dot(B_hat_)).T.dot(Y-X.dot(B_hat_))

  return B_hat_, var_hat_, Se_2_

"""###Mult LR partition TSS

This function takes input and output samples as a numpy array and returns TSS, ResidSs, RegSS by using Corollary 3.5
"""

def Mult_LR_partition_TSS (X, Y):

  B_hat_, var_hat_, Se_2_ = Mult_LR_Least_squares(X, Y)
  Y_hat_ = X.dot(B_hat_)
  Y_bar_ = Y.mean()

  XTX = X.T.dot(X)

  TSS_ = sum([(y-Y_bar_)**2 for y in Y_hat_])
  ResidSS_ = sum([(Y.tolist()[i][0]-Y_hat_[i])**2 for i in range(len(Y_hat_))])
  RegSS_ = TSS_-ResidSS_

  return TSS_, ResidSS_, RegSS_

"""###Mult norm LR simul CI

This function takes input and output samples as a numpy array and alpha as a float variable and returns a vector of simultaneous confidence interval. FWER is preserved by using Bonferroni Correction. Equation 3.47 is used to create confidence intervals
"""

def Mult_norm_LR_simul_CI (X, Y, alpha):

  k = X.shape[1]-1
  n = X.shape[0]
  Fwer_alp = Bonferroni_correction(alpha, k)

  t=stats.t.ppf(1-(Fwer_alp/(2)),n-k-1)
  B_hat_, var_, Se_2_ = Mult_LR_Least_squares(X, Y)
  Se_ = sqrt(Se_2_)

  CI = []

  for i in range(k):

    Confidence_interval = ( (B_hat_[i] - (t*Se_*sqrt(np.linalg.inv(X.T.dot(X))[i , i])))[0], (B_hat_[i] + t*Se_*sqrt(np.linalg.inv(X.T.dot(X))[i , i]))[0] )
    CI.append(Confidence_interval)

  return CI

"""###Mult norm LR CR

This function takes input and output samples as a numpy array and C as a full rank matrix and alpha as a float value and returns parameters to create a Confidence Region. Equation 3.61 is used to implement this function.
"""

def Mult_norm_LR_CR (X,Y,C,alpha):

  r = np.linalg.matrix_rank(C)
  k = X.shape[1]-1
  n = X.shape[0]
  
  f = stats.f.ppf(1-alpha, r, n-k-1)

  XTX = X.T.dot(X)
  B_hat_, var_, Se_2_ = Mult_LR_Least_squares(X, Y)
  parameter1 = (C.dot(B_hat_)).T
  parameter2 = np.linalg.inv(C.dot(np.linalg.inv((XTX))).dot(C.T))
  parameter3 = C.dot(B_hat_)
  rhs = r*Se_2_*f

  return parameter1, parameter2, parameter3, rhs

"""###Mult norm LR is in CR

This function takes input and output values as a numpy array and C as a full rank matrix and $c_0$ vector and alpha as a float variable and returns whether given $c_0$ is inside of Confidence Region or not. Returns True if it is in Critical Region, otherwise False. Equation 3.61 is used to implement this function.
"""

def Multnorm_LR_is_in_CR (X,Y,C,c0,alpha):
  
  r = np.linalg.matrix_rank(C)
  k = X.shape[1]-1
  n = X.shape[0]

  f = stats.f.ppf(1-alpha, r, n-k-1)

  XTX = X.T.dot(X)
  B_hat_, var_, Se_2_ = Mult_LR_Least_squares(X, Y)
  parameter1 = (C.dot(B_hat_)-c0).T
  parameter2 = np.linalg.inv(C.dot(np.linalg.inv((XTX))).dot(C.T))
  parameter3 = C.dot(B_hat_)-c0
  rhs = r*Se_2_*f
  test = parameter1.dot(parameter2).dot(parameter3)

  if test[0] > rhs:

    return False

  else:

    return True

"""###Mult norm LR test general

This function takes input and output values as a numpy array and C as a full rank matrix and $c_0$ vector and alpha as a float variable and returns the test result. Equation 3.62 is used to implement this function.
"""

def Mult_norm_LR_test_general (X, Y, C, c0, alpha):

  r = np.linalg.matrix_rank(C)
  k = X.shape[1]-1
  n = X.shape[0]

  f = stats.f.ppf(1-alpha, r, n-k-1)

  XTX = X.T.dot(X)
  B_hat_, var_, Se_2_ = Mult_LR_Least_squares(X, Y)
  parameter1 = (C.dot(B_hat_)-c0).T
  parameter2 = np.linalg.inv(C.dot(np.linalg.inv((XTX))).dot(C.T))
  parameter3 = C.dot(B_hat_)-c0
  rhs = r*Se_2_*f
  test = parameter1.dot(parameter2).dot(parameter3)

  if test[0] > r*Se_2_*f:

    return "Reject"

  else:

    return "Don't Reject"

"""###Mult norm LR test comp

This function takes input and output as a numpy array, alpha and j as a float variable and returns hypothesis test result. j represents the part of $\beta$ where that part ($\beta_1, \beta_2,...\beta_j,$) is the one to be considered ($\beta_0$ is not considered). To implement this function, Equation 3.63 is used.
"""

def Mult_norm_LR_test_comp(X, Y, alpha, j):

  ghost = np.zeros((j, 1))
  ghost1 = np.diag([1 for _ in range(j)])
  
  C = np.concatenate((ghost,ghost1),axis=1)
  c0 = np.array([[0] for _ in range(C.shape[0])])

  return Mult_norm_LR_test_general (X, Y, C, c0, alpha)

"""###Mult norm LR test linear reg

This function takes input and output as a numpy array, alpha as a float variable and returns hypothesis test that is considerinf all $\beta s$. To implement this function, Equation 3.63 is used.
"""

def Mult_norm_LR_test_linear_reg (X,Y,alpha):

  return Mult_norm_LR_test_comp(X, Y, alpha, X.shape[1]-1)

"""###Mult_norm_LR_pred_CI

This function takes input and output as a numpy array, matrix D, alpha with a float variable and a metod. It returns a list of simultaneos confidence intervals. To implement this task, Theorem 3.19 is used. You need to choose either Bonferroni or Scheffe
"""

def Mult_norm_LR_pred_CI(X,Y,D,alpha,method):

    beta_hat,var,Se_2 = Mult_LR_Least_squares(X,Y)

    k = X.shape[1]
    n = X.shape[0] 
    m = D.shape[1]

    XTX = np.linalg.inv(X.T.dot(X))
    CI = []

    if method =="Bonferroni":
      fwer_alp = Bonferroni_correction(alpha, m)
      for x0 in D:
          t = stats.t.ppf(1- (fwer_alp/(2*m)), n-2) * sqrt(Se_2 *(x0.T.dot(XTX).dot(x0)))
          pred = x0.T.dot(beta_hat)
          CI.append((pred-t,pred+t))

      return CI

    elif method =="Scheffe":
      fwer_alp = Sidak_correction(alpha, m)
      for x0 in D:
        f = sqrt(k*stats.f.ppf(1-fwer_alp,k,n-k)) * sqrt(Se_2) * sqrt((x0.T.dot(XTX).dot(x0))) 
        pred = x0.T.dot(beta_hat)
        CI.append((pred-f,pred+f))

      return CI

"""##Main Script
---

###Import Libraries
"""

import numpy as np
import pandas as pd
from scipy import stats
from math import sqrt
from statsmodels.stats.libqsturng import qsturng

"""###Scirpt

####One Way Anova
"""

#This dataset contains equal sized samples

A = [25, 30, 28, 36, 29]
B = [45, 55, 29, 56, 40]
C = [30, 29, 33, 37, 27]
D = [54, 60, 51, 62, 73]

data = {"A": A, "B": B, "C": C, "D": D}

# This dataset contains non-equal sized samples

A = [42, 53, 49, 53, 43, 44, 45, 52, 54]
B = [69, 54, 58, 64, 64, 55, 56, np.nan, np.nan]
C = [35, 40, 53, 42, 50, 39, 55, 39, 40]

data = {"A": A, "B": B, "C": C}

# Create a dataset for experiments (ANOVA)
#------------------------------------------#

df_main = pd.DataFrame(data, columns = data.keys(), )

df_main

# ANOVA1_partition_TSS

SStotal, SSw, SSb = ANOVA1_partition_TSS (df_main)
print("SStotal: {}\nSSw: {}\nSSb: {}".format(SStotal, SSw, SSb))

# Create One-Way ANOVA Table

alph = 0.1
ANOVA1_test_equality(df_main,alph)

c1 = [-1, 0, 1]
c2 = [0, 1, -1]
n = [len(df_main[i])-sum(df_main[i].isna()) for i in df_main.columns]

# Check Contrasts of Linear Combinations

c1_ = ANOVA1_is_contrast(c1)
c2_ = ANOVA1_is_contrast(c2)

print ("c1 is contrast: {}\nc2 is contrast: {}".format(c1_, c2_))

# Check orthogonality

ANOVA1_is_orthogonal(n, c1, c2)

# Find m

m = len([c1, c2])

# Apply Bonferroni Correction

corr_alph = Bonferroni_correction(alph, m)

print("Corrected alpha is:", corr_alph, "\nFWER is:", 1-(1-corr_alph)**m)

# Apply Sidak Correction

corr_alph = Sidak_correction(alph, m)

print("Corrected alpha is:", corr_alph, "\nFWER is:", 1-(1-corr_alph)**m)

C_ = np.array([c1, c2])

print(C_)

Confi = ANOVA1_CI_linear_combs(df_main, alph, C_, "Scheffe")

print(Confi)

print(ANOVA1_test_linear_combs(df_main, alph, C_, "Scheffe"))

"""####Two Way Anova"""

#Dataset for Two Way ANOVA

dict_ = {"Compost 1":{"soil A": [8,7,7,6], "soil B": [13, 15, 14, 15]}, "Compost 2":{"soil A": [5,6,6,4], "soil B": [12, 14, 13, 14]}}

df_anova2 = pd.DataFrame(dict_)

df_anova2.head()

SStotal, SSA, SSB, SSAB, SSE = ANOVA2_partition_TSS(df_anova2.values)

print("SStotal: {}\nSSA: {}\nSSB: {}\nSSAB: {}\nSSE: {}".format(SStotal, SSA, SSB, SSAB, SSE))

mu_hat, ai_hat, bj_hat, gammaij_hat =ANOVA2_MLE(df_anova2.values)

print("mu_hat: {}\nai_hat: {}\nbj_hat: {}\ngammaij_hat: {}".format(mu_hat, ai_hat, bj_hat, gammaij_hat))

Two_Way_ANOVA_Table = ANOVA2_test_equality(df_anova2.values, 0.1, "AB")

Two_Way_ANOVA_Table.head()

"""####Linear Regression"""

# Dataset for Linear Regression

attend_lectures = [1, 0.5, 0.2, 0.4, 0.5, 0.7, 0.8, 0.9, 0.6, 0.1, 0, 0, 0.7, 0.8, 1]
homework = [0.25, 1, 0.5, 1, 1, 0.75, 1, 0.25, 0, 0, 1, 0.5, 0.25, 0.75, 1]
overall_grade = [60, 65, 40, 70, 65, 70, 85, 70, 44, 20, 40, 30, 50, 77, 90]

X_ = np.array([[1 for a in range(len(attend_lectures))], attend_lectures, homework]).T
Y_ = np.array([overall_grade]).T

Beta, var, Se_2 = Mult_LR_Least_squares(X_, Y_)

print("Beta: \n{}\nVariance: {}\nUnbiased-Variance: {}".format(Beta, var, Se_2))

TSS, ResidSS, RegSS = Mult_LR_partition_TSS(X_, Y_)

print("TSS: {}\nResidSS: {}\nRegSS: {}".format(TSS, ResidSS, RegSS))

print(Mult_norm_LR_simul_CI(X_, Y_, 0.1))

C_ = [[0,1,0], [0,0,1]]
C_ = np.array(C_)
c0_ = np.array([[0], [0]])

print(C_)

print(c0_)

print(Mult_norm_LR_CR(X_, Y_, C_, 0.1))

print(Multnorm_LR_is_in_CR (X_,Y_,C_,c0_,0.1))

print(Mult_norm_LR_test_general (X_,Y_,C_,c0_,0.1))

print(Mult_norm_LR_test_comp(X_,Y_,0.1, 2))

print(Mult_norm_LR_test_linear_reg(X_, Y_, 0.1))

D_ = np.array([[0,1,0],[0,0,1]])

print(Mult_norm_LR_pred_CI(X_,Y_,D_,0.1, "Scheffe"))