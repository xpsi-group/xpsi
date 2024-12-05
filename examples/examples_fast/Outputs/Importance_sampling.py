#!/usr/bin/env python
# coding: utf-8

#    # <center> This notebook will walk you through how to perform importance sampling using XPSI

#
# ```markdown
# # Performing Importance Sampling with X-PSI
#
# Importance sampling can be performed using X-PSI under the following circumstances:
#
# 1. # Model Update:
#    - The model itself has changed, which implies a modification of the likelihood function.
#
# 2. # Change in Parameter Priors:
#    - The prior for one or more parameters has been modified.
#    - However, the new priors must still lie within the hypervolume of the original prior distribution.
#
# 3. # Both Changes:
#    - Both the prior and the model have been updated simultaneously.
# ```

# # <center>  1. Model Update

# In the Module forder, we made a file called : main_IS_likelihood.py where we change the num_energies  in the likelihood class from 64 to 256
# This implies that our model has changed and the likelihood computation is more accurate

# In[8]:


import xpsi
from xpsi.Sample import importance


# In[1]:


# Let's import the old and the new models
from Modules import main as old_model
from Modules import main_IS_likelihood as new_model


# In[6]:


# Just me doing some lazy things with Serena script
def derive_names(model):
    out = str(model.likelihood)
    outA = out.split('\n')
    outA = outA[2:-1]
    names = []
    for i in range(len(outA)):
        a = outA[i]
        ind_end = a.find(':')
        if a[:ind_end].find(' ')<0:
            names.append(a[:ind_end])
    return names

names = derive_names(old_model)


# In[7]:


names


# In[13]:

import time

start = time.time()
# Doing importance Sampling
# Doing importance Sampling
importance(new_model.likelihood,
                  old_model.likelihood,
                  'Outputs/ST_live_1000_eff_0.3_seed42',
                  names = names,
                  likelihood_change = True,
                  prior_change = False,
                  weight_threshold=1.0e-30,
                  overwrite=True)

print(time.time()-start)
