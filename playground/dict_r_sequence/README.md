# Motivation:

According to Knuth, using linear Probing with filling rate α in a dictionary with M slots

On average:
- 1/2 ( 1 + 1/(1-α)) probes, if the element exists in the dictionary
- 1/2 ( 1 + 1/(1-α)**2 ) probes, if the element doesnt't exist in the dictionary

In our case, we have the liberty to add an element to a dictionary or discard it since we can generate a new one.


One might consider the following optimization:
- To search for an element inside the dicitonary, make at most O( 1/2 ( 1 + 1/(1-α)) ) propes.
- It is possible to exploit the discarding elements by excluding all elements that requires more than c0 * 1 + 1/(1-α))


The goal of this experiment is to determine c0 such that 0.95α elements are added to the dictionary. 


# Experiment Steps

For a given filling rate α, we would like to find the distribution of (nprobes, nelements that need nprobes to be inserted). Our initial guess this follows a Poisson Distribution. 


1. Determine the mean lamba of poisson distribution for given α
2. Plot the experimental distributions against the estimated poisson dist





