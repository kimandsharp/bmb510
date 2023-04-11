# bmb510
Code for bmb510 data analysis class, taught to 1st year Ph.D students
in the Biochemistry and Molecular Biophysics program at the 
University of Pennsylvania.
Requires math, matplotlib, and numpy packages
Bayesian versions of most basics stats operations used in conventional statistics:
estimation of mean, std. dev of populations, difference in means (T-tests), Poisson, Binomial, Rank test
Linear Regression, etc, with Bayesian versions of a few more advanced cases: periodic data, survival data.
Written as a set of standalone Python programs, each designed to do one thing well (I hope!),
implementing 'exact' numerical integration for the posterior probability distributions, and providing
graphs as well as numbers.
Goals are: 
Provide a better alternative to usual null-hypothesis/significance testing way of teaching stats.
Act as a prequel to more computationally advanced Bayesian stats described by A. Gelman et al., 
"Bayesian Data Analysis", 3rd ed., John Kruschke "Doing Bayesian Data Analysis" as 
implemented in "PyMC3"
