The core definitions of this library are that of statistical manifolds. A
statistical manifold is a manifold of probability distributions, such that
each point on the manifold is a probability distribution. In Goal, we can
create a normal distribution, for example, by writing

    x :: Standard :#: Normal
    x = fromList Normal [0,1]

where 0 is the mean and 1 is the variance. We can create a binomial distribution
with

    x :: Standard :#: Binomial
    x = fromList (Binomial 5) [0.5]

which corresponds to the distribution over 5 fair coin tosses. Where
possible, this library provides implementations for evaluating the
probability density function of a given distribution, generating samples,
and computing the maximum likelihood estimator.

Exponential families are also a core part of this library. An exponential
family is a kind of statistical manifold with a particular convex structure,
which comes with a pair of coordinate systems known as the mixture and
natural coordinates which are dual to each other. These coordinate systems
are not usually as intuitive as the standard coordinates, and so to work
with a distribution in natural coordinates for example, we would usually
write

    x :: Natural :#: Normal
    x = transition $ fromList Normal [0,1]

With these exponential family structures, it is also possible to define
exponential family harmoniums (e.g. Boltzmann machines) and multilayer
perceptrons. A restricted boltzmann machine for example, can be created with

    d :: Int
    d = 10

    brnls :: Replicated Bernoulli
    brnls = Replicated Bernoulli d

    hrm :: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    hrm = Harmonium brnls brnls

    x :: Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    x = fromList hrm $ replicate (dimension hrm) 0.1

which is a restricted Boltzmann machine with 10 hidden and 10 visible units,
all the weights of which have been initialized to 0.1.

*Scripts*

**backpropagation**: Applies a multilayer perceptron to a simple regression
problem.

**cross-entropy-descent**: Gradient descent on the negative log-likelihood
in a simple experiment, comparing natural and vanilla gradient descent.

**divergence**: Visualizes the KL-divergence with contour lines. In
particular, we compute the contour lines of the divergence between two
points in mixture and natural coordinates, for both the Bernoulli and
Poisson exponential families.

**multivariate**: Vizualises the multivariate normal distribution and a
corresponding MLE problem.

**poisson-binomial**: Shows how Poisson and binomial distributions can
approximate each other.

**population-code**: Visualizes population codes which are a common model in
computational neuroscience, which can be defined as a kind of harmonium.

**univariate**: Visualizes functionality of the Bernoulli, Categorical,
Poisson, and Normal distributions.

**von-mises**: Compares samples generated from the von-Mises distribution to
the unnormalized density.
