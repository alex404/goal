The core definition of this library is that of a statistical manifolds. A
statistical manifold is a manifold of probability distributions, such that
each point on the manifold is a probability distribution. In Goal, we can
create a normal distribution, for example, by writing

```haskell
nrm :: Source # Normal
nrm = Point $ S.fromTuple (0,1)
```
where 0 is the mean and 1 is the variance. For each `Statistical` `Manifold`,
the `Source` coordinate system represents some standard parameterization, e.g.
as one typically finds on Wikipedia. We can create a binomial distribution with
```haskell
bnm :: Source # Binomial 5
bnm = Point $ S.singleton 0.5
```
which is a binomial distribution over 5 fair coin tosses. Where possible, this
library provides implementations for evaluating the probability density function
of a given distribution, generating samples, and computing the maximum
likelihood estimator.

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

For more in-depth tutorials visit my
[blog](https://sacha-sokoloski.gitlab.io/website/pages/blog.html).

