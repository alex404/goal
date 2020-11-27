provides tools for implementing and applying statistical and machine learning
algorithms. The core concept of goal-probability is that of a statistical
manifold, i.e. manifold of probability distributions, with a focus on
exponential family distributions. Various graphical models are also defined
here, e.g. mixture models and restricted Boltzmann machines, as well as
algorithms for fitting them e.g. expectation maximization and contrastive
divergence minimization. What follows is brief introduction to how we define and
work with statistical manifolds in Goal.

The core definition of this library is that of a `Statistical` `Manifold`.
```haskell
class Manifold x => Statistical x where
    type SamplePoint x :: Type
```
A `Statistical` `Manifold` is a `Manifold` of probability distributions, such
that each point on the manifold is a probability distribution over the specified
space of `SamplePoint`s. We may evaluate the probability mass/density of a `SamplePoint` under a given distribution as long as the distribution is `AbsolutelyContinous`.
```haskell
class Statistical x => AbsolutelyContinuous c x where
    density :: Point c x -> SamplePoint x -> Double
    densities :: Point c x -> Sample x -> [Double]
```
Similarly, we may generate a `Sample` from a given distribution as long as it is `Generative`.
```haskell
type Sample x = [SamplePoint x]

class Statistical x => Generative c x where
    samplePoint :: Point c x -> Random r (SamplePoint x)
    sample :: Int -> Point c x -> Random r (Sample x)
```
In both these cases, class methods are defined both both single and bulk
evaluation, to potentially take advantage of bulk linear algebra operations.

Let us now look at some example distributions that we may define; for the sake
of brevity, I will not introduce every bit of necessary code. In
Goal we create a normal distribution by writing
```haskell
nrm :: Source # Normal
nrm = fromTuple (0,1)
```
where 0 is the mean and 1 is the variance. For each `Statistical` `Manifold`,
the `Source` coordinate system represents some standard parameterization, e.g.
as one typically finds on Wikipedia. Similarly, we can create a binomial
distribution with
```haskell
bnm :: Source # Binomial 5
bnm = Point $ S.singleton 0.5
```
which is a binomial distribution over 5 fair coin tosses.

Exponential families are also a core part of this library. An `ExponentiaFamily`
is a kind of `Statistical` `Manifold` defined in terms of a
`sufficientStatistic` and a `baseMeasure`.
```haskell
class Statistical x => ExponentialFamily x where
    sufficientStatistic :: SamplePoint x -> Mean # x
    baseMeasure :: Proxy x -> SamplePoint x -> Double
```

Exponential families may always be parameterized in terms of the so-called
`Natural` and `Mean` parameters. Mean coordinates are equal to the average value
of the `sufficientStatistic` under the given distribution. The `Natural`
coordinates are arguably less intuitive, but they are critical for implementing
evaluating exponential family distributions numerically. For example, the
unnormalized density function of an `ExponentialFamily` distribution is
given by
```haskell
unnormalizedDensity :: forall x . ExponentialFamily x => Natural # x -> SamplePoint x -> Double
unnormalizedDensity p x =
    exp (p <.> sufficientStatistic x) * baseMeasure (Proxy @ x) x
```

For in-depth tutorials visit my
[blog](https://sacha-sokoloski.gitlab.io/website/pages/blog.html).

