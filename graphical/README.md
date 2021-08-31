This library provides definitions and algorithms for various graphical models
such as mixture models, kalman Filters, and restricted Boltzmann machines, as
well as algorithms for fitting them e.g.  expectation maximization and
contrastive divergence minimization. Underlying all of these models is is a
generalized linear object known as a Harmonium, and in the following I will
briefly introduce them.

The core definition of this library is a `Manifold` of joint distributions I
call of an `AffineHarmonium`
```haskell
newtype AffineHarmonium f y x z w = AffineHarmonium (Affine f y z x, w)
```
which is a product `Manifold` composed of a of a `Manifold` of likelihood
functions `Affine f y z x`, and a `Manifold` of distributions `w` that partially
define the space of priors. `AffineHarmoniums` provide a bit more flexibility
than what I call a `Harmonium`
```haskell
type Harmonium f z w = AffineHarmonium f z w z w
```
which is a simpler object. Nevertheless, from a theoretical point of view, an
`AffineHarmonium` is a special case of a `Harmonium`, and we may think of them
more or less equivalently.

A `Harmonium` is a model over a observable variables and latent variables, and represents a sort of generalized linear joint distribution over the two of them. The theory for `Harmonium`s is summarized well by [this paper](https://papers.nips.cc/paper/2004/hash/0e900ad84f63618452210ab8baae0218-Abstract.html) and [this paper](https://proceedings.neurips.cc/paper/2013/hash/28f0b864598a1291557bed248a998d4e-Abstract.html). Although `Harmonium`s might seem like a little-studied, and esoteric object, various well known models, such as mixture models and restricted Boltzmann machines, are in fact `Harmonium`s, and various other models, such as factor analysis, Kalman filters, or hidden Markov models, can be expressed in terms of them.

All of the aforementioned models can be fit with Expectation-Maximization (EM), and EM can be expressed in an entirely general manner for `Harmonium`s. Firstly, the expectation step of a `Harmonium` is implemented by
```haskell
expectationStep
    :: ( ExponentialFamily z, Map Natural f x y, Bilinear f y x
       , Translation z y, Translation w x, LegendreExponentialFamily w )
    => Sample z -- ^ Model Samples
    -> Natural # AffineHarmonium f y x z w -- ^ Harmonium
    -> Mean # AffineHarmonium f y x z w -- ^ Harmonium expected sufficient statistics
expectationStep zs hrm =
    let mzs = sufficientStatistic <$> zs
        mys = anchor <$> mzs
        pstr = fst . split $ transposeHarmonium hrm
        mws = transition <$> pstr >$> mys
        mxs = anchor <$> mws
        myx = (>$<) mys mxs
     in joinHarmonium (average mzs) myx $ average mws
```
In summary, what we do is

- take some observations,
- compute their `sufficientStatistics`,
- map these statistics into the predictions of the latent variables,
- transition these latent predictions from `Natural` coordinates to `Mean` coordinates,
- and assemble the results into the `Mean` `sufficientStatistics` of the joint distribution.

The maximization step is then simply mapping the whole joint distribution from `Mean` back to `Natural` coordinates, such that `expectationMaximization` may be expressed as
```haskell
expectationMaximization
    :: ( DuallyFlatExponentialFamily (AffineHarmonium f y x z w)
       , ExponentialFamily z, Map Natural f x y, Bilinear f y x
       , Translation z y, Translation w x, LegendreExponentialFamily w )
    => Sample z
    -> Natural # AffineHarmonium f y x z w
    -> Natural # AffineHarmonium f y x z w
expectationMaximization zs hrm = transition $ expectationStep zs hrm
```

As such, for a wide variety of models, we may reduce implementing expectation maximization to instantating the class requirements of the `expectationMaximization` function. This is rarely trivial, but in some sense, much more straight-forward and well-defined that deriving EM algorithms from scratch.

For in-depth tutorials visit my
[blog](https://sacha-sokoloski.gitlab.io/website/pages/blog.html).
