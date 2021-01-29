{-# LANGUAGE Arrows #-}
-- | A collection of algorithms for optimizing harmoniums.

module Goal.Graphical.Learning
    ( -- * Expectation Maximization
      expectationMaximization
    , expectationMaximizationAscent
    , gibbsExpectationMaximization
    , latentProcessExpectationMaximization
    , latentProcessExpectationMaximizationAscent
    , timeDependentLatentProcessExpectationMaximizationAscent
    -- * Differentials
    , harmoniumInformationProjectionDifferential
    , contrastiveDivergence
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Generative
import Goal.Graphical.Generative.Harmonium
import Goal.Graphical.Generative.Dynamic
import Goal.Graphical.Inference

import qualified Data.List as L


--- Differentials ---


-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: ( Map Natural f y x, LegendreExponentialFamily z
       , SamplePoint w ~ SamplePoint x, Translation z y
       , ExponentialFamily x, ExponentialFamily w, Generative Natural w )
    => Int
    -> Natural # Harmonium f y x z w -- ^ Harmonium
    -> Natural # w -- ^ Model Distribution
    -> Random r (Mean # w) -- ^ Differential Estimate
harmoniumInformationProjectionDifferential n hrm px = do
    xs <- sample n px
    let (lkl,nw) = split hrm
        mys0 = lkl >$>* xs
        mws = sufficientStatistic <$> xs
        mys = zipWith (\mw my0 -> mw <.> (px - nw) - potential my0) mws mys0
        ln = fromIntegral $ length xs
        mwht = average mws
        myht = sum mys / ln
        foldfun (mw,my) (k,z0) = (k+1,z0 + ((my - myht) .> (mw - mwht)))
    return . uncurry (/>) . foldr foldfun (-1,0) $ zip mws mys

-- | Contrastive divergence on harmoniums (<https://www.mitpressjournals.org/doi/abs/10.1162/089976602760128018?casa_token=x_Twj1HaXcMAAAAA:7-Oq181aubCFwpG-f8Lo1wRKvGnmujzl8zjn9XbeO5nGhfvKCCQjsu4K4pJCkMNYUYWqc2qG7TRXBg Hinton, 2019>).
contrastiveDivergence
    :: ( Generative Natural z, ExponentialFamily z, Translation w x
       , Generative Natural w, ExponentialFamily y, Translation z y
       , LegendreExponentialFamily w, Bilinear f y x, Map Natural f x y
       , Map Natural f y x, SamplePoint y ~ SamplePoint z
       , SamplePoint x ~ SamplePoint w, ExponentialFamily x )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample z -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium f y x z w -- ^ The harmonium
      -> Random s (Mean # Harmonium f y x z w) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
    return $ stochasticRelativeEntropyDifferential xzs0 xzs1


--- Expectation Maximization ---


-- | EM for latent variable odels.
expectationMaximization
    :: ( DuallyFlatExponentialFamily f, ExpectationMaximization f )
    => Observations f
    -> Natural # f
    -> Natural # f
expectationMaximization zs hrm = transition $ expectationStep zs hrm

-- | Ascent of the EM objective on harmoniums for when the expectation
-- step can't be computed in closed-form. The convergent harmonium distribution
-- of the output harmonium-list is the result of 1 iteration of the EM
-- algorithm.
expectationMaximizationAscent
    :: (LegendreExponentialFamily f, ExpectationMaximization f)
    => Double -> GradientPursuit -> Observations f -> Natural # f -> [Natural # f]
expectationMaximizationAscent eps gp zs nhrm =
    let mhrm' = expectationStep zs nhrm
     in vanillaGradientSequence (relativeEntropyDifferential mhrm') (-eps) gp nhrm

-- | Ascent of the EM objective on harmoniums for when the expectation
-- step can't be computed in closed-form. The convergent harmonium distribution
-- of the output harmonium-list is the result of 1 iteration of the EM
-- algorithm.
gibbsExpectationMaximization
    :: ( ExponentialFamily z, Map Natural f x y, Manifold w, Map Natural f y x
       , Translation z y, Translation w x, SamplePoint y ~ SamplePoint z
       , SamplePoint w ~ SamplePoint x, ExpectationMaximization (Harmonium f y x z w)
       , ExponentialFamily y, Generative Natural w, ExponentialFamily x
       , Generative Natural z, Manifold (Harmonium f y x z w)
       , Bilinear f y x, LegendreExponentialFamily w )
    => Double
    -> Int
    -> Int
    -> GradientPursuit
    -> Sample z -- ^ Observations
    -> Natural # Harmonium f y x z w -- ^ Current Harmonium
    -> Chain (Random r) (Natural # Harmonium f y x z w) -- ^ Harmonium Chain
gibbsExpectationMaximization eps cdn nbtch gp zs0 nhrm0 =
    let mhrm0 = expectationStep zs0 nhrm0
     in chainCircuit nhrm0 $ proc nhrm -> do
         zs <- minibatcher nbtch zs0 -< ()
         xzs0 <- arrM (uncurry initialPass) -< (nhrm,zs)
         xzs1 <- arrM (\(x,y) -> iterateM' cdn (gibbsPass x) y) -< (nhrm,xzs0)
         let dff = mhrm0 - averageSufficientStatistic xzs1
         gradientCircuit eps gp -< (nhrm,vanillaGradient dff)

latentProcessExpectationStep
    :: ( ConjugatedLikelihood g x x w w, ConjugatedLikelihood f y x z w
       , Transition Natural Mean w, Transition Natural Mean (Harmonium g x x w w)
       , Manifold (Harmonium g x x w w)
       , Bilinear g x x, Map Natural f x y, Bilinear f y x
       , SamplePoint y ~ SamplePoint z )
    => Observations (LatentProcess f g y x z w)
    -> Natural # LatentProcess f g y x z w
    -> (Mean # w, Mean # Harmonium f y x z w, Mean # Harmonium g x x w w)
latentProcessExpectationStep zss ltnt =
    let (prr,emsn,trns) = splitLatentProcess ltnt
        (smthss,hrmss) = unzip $ conjugatedSmoothing trns emsn prr <$> zss
        mprr = average $ toMean . head <$> smthss
        mtrns = average $ toMean <$> concat hrmss
        mws = toMean <$> concat smthss
        mzs = sufficientStatistic <$> concat zss
        mys = anchor <$> mzs
        mxs = anchor <$> mws
        memsn = joinHarmonium (average mzs) (mys >$< mxs) (average mws)
     in (mprr,memsn,mtrns)

latentProcessExpectationMaximization
    :: ( ConjugatedLikelihood g x x w w, ConjugatedLikelihood f y x z w
       , Transition Natural Mean w, Transition Natural Mean (Harmonium g x x w w)
       , Transition Mean Natural w
       , Transition Mean Natural (Harmonium f y x z w)
       , Transition Mean Natural (Harmonium g x x w w)
       , Manifold (Harmonium g x x w w)
       , Bilinear g x x, Map Natural f x y, Bilinear f y x
       , SamplePoint y ~ SamplePoint z )
    => Observations (LatentProcess f g y x z w)
    -> Natural # LatentProcess f g y x z w
    -> Natural # LatentProcess f g y x z w
latentProcessExpectationMaximization zss ltnt =
    let (mprr,memsn,mtrns) = latentProcessExpectationStep zss ltnt
        prr' = toNatural mprr
        emsn' = fst . split $ toNatural memsn
        trns' = fst . split $ toNatural mtrns
     in joinLatentProcess prr' emsn' trns'

latentProcessExpectationMaximizationAscent
    :: ( ConjugatedLikelihood g x x w w, ConjugatedLikelihood f y x z w
       , DuallyFlatExponentialFamily w
       , LegendreExponentialFamily (Harmonium f y x z w)
       , LegendreExponentialFamily (Harmonium g x x w w)
       , Bilinear g x x, Map Natural f x y, Bilinear f y x
       , SamplePoint y ~ SamplePoint z )
    => Double
    -> Int
    -> GradientPursuit
    -> [Sample z]
    -> Natural # LatentProcess f g y x z w
    -> Natural # LatentProcess f g y x z w
latentProcessExpectationMaximizationAscent eps nstps gp zss ltnt =
    let (mprr,mehrm,mthrm) = latentProcessExpectationStep zss ltnt
        (nprr,nemsn,ntrns) = splitLatentProcess ltnt
        neql0 = toNatural . snd $ split mehrm
        neql1 = toNatural . snd $ split mthrm
        nehrm = joinConjugatedHarmonium nemsn neql0
        nthrm = joinConjugatedHarmonium ntrns neql1
        nprr' = (!! nstps)
            $ vanillaGradientSequence (relativeEntropyDifferential mprr) (-eps) gp nprr
        nemsn' = fst . split . (!! nstps)
            $ vanillaGradientSequence (relativeEntropyDifferential mehrm) (-eps) gp nehrm
        ntrns' = fst . split . (!! nstps)
            $ vanillaGradientSequence (relativeEntropyDifferential mthrm) (-eps) gp nthrm
     in joinLatentProcess nprr' nemsn' ntrns'

timeDependentLatentProcessExpectationStep
    :: ( ConjugatedLikelihood g x x w w, ConjugatedLikelihood f y x z w
       , Transition Natural Mean w, Transition Natural Mean (Harmonium g x x w w)
       , Manifold (Harmonium g x x w w), KnownNat n
       , Bilinear g x x, Map Natural f x y, Bilinear f y x
       , SamplePoint y ~ SamplePoint z )
    => Observations (LatentProcess f g y x z w)
    -> Natural # Affine Tensor y (LatentProcess f g y x z w) (Categorical n)
    -> (Mean # w, [Mean # Harmonium f y x z w], Mean # Harmonium g x x w w)
timeDependentLatentProcessExpectationStep zss cltnt =
    let (prr,emsn,trns) = splitLatentProcess . fst $ split cltnt
        (smthss,hrmss) = unzip $ conjugatedSmoothing trns emsn prr <$> zss
        mprr = average $ toMean . head <$> smthss
        mtrns = average $ toMean <$> concat hrmss
        memsns = do
            (smths,zs) <- zip (L.transpose smthss) (L.transpose zss)
            let mws = toMean <$> smths
                mzs = sufficientStatistic <$> zs
                mys = anchor <$> mzs
                mxs = anchor <$> mws
            return $ joinHarmonium (average mzs) (mys >$< mxs) (average mws)
     in (mprr,memsns,mtrns)

timeDependentLatentProcessExpectationMaximizationAscent
    :: forall f g y x z w n
    . ( ConjugatedLikelihood g x x w w, ConjugatedLikelihood f y x z w
      , DuallyFlatExponentialFamily w, KnownNat n
      , LegendreExponentialFamily (Harmonium f y x z w)
      , LegendreExponentialFamily (Harmonium g x x w w)
      , Bilinear g x x, Map Natural f x y, Bilinear f y x
      , SamplePoint y ~ SamplePoint z )
    => Double
    -> Int
    -> GradientPursuit
    -> [Sample z]
    -> Natural # Affine Tensor y (LatentProcess f g y x z w) (Categorical n)
    -> Natural # Affine Tensor y (LatentProcess f g y x z w) (Categorical n)
timeDependentLatentProcessExpectationMaximizationAscent eps nstps gp zss cltnt =
    let (mprr,mehrms,mthrm) = timeDependentLatentProcessExpectationStep zss cltnt
        (ltnt,nyk) = split cltnt
        (nprr,nemsn,ntrns) = splitLatentProcess ltnt
        nsmth = toNatural . average $ snd . split <$> mehrms
        nehrm = joinConjugatedHarmonium nemsn nsmth
        nthrm = joinConjugatedHarmonium ntrns nsmth
        nprr' = (!! nstps)
            $ vanillaGradientSequence (relativeEntropyDifferential mprr) (-eps) gp nprr
        ntrns' = fst . split . (!! nstps)
            $ vanillaGradientSequence (relativeEntropyDifferential mthrm) (-eps) gp nthrm
        cnehrm :: Natural # Affine Tensor y (Harmonium f y x z w) (Categorical n)
        cnehrm = join nehrm nyk
        myxs = zip mehrms $ sufficientStatistic <$> [0..]
        backpropper = backpropagation relativeEntropyDifferential myxs
        cnehrm' = vanillaGradientSequence backpropper (-eps) gp cnehrm !! nstps
        (nehrm',nyk') = split cnehrm'
        nemsn' = fst $ split nehrm'
     in join (joinLatentProcess nprr' nemsn' ntrns') nyk'
