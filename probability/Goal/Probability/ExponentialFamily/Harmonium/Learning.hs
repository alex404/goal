{-# LANGUAGE
    RankNTypes,
    TypeOperators,
    DataKinds,
    GADTs,
    FlexibleContexts,
    ScopedTypeVariables
#-}

-- | A collection of algorithms for optimizing harmoniums.

module Goal.Probability.ExponentialFamily.Harmonium.Learning
    ( -- * Differentials
      stochasticConjugatedHarmoniumDifferential
    , harmoniumInformationProjectionDifferential
    , stochasticMixtureModelDifferential
    , contrastiveDivergence
      -- ** Conditional
    , mixtureStochasticConditionalCrossEntropyDifferential
    , conditionalHarmoniumConjugationDifferential
    -- * Algorithms
    -- ** Expectation Maximization
    , mixtureExpectationMaximization
    , deepMixtureModelExpectationStep
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Conjugation
import Goal.Probability.ExponentialFamily.Harmonium.Conditional


--- Differentials ---


-- | Estimates the stochastic cross entropy differential of a conjugated harmonium with
-- respect to the relative entropy, and given an observation.
stochasticConjugatedHarmoniumDifferential
    :: ( Map Mean Natural f z x, Bilinear f z x, ExponentialFamily z
       , ExponentialFamily x, Generative Natural z, Generative Natural x )
       => Sample z -- ^ Observations
       -> Natural # x -- ^ Conjugation Parameters
       -> Natural # Harmonium f z x -- ^ Harmonium
       -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ Differential
{-# INLINE stochasticConjugatedHarmoniumDifferential #-}
stochasticConjugatedHarmoniumDifferential zs rprms hrm = do
    pzxs <- initialPass hrm zs
    qzxs <- sampleConjugatedHarmonium (length zs) (toSingletonSum rprms) hrm
    return $ stochasticCrossEntropyDifferential' pzxs qzxs

-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: (Generative Natural x, ExponentialFamily x, Map Mean Natural f z x, Legendre Natural z)
    => Int
    -> Natural # Harmonium f z x -- ^ Harmonium
    -> Natural # x -- ^ Model Distribution
    -> Random s (CotangentVector Natural x) -- ^ Differential Estimate
{-# INLINE harmoniumInformationProjectionDifferential #-}
harmoniumInformationProjectionDifferential n hrm px = do
    xs <- sample n px
    let (affmn,nm0) = splitBottomHarmonium hrm
        (nn,nmn) = splitAffine affmn
        nm = fromOneHarmonium nm0
        mxs = sufficientStatistic <$> xs
        mys0 = nmn >$> mxs
        mys = zipWith (\mx my0 -> mx <.> (px <-> nm) - potential (nn <+> my0)) mxs mys0
        ln = fromIntegral $ length xs
        mxht = averagePoint mxs
        myht = sum mys / ln
        foldfun (mx,my) (k,z0) = (k+1,z0 <+> ((my - myht) .> (mx <-> mxht)))
        cvr = uncurry (/>) . foldr foldfun (-1,zero) $ zip mxs mys
    return $ primalIsomorphism cvr

-- | The stochastic cross entropy differential of a mixture model.
stochasticMixtureModelDifferential
    :: ( Enum e, Legendre Natural z, KnownNat n, ExponentialFamily z )
      => Sample z -- ^ Observations
      -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Categorical harmonium
      -> CotangentVector Natural (Harmonium Tensor z (Categorical e n)) -- ^ Differential
{-# INLINE stochasticMixtureModelDifferential #-}
stochasticMixtureModelDifferential zs hrm =
    let pxs = harmoniumEmpiricalExpectations zs hrm
        qxs = dualTransition hrm
     in primalIsomorphism $ qxs <-> pxs

contrastiveDivergence
    :: ( Generative Natural z, ExponentialFamily z, ExponentialFamily x
       , Generative Natural x, Map Mean Natural f z x, Bilinear f z x )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample z -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium f z x -- ^ The harmonium
      -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
    return $ stochasticCrossEntropyDifferential' xzs0 xzs1

-- | The stochastic conditional cross-entropy differential, based on target
-- inputs and outputs expressed as distributions in mean coordinates.
mixtureStochasticConditionalCrossEntropyDifferential
    :: ( ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k )
    => Sample x -- ^ Input mean distributions
    -> Sample z -- ^ Output mean distributions
    -> Mean #> Natural # MixtureGLM z k x -- ^ Function
    -> CotangentVector (Mean #> Natural) (MixtureGLM z k x) -- ^ Differential
{-# INLINE mixtureStochasticConditionalCrossEntropyDifferential #-}
mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm =
    -- This could be better optimized but not throwing out the second result of propagate
    let dmglms = dualIsomorphism
            <$> zipWith stochasticMixtureModelDifferential ((:[]) <$> zs) (mglm >$>* xs)
        dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmglm | dmglm <- dmglms ]
        f = snd $ splitBottomSubLinear mglm
        df = fst $ propagate dzs (sufficientStatistic <$> xs) f
     in primalIsomorphism $ joinBottomSubLinear (averagePoint dmglms) df

-- | A gradient for conjugateing gains which won't allow them to be negative.
conditionalHarmoniumConjugationDifferential
    :: ( ExponentialFamily x, Map Mean Natural f z x
       , Legendre Natural (DeepHarmonium gs (z : zs)) )
    => Double -- ^ Conjugation shift
    -> Natural # x -- ^ Conjugation parameters
    -> Sample x -- ^ Sample points
    -> Mean #> Natural # f z x -- ^ linear part of ppc
    -> Natural # DeepHarmonium gs (z : zs) -- ^ Gains
    -> CotangentPair Natural (DeepHarmonium gs (z : zs)) -- ^ Conjugated PPC
{-# INLINE conditionalHarmoniumConjugationDifferential #-}
conditionalHarmoniumConjugationDifferential rho0 rprms xsmps tns dhrm =
    let lkl = joinBottomSubLinear dhrm tns
        rcts = conjugationCurve rho0 rprms xsmps
        ndhrmlkls = lkl >$>* xsmps
        mdhrmlkls = dualTransition <$> ndhrmlkls
        ptns = potential <$> ndhrmlkls
     in joinTangentPair dhrm . averagePoint
         $ [ primalIsomorphism $ (ptn - rct) .> mdhrmlkl | (rct,mdhrmlkl,ptn) <- zip3 rcts mdhrmlkls ptns ]

-- | EM implementation for mixture models/categorical harmoniums.
mixtureExpectationMaximization
    :: ( Enum e, Legendre Natural z, KnownNat n, ExponentialFamily z, Transition Mean Natural z )
    => Sample z -- ^ Observations
    -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Current Harmonium
    -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Updated Harmonium
{-# INLINE mixtureExpectationMaximization #-}
mixtureExpectationMaximization zs hrm =
    let zs' = hSingleton <$> zs
        (cats,mzs) = deepMixtureModelExpectationStep zs' $ transposeHarmonium hrm
     in joinMixtureModel (S.map (toNatural . fromOneHarmonium) mzs) cats

-- | E-step implementation for deep mixture models/categorical harmoniums. Note
-- that for the sake of type signatures, this acts on transposed harmoniums
-- (i.e. the categorical variables are at the bottom of the hierarchy).
deepMixtureModelExpectationStep
    :: ( KnownNat n, Enum e, ExponentialFamily x, ExponentialFamily (DeepHarmonium fs (x : zs)) )
    => Sample (DeepHarmonium fs (x ': zs)) -- ^ Observations
    -> Natural # DeepHarmonium (Tensor ': fs) (Categorical e n ': x ': zs) -- ^ Current Harmonium
    -> (Natural # Categorical e n, S.Vector (n+1) (Mean # DeepHarmonium fs (x ': zs)))
{-# INLINE deepMixtureModelExpectationStep #-}
deepMixtureModelExpectationStep xzs dhrm =
    let aff = fst $ splitBottomHarmonium dhrm
        muss = toMean <$> aff >$>* fmap hHead xzs
        sxzs = sufficientStatistic <$> xzs
        (cmpnts0,nrms) = foldr folder (S.replicate zero, S.replicate 0) $ zip muss sxzs
     in (toNatural $ averagePoint muss, S.zipWith (/>) nrms cmpnts0)
    where folder (Point cs,sxz) (cmpnts,nrms) =
              let ws = S.cons (1 - S.sum cs) cs
                  cmpnts' = S.map (.> sxz) ws
               in (S.zipWith (<+>) cmpnts cmpnts', S.add nrms ws)

--iterativeMixtureModelOptimization
--    :: forall e n z . ( Enum e, KnownNat n , Legendre Natural z, ExponentialFamily z )
--    => Int -- ^ Number of gradient steps per iteration
--    -> Double -- ^ Step size
--    -> GradientPursuit -- ^ Gradient Pursuit Algorithm
--    -> Maybe Int -- ^ Minibatch size (or just full batch)
--    -> Double -- ^ New component ratio (must be between 0 and 1)
--    -> Sample z -- ^ Observations
--    -> Natural # z -- ^ Initial Distribution
--    -> (Natural # Harmonium Tensor z (Categorical e n),[[Double]]) -- ^ (Final Mixture and gradient descent)
--iterativeMixtureModelOptimization nstps eps grd mnbtch rto zss nz =
--    let nbtch = fromMaybe (length zss) mnbtch
--     in iterativeMixtureModelOptimization0 nstps eps grd nbtch rto zss (toNullMixture nz,[])
--
---- | An iterative algorithm for training a mixture model.
--iterativeMixtureModelOptimization0
--    :: forall e n k z . ( Enum e, KnownNat n, KnownNat k, Legendre Natural z, ExponentialFamily z )
--    => Int -- ^ Number of gradient steps per iteration
--    -> Double -- ^ Step size
--    -> GradientPursuit -- ^ Gradient Pursuit Algorithm
--    -> Int -- ^ Minibatch size
--    -> Double -- ^ New component ratio (must be between 0 and 1)
--    -> Sample z -- ^ Observations
--    -> (Natural # Harmonium Tensor z (Categorical e k),[[Double]]) -- ^ Initial Mixture and LLs
--    -> (Natural # Harmonium Tensor z (Categorical e n),[[Double]]) -- ^ Final Mixture and LLs
--iterativeMixtureModelOptimization0 nstps eps grd nbtch rto zss (hrm0,llss) =
--    let trncrc :: Circuit Identity [SamplePoint z] (Natural # Harmonium Tensor z (Categorical e k),Double)
--        trncrc = accumulateCircuit hrm0 $ proc (zs,hrm) -> do
--            let ll = average $ log . mixtureDensity hrm <$> zs
--                dhrm = stochasticMixtureModelDifferential zs hrm
--            hrm' <- gradientCircuit eps grd -< joinTangentPair hrm $ vanillaGradient dhrm
--            returnA -< ((hrm,ll),hrm')
--        (hrms,lls) = unzip . runIdentity . streamCircuit trncrc . take nstps . breakEvery nbtch $ cycle zss
--        llss' = lls : llss
--        hrm1 = last hrms
--     in case sameNat (Proxy @ k) (Proxy @ n) of
--          Just Refl -> (hrm1, reverse llss')
--          Nothing -> iterativeMixtureModelOptimization0 nstps eps grd nbtch rto zss
--              (expandMixtureModel rto hrm1,llss')
--
--expandMixtureModel
--    :: forall e z n
--     . ( Enum e, Legendre Natural z, KnownNat n, ExponentialFamily z )
--    => Double -- ^ Weight fraction for new component (0 < x < 1)
--    -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Current Harmonium
--    -> Natural # Harmonium Tensor z (Categorical e (n+1)) -- ^ Updated Harmonium
--expandMixtureModel rto hrm =
--    let (nzs,nx) = splitMixtureModel hrm
--        sx = toSource nx
--        nwght = density nx (toEnum 0)
--        (mxwght,mxidx) = maximumBy (comparing fst) $ zip (density nx (toEnum 0) : listCoordinates sx) [0..]
--     in if mxidx == nidx
--           then let sx' :: Source # Categorical e (n+1)
--                    sx' = Point . S.cons (rto * nwght) $ coordinates sx
--                    nzs' = S.cons (S.last nzs) nzs
--                 in buildMixtureModel nzs' $ toNatural sx'
--           else let sx' :: Source # Categorical e n
--                    sx' = Point $ S.unsafeUpd (coordinates sx) [(mxidx,(1-rto)*mxwght)]
--                    sx'' :: Source # Categorical e (n+1)
--                    sx'' = Point . S.cons (rto * nwght) $ coordinates sx'
--                    nzs' = S.cons (S.last nzs) nzs
--                 in buildMixtureModel nzs' $ toNatural sx''
--

--dualContrastiveDivergence
--    :: forall s f z x
--    . ( Generative Natural z, ExponentialFamily z, ExponentialFamily x, Generative Natural x
--      , Map Mean Natural f x z, Bilinear f z x, Bilinear f x z )
--      => Int -- ^ The number of contrastive divergence steps
--      -> Int -- ^ The number of samples
--      -> Natural # x -- ^ Target marginal
--      -> Natural # Harmonium f z x -- ^ The harmonium
--      -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ The gradient estimate
--dualContrastiveDivergence cdn nsmps prr hrm = do
--    xs <- sample nsmps prr
--    dhrm' <- contrastiveDivergence cdn xs $ transposeHarmonium hrm
--    return $ primalIsomorphism . transposeHarmonium $ dualIsomorphism dhrm'
--
----class FitConjugationParameters (fs :: [* -> * -> *]) (ms :: [*]) where
----    fitConjugationParameters
----        :: Double
----        -> Maybe Int
----        -> Natural # DeepHarmonium fs ms
----        -> Natural # Sum (Tail ms)
----        -> Random s (Natural # Sum (Tail ms))
----
----instance FitConjugationParameters '[] '[m] where
----    {-# INLINE fitConjugationParameters #-}
----    fitConjugationParameters _ _ _ _ = zero
----
----instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f z x, Manifold (Sum ms)
----         , ExponentialFamily n, SampleConjugated fs (n : ms), Generative Natural m
----         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
----  => SampleConjugated (f : fs) (m : n : ms) where
----    {-# INLINE sampleConjugated #-}
----    sampleConjugated rprms dhrm = do
----        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
----            (rprm,rprms') = splitSum rprms
----        (ys,xs) <- fmap hUnzip . sampleConjugated rprms' $ biasBottom rprm dhrm'
----        zs <- samplePoint $ mapReplicatedPoint (pn <+>) (pf >$>* ys)
----        return . hZip zs $ hZip ys xs
----
----
