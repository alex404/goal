-- | A collection of algorithms for optimizing harmoniums.

module Goal.Probability.ExponentialFamily.Harmonium.Learning
    ( -- * Differentials
      harmoniumInformationProjectionDifferential
    , contrastiveDivergence
      -- ** Conditional
    --, mixtureStochasticConditionalCrossEntropyDifferential
    --, mixtureStochasticConditionalEMAscent
    -- * Expectation Maximization
    , expectationMaximization
    , expectationMaximizationPursuit
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Conditional

import qualified Data.Vector as V
import System.Random.MWC.Probability hiding (initialize,sample)
import System.Random.MWC.Distributions (uniformShuffle)



--- Differentials ---


-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: ( Map Mean Natural f z x, LegendreExponentialFamily z
       , ExponentialFamily x, Generative Natural x)
    => Int
    -> Natural # Harmonium f z x -- ^ Harmonium
    -> Natural # x -- ^ Model Distribution
    -> Random r (Mean # x) -- ^ Differential Estimate
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
    return . uncurry (/>) . foldr foldfun (-1,zero) $ zip mxs mys

contrastiveDivergence
    :: ( Generative Natural z, ExponentialFamily z, Generative Natural x
       , ExponentialFamily x, Map Mean Natural f z x, Bilinear f z x )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample z -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium f z x -- ^ The harmonium
      -> Random s (Mean # Harmonium f z x) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
    return $ stochasticCrossEntropyDifferential xzs0 xzs1


--- Expectation Maximization ---


-- | EM implementation for mixture models/categorical harmoniums.
expectationMaximization
    :: ( DuallyFlatExponentialFamily (Harmonium f z x), LegendreExponentialFamily x
       , ExponentialFamily z, Bilinear f z x, Map Mean Natural f x z )
    => Sample z -- ^ Observations
    -> Natural # Harmonium f z x -- ^ Current Harmonium
    -> Natural # Harmonium f z x -- ^ Updated Harmonium
{-# INLINE expectationMaximization #-}
expectationMaximization zs hrm = transition $ harmoniumEmpiricalExpectations zs hrm

expectationMaximizationPursuit
    :: ( LegendreExponentialFamily (Harmonium f z x), LegendreExponentialFamily x
       , ExponentialFamily z, Bilinear f z x, Map Mean Natural f x z )
    => Double
    -> GradientPursuit
    -> Sample z -- ^ Observations
    -> Natural # Harmonium f z x -- ^ Current Harmonium
    -> [Natural # Harmonium f z x] -- ^ Updated Harmonium
expectationMaximizationPursuit eps gp zs nhrm =
    let mhrm' = harmoniumEmpiricalExpectations zs nhrm
     in vanillaGradientSequence (crossEntropyDifferential mhrm') (-eps) gp nhrm


---- | Estimates the stochastic cross entropy differential of a conjugated harmonium with
---- respect to the relative entropy, and given an observation.
--stochasticConjugatedHarmoniumDifferential
--    :: ( Map Mean Natural f z x, Bilinear f z x, ExponentialFamily z
--       , ExponentialFamily x, Generative Natural z, Generative Natural x )
--       => Sample z -- ^ Observations
--       -> Natural # x -- ^ Conjugation Parameters
--       -> Natural # Harmonium f z x -- ^ Harmonium
--       -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ Differential
--{-# INLINE stochasticConjugatedHarmoniumDifferential #-}
--stochasticConjugatedHarmoniumDifferential zs rprms hrm = do
--    pzxs <- initialPass hrm zs
--    qzxs <- sampleConjugatedHarmonium (length zs) (toSingletonSum rprms) hrm
--    return $ stochasticCrossEntropyDifferential' pzxs qzxs
--
---- | The stochastic conditional cross-entropy differential, based on target
---- inputs and outputs expressed as distributions in mean coordinates.
--mixtureStochasticConditionalCrossEntropyDifferential
--    :: ( ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k )
--    => Sample x -- ^ Input mean distributions
--    -> Sample z -- ^ Output mean distributions
--    -> Mean #> Natural # MixtureGLM z k x -- ^ Function
--    -> CotangentVector (Mean #> Natural) (MixtureGLM z k x) -- ^ Differential
--{-# INLINE mixtureStochasticConditionalCrossEntropyDifferential #-}
--mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm =
--    -- This could be better optimized but not throwing out the second result of propagate
--    let dmglms = dualIsomorphism
--            <$> zipWith stochasticMixtureDifferential ((:[]) <$> zs) (mglm >$>* xs)
--        dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmglm | dmglm <- dmglms ]
--        f = snd $ splitBottomSubLinear mglm
--        df = fst $ propagate dzs (sufficientStatistic <$> xs) f
--     in primalIsomorphism $ joinBottomSubLinear (averagePoint dmglms) df
--
---- | NB: Write now this never shuffles the training data, as it probably should.
--mixtureStochasticConditionalEMAscent
--    :: ( ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k )
--    => Double -- ^ Learning Rate
--    -> GradientPursuit -- ^ Gradient Pursuit Algorithm
--    -> Int -- ^ Batch size
--    -> Int -- ^ Number of iterations
--    -> Sample x -- ^ Conditioning sample
--    -> Sample z -- ^ Observable Sample
--    -> Mean #> Natural # MixtureGLM z k x
--    -> Random r (Mean #> Natural # MixtureGLM z k x)
--{-# INLINE mixtureStochasticConditionalEMAscent #-}
--mixtureStochasticConditionalEMAscent eps gp nbtch nstps xs0 zs0 mglm0 = do
--    let mxtrs0 = mglm0 >$>* xs0
--        szcs0 = zipWith harmoniumEmpiricalExpectations ((:[]) <$> zs0) mxtrs0
--        dmglmcrc = loopCircuit' mglm0 $ proc (xsczs,mglm) -> do
--            let (xs,sczs) = unzip xsczs
--            let mxtrs = dualTransition <$> mglm >$>* xs
--                f = snd $ splitBottomSubLinear mglm
--                dmxtrs = zipWith (<->) sczs mxtrs
--                dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmxtr | dmxtr <- dmxtrs ]
--                df = fst $ propagate dzs (sufficientStatistic <$> xs) f
--                dmglm = primalIsomorphism $ joinBottomSubLinear (averagePoint dmxtrs) df
--            gradientCircuit eps gp -< joinTangentPair mglm $ vanillaGradient dmglm
--        ncycs = 1 + div (length xs0 - 1) (nstps * nbtch)
--    smpss <- replicateM ncycs (shuffleList $ zip xs0 szcs0)
--    let smps = take nstps . breakEvery nbtch $ concat smpss
--    iterateCircuit dmglmcrc smps
--    ---- This could be better optimized but not throwing out the second result of propagate
--    --let dmglms = dualIsomorphism
--    --        <$> zipWith stochasticMixtureDifferential ((:[]) <$> zs) (mglm >$>* xs)
--    --    dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmglm | dmglm <- dmglms ]
--    --    f = snd $ splitBottomSubLinear mglm
--    --    df = fst $ propagate dzs (sufficientStatistic <$> xs) f
--    -- in primalIsomorphism $ joinBottomSubLinear (averagePoint dmglms) df
--
--shuffleList :: [a] -> Random r [a]
--shuffleList xs = fmap V.toList . Prob $ uniformShuffle (V.fromList xs)
--
---- | A gradient for conjugateing gains which won't allow them to be negative.
--conditionalHarmoniumConjugationDifferential
--    :: ( ExponentialFamily x, Map Mean Natural f z x
--       , Legendre Natural (DeepHarmonium gs (z : zs)) )
--    => Double -- ^ Conjugation shift
--    -> Natural # x -- ^ Conjugation parameters
--    -> Sample x -- ^ Sample points
--    -> Mean #> Natural # f z x -- ^ linear part of ppc
--    -> Natural # DeepHarmonium gs (z : zs) -- ^ Gains
--    -> CotangentPair Natural (DeepHarmonium gs (z : zs)) -- ^ Conjugated PPC
--{-# INLINE conditionalHarmoniumConjugationDifferential #-}
--conditionalHarmoniumConjugationDifferential rho0 rprms xsmps tns dhrm =
--    let lkl = joinBottomSubLinear dhrm tns
--        rcts = conjugationCurve rho0 rprms xsmps
--        ndhrmlkls = lkl >$>* xsmps
--        mdhrmlkls = dualTransition <$> ndhrmlkls
--        ptns = potential <$> ndhrmlkls
--     in joinTangentPair dhrm . averagePoint
--         $ [ primalIsomorphism $ (ptn - rct) .> mdhrmlkl | (rct,mdhrmlkl,ptn) <- zip3 rcts mdhrmlkls ptns ]
--
---- | EM implementation for mixture models/categorical harmoniums.
--mixtureExpectationMaximization
--    :: ( ClosedFormExponentialFamily z, KnownNat n )
--    => Sample z -- ^ Observations
--    -> Natural # Mixture z n -- ^ Current Harmonium
--    -> Natural # Mixture z n -- ^ Updated Harmonium
--{-# INLINE mixtureExpectationMaximization #-}
--mixtureExpectationMaximization zs hrm = mixtureParameters $ harmoniumEmpiricalExpectations zs hrm
--
--mixtureEMPursuit
--    :: ( ExponentialFamily z, Legendre Natural z, KnownNat n )
--    => Double
--    -> GradientPursuit
--    -> Sample z -- ^ Observations
--    -> Natural # Mixture z n -- ^ Current Harmonium
--    -> [Natural # Mixture z n] -- ^ Updated Harmonium
--mixtureEMPursuit eps gp zs nhrm =
--    let mhrm' = harmoniumEmpiricalExpectations zs nhrm
--     in vanillaGradientSequence (pairTangentFunction $ crossEntropyDifferential mhrm') eps gp nhrm
--
--
------ | E-step implementation for deep mixture models/categorical harmoniums. Note
------ that for the sake of type signatures, this acts on transposed harmoniums
------ (i.e. the categorical variables are at the bottom of the hierarchy).
----deepMixtureExpectationStep
----    :: ( KnownNat n, Enum e, ExponentialFamily x, ExponentialFamily (DeepHarmonium fs (x : zs)) )
----    => Sample (DeepHarmonium fs (x ': zs)) -- ^ Observations
----    -> Natural # DeepHarmonium (Tensor ': fs) (Categorical e n ': x ': zs) -- ^ Current Harmonium
----    -> (Natural # Categorical e n, S.Vector (n+1) (Mean # DeepHarmonium fs (x ': zs)))
----{-# INLINE deepMixtureExpectationStep #-}
----deepMixtureExpectationStep xzs dhrm =
----    let aff = fst $ splitBottomHarmonium dhrm
----        muss = toMean <$> aff >$>* fmap hHead xzs
----        sxzs = sufficientStatistic <$> xzs
----        (cmpnts0,nrms) = foldr folder (S.replicate zero, S.replicate 0) $ zip muss sxzs
----     in (toNatural $ averagePoint muss, S.zipWith (/>) nrms cmpnts0)
----    where folder (Point cs,sxz) (cmpnts,nrms) =
----              let ws = S.cons (1 - S.sum cs) cs
----                  cmpnts' = S.map (.> sxz) ws
----               in (S.zipWith (<+>) cmpnts cmpnts', S.add nrms ws)
--
----iterativeMixtureOptimization
----    :: forall e n z . ( Enum e, KnownNat n , Legendre Natural z, ExponentialFamily z )
----    => Int -- ^ Number of gradient steps per iteration
----    -> Double -- ^ Step size
----    -> GradientPursuit -- ^ Gradient Pursuit Algorithm
----    -> Maybe Int -- ^ Minibatch size (or just full batch)
----    -> Double -- ^ New component ratio (must be between 0 and 1)
----    -> Sample z -- ^ Observations
----    -> Natural # z -- ^ Initial Distribution
----    -> (Natural # Harmonium Tensor z (Categorical e n),[[Double]]) -- ^ (Final Mixture and gradient descent)
----iterativeMixtureOptimization nstps eps grd mnbtch rto zss nz =
----    let nbtch = fromMaybe (length zss) mnbtch
----     in iterativeMixtureOptimization0 nstps eps grd nbtch rto zss (toNullMixture nz,[])
----
------ | An iterative algorithm for training a mixture model.
----iterativeMixtureOptimization0
----    :: forall e n k z . ( Enum e, KnownNat n, KnownNat k, Legendre Natural z, ExponentialFamily z )
----    => Int -- ^ Number of gradient steps per iteration
----    -> Double -- ^ Step size
----    -> GradientPursuit -- ^ Gradient Pursuit Algorithm
----    -> Int -- ^ Minibatch size
----    -> Double -- ^ New component ratio (must be between 0 and 1)
----    -> Sample z -- ^ Observations
----    -> (Natural # Harmonium Tensor z (Categorical e k),[[Double]]) -- ^ Initial Mixture and LLs
----    -> (Natural # Harmonium Tensor z (Categorical e n),[[Double]]) -- ^ Final Mixture and LLs
----iterativeMixtureOptimization0 nstps eps grd nbtch rto zss (hrm0,llss) =
----    let trncrc :: Circuit Identity [SamplePoint z] (Natural # Harmonium Tensor z (Categorical e k),Double)
----        trncrc = accumulateCircuit hrm0 $ proc (zs,hrm) -> do
----            let ll = average $ log . mixtureDensity hrm <$> zs
----                dhrm = stochasticMixtureDifferential zs hrm
----            hrm' <- gradientCircuit eps grd -< joinTangentPair hrm $ vanillaGradient dhrm
----            returnA -< ((hrm,ll),hrm')
----        (hrms,lls) = unzip . runIdentity . streamCircuit trncrc . take nstps . breakEvery nbtch $ cycle zss
----        llss' = lls : llss
----        hrm1 = last hrms
----     in case sameNat (Proxy @ k) (Proxy @ n) of
----          Just Refl -> (hrm1, reverse llss')
----          Nothing -> iterativeMixtureOptimization0 nstps eps grd nbtch rto zss
----              (expandMixture rto hrm1,llss')
----
----expandMixture
----    :: forall e z n
----     . ( Enum e, Legendre Natural z, KnownNat n, ExponentialFamily z )
----    => Double -- ^ Weight fraction for new component (0 < x < 1)
----    -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Current Harmonium
----    -> Natural # Harmonium Tensor z (Categorical e (n+1)) -- ^ Updated Harmonium
----expandMixture rto hrm =
----    let (nzs,nx) = splitMixture hrm
----        sx = toSource nx
----        nwght = density nx (toEnum 0)
----        (mxwght,mxidx) = maximumBy (comparing fst) $ zip (density nx (toEnum 0) : listCoordinates sx) [0..]
----     in if mxidx == nidx
----           then let sx' :: Source # Categorical e (n+1)
----                    sx' = Point . S.cons (rto * nwght) $ coordinates sx
----                    nzs' = S.cons (S.last nzs) nzs
----                 in buildMixture nzs' $ toNatural sx'
----           else let sx' :: Source # Categorical e n
----                    sx' = Point $ S.unsafeUpd (coordinates sx) [(mxidx,(1-rto)*mxwght)]
----                    sx'' :: Source # Categorical e (n+1)
----                    sx'' = Point . S.cons (rto * nwght) $ coordinates sx'
----                    nzs' = S.cons (S.last nzs) nzs
----                 in buildMixture nzs' $ toNatural sx''
----
--
----dualContrastiveDivergence
----    :: forall s f z x
----    . ( Generative Natural z, ExponentialFamily z, ExponentialFamily x, Generative Natural x
----      , Map Mean Natural f x z, Bilinear f z x, Bilinear f x z )
----      => Int -- ^ The number of contrastive divergence steps
----      -> Int -- ^ The number of samples
----      -> Natural # x -- ^ Target marginal
----      -> Natural # Harmonium f z x -- ^ The harmonium
----      -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ The gradient estimate
----dualContrastiveDivergence cdn nsmps prr hrm = do
----    xs <- sample nsmps prr
----    dhrm' <- contrastiveDivergence cdn xs $ transposeHarmonium hrm
----    return $ primalIsomorphism . transposeHarmonium $ dualIsomorphism dhrm'
----
------class FitConjugationParameters (fs :: [* -> * -> *]) (ms :: [*]) where
------    fitConjugationParameters
------        :: Double
------        -> Maybe Int
------        -> Natural # DeepHarmonium fs ms
------        -> Natural # Sum (Tail ms)
------        -> Random s (Natural # Sum (Tail ms))
------
------instance FitConjugationParameters '[] '[m] where
------    {-# INLINE fitConjugationParameters #-}
------    fitConjugationParameters _ _ _ _ = zero
------
------instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f z x, Manifold (Sum ms)
------         , ExponentialFamily n, SampleConjugated fs (n : ms), Generative Natural m
------         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
------  => SampleConjugated (f : fs) (m : n : ms) where
------    {-# INLINE sampleConjugated #-}
------    sampleConjugated rprms dhrm = do
------        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
------            (rprm,rprms') = splitSum rprms
------        (ys,xs) <- fmap hUnzip . sampleConjugated rprms' $ biasBottom rprm dhrm'
------        zs <- samplePoint $ mapReplicatedPoint (pn <+>) (pf >$>* ys)
------        return . hZip zs $ hZip ys xs
------
------
