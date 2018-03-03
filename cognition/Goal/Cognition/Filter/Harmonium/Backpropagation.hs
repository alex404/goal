{-# LANGUAGE Arrows #-}

-- | Training multilayer perceptrons to compute predictions for filtering by
-- contrastive divergence minimization.

module Goal.Cognition.Filter.Harmonium.Backpropagation
    ( -- * Backpropagation
      harmoniumFilterChainBackpropagation
    , harmoniumFilterChainBackpropagation'
    , harmoniumFilterChainRectifiedBackpropagation
    , harmoniumFilterChainRectifiedBackpropagation'
    -- * HarmoniumCircuits
    , harmoniumCircuitTrainer
    , harmoniumCircuitTrainer'
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import Goal.Cognition.Filter.Harmonium


--- Backpropagation ---

harmoniumFilterChainBackpropagation
    :: ( NeuralNetwork f, ExponentialFamily x, StandardGenerative Natural x
       , ExponentialFamily n, StandardGenerative Natural n )
    => Int -- ^ Contrastive Divergence-N
    -> Function Mixture Natural :#: Tensor x (Codomain f) -- ^ Prediction Decoder
    -> Function Mixture Natural :#: Affine n x -- ^ Emission Distribution
    -> Function Mixture Mixture :#: f -- ^ Neural network
    -> Mixture :#: Domain f -- ^ Current Belief (and possibly Action)
    -> Sample n -- ^ Subsequent Observation
    -> RandST r (Mixture :#: Codomain f, Differentials :#: Tangent (Function Mixture Mixture) f)
harmoniumFilterChainBackpropagation cdn dcdy ems g z n' = do
    let (prps,[my']) = feedForward g [z]
    err <- harmoniumFilterContrastiveDivergence cdn dcdy ems my' n'
    return (my',feedBackward g [z] prps [err])

harmoniumFilterChainBackpropagation'
    :: ( NeuralNetwork f, ClosedFormExponentialFamily x, StandardGenerative Natural x
       , ExponentialFamily n, StandardGenerative Natural n )
    => Int -- ^ Contrastive Divergence-N
    -> Function Mixture Natural :#: Tensor x (Codomain f) -- ^ Prediction Decoder
    -> Function Mixture Natural :#: Affine n x -- ^ Emission Distribution
    -> Function Mixture Mixture :#: f -- ^ Neural network
    -> Mixture :#: Domain f -- ^ Current Belief (and possibly Action)
    -> Sample n -- ^ Subsequent Observation
    -> RandST r (Mixture :#: Codomain f, Differentials :#: Tangent (Function Mixture Mixture) f)
harmoniumFilterChainBackpropagation' cdn dcdy ems g z n' = do
    let (prps,[my']) = feedForward g [z]
    err <- harmoniumFilterContrastiveDivergence' cdn dcdy ems my' n'
    return (my',feedBackward g [z] prps [err])

harmoniumFilterChainRectifiedBackpropagation
    :: ( NeuralNetwork f, ExponentialFamily x, StandardGenerative Natural x, ExponentialFamily n )
    => Function Mixture Natural :#: Tensor x (Codomain f) -- ^ Prediction Decoder
    -> Function Mixture Natural :#: Affine n x -- ^ Emission Distribution
    -> Function Mixture Mixture :#: f -- ^ Neural network
    -> Mixture :#: Domain f -- ^ Current Belief (and possibly Action)
    -> Sample n -- ^ Subsequent Observation
    -> RandST r (Mixture :#: Codomain f, Differentials :#: Tangent (Function Mixture Mixture) f)
harmoniumFilterChainRectifiedBackpropagation dcdy ems g z n' = do
    let (prps,[my']) = feedForward g [z]
    err <- harmoniumFilterRectifiedDifferentials dcdy ems my' n'
    return (my',feedBackward g [z] prps [err])

harmoniumFilterChainRectifiedBackpropagation'
    :: ( NeuralNetwork f, ClosedFormExponentialFamily x, ExponentialFamily n )
    => Function Mixture Natural :#: Tensor x (Codomain f) -- ^ Prediction Decoder
    -> Function Mixture Natural :#: Affine n x -- ^ Emission Distribution
    -> Function Mixture Mixture :#: f -- ^ Neural network
    -> Mixture :#: Domain f -- ^ Current Belief (and possibly Action)
    -> Sample n -- ^ Subsequent Observation
    -> RandST r (Mixture :#: Codomain f, Differentials :#: Tangent (Function Mixture Mixture) f)
harmoniumFilterChainRectifiedBackpropagation' dcdy ems g z n' = do
    let (prps,[my']) = feedForward g [z]
        err = harmoniumFilterRectifiedDifferentials' dcdy ems my' n'
    return (my',feedBackward g [z] prps [err])

harmoniumCircuitTrainer
    :: ( NeuralNetwork g, ClosedFormExponentialFamily n, StandardGenerative Natural n
       , ExponentialFamily x, StandardGenerative Natural x )
    => Double
    -> Double
    -> Double
    -> Double
    -> Int
    -> HarmoniumCircuit g n x
    -> RandST s (Mealy (Sample n) (HarmoniumCircuit g n x))
harmoniumCircuitTrainer eps bt1 bt2 rg nstpz (HarmoniumCircuit emsn g0 amtx bmtx dcdy dcdz y0 (Just cdn)) = do
    tflt <- harmoniumTrainingFilter y0 amtx bmtx nstpz (harmoniumFilterChainBackpropagation cdn dcdy emsn)
    return . accumulateMealy g0 $ proc (n,g) -> do
        --g' <- snd ^<< second (adamDescent eps bt1 bt2 rg) <<< tflt -< (n,g)
        dg <- snd ^<< tflt -< (n,g)
        let tst = maximum $ abs <$> listCoordinates dg
        g' <- adamDescent eps bt1 bt2 rg -< trace ("foo" ++ show tst) dg
        returnA -< (HarmoniumCircuit emsn g' amtx bmtx dcdy dcdz y0 (Just cdn),g')
harmoniumCircuitTrainer eps bt1 bt2 rg nstpz (HarmoniumCircuit emsn g0 amtx bmtx dcdy dcdz y0 Nothing) = do
    tflt <- harmoniumTrainingFilter y0 amtx bmtx nstpz (harmoniumFilterChainRectifiedBackpropagation dcdy emsn)
    return . accumulateMealy g0 $ proc (n,g) -> do
        --g' <- snd ^<< second (adamDescent eps bt1 bt2 rg) <<< tflt -< (n,g)
        dg <- snd ^<< tflt -< (n,g)
        let tst = maximum $ abs <$> listCoordinates dg
        g' <- adamDescent eps bt1 bt2 rg -< trace ("bar" ++ show tst) dg
        returnA -< (HarmoniumCircuit emsn g' amtx bmtx dcdy dcdz y0 Nothing,g')

harmoniumCircuitTrainer'
    :: ( NeuralNetwork g, ClosedFormExponentialFamily n, StandardGenerative Natural n
       , ClosedFormExponentialFamily x, StandardGenerative Natural x )
    => Double
    -> Double
    -> Double
    -> Double
    -> Int
    -> HarmoniumCircuit g n x
    -> RandST s (Mealy (Sample n) (HarmoniumCircuit g n x))
harmoniumCircuitTrainer' eps bt1 bt2 rg nstpz (HarmoniumCircuit emsn g0 amtx bmtx dcdy dcdz y0 (Just cdn)) = do
    tflt <- harmoniumTrainingFilter y0 amtx bmtx nstpz (harmoniumFilterChainBackpropagation' cdn dcdy emsn)
    return . accumulateMealy g0 $ proc (n,g) -> do
        g' <- snd ^<< second (adamDescent eps bt1 bt2 rg) <<< tflt -< (n,g)
        returnA -< (HarmoniumCircuit emsn g' amtx bmtx dcdy dcdz y0 (Just cdn),g')
harmoniumCircuitTrainer' eps bt1 bt2 rg nstpz (HarmoniumCircuit emsn g0 amtx bmtx dcdy dcdz y0 Nothing) = do
    tflt <- harmoniumTrainingFilter y0 amtx bmtx nstpz (harmoniumFilterChainRectifiedBackpropagation' dcdy emsn)
    return . accumulateMealy g0 $ proc (n,g) -> do
        g' <- snd ^<< second (adamDescent eps bt1 bt2 rg) <<< tflt -< (n,g)
        returnA -< (HarmoniumCircuit emsn g' amtx bmtx dcdy dcdz y0 Nothing,g')


--- Internal ---


-- | If learning a differential function, multiply the result by 'dt'.
harmoniumFilterContrastiveDivergence
    :: ( ExponentialFamily x, StandardGenerative Natural x, ExponentialFamily n
       , StandardGenerative Natural n, Manifold y )
    => Int
    -> Function Mixture Natural :#: Tensor x y
    -> Function Mixture Natural :#: Affine n x
    -> Mixture :#: y
    -> Sample n
    -> Rand (ST s) (Differentials :#: Tangent Mixture y)
harmoniumFilterContrastiveDivergence cdn dcdy ems y n = do
    let hrm = affineToHarmonium (dcdy >.> y) ems
    gbschn <- gibbsChain hrm n
    let gbs = take (cdn + 1) $ streamChain gbschn
        x0 = sufficientStatistic (domain $ manifold ems) . fst $ head gbs
        xN = sufficientStatistic (domain $ manifold ems) . fst $ last gbs
    return . fromCoordinates (Tangent y) . coordinates $ matrixTranspose dcdy >.> (xN <-> x0)

harmoniumFilterContrastiveDivergence'
    :: ( ClosedFormExponentialFamily x, StandardGenerative Natural x, ExponentialFamily n
       , StandardGenerative Natural n, Manifold y )
    => Int
    -> Function Mixture Natural :#: Tensor x y
    -> Function Mixture Natural :#: Affine n x
    -> Mixture :#: y
    -> Sample n
    -> Rand (ST s) (Differentials :#: Tangent Mixture y)
harmoniumFilterContrastiveDivergence' cdn dcdy ems y n = do
    let hrm = affineToHarmonium (dcdy >.> y) ems
        sn = sufficientStatistic (codomain $ manifold ems) n
    gbschn <- gibbsChain' hrm sn
    let gbs = take (cdn + 1) $ streamChain gbschn
        snN = snd $ last gbs
        project sn' = dualTransition
            $ dcdy >.> y <+> (matrixTranspose . snd $ splitAffine ems) >.> sn'
    return . fromCoordinates (Tangent y) . coordinates $ matrixTranspose dcdy >.> (project snN <-> project sn)

harmoniumFilterRectifiedDifferentials
    :: ( ExponentialFamily x, StandardGenerative Natural x, ExponentialFamily n, Manifold y )
    => Function Mixture Natural :#: Tensor x y
    -> Function Mixture Natural :#: Affine n x
    -> Mixture :#: y
    -> Sample n
    -> RandST s (Differentials :#: Tangent Mixture y)
harmoniumFilterRectifiedDifferentials dcdy ems y n = do
    x0 <- standardGenerate $ dcdy >.> y <+> (matrixTranspose . snd $ splitAffine ems) >.>* n
    xN <- standardGenerate $ dcdy >.> y
    let xm = codomain $ manifold dcdy
    return . fromCoordinates (Tangent y) . coordinates $ matrixTranspose dcdy >.> (sufficientStatistic xm xN <-> sufficientStatistic xm x0)

harmoniumFilterRectifiedDifferentials'
    :: ( ClosedFormExponentialFamily x, ExponentialFamily n, Manifold y )
    => Function Mixture Natural :#: Tensor x y
    -> Function Mixture Natural :#: Affine n x
    -> Mixture :#: y
    -> Sample n
    -> (Differentials :#: Tangent Mixture y)
harmoniumFilterRectifiedDifferentials' dcdy ems y n =
    let t = dualTransition $ dcdy >.> y <+> (matrixTranspose . snd $ splitAffine ems) >.>* n
     in fromCoordinates (Tangent y) . coordinates $ matrixTranspose dcdy >.> (dualTransition (dcdy >.> y) <-> t)



--- Graveyard ---

{-
-- | Given a sample path from the belief parameters and observations, applies
-- 'harmoniumFilterBackpropagation'' to train the given 'NeuralNetwork' on the belief
-- parameter-observation transitions in the path.
harmoniumFilterTrainingEpoch'
    :: ( ClosedFormExponentialFamily x, Generative Natural x, ExponentialFamily n
       , Generative Natural n, ClosedFormExponentialFamily z, Riemannian Natural z
       , ClosedFormExponentialFamily y, Riemannian Natural y )
    => [(Sample n, Mixture :#: z)] -- ^ The sample path
    -> Int -- ^ Contrastive Divergence-N
    -> Function Mixture Mixture :#: Tensor z n -- ^ Observation-Encoding mapping
    -> Function Mixture Natural :#: Tensor x z -- ^ Decoding Distribution
    -> Function Mixture Natural :#: Affine n x -- ^ Emission Distribution
    -> Function Mixture Mixture :#: NeuralNetwork z y z -- ^ Neural network
    -> Rand (ST r) (Differentials :#: Tangent (Function Mixture Mixture) (NeuralNetwork z y z))
harmoniumFilterTrainingEpoch' stps cdn mmtx dcd ems g = do
    let zns' = zip (snd <$> stps) (fst <$> tail stps)
    dlppcfs <- mapM (uncurry $ harmoniumFilterChainBackpropagation' cdn mmtx dcd ems g) zns'
    return $ averagePoint dlppcfs

-- | This function combines backpropagation and contrastive divergence
-- minimization. Contrastive divergence is driven by a 'Harmonium' determined by
-- a prior computed with the given 'NeuralNetwork' and the emission
-- distribution. Backpropagation is then applied to the approximate negative-log
-- likelihood gradient computed with contrastive divergence.
harmoniumFilterFlowBackpropagation
    :: ( ExponentialFamily x, Generative Natural x, ExponentialFamily n, Generative Natural n
       , Backpropagation Mixture Partials f, ClosedFormExponentialFamily z
       , Codomain f ~ Tangent Mixture z )
    => Int -- ^ Contrastive Divergence-N
    -> Function Mixture Natural :#: Tensor x z -- ^ Decoding Distribution
    -> Function Mixture Natural :#: Affine n x -- ^ Emission Distribution
    -> Function Mixture Partials :#: f
    -> Double -- ^ Time Step
    -> Mixture :#: Domain f -- ^ Current Belief (and possibly Action)
    -> Sample n -- ^ Subsequent Observation
    -> RandST r (Differentials :#: Tangent (Function Mixture Partials) f)
harmoniumFilterFlowBackpropagation cdn dcd ems g dt zu n' = do
    let gz = g >.> zu
        z' = gradientStep dt gz
    err <- harmoniumFilterContrastiveDivergence cdn dcd ems z' n'
    return $ regression g [zu] [(dt .>) . fromCoordinates (Tangent gz) $ coordinates err]

-- | In this version we rely on the latent 'Manifold' being a
-- 'ClosedFormExponentialFamily'. In this case the gradient can be approximated
-- slightly more efficiently.
harmoniumFilterFlowBackpropagation'
    :: ( ClosedFormExponentialFamily x, Generative Natural x, ExponentialFamily n
       , Generative Natural n, Backpropagation Mixture Partials f, ClosedFormExponentialFamily z
       , Codomain f ~ Tangent Mixture z )
    => Int -- ^ Contrastive Divergence-N
    -> Function Mixture Natural :#: Tensor x z -- ^ Decoding Distribution
    -> Function Mixture Natural :#: Affine n x -- ^ Emission Distribution
    -> Function Mixture Partials :#: f -- ^ Neural network
    -> Double -- ^ Time step
    -> Mixture :#: Domain f -- ^ Current Belief (and possibly Action)
    -> Sample n -- ^ Subsequent Observation
    -> RandST r (Differentials :#: Tangent (Function Mixture Partials) f)
harmoniumFilterFlowBackpropagation' cdn dcd ems g dt zu n' = do
    let gz = g >.> zu
        z' = gradientStep dt gz
    err <- harmoniumFilterContrastiveDivergence' cdn dcd ems z' n'
    return $ regression g [zu] [(dt .>) . fromCoordinates (Tangent gz) $ coordinates err]


-}
