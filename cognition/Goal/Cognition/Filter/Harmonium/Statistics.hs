-- | Training multilayer perceptrons to compute predictions for filtering by the method of maximum likelihood.

module Goal.Cognition.Filter.Harmonium.Statistics
    ( -- * Learning Statistics
      predictedVectorField
    , initialPredictedTuningCurves
    , predictionTuningCurves
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability


--- Backpropagation ---


--- Learning Statistics ---

initialPredictedTuningCurves
    :: ( NeuralNetwork g, Propagation g ~ (a,c :#: h,b)
       , ExponentialFamily n, StandardGenerative Natural n, ExponentialFamily x )
    => (Function Mixture Natural :#: Affine n x) -- ^ Emission distribution
    -> (Function Mixture Mixture :#: Tensor (Domain g) n) -- ^ Observation-Encoding mapping
    -> [Sample x]
    -> Int
    -> (Function Mixture Mixture :#: g) -- ^ The neural network
    -> RandST s [[Double]] -- ^ Every sublist is a tuning curve
initialPredictedTuningCurves emsn amtx smps k g = do
    nss <- mapM standardGenerate $ emsn >$>* concat (replicate k <$> smps)
    let (prps,_) = feedForward g $ amtx >$>* nss
        (_,hss,_) = unzip3 prps
        hs = averagePoint <$> breakEvery k hss
    return . transpose $ listCoordinates <$> hs

-- | Screws up when the bins are empty!!!
predictionTuningCurves
    :: ( NeuralNetwork g,  Propagation g ~ (a,c :#: h,b) )
    => [(x, Mixture :#: Domain g)] -- ^ Latent states and beliefs
    -> [x -> Bool] -- ^ Latent state bins
    -> (Function Mixture Mixture :#: g) -- ^ The neural network
    -> [[Double]] -- ^ Average hidden layer at each bin
predictionTuningCurves xzs0 xbls g =
    transpose . map listCoordinates . reverse . fst $ foldl foldFun ([],xzs0) xbls
        where foldFun (ahs,xzs) xbl =
                  let (txzs,fxzs) = partition (xbl . fst) xzs
                      (prps,_) = feedForward g $ snd <$> txzs
                      (_,hs,_) = unzip3 prps
                   in (averagePoint hs:ahs,fxzs)


-- | This function tries to approximate the vector field over the latent states
-- encoded by the 'ThreeLayerPerceptron. This is primarily for plotting purposes, and
-- can help visualize the learning of the 'ThreeLayerPerceptron.
predictedVectorField
    :: ( ClosedFormExponentialFamily n, RiemannianExponentialFamily y, ExponentialFamily x
       , RiemannianExponentialFamily z, Transition Natural Standard x, Sample x ~ [Double])
    => (Function Mixture Natural :#: Tensor x z) -- ^ Decoding distribution
    -> (Function Mixture Natural :#: Affine n x) -- ^ Emission distribution
    -> (Function Mixture Mixture :#: Tensor z n) -- ^ Observation-Encoding mapping
    -> (Function Mixture Mixture :#: ThreeLayerPerceptron z y z) -- ^ The neural network
    -> Double -- ^ Time step
    -> Double -- ^ Vector scale
    -> Int -- ^ x coordinate index
    -> Int -- ^ y coordinate index
    -> [Double] -- ^ Missing coordinates
    -> (Double,Double) -- ^ Input (x,y) coordinates
    -> (Double,Double) -- ^ Means of transitioned (x,y)
predictedVectorField dcd ems mtx g dt scl ix iy xs (x,y) =
    let xs' = take ix xs ++ x : drop (ix+1) xs
        xs'' = take iy xs' ++ y : drop (iy+1) xs'
        sp = transitionTo Standard . (dcd >.>)
            . (g >.>) . (mtx >.>) . dualTransition $ ems >.>* xs''
     in (scl/dt * (coordinate (2*ix) sp - (xs'' !! ix)), scl/dt * (coordinate (2*iy) sp - (xs'' !! iy)))
