-- | Training multilayer perceptrons to compute predictions for filtering by
-- contrastive divergence minimization.

module Goal.Cognition.Filter.Approximate where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

--- Backpropagation ---

extendedKalmanPrediction
    :: (Coordinates -> Coordinates) -- ^ System Dynamics
    -> (Coordinates -> H.Matrix Double) -- ^ System Dynamics Jacobian
    -> H.Matrix Double -- ^ Step-Dependent Process Noise
    -> (Standard :#: MultivariateNormal) -- ^ Prior Distribution
    -> (Standard :#: MultivariateNormal) -- ^ Prediction Distribution
extendedKalmanPrediction xdnm xjcb xns p =
    let (xmu,xsgma) = splitMultivariateNormal p
        xmu0 = xdnm xmu
        xjcb' = xjcb xmu
        xsgma0 = (xjcb' <> xsgma <> H.tr xjcb')+ xns
    in joinMultivariateNormal xmu0 xsgma0

extendedKalmanPrediction1D
    :: (Double -> Double) -- ^ System Dynamics
    -> (Double -> Double) -- ^ System Dynamics Jacobian
    -> Double -- ^ Step-Dependent Process Noise
    -> (Standard :#: Normal) -- ^ Prior Distribution
    -> (Standard :#: Normal) -- ^ Prediction Distribution
extendedKalmanPrediction1D xdnm xjcb xns p =
    let p' = fromCoordinates (MultivariateNormal 1) $ coordinates p
        xdnm' = C.singleton . xdnm . (C.! 0)
        xjcb' = H.fromLists . (:[]) . (:[]) . xjcb . (C.! 0)
        xns' = H.fromLists [[xns]]
     in fromCoordinates Normal . coordinates $ extendedKalmanPrediction xdnm' xjcb' xns' p'

extendedKalmanInference
    :: (Coordinates -> Coordinates) -- ^ Observation Model
    -> (Coordinates -> H.Matrix Double) -- ^ Observation Model Jacobian
    -> H.Matrix Double -- ^ Observation Noise
    -> Coordinates -- ^ Observation
    -> (Standard :#: MultivariateNormal) -- ^ Prior Distribution
    -> (Standard :#: MultivariateNormal) -- ^ Prediction Distribution
extendedKalmanInference ydnm yjcb yns y p =
    let (xmu0,xsgma0) = splitMultivariateNormal p
        y' = y - ydnm xmu0
        yjcb' = yjcb xmu0
        smtx' = (yjcb' <> xsgma0 <> H.tr yjcb') + yns
        kmtx' = xsgma0 <> H.tr yjcb' <> H.pinv smtx'
        xmu' = xmu0 + (kmtx' H.#> y')
        xsgma' = (H.ident (H.rows kmtx') - (kmtx' <> yjcb')) <> xsgma0
     in joinMultivariateNormal xmu' xsgma'

extendedKalmanInference1D
    :: (Double -> Double) -- ^ Observation Model
    -> (Double -> Double) -- ^ Observation Model Jacobian
    -> Double -- ^ Observation Noise
    -> Double -- ^ Observation
    -> (Standard :#: Normal) -- ^ Prediction Distribution
    -> (Standard :#: Normal) -- ^ Updated Beliefs
extendedKalmanInference1D ydnm yjcb yns y p =
    let p' = fromCoordinates (MultivariateNormal 1) $ coordinates p
        ydnm' = C.singleton . ydnm . (C.! 0)
        yjcb' = H.fromLists . (:[]) . (:[]) . yjcb . (C.! 0)
        yns' = H.fromLists [[yns]]
     in fromCoordinates Normal . coordinates $ extendedKalmanInference ydnm' yjcb' yns' (C.singleton y) p'


{-
informationFilter1D
    :: (Natural :#: Normal) -- ^ Prior Distribution
    -> Double -- ^ Linear Dynamics
    -> Double -- ^ Process Noise
    -> Double -- ^ Observation Model
    -> (Int -> Double) -- ^ Observation Noise
    -> Mealy Double (Natural :#: Normal) -- ^ Belief dynamics
informationFilter1D p0 a b ydnm yns = accumulateFunction (p0,0) $ \y (p,k) ->
    let [tht1,tht2] = listCoordinates p
        tht2str = -2*tht2
        tht1br = a*tht1/(a^2 + tht2str*b)
        tht2br = tht2str/(a^2 + tht2str*b)
        qk = yns k
        tht1' = ydnm*y/qk + tht1br
        tht2' = (ydnm^2)/qk + tht2br
        p' = fromList (manifold p0) [tht1',-0.5*tht2']
     in (p',(p',k+1))

kalmanFilter1D
    :: (Standard :#: Normal) -- ^ Prior Distribution
    -> Double -- ^ Linear Dynamics
    -> Double -- ^ Process Noise
    -> Double -- ^ Observation Model
    -> Double -- ^ Observation Noise
    -> Mealy Double (Standard :#: Normal) -- ^ Belief dynamics
kalmanFilter1D z0 xdnm xns ydnm yns =
    let z0' = fromCoordinates (MultivariateNormal 1) $ coordinates z0
        xdnm' = H.fromLists [[xdnm]]
        xns' = H.fromLists [[xns]]
        ydnm' = H.fromLists [[ydnm]]
        yns' = H.fromLists [[yns]]
     in arr realToFrac >>> kalmanFilter z0' xdnm' xns' ydnm' yns' >>> arr (fromCoordinates Normal . coordinates)

kalmanFilter
    :: (Standard :#: MultivariateNormal) -- ^ Prior Distribution
    -> H.Matrix Double -- ^ Linear Dynamics
    -> H.Matrix Double -- ^ Process Noise
    -> H.Matrix Double -- ^ Observation Model
    -> H.Matrix Double -- ^ Observation Noise
    -> Mealy Coordinates (Standard :#: MultivariateNormal) -- ^ Belief dynamics
kalmanFilter z0 xdnm xns ydnm yns =
    let xdnm' = (xdnm H.#>)
        xjcb' = const xdnm
        xns' = const xns
        ydnm' = (ydnm H.#>)
        yjcb' = const ydnm
        yns' = const yns
     in extendedKalmanFilter z0 xdnm' xjcb' xns' ydnm' yjcb' yns'

extendedKalmanFilter
    :: (Standard :#: MultivariateNormal) -- ^ Prior Distribution
    -> (Coordinates -> Coordinates) -- ^ System Dynamics
    -> (Coordinates -> H.Matrix Double) -- ^ System Dynamics Jacobian
    -> (Int -> H.Matrix Double) -- ^ Step-Dependent Process Noise
    -> (Coordinates -> Coordinates) -- ^ Observation Model
    -> (Coordinates -> H.Matrix Double) -- ^ Observation Model Jacobian
    -> (Int -> H.Matrix Double) -- ^ Step-Dependent Observation Noise
    -> Mealy Coordinates (Standard :#: MultivariateNormal) -- ^ Belief dynamics
extendedKalmanFilter z0 xdnm xjcb xns ydnm yjcb yns = accumulateFunction (z0,0) $ \y (z,stp) ->
    let (xmu,xsgma) = splitMultivariateNormal z
        xmu0 = xdnm xmu
        xjcb' = xjcb xmu
        xsgma0 = (xjcb' <> xsgma <> H.tr xjcb')+ xns stp
        y' = y - ydnm xmu0
        yjcb' = yjcb xmu0
        smtx' = (yjcb' <> xsgma0 <> H.tr yjcb') + yns stp
        kmtx' = xsgma0 <> H.tr xjcb' <> H.pinv smtx'
        (xmu',xsgma') = (xmu0 + (kmtx' H.#> y'), (H.ident (H.rows kmtx') - (kmtx' <> yjcb')) <> xsgma0)
        z' = joinMultivariateNormal xmu' xsgma'
    in (z', (z',stp+1))

extendedInformationFilter1D
    :: (Natural :#: Normal) -- ^ Prior Distribution
    -> (Double -> Double) -- ^ Dynamics
    -> (Double -> Double) -- ^ System Dynamics Derivative
    -> (Int -> Double) -- ^ Step-Dependent Process Noise
    -> (Double -> Double) -- ^ Observation Model
    -> (Double -> Double) -- ^ Observation Model Jacobian
    -> (Int -> Double) -- ^ Step-Dependent Observation Noise (in Natural Coordinates)
    -> Mealy Double (Natural :#: Normal) -- ^ Belief dynamics
extendedInformationFilter1D z0 xdnm xjcb xns ydnm yjcb yns =
    let z0' = fromCoordinates (MultivariateNormal 1) $ coordinates z0
        xdnm' = C.singleton . xdnm . (C.! 0)
        xjcb' = H.fromLists . (:[]) . (:[]) . xjcb . (C.! 0)
        xns' stp = H.fromLists [[xns stp]]
        ydnm' = C.singleton . ydnm . (C.! 0)
        yjcb' = H.fromLists . (:[]) . (:[]) . yjcb . (C.! 0)
        yns' stp = H.fromLists [[yns stp]]
     in arr realToFrac >>> extendedInformationFilter z0' xdnm' xjcb' xns' ydnm' yjcb' yns' >>> arr (fromCoordinates Normal . coordinates)

informationFilter
    :: (Natural :#: MultivariateNormal) -- ^ Prior Distribution
    -> H.Matrix Double -- ^ Linear Dynamics
    -> H.Matrix Double -- ^ Process Noise
    -> H.Matrix Double -- ^ Observation Model
    -> H.Matrix Double -- ^ Observation Noise
    -> Mealy Coordinates (Natural :#: MultivariateNormal) -- ^ Belief dynamics
informationFilter z0 xdnm xns ydnm yns =
    let xdnm' = (xdnm H.#>)
        xjcb' = const xdnm
        xns' = const xns
        ydnm' = (ydnm H.#>)
        yjcb' = const ydnm
        yns' = const yns
     in extendedInformationFilter z0 xdnm' xjcb' xns' ydnm' yjcb' yns'

extendedInformationFilter
    :: (Natural :#: MultivariateNormal) -- ^ Prior Distribution
    -> (Coordinates -> Coordinates) -- ^ System Dynamics
    -> (Coordinates -> H.Matrix Double) -- ^ System Dynamics Jacobian
    -> (Int -> H.Matrix Double) -- ^ Step-Dependent Process Noise
    -> (Coordinates -> Coordinates) -- ^ Observation Model
    -> (Coordinates -> H.Matrix Double) -- ^ Observation Model Jacobian
    -> (Int -> H.Matrix Double) -- ^ Step-Dependent Observation Noise (in Natural Coordinates)
    -> Mealy Coordinates (Natural :#: MultivariateNormal) -- ^ Belief dynamics
extendedInformationFilter z0 xdnm xjcb xns ydnm yjcb yns = accumulateFunction (z0,0) $ \y (z,stp) ->
    let (xmu,xsgma) = splitMultivariateNormal . chart Standard $ transition z
        xmu0 = xdnm xmu
        xjcb' = xjcb xmu
        xsgma0 = (xjcb' <> xsgma <> H.tr xjcb') + xns stp
        (eta0,omg0) = splitMultivariateNormal . chart Natural . transition
            . chart Standard $ joinMultivariateNormal xmu0 xsgma0
        yjcb' = yjcb xmu0
        yns' = yns stp
        omg' = omg0 + (H.tr yjcb' <> yns' <> yjcb')
        eta' = eta0 + H.scale (-2) (H.tr yjcb' <> yns') H.#> (y - ydnm xmu0 + (yjcb' H.#> xmu0))
        z' = joinMultivariateNormal eta' omg'
    in (z', (z',stp+1))
-}
