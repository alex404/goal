{-# LANGUAGE FlexibleContexts,DataKinds,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import Goal.Geometry
import Goal.Probability

import Data.List

--- Program ---


-- Globals --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

prjdr :: String
prjdr = "extra/recovering-ppcs"

-- PPC --

type NNeurons = 10

mus :: S.Vector NNeurons Double
mus = S.init $ S.range mnx mxx

randomTuningCurves :: Double -> Random s (S.Vector NNeurons (Source # VonMises))
randomTuningCurves kp0 = do
    let knrm :: Source # Normal
        knrm = Point $ S.fromTuple (kp0,kp0/4)
    kps <- sample knrm
    return . S.zipWith (\x y -> Point $ S.doubleton x y) mus $ G.convert kps

randomGains :: Double -> Random s (S.Vector NNeurons Double)
randomGains gn0 = do
    let gnrm :: Source # Normal
        gnrm = Point $ S.fromTuple (gn0,gn0/4)
    G.convert <$> sample gnrm

randomPPC :: Double -> Double -> Random s (Mean ~> Natural # R NNeurons Poisson <* VonMises)
randomPPC kp0 gn0 = do
    sps <- randomTuningCurves kp0
    gns <- randomGains gn0
    return $ vonMisesPopulationEncoder True (Right gns) sps

-- Training --

prr :: Source # Categorical Int 8
prr = Point $ S.replicate (1/8)

samplePPC :: Mean ~> Natural # R NNeurons Poisson <* VonMises
                 -> Random s (Sample NSamples VonMises, Sample NSamples (R NNeurons Poisson))
samplePPC lkl = do
    ns <- sample prr
    let scl = 1/8 * 2 * pi
        xs = (*scl) . fromIntegral <$> ns
    zs <- samplePoint $ lkl >$>* xs
    return (xs,zs)

type NSamples = 50

nepchs :: Int
nepchs = 2000

eps :: Double
eps = -0.02

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8

-- Plot --

pltsmps :: B.Vector 200 Double
pltsmps = B.range mnx mxx

tuningCurveLayout
    :: Mean ~> Natural # R NNeurons Poisson <* VonMises
    -> Maybe (Mean ~> Natural # R NNeurons Poisson <* VonMises)
    -> Maybe [Double]
    -> Layout Double Double
tuningCurveLayout ppc1 mtruppc mnubxs = execEC $ do

    goalLayout

    radiansAbscissa

    layout_x_axis . laxis_title .= "(Preferred) Stimulus"

    layout_y_axis . laxis_title .= "Rate"

    plot . liftEC $ do
        plot_lines_title .= "Model"
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= toList (toList <$> tuningCurves pltsmps ppc1)

    case mtruppc of
      Nothing -> return ()
      (Just truppc) -> plot . liftEC $ do
                         plot_lines_title .= "Target"
                         plot_lines_style .= solidLine 5 (red `withOpacity` 0.5)
                         plot_lines_values .= toList (toList <$> tuningCurves pltsmps truppc)

    case mnubxs of
      Nothing -> return ()
      (Just nubxs) ->
          sequence_ [plot . return $ vlinePlot "" (solidLine 3 $ opaque black) x | x <- nubxs]


-- Main --


main :: IO ()
main = do

    truppc <- realize $ randomPPC 10 1
    ppc0 <- realize $ randomPPC 1 1

    (xs,zs) <- realize $ samplePPC truppc

    let nubxs = nub $ B.toList xs

    let backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs zs p
        ppcs = take nepchs $ vanillaAdamSequence eps b1 b2 rg backprop ppc0
        tclyt0 = tuningCurveLayout ppc0 Nothing Nothing
        tclyt1 = tuningCurveLayout (last ppcs) (Just truppc) (Just nubxs)

    void . goalRenderableToSVG prjdr "initial-ppc" 600 300 $ toRenderable tclyt0
    void . goalRenderableToSVG prjdr "trained-ppc" 600 300 $ toRenderable tclyt1
