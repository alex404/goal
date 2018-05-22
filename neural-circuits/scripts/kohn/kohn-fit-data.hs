{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.NeuralCircuits

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S

import qualified Data.Map as M
import Data.List

--- Globals ---

preflnm,pstflnm,presflnm,pstsflnm :: String
preflnm = "prestms"
pstflnm = "pststms"
presflnm = "prestrm"
pstsflnm = "pststrm"

-- Training --

nepchs :: Int
nepchs = 10000

eps :: Double
eps = -0.005

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8


--- Functions ---


unsafeFromListB :: KnownNat k => [x] -> B.Vector k x
unsafeFromListB = fromJust . B.fromList

mapToSample
    :: KnownNat nn
    => Double
    -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
    -> (B.Vector NStimuli (SamplePoint VonMises), B.Vector NStimuli (SamplePoint (R nn Poisson)))
mapToSample nrm stmmp =
    let xs = unsafeFromListB $ M.keys stmmp
        ys = unsafeFromListB
            [ unsafeFromListB . map snd . sort . M.toList $ round . (/nrm) . genericLength <$> nmp
              | nmp <- M.elems stmmp ]
     in (xs, ys)

streamToRawSum
    :: KnownNat t
    => [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> (B.Vector t Stimulus, B.Vector t Int)
streamToRawSum sstrm =
    let (xs,ns0) = B.unzip $ unsafeFromListB sstrm
     in ( xs, length . M.foldr (++) [] <$> ns0)


--- Main ---


fitData
    :: forall nn t1 t2
    . (KnownNat nn, KnownNat t1, KnownNat t2, 1 <= nn, 1 <= t1, 1 <= t2)
    => KohnExperiment nn t1 t2
    -> IO ()
fitData kxp = do

    let sbdr = kohnProjectPath kxp

    adpt <- getAdaptor kxp

    (prestms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile sbdr preflnm
    (pststms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile sbdr pstflnm
    (prestrm :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile sbdr presflnm
    (pststrm :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile sbdr pstsflnm

    let sps :: S.Vector nn (Source # VonMises)
        sps = fromJust $ S.fromList
            [ Point $ S.doubleton mu 1 | mu <- tail $ range 0 (2*pi) (1 + natValInt (Proxy :: Proxy nn)) ]

    let ppc0 :: Mean ~> Natural # R nn Poisson <* VonMises
        ppc0 = vonMisesPopulationEncoder sps 1

    let (xs,prens) = mapToSample 25 prestms
        pstns = snd $ mapToSample 20 pststms

        prerwsm :: (B.Vector t1 Stimulus, B.Vector t1 Int)
        prerwsm = streamToRawSum prestrm
        pstrwsm :: (B.Vector t2 Stimulus, B.Vector t2 Int)
        pstrwsm = streamToRawSum pststrm

    let cost = stochasticConditionalCrossEntropy xs

    let backprop ns p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs ns p

        preppcs = take nepchs $ vanillaAdamSequence eps b1 b2 rg (backprop prens) ppc0
        pstppcs = take nepchs $ vanillaAdamSequence eps b1 b2 rg (backprop pstns) ppc0
        preppc = last preppcs
        pstppc = last pstppcs

        nlllyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Negative Log-Likelihood"

            plot . liftEC $ do

                plot_lines_title .= "pre"
                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip [(0 :: Int)..] $ cost prens <$> preppcs ]

            plot . liftEC $ do

                plot_lines_title .= "post"
                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ zip [0..] $ cost pstns <$> pstppcs ]

    goalRenderableToSVG sbdr "pre-population" 1200 600 . toRenderable $ vonMisesFits preppc adpt
    goalRenderableToSVG sbdr "post-population" 1200 600 . toRenderable $ vonMisesFits pstppc adpt
    goalRenderableToSVG sbdr "raw-pre-population" 1200 600 . toRenderable $ rawDataFits prerwsm
    goalRenderableToSVG sbdr "raw-post-population" 1200 600 . toRenderable $ rawDataFits pstrwsm
    goalRenderableToSVG sbdr "negative-log-likelihood" 1200 600 . toRenderable $ nlllyt

    sequence_ $  do
        (k,hst) <- zip [(0 :: Int)..] $ rawDataHistograms prerwsm
        let sbdr' = sbdr ++ "/histograms"
            ttl' = "histogram-pre-population-" ++ show k
        return $ goalRenderableToSVG sbdr' ttl'  1200 600 $ toRenderable hst

    sequence_ $  do
        (k,hst) <- zip [(0 :: Int)..] $ rawDataHistograms pstrwsm
        let sbdr' = sbdr ++ "/histograms"
            ttl' = "histogram-post-population-" ++ show k
        return $ goalRenderableToSVG sbdr' ttl'  1200 600 $ toRenderable hst

main :: IO ()
main = do
    fitData experiment112l44
    fitData experiment112l45
    fitData experiment112r35
    fitData experiment112r36
    fitData experiment105r62
    fitData experiment107l114
    fitData experiment112l16
    fitData experiment112r29
    fitData experiment112r32

