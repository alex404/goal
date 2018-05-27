{-# LANGUAGE TypeOperators,FlexibleContexts #-}

--- Imports ---


import MNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Data.Vector as V

--- Globals ---


gridify n lyts =
    let rnblss = breakEvery n $ weights (1,1) . tval . toRenderable <$> lyts
     in gridToRenderable . aboveN $ besideN <$> rnblss

pltdr = mnstdr ++ "/deep/layer1/plot"

--- Functions ---



main = do

    rfl <- goalReadFile mnstlyr1dr "rectified/ryhrm20"
    sfl <- goalReadFile mnstlyr1dr "standard/ryhrm20"
    rnllsfl <- goalReadFile mnstlyr1dr "rectified/rnlls"

    let rry0 :: Natural :#: Replicated Bernoulli
        rnyz0 :: Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
        (rry0,rnyz0) = read rfl
    let sry0 :: Natural :#: Replicated Bernoulli
        snyz0 :: Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
        (sry0,snyz0) = read sfl
        anlls = read rnllsfl

    mflds <- runWithSystemRandom $ marginalReceptiveFields 100 rry0 rnyz0

    let rrflds = receptiveFields rnyz0
        srflds = receptiveFields snyz0

        lrnbl = toRenderable . execEC $ do

           goalLayout

           layout_x_axis . laxis_title .= "Epoch"
           layout_y_axis . laxis_title .= "-Log-Likelihood"

           plot . liftEC $ do

               plot_lines_style .= solidLine 3 (opaque black)
               plot_lines_values .= [zip [(0 :: Int)..] (anlls :: [Double])]

    goalRenderableToPDF pltdr "ll-descent" 500 150 lrnbl
    goalRenderableToPDF pltdr "rrflds" 5000 5000 $ gridify 8 rrflds
    goalRenderableToPDF pltdr "srflds" 5000 5000 $ gridify 8 srflds
    goalRenderableToPDF pltdr "mflds" 5000 5000 $ gridify 10 mflds

{-
    sequence_ [goalRenderableToSVG (pltdr ++ "/rectified-receptive-fields") ("rreceptive-"++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] rrflds]
    sequence_ [goalRenderableToSVG (pltdr ++ "/standard-receptive-fields") ("sreceptive-"++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] srflds]
    sequence_ [goalRenderableToSVG (pltdr ++ "/marginal-fields") ("marginal-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] mflds]
-}

