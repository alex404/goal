{-# LANGUAGE Arrows #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Simulation

import Goal.Probability

--- Main ---


main :: IO ()
main = do
    let sz = 5
        ntrl = 20
        nstps = 200
        fps = 10
        x0 = 1
        xs = [1..sz]
        smtx = fromList (markovTensor xs) . concat . replicate sz $ replicate sz (1/fromIntegral sz)

    chn <- runWithSystemRandom $ chain x0 (markovTransition smtx)

    let chainToRenderable ln = toRenderable . execEC $ do

           layout_y_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,sz + 1)

           plot . liftEC $ do
                   plot_lines_style .= solidLine 3 (opaque red)
                   plot_lines_title .= "Markov Chain"
                   plot_lines_values .= [ln]

    let mly = accumulateMealy (0 :: Int) $ proc ((),t) -> do
            x <- chn -< ()
            stps <- chainWindow ntrl -< (t,x)
            returnA -< (chainToRenderable stps,t+1)

    goalRenderablesToAnimation "simulation" "markov-chain"  fps 800 600 . take nstps $ streamChain mly

