{-# LANGUAGE TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core hiding (Down)
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition


--- Globals ---


nstps = 100


--- Types ---


data Colours = Red | Green | Blue deriving (Eq, Show)

data Sides = Down | Up deriving (Eq, Show)


--- HMM ---

xcat = [Red, Green, Blue]
xcatm = Categorical xcat

trns :: Colours -> Standard :#: Categorical [Colours]
trns Red = fromList xcatm [0.7,0.25]
trns Green = fromList xcatm [0.3,0.4]
trns Blue = fromList xcatm [0.05,0.25]

ycat = [Down, Up]
ycatm = Categorical ycat

emsn :: Colours -> Standard :#: Categorical [Sides]
emsn Red = fromList ycatm [0.95]
emsn Green = fromList ycatm [0.5]
emsn Blue = fromList ycatm [0.05]

clrs :: Colours -> Colour Double
clrs Red = red
clrs Green = green
clrs Blue = blue

sds :: Sides -> Int
sds Down = -1
sds Up = 1

estpos :: Colours -> Int
estpos Red = 2
estpos Green = 3
estpos Blue = 4


main = do

    let x0 = Green
        p0 = fromList xcatm [0.33,0.33]
    xchn <- runWithSystemRandom $ chain x0 (generate . trns)
    ychn <- runWithSystemRandom $ accumulateRandomFunction0 (generate . emsn)
    let zflt = parametricFilter (discretePrediction trns) (discreteInference emsn) p0

    let xyzchn = filterChain xchn ychn zflt
        (xs,ys,zs) = unzip3 . take nstps $ streamChain xyzchn

    let lyt = execEC $ do

                layout_y_axis . laxis_generate .= scaledIntAxis defaultIntAxis (-2,5)
                let ns = [0..]

                plot . liftEC $ do
                    plot_points_style .= filledPolygon 7 4 False (opaque black)
                    plot_points_values .= zip ns (sds <$> ys)

                sequence_ $ do

                    x <- xcat
                    let ns' = fst <$> filter ((== x) . snd) (zip ns xs)

                    return . plot . liftEC $ do
                        plot_points_style .= filledCircles 5 (opaque $ clrs x)
                        plot_points_values .= zip ns' (repeat 0)

                sequence_ $ do

                    (n,z) <- zip ns zs
                    x <- xcat

                    return . plot . liftEC $ do
                        plot_points_style .= filledPolygon 7 4 False (clrs x `withOpacity` density z x)
                        plot_points_values .= [(n,estpos x)]

    goalRenderableToSVG "cognition" "hidden-markov-model" 1200 200 $ toRenderable (lyt :: Layout Int Int)
