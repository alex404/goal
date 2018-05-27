{-# LANGUAGE Arrows,TypeOperators,FlexibleContexts #-}

--- Imports ---


import MNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation


--- Globals ---


-- Initialization --

hrmxs = 256
zm = Replicated Bernoulli $ dghght * dgwdth
xm = Replicated Bernoulli hrmxs
xcat = Categorical [0..9]
caffm = Affine xcat zm
hrmaffm = Affine xcat xm
w0 = Standard # fromList Normal [0,0.0001]

-- Training --

eps = 0.001
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
nbtch = 10
trnepchn = 1000
trnbrnn = 100
nepchs = 20

-- Testing --

tstepchn = 1000

-- Functions --

validate
    :: [(Int,[Bool])]
    -> Function Mixture Natural :#: Affine (Categorical [Int]) (Replicated Bernoulli)
    -> Double
validate lblss aff =
    let classify (l,bls) =
            if l == (fst . maximumBy (comparing snd) . zip [0..] . listCoordinates . dualTransition $ aff >.>* bls)
                   then 1
                   else 0
     in average $ classify <$> lblss

--- Main ---


main = do

    str <- goalReadFile "mnist" "rxrhrmsxshrm"
    let ((_,rhrm),(_,shrm)) = read str :: ((Natural :#: Replicated Bernoulli, Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)),(Natural :#: Replicated Bernoulli, Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)))

    ldgs <- mnistTrainingData
    caff0 <- runWithSystemRandom $ initialize w0 caffm
    raff0 <- runWithSystemRandom $ initialize w0 hrmaffm
    saff0 <- runWithSystemRandom $ initialize w0 hrmaffm

    let cmly aff0 = accumulateMealy0 aff0 $ proc (lblss,aff) -> do

            let daff = glmTrain lblss aff
            aff' <- adamDescent eps bt1 bt2 rg -< daff
            returnA -< aff'


    cldgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (replicateM nbtch $ sampleRandomLabelledMNIST ldgs))
    rldgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (sampleHarmoniumLabelledMNIST nbtch ldgs rhrm))
    sldgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (sampleHarmoniumLabelledMNIST nbtch ldgs shrm))

    let caffs = take nepchs . takeEvery trnepchn . drop trnbrnn . streamChain $ cmly caff0 <<< cldgchn
    let raffs = take nepchs . takeEvery trnepchn . drop trnbrnn . streamChain $ cmly raff0 <<< rldgchn
    let saffs = take nepchs . takeEvery trnepchn . drop trnbrnn . streamChain $ cmly saff0 <<< sldgchn

    ldgs' <- mnistTestData

    clblss <- runWithSystemRandom (replicateM tstepchn (sampleRandomLabelledMNIST ldgs'))
    rlblss <- runWithSystemRandom (sampleHarmoniumLabelledMNIST tstepchn ldgs' rhrm)
    slblss <- runWithSystemRandom (sampleHarmoniumLabelledMNIST tstepchn ldgs' shrm)

    sequence_ $ print <$> transpose [validate clblss <$> caffs, validate rlblss <$> raffs, validate slblss <$> saffs]



{-
   GRAVEYARD


correlations :: [[Bool]] -> [Double]
correlations blss =
    let f bl = if bl then 1 else 0
        smpss = map f <$> blss
        avgs = average <$> transpose smpss
        sds = trace (show $ head <$> smpss) . map sqrt $ zipWith (-) avgs $ (^2) <$> avgs
        correlation i j = (average [ (smps !! i) * (smps !! j) | smps <- smpss] - (avgs !! i) * (avgs !! j))/((sds !! i)*(sds !!j))
        idxs = [0..length (head blss) - 1]
     in [ correlation i j | i <- idxs,j <- idxs, j > i]

harmoniumCorellationsLayout
    :: [Bool]
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> RandST s (Layout Double Int)
harmoniumCorellationsLayout z rhrm shrm = do

    rgbbschn <- gibbsChain rhrm z
    sgbbschn <- gibbsChain shrm z

    let rsmps = take gbbssmpn . drop gbbsbrnn $ fst <$> streamChain rgbbschn
        ssmps = take gbbssmpn . drop gbbsbrnn $ fst <$> streamChain sgbbschn
        rcrs = take crlsn $ correlations rsmps
        scrs = take crlsn $ correlations ssmps
        mn = -1
        mx = 1

        lyt = execEC $ do

            goalLayout
            histogramLayout 20 mn mx

            plot . fmap plotBars . liftEC $ do

                plot_bars_titles .= ["rectified","standard"]
                histogramPlot 20 mn mx [rcrs,scrs]

    return lyt
    -}
