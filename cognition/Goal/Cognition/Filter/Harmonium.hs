{-# LANGUAGE RankNTypes,Arrows #-}

-- | Training multilayer perceptrons to compute predictions for filtering by
-- contrastive divergence minimization.

module Goal.Cognition.Filter.Harmonium where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation


--- Mealys ---


harmoniumEncodedFilter
    :: ( Apply Mixture Mixture g, ExponentialFamily n )
    => Mixture :#: Codomain g -- ^ Prior parameters
    -> Function Mixture Mixture :#: Tensor (Domain g) n -- ^ Observation-Posterior mapping
    -> Function Mixture Mixture :#: Tensor (Domain g) (Codomain g) -- ^ Prior-Posterior mapping
    -> (Function Mixture Mixture :#: g) -- ^ Prediction function
    -> Mealy (Sample n) (Mixture :#: Domain g) -- ^ Filter
harmoniumEncodedFilter y0 amtx bmtx g = accumulateFunction y0 $ \n y ->
    let z = amtx >.>* n <+> bmtx >.> y
     in (z, g >.> z)

harmoniumTrainingFilter
    :: (Manifold z, ExponentialFamily n, Manifold y, Manifold g)
    => Mixture :#: y -- ^ Prior parameters
    -> Function Mixture Mixture :#: Tensor z n -- ^ Observation-Posterior mapping
    -> Function Mixture Mixture :#: Tensor z y -- ^ Prior-Posterior mapping
    -> Int -- ^ Number of steps before resetting beliefs (0 means never)
    -> ( forall s . Function Mixture Mixture :#: g -> Mixture :#: z -> Sample n
       -> RandST s (Mixture :#: y, Differentials :#: Tangent (Function Mixture Mixture) g) )
    -> RandST r ( Mealy (Sample n, Function Mixture Mixture :#: g)
        (Mixture :#: z, Differentials :#: Tangent (Function Mixture Mixture) g) ) -- ^ Filter
harmoniumTrainingFilter y0 amtx bmtx kmx rf = do
    f <- accumulateRandomFunction0 (uncurry (uncurry rf))
    return . accumulateMealy (1,Nothing) $ proc ((n',g'),(k,mz)) ->
        if isNothing mz
        then do
            let dg' = zero $ Tangent g'
                z' = amtx >.>* n' <+> bmtx >.> y0
            returnA -< ((z',dg'),(k,Just z'))
        else do
            let z = fromJust mz
            (y1,dg') <- f -< ((g',z),n')
            let (k',y') = if k == kmx
                    then (1,y0)
                    else (k+1,y1)
                z' = amtx >.>* n' <+>  bmtx >.> y'
            returnA -< ((z',dg'),(k',Just z'))

harmoniumCircuitFilter
    :: (NeuralNetwork g, ExponentialFamily n, Manifold x)
    => HarmoniumCircuit g n x
    -> Mealy (Sample n) (Natural :#: x)
harmoniumCircuitFilter (HarmoniumCircuit _ g amtx bmtx _ dcdz y0 _) =
    harmoniumEncodedFilter y0 amtx bmtx g >>^ (dcdz >.>)

data HarmoniumCircuit g n x = HarmoniumCircuit
    { emissionDistribution :: Function Mixture Natural :#: Affine n x
    , neuralNetwork :: Function Mixture Mixture :#: g
    , responseRecoder :: Function Mixture Mixture :#: Tensor (Domain g) n
    , predictionRecoder :: Function Mixture Mixture :#: Tensor (Domain g) (Codomain g)
    , predictionDecoder :: Function Mixture Natural :#: Tensor x (Codomain g)
    , beliefDecoder :: Function Mixture Natural :#: Tensor x (Domain g)
    , priorRates :: Mixture :#: Codomain g
    , maybeCDN :: Maybe Int
    }


-- | Creates a population code which encodes the various parameters of the
-- encoded distribution with orthogonal basis vectors, and where the basis
-- vectors are also orthogonal to the vector of 1st.
orthogonalCode
    :: (Manifold x, Manifold n, Manifold z)
    => z -- ^ Encoding space
    -> Double -- ^ Spread parameter
    -> Function Mixture Natural :#: Affine n x -- ^ Emission distribution
    -> ( Function Mixture Natural :#: Tensor x z
       , Function Mixture Mixture :#: Tensor z n ) -- ^ Belief decoder and observation recoder
orthogonalCode mz tht emsn =
    let (d0,rmnd) = flip divMod 2 . dimension . domain $ manifold emsn
        (d,tailfun) = if rmnd == 1 then (d0+1,tail) else (d0,id)
        n = div (dimension mz) d
        nf = fromIntegral n
        nfrng = fromIntegral <$> [1..n]
        i0 = (nf + 1)/2
        zro = replicate n 0
        ztz0 = [ 1/(nf*tht) * sin (2*pi*(i-i0)/nf) | i <- nfrng ]
        omz0 = [ 1/(nf*tht) * cos (2*pi*(i-i0)/nf) | i <- nfrng ]
        bufferedLists zs0 = [ concat $ replicate k zro ++ [zs0] ++ replicate (d-1-k) zro | k <- [0..d-1] ]
        ztzs = (Natural #) . fromList mz <$> bufferedLists ztz0
        omzs = tailfun $ (Natural #) . fromList mz <$> bufferedLists omz0
        dcdz = Function Mixture Natural # fromRows (domain $ manifold emsn) (ztzs ++ omzs)
        ztzstr0 = [ tht*2 * sin (2*pi*(i-i0)/nf) | i <- nfrng ]
        omzstr0 = [ tht*2 * cos (2*pi*(i-i0)/nf) | i <- nfrng ]
        ztzstrs = (Mixture #) . fromList mz <$> bufferedLists ztzstr0
        omzstrs = tailfun $ (Mixture #) . fromList mz <$> bufferedLists omzstr0
        dcdn = matrixTranspose (snd $ splitAffine emsn)
        (ztns,omns) = splitAt 2 $ toRows dcdn
        amtx = foldl1' (<+>) [str >.< nrw | (str,nrw) <- zip (ztzstrs ++ omzstrs) (ztns ++ omns)]
     in (dcdz,amtx)


{-
averageRateOrthogonalCode2D
    :: (Manifold x, Manifold n, Manifold z)
    => z
    -> Double
    -> Function Mixture Natural :#: Affine n x
    -> ( Function Mixture Natural :#: Tensor x z
       , Function Mixture Mixture :#: Tensor z n )
averageRateOrthogonalCode2D mz tht emsn =
    let n = dimension mz
        nf = fromIntegral n
        nfrng = fromIntegral <$> [1..n]
        i0 = (nf + 1)/2
        ztz = Natural # fromList mz [ 1/(nf*tht) * sin (2*pi*(i-i0)/nf) | i <- nfrng ]
        omz = Natural # fromList mz [ 1/(nf*tht) * cos (2*pi*(i-i0)/nf) | i <- nfrng ]
        dcdz = Function Mixture Natural # fromRows (domain $ manifold emsn) [ztz,omz]
        ztzstr = Mixture # fromList mz [ tht*2 * sin (2*pi*(i-i0)/nf) | i <- nfrng ]
        omzstr = Mixture # fromList mz [ tht*2 * cos (2*pi*(i-i0)/nf) | i <- nfrng ]
        [ztn,omn] = toRows . matrixTranspose . snd $ splitAffine emsn
     in (dcdz,(ztzstr >.< ztn) <+> (omzstr >.< omn))

averageRateOrthogonalCode
    :: (Manifold x, Manifold n, Manifold z)
    => z
    -> Double
    -> Function Mixture Natural :#: Affine n x
    -> ( Function Mixture Natural :#: Tensor x z
       , Function Mixture Mixture :#: Tensor z n )
averageRateOrthogonalCode mz tht emsn =
    let d = dimension . domain $ manifold emsn
        n = div (dimension mz) d
        nf = fromIntegral n
        nfrng = fromIntegral <$> [1..n]
        i0 = (nf + 1)/2
        zro = replicate n 0
        zs0 = [ 1/(nf*tht) * sin (2*pi*(i-i0)/nf) | i <- nfrng ]
        zss = (Natural #) . fromList mz <$> bufferedLists d zro zs0
        dcdz = fromRows (domain $ manifold emsn) zss
        zstrs0 = [ tht*2 * sin (2*pi*(i-i0)/nf) | i <- nfrng ]
        zstrss = (Mixture #). fromList mz <$> bufferedLists d zro zstrs0
        dcdn = matrixTranspose (snd $ splitAffine emsn)
        amtx = foldl1' (<+>) [zstr >.< nrw | (zstr,nrw) <- zip zstrss $ toRows dcdn]
     in (dcdz,amtx)
     -}
