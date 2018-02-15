{-# LANGUAGE TypeFamilies,MultiParamTypeClasses,TypeOperators,ExplicitNamespaces,UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation.
module Goal.Probability.ExponentialFamily.NeuralNetwork
    ( -- * Neural Networks
      InterLayer
    , type (<*<)
    , type (:+:)
    , splitInterLayer
    , joinInterLayer
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability.ExponentialFamily

--- Multilayer ---

data InterLayer f g

type (f :+: g) = InterLayer f g
infixr 3 :+:

type (m <*< g) = m <* Codomain g :+: g
infixr 3 <*<


splitInterLayer :: (Manifold m, Manifold n) => Point c (InterLayer m n) x -> (Point c m x, Point c n x)
{-# INLINE splitInterLayer #-}
splitInterLayer (Point xs) =
    let (xms,xns) = splitV xs
     in (Point xms, Point xns)

joinInterLayer :: (Manifold m, Manifold n) => Point c m x -> Point c n x -> Point c (InterLayer m n) x
{-# INLINE joinInterLayer #-}
joinInterLayer (Point xms) (Point xns) =
    Point $ joinV xms xns

instance (Manifold f, Manifold g) => Manifold (InterLayer f g) where
    type Dimension (InterLayer f g) = Dimension f + Dimension g

instance (Map f, Map g, Codomain g ~ Domain f) => Map (InterLayer f g) where
    type Domain (InterLayer f g) = Domain g
    type Codomain (InterLayer f g) = Codomain f

instance (Apply Mean Natural f, Apply Mean Natural g, Transition Natural Mean (Codomain g), Codomain g ~ Domain f)
  => Apply Mean Natural (InterLayer f g) where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = splitInterLayer fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitInterLayer fg
         in f >$> (transition <$> g >$> xs)



{-
import qualified Data.Vector.Storable as C


--- Neural Networks ---

class Apply Mixture Mixture f => EFMLP f where
    type Propagation f :: *
    feedForward
        :: Function Mixture Mixture :#: f
        -> [Mixture :#: Domain f] -- ^ The list of inputs
        -> ([Propagation f], [Mixture :#: Codomain f])
    feedBackward
        :: Function Mixture Mixture :#: f -- ^ The neural network
        -> [Mixture :#: Domain f] -- ^ The inputs
        -> [Propagation f] -- ^ The output layer natural coordinates
        -> [Differentials :#: Tangent Mixture (Codomain f)] -- ^ The computed errors
        -> Differentials :#: Tangent (Function Mixture Mixture) f -- ^ The differentials

{-
-- | Backpropagation algorithm with the mean squared error function.
backpropagation
    :: NeuralNetwork f
    => (Function Mixture Mixture :#: f) -- ^ The neural network
    -> [Sample (Domain f)] -- ^ The input
    -> [Sample (Codomain f)] -- ^ The target outputs
    -> Differentials :#: Tangent (Function Mixture Mixture) f -- ^ The neural network gradient
backpropagation f xs ys =
    let mxs = sufficientStatistic (domain $ manifold f) <$> xs
        mys = sufficientStatistic (codomain $ manifold f) <$> ys
    let (prps,mys') = feedForward f xs
        errs = [ (-1) .> (potentialDifferentials q <-> potentialDifferentials (transitionTo Mixture sp1)) | (my,my') <- zip mys mys' ]
     in feedBackward f xs prps errs
     -}

-- | Backpropagation algorithm with the mean squared error function.
meanSquaredBackpropagation
    :: NeuralNetwork f
    => (Function Mixture Mixture :#: f) -- ^ The neural network
    -> [Mixture :#: Domain f] -- ^ The inputs
    -> [Mixture :#: Codomain f] -- ^ The target outputs
    -> Differentials :#: Tangent (Function Mixture Mixture) f -- ^ The neural network gradient
meanSquaredBackpropagation f xs ts =
    let (prps,ys) = feedForward f xs
        errs = [ fromCoordinates (Tangent z) . coordinates $ z <-> t | (t,z) <- zip ts ys ]
     in feedBackward f xs prps errs



-- | A multilayer perceptron with three layers.
data ThreeLayerPerceptron m n o = ThreeLayerPerceptron m n o deriving (Eq, Read, Show)


--- Functions ---

-- | Splits the 'ThreeLayerPerceptron' into its component affine transformations.
splitThreeLayerPerceptron
    :: (Manifold m, Manifold n)
    => (Function Mixture Mixture :#: ThreeLayerPerceptron m n o) -- ^ Neural Network
    -> (Function Mixture Natural :#: Affine m n, Function Mixture Natural :#: Affine n o) -- ^ Resulting Affine Components
splitThreeLayerPerceptron nnp =
    let (ThreeLayerPerceptron m n o) = manifold nnp
        (cs1,cs2) = C.splitAt (dimension m * dimension n + dimension m) $ coordinates nnp
     in (fromCoordinates (Affine m n) cs1, fromCoordinates (Affine n o) cs2)

-- | Construct a 'ThreeLayerPerceptron' from component affine transformations.
joinThreeLayerPerceptron
    :: (Function Mixture Natural :#: Affine m n) -- ^ Affine function on the hidden layer
    -> (Function Mixture Natural :#: Affine n o) -- ^ Affine function on the input
    -> Function Mixture Mixture :#: ThreeLayerPerceptron m n o -- ^ Neural Network
joinThreeLayerPerceptron aff1 aff2 =
    let (Affine m n) = manifold aff1
        (Affine _ o) = manifold aff2
     in fromCoordinates (ThreeLayerPerceptron m n o) $ coordinates aff1 C.++ coordinates aff2

-- | Feeds an input forward through the network, and returns every step of
-- the computation.
feedForward3LP
    :: (ClosedFormExponentialFamily m, ClosedFormExponentialFamily n, Manifold o)
    => (Function Mixture Mixture :#: ThreeLayerPerceptron m n o) -- ^ The neural network
    -> [Mixture :#: o] -- ^ The list of inputs
    -> ([(Natural :#: n, Mixture :#: n, Natural :#: m)], [Mixture :#: m]) -- ^ The outputs and evaluation steps
feedForward3LP nnp xps =
    let (aff1,aff2) = splitThreeLayerPerceptron nnp
        nyps = aff2 >$> xps
        yps = dualTransition <$> nyps
        nzps = aff1 >$> yps
        zps = dualTransition <$> nzps
     in (zip3 nyps yps nzps, zps)

-- | Given the results of a feed forward application, back propagates a
-- given error (last input) through the network.
feedBackward3LP
    :: (RiemannianExponentialFamily m, RiemannianExponentialFamily n, Manifold o)
    => (Function Mixture Mixture :#: ThreeLayerPerceptron m n o) -- ^ The neural network
    -> [Mixture :#: o] -- ^ The inputs
    -> [(Natural :#: n, Mixture :#: n, Natural :#: m)]
    -> [Differentials :#: Tangent Mixture m] -- ^ The computed errors
    -> Differentials :#: Tangent (Function Mixture Mixture) (ThreeLayerPerceptron m n o) -- ^ The neural network Gradient
feedBackward3LP nnp xps nyynzps errs0 =
    let (nyps,yps,nzps) = unzip3 nyynzps
        errs1 = fromCoordinates (codomain $ manifold nnp) . coordinates <$> errs0
        (aff1,_) = splitThreeLayerPerceptron nnp
        (_,mtx1) = splitAffine aff1
        dmps = zipWith legendreFlat nzps errs1
        dmtx1s = [ dmp >.< yp | (dmp,yp) <- zip dmps yps ]
        errs2 = matrixTranspose mtx1 >$> dmps
        dnps = zipWith legendreFlat nyps errs2
        dmtx2s = [ dnp >.< xp | (dnp,xp) <- zip dnps xps ]
     in fromCoordinates (Tangent nnp) $ coordinates (averagePoint dmps) C.++ coordinates (averagePoint dmtx1s)
            C.++ coordinates (averagePoint dnps) C.++ coordinates (averagePoint dmtx2s)


--- Instances ---


instance (Manifold m, Manifold n, Manifold o) => Manifold (ThreeLayerPerceptron m n o) where
    dimension (ThreeLayerPerceptron m n o) =
        dimension m + dimension m * dimension n + dimension n + dimension n * dimension o

instance (ClosedFormExponentialFamily m, ClosedFormExponentialFamily n, Manifold o)
    => Map (ThreeLayerPerceptron m n o) where
    type Domain (ThreeLayerPerceptron m n o) = o
    domain (ThreeLayerPerceptron _ _ o) = o
    type Codomain (ThreeLayerPerceptron m n o) = m
    codomain (ThreeLayerPerceptron m _ _) = m

instance (RiemannianExponentialFamily m, RiemannianExponentialFamily n, Manifold o)
    => NeuralNetwork (ThreeLayerPerceptron m n o) where
    type Propagation (ThreeLayerPerceptron m n o) = (Natural :#: n, Mixture :#: n, Natural :#: m)
    feedForward = feedForward3LP
    feedBackward = feedBackward3LP

instance (RiemannianExponentialFamily m, RiemannianExponentialFamily n, Manifold o)
    => Apply Mixture Mixture (ThreeLayerPerceptron m n o) where
    (>$>) f xs = snd $ feedForward f xs

-}
{-
--backpropagation :: NeuralNetwork (m ': ms) -> (Mixture :#: m -> Mixture :#: m) -> Differential :#:
backpropagate :: NeuralNetwork (m ': ms) -> Mixture :#: m -> Differential :#: NeuralNetwork (m ': ms)
backpropagate nnp dp =



--- Internal ---


popManifold :: NeuralNetwork (m ': ms) -> (m, NeuralNetwork ms)
popManifold (Layer m ms) = (m,ms)

popNeuralNetwork
    :: (Manifold m, Manifold n, Manifold (NeuralNetwork (n ': ms)))
    => Function Mixture Mixture :#: NeuralNetwork (m ': n ': ms)
    -> (Natural :#: m, Function Mixture Natural :#: Tensor m n, Function Mixture Mixture :#: NeuralNetwork (n ': ms))
popNeuralNetwork nnp =
    let (m,nn') = popManifold $ manifold nnp
        (n,_) = popManifold nn'
        tns = Tensor m n
        css = coordinates nnp
        (mcs,css') = C.splitAt (dimension m) css
        (mtxcs,nncs') = C.splitAt (dimension tns) css'
        mp = fromCoordinates m mcs
        mtx = fromCoordinates tns mtxcs
        nnp' = fromCoordinates nn' nncs'
     in (mp,mtx,nnp')

feedForward
    :: Function Mixture Mixture :#: NeuralNetwork ms
    -> [Mixture :#: Domain (NeuralNetwork ms)]
    -> [Mixture :#: Responses ms]
feedForward nnp0 xps0 =
    recurse nnp0 xps0 [ chart Mixture . fromCoordinates (Responses $ Layer (manifold xp) Nub) | xp <- xps ]
        where recurse nnp xps rss =
                  let (b,mtx,nnp') = popNeuralNetwork nnp
                      yps = nnp' >$> xps
                   in map (dualTransition . (<+> b)) $ mtx >$> ys


feedBackward
    :: [Mixture :#: Codomain (NeuralNetwork ms)]
    -> [Mixture :#: Responses ms]
    -> Differential :#: Tangent (Function Mixture Mixture) (NeuralNetwork ms)
feedBackward = undefined

--- Instances ---


-- Responses --

instance Eq (Responses '[]) where
    (==) _ _ = True

instance (Eq m, Eq (NeuralNetwork ms)) => Eq (Responses (m ': ms)) where
    (==) (Responses (Layer m ms)) (Responses (Layer m' ms'))
        | m == m' = ms == ms'
        | otherwise = False

instance Manifold (Responses '[]) where
    dimension _ = 0


instance (Manifold m, Manifold (NeuralNetwork ms)) => Manifold (Responses (m ': ms)) where
    dimension (Responses (Layer m ms)) =  dimension m + dimension ms


-- NeuralNetwork --

instance Eq (NeuralNetwork '[]) where
    (==) _ _ = True

instance (Eq m, Eq (NeuralNetwork ms)) => Eq (NeuralNetwork (m ': ms)) where
    (==) (Layer m ms) (Layer m' ms')
        | m == m' = ms == ms'
        | otherwise = False

instance Manifold (NeuralNetwork '[]) where
    dimension _ = 0

instance Manifold m => Manifold (NeuralNetwork '[m]) where
    dimension _ = 0

instance (Manifold m, Manifold n, Manifold (NeuralNetwork (n ': ms))) => Manifold (NeuralNetwork (m ': n ': ms)) where
    dimension (Layer m (Layer n ms)) =  dimension m + dimension m * dimension n + dimension (Layer n ms)

instance Manifold m => Map (NeuralNetwork '[m]) where
    type Domain (NeuralNetwork '[m]) = m
    domain (Layer m _) = m
    type Codomain (NeuralNetwork '[m]) = m
    codomain (Layer m _) = m

instance (ExponentialFamily m, Manifold n) => Apply Mixture Mixture (NeuralNetwork '[m,n]) where
    (>$>) p xs =
        let (b,mtx,_) = popNeuralNetwork p
         in map (dualTransition . (<+> b)) $ mtx >$> xs

instance (ExponentialFamily m, Manifold n, Map (NeuralNetwork (n ': ms)))
    => Map (NeuralNetwork (m ': n ': ms)) where
    type Domain (NeuralNetwork (m ': n ': ms)) = Domain (NeuralNetwork (n ': ms))
    domain (Layer _ nn) = domain nn
    type Codomain (NeuralNetwork (m ': n ': ms)) = m
    codomain (Layer m _) = m

instance (ExponentialFamily m, Manifold n, Apply Mixture Mixture (NeuralNetwork (n ': o ': ms)))
    => Apply Mixture Mixture (NeuralNetwork (m ': n ': o ': ms)) where
    (>$>) p xs =
        let (b,mtx,p') = popNeuralNetwork p
            ys = p' >$> xs
         in map (dualTransition . (<+> b)) $ mtx >$> ys
    -}
