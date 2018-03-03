{-# LANGUAGE UndecidableInstances #-}

-- | Provides additional simulation tools.

module Goal.Cognition.Filter.Harmonium.NeuralField where

--- Imports ---


-- Goal --

import Goal.Geometry
import Goal.Probability


--- Neural Fields ---

data NeuralField h z = NeuralField h z Double deriving (Eq, Read, Show)

instance (Manifold h, Manifold z) => Manifold (NeuralField h z) where
    dimension (NeuralField h z _) =
        dimension z + 2 * dimension z * dimension h + dimension h

instance (Manifold h, Manifold z) => Map (NeuralField h z) where
    type Domain (NeuralField h z) = z
    domain (NeuralField _ z _) = z
    type Codomain (NeuralField h z) = z
    codomain (NeuralField _ z _) = z

instance (RiemannianExponentialFamily h, RiemannianExponentialFamily z) => Apply Mixture Mixture (NeuralField h z) where
    (>$>) g zs = snd $ feedForward g zs

instance (RiemannianExponentialFamily h, RiemannianExponentialFamily z) => NeuralNetwork (NeuralField h z) where
    type Propagation (NeuralField h z) = (Natural :#: h, Mixture :#: h, Natural :#: z)
    feedForward = feedForwardNF
    feedBackward g zs prps errs =
        let (NeuralField _ _ dt) = manifold g
            g' = neuralFieldToThreeLayerPerceptron g
            dg' = feedBackward g' zs prps errs
         in dt .> fromCoordinates (Tangent g) (coordinates dg')

feedForwardNF
    :: (ClosedFormExponentialFamily z, ClosedFormExponentialFamily h)
    => (Function Mixture Mixture :#: NeuralField h z) -- ^ The neural network
    -> [Mixture :#: z] -- ^ The list of inputs
    -> ([(Natural :#: h, Mixture :#: h, Natural :#: z)], [Mixture :#: z]) -- ^ The outputs and evaluation steps
feedForwardNF g zs =
    let (NeuralField _ _ dt) = manifold g
        (aff1,aff2) = splitThreeLayerPerceptron $ neuralFieldToThreeLayerPerceptron g
        nhs = aff2 >$> zs
        mhs = dualTransition <$> nhs
        nzs' = [ dt .> nz0' <+> dualTransition z | (z,nz0') <- zip zs $ aff1 >$> mhs ]
        mzs' = dualTransition <$> nzs'
     in (zip3 nhs mhs nzs', mzs')


--- Internal ---


neuralFieldToThreeLayerPerceptron
    :: Function Mixture Mixture :#: NeuralField h z
    -> Function Mixture Mixture :#: ThreeLayerPerceptron z h z
neuralFieldToThreeLayerPerceptron g =
    let (NeuralField y z _) = manifold g
     in fromCoordinates (ThreeLayerPerceptron z y z) $ coordinates g


