module NeuralData.Conditional.VonMises
    ( -- * IO
    getMixtureLikelihood
    , strengthenMixtureLikelihood
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

import NeuralData


--- IO ---


getMixtureLikelihood
    :: String -- ^ Experiment name
    -> String -- ^ Dataset
    -> IO (Maybe (NatNumber,NatNumber,[Double]))
getMixtureLikelihood expnm dst =
    fmap read <$> goalReadDataset (Experiment prjnm expnm) (dst ++ "-parameters")

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => [Double]
    -> Natural #> ConditionalMixture (Neurons k) n VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs
