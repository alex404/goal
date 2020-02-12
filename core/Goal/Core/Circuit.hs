{-# LANGUAGE Arrows,LambdaCase #-}
-- | A set of functions for working with the 'Arrow' known as a Mealy automata,
-- here referred to as 'Circuit's. Circuits are essentialy a way of building
-- composable fold and iterator operations, where some of the values being
-- processed can be hidden.
module Goal.Core.Circuit
    ( -- * Circuits
    Circuit (Circuit)
    , accumulateFunction
    , accumulateCircuit
    , streamCircuit
    , iterateCircuit
    , loopCircuit
    , loopCircuit'
    , arrM
    -- * Chains
    , Chain
    , chain
    , chainCircuit
    , streamChain
    , iterateChain
    , skipChain
    , skipChain0
    , sortChains
    -- ** Recursive Computations
    , iterateM
    , iterateM'
    ) where

--- Imports ---


-- Unqualified --

import Control.Arrow
import Control.Monad

-- Qualified --

import qualified Control.Category as C
import qualified Data.List as L
import Data.Ord

--- Circuits ---

-- | An arrow which takes an input, monadically produces an output, and updates
-- an (inaccessable) internal state.
newtype Circuit m a b = Circuit (a -> m (b, Circuit m a b))

-- | Takes a function from a value and an accumulator (e.g. just a sum value or
-- an evolving set of parameters for some model) to a value and an accumulator.
-- The accumulator is then looped back into the function, returning a Circuit
-- from a to b, which updates the accumulator every step.
accumulateFunction :: Monad m => acc -> (a -> acc -> m (b,acc)) -> Circuit m a b
accumulateFunction acc f = Circuit $ \a -> do
    (b,acc') <- f a acc
    return (b,accumulateFunction acc' f)

-- | accumulateCircuit takes a Circuit with an accumulating parameter and loops it.
accumulateCircuit :: Monad m => acc -> Circuit m (a,acc) (b,acc) -> Circuit m a b
accumulateCircuit acc0 mly0 = accumulateFunction (acc0,mly0) $ \a (acc,Circuit crcf) -> do
    ((b,acc'),mly') <- crcf (a,acc)
    return (b,(acc',mly'))

-- | Takes a Circuit with an accumulating parameter and loops it, but continues
-- to return the output and calculated parameter.
loopCircuit :: Monad m => acc -> Circuit m (a,acc) (b,acc) -> Circuit m a (b,acc)
loopCircuit acc0 mly0 = accumulateFunction (acc0,mly0) $ \a (acc,Circuit crcf) -> do
    ((b,acc'),mly') <- crcf (a,acc)
    return ((b,acc'),(acc',mly'))

-- | Takes a Circuit over an accumulator and no output and loops it.
loopCircuit' :: Monad m => acc -> Circuit m (a,acc) acc -> Circuit m a acc
loopCircuit' acc0 mly0 = accumulateFunction (acc0,mly0) $ \a (acc,Circuit crcf) -> do
    (acc',mly') <- crcf (a,acc)
    return (acc',(acc',mly'))


-- | Feeds a list of inputs into a Circuit automata and gathers a list of outputs.
streamCircuit :: Monad m => Circuit m a b -> [a] -> m [b]
streamCircuit _ [] = return []
streamCircuit (Circuit mf) (a:as) = do
    (b,crc') <- mf a
    (b :) <$> streamCircuit crc' as

-- | Feeds a list of inputs into a Circuit automata and returns the final
-- output. Throws an error on the empty list.
iterateCircuit :: Monad m => Circuit m a b -> [a] -> m b
iterateCircuit _ [] = error "Empty list fed to iterateCircuit"
iterateCircuit (Circuit mf) [a] = fst <$> mf a
iterateCircuit (Circuit mf) (a:as) = do
    (_,crc') <- mf a
    iterateCircuit crc' as

-- | Turn a monadic function into a circuit.
arrM :: Monad m => (a -> m b) -> Circuit m a b
arrM mf = Circuit $ \a -> do
    b <- mf a
    return (b, arrM mf)


--- Chains ---


-- | A 'Chain' is a form of iterator built on a 'Circuit'.
type Chain m x = Circuit m () x

-- | Creates a 'Chain' from an initial state and a transition function. The
-- first step of the chain returns the initial state, and then continues with
-- generated states.
chain
    :: Monad m
    => x -- ^ The initial state
    -> (x -> m x) -- ^ The transition function
    -> Chain m x -- ^ The resulting 'Chain'
chain x0 mf = accumulateFunction x0 $ \() x -> do
    x' <- mf x
    return (x,x')

-- | Creates a 'Chain' from an initial state and a transition circuit. The
-- first step of the chain returns the initial state, and then continues with
-- generated states.
chainCircuit
    :: Monad m
    => x -- ^ The initial state
    -> Circuit m x x -- ^ The transition circuit
    -> Chain m x -- ^ The resulting 'Chain'
chainCircuit x0 crc = accumulateCircuit x0 $ proc ((),x) -> do
    x' <- crc -< x
    returnA -< (x,x')

-- | Returns a list of the given size of the given 'Chain's output.
streamChain :: Monad m => Int -> Chain m x -> m [x]
streamChain n chn = streamCircuit chn $ replicate (n+1) ()

-- | Returns the given index of the given 'Chain's output.
iterateChain :: Monad m => Int -> Chain m x -> m x
iterateChain 0 (Circuit mf) = fst <$> mf ()
iterateChain k (Circuit mf) = mf () >>= iterateChain (k-1) . snd

-- | Skip every 'n' outputs between each output of the given 'Chain' aftert the first.
skipChain :: Monad m => Int -> Chain m x -> Chain m x
skipChain n (Circuit mf) = Circuit $ \() -> do
    (x',crc') <- mf ()
    return (x', skipChain0 n crc')

-- | Skip every 'n' outputs between each output of the given 'Chain'.
skipChain0 :: Monad m => Int -> Chain m x -> Chain m x
skipChain0 n crc = Circuit $ \() -> do
    (Circuit mf) <- iterateM' n iterator crc
    (x',crc') <- mf ()
    return (x', skipChain0 n crc')
        where iterator (Circuit mf') = snd <$> mf' ()


-- | Iterate a monadic action the given number of times, returning the complete
-- sequence of computations.
iterateM :: Monad m => Int -> (x -> m x) -> x -> m [x]
iterateM n mf x0 = streamChain n $ chain x0 mf

-- | Iterate a monadic action the given number of times.
iterateM' :: Monad m => Int -> (x -> m x) -> x -> m x
iterateM' n mf x0 = iterateChain n $ chain x0 mf

-- | A convenience function for numerically comparing the result of chains.
sortChains
    :: (Monad m, RealFloat x, Ord x)
    => Int -- ^ Number of chain steps
    -> (a -> x) -- ^ Objective function
    -> [Chain m a] -- ^ Chains to test
    -> m ([(x,[a])], [(x,[a])]) -- ^ (Sorted (objective,stream), Infinite/NaN strms)
sortChains nstps objective chns = do
    strms <- mapM (streamChain nstps) chns
    let xstrms = [ (objective $ last strm, strm) | strm <- strms ]
        (dvgs,cvgs) = L.partition (\(obj,_) -> isNaN obj || isInfinite obj) xstrms
        cvgs' = L.sortBy (comparing fst) cvgs
    return (cvgs',dvgs)


--- Instances ---


instance Monad m => C.Category (Circuit m) where
    --id :: Circuit a a
    id = Circuit $ \a -> return (a,C.id)
    --(.) :: Circuit b c -> Circuit a b -> Circuit a c
    (.) = dot
        where dot (Circuit crc1) (Circuit crc2) = Circuit $ \a -> do
                  (b, crcA2') <- crc2 a
                  (c, crcA1') <- crc1 b
                  return (c, crcA1' `dot` crcA2')

instance Monad m => Arrow (Circuit m) where
    --arr :: (a -> b) -> Circuit a b
    arr f = Circuit $ \a -> return (f a, arr f)
    --first :: Circuit a b -> Circuit (a,c) (b,c)
    first (Circuit crcf) = Circuit $ \(a,c) -> do
        (b, crcA') <- crcf a
        return ((b,c), first crcA')

instance Monad m => ArrowChoice (Circuit m) where
    --left :: Circuit a b -> Circuit (Either a c) (Either b c)
    left crcA@(Circuit crcf) = Circuit $
        \case
          Left a -> do
              (b,crcA') <- crcf a
              return (Left b,left crcA')
          Right c -> return (Right c,left crcA)
