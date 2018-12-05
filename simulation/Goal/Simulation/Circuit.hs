{-# LANGUAGE
    ExplicitNamespaces,
    TypeOperators,
    RankNTypes,
    LambdaCase
    #-}

-- | A set of functions for working with the 'Arrow' known as Mealy automata,
-- here referred to as 'Circuit's.
module Goal.Simulation.Circuit
    ( -- * Circuits
    Circuit (Circuit)
    , type (>>>)
    -- ** Accumulation
    , accumulateFunction
    , accumulateFunction'
    , accumulateFunction0
    -- ** Re-accumulation
    , accumulateCircuit
    , accumulateCircuit'
    , accumulateCircuit0
    -- ** Stochastic Accumulation
    , accumulateRandomFunction
    , accumulateRandomFunction'
    , accumulateRandomFunction0
    -- * Execution
    , stream
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Probability

import qualified Control.Category as C

-- Reexporting --

import qualified Control.Monad.ST as ST


--- Circuits ---

-- | An arrow which takes an input, produces an output, and updates an
-- (inaccessable) internal state.
newtype Circuit a b = Circuit (a -> (b, Circuit a b))

-- | Infix synonym for 'Circuit'.
type (>>>) = Circuit
infixl 2 >>>

-- | accumulateFunction takes a function from a value and an accumulator (e.g. just a sum
-- value or an evolving set of parameters for some model) to a value and an accumulator.
-- The accumulator is then looped back into the function, returning a Circuit from a to
-- b, which updates the accumulator every step.
accumulateFunction :: acc -> (a -> acc -> (b,acc)) -> Circuit a b
{-# INLINE accumulateFunction #-}
accumulateFunction acc f = Circuit $ \a ->
    let (b,acc') = f a acc
     in (b,accumulateFunction acc' f)

-- | accumulateFunction' acts like accumulateFunction but the Circuit automata will
-- continue to return the accumulator as it generates it.
accumulateFunction' :: acc -> (a -> acc -> (b,acc)) -> Circuit a (b,acc)
{-# INLINE accumulateFunction' #-}
accumulateFunction' acc0 f =
    accumulateFunction acc0 f'
    where f' a acc =
              let (b,acc') = f a acc
               in ((b,acc'),acc')

-- | Similar to 'accumulateFunction'', except without only an output accumulator.
accumulateFunction0 :: acc -> (a -> acc -> acc) -> Circuit a acc
{-# INLINE accumulateFunction0 #-}
accumulateFunction0 acc0 f =
    accumulateFunction acc0 f'
    where f' a acc =
              let acc' = f a acc
               in (acc',acc')

-- | accumulateRandomFunction is analogous to accumulateFunction, but takes as an
-- argument a function which returns a random variable.
accumulateRandomFunction :: acc -> (a -> acc -> forall s . Random s (b,acc)) -> Random s' (Circuit a b)
{-# INLINE accumulateRandomFunction #-}
accumulateRandomFunction acc0 rf = do
    rf' <-  accumulateRandomFunction0 (uncurry rf)
    return $ accumulateCircuit acc0 rf'

-- | accumulateRandomFunction' is analogous to accumulateRandomFunction, except
-- it returns the internal state along with the output.
accumulateRandomFunction' :: acc -> (a -> acc -> forall s . Random s (b,acc)) -> Random s' (Circuit a (b,acc))
{-# INLINE accumulateRandomFunction' #-}
accumulateRandomFunction' acc0 rf = do
    rf' <- accumulateRandomFunction0 (uncurry rf)
    return $ accumulateCircuit' acc0 rf'

-- | accumulateRandomFunction' Arrow-izes stateless random functions.
accumulateRandomFunction0 :: (a -> forall s . Random s b) -> Random s' (Circuit a b)
{-# INLINE accumulateRandomFunction0 #-}
accumulateRandomFunction0 rf = do
    sd <- seed
    return $ accumulateFunction sd f
    where f a sd = ST.runST $ do
              gn <- restore sd
              let (Prob sampler) = rf a
              b <- sampler gn
              sd' <- save gn
              return (b,sd')

-- | accumulateCircuit takes a Circuit with an accumulating parameter and loops it.
accumulateCircuit :: acc -> Circuit (a,acc) (b,acc) -> Circuit a b
{-# INLINE accumulateCircuit #-}
accumulateCircuit acc0 mly0 =
    accumulateFunction (acc0,mly0) f
    where f a (acc,Circuit cf) =
              let ((b,acc'),mly') = cf (a,acc)
               in (b,(acc',mly'))

-- | accumulateCircuit except with a returned accumulator.
accumulateCircuit' :: acc -> Circuit (a,acc) (b,acc) -> Circuit a (b,acc)
{-# INLINE accumulateCircuit' #-}
accumulateCircuit' acc0 mly0 =
    accumulateFunction (acc0,mly0) f
    where f a (acc,Circuit cf) =
              let ((b,acc'),mly') = cf (a,acc)
               in ((b,acc'),(acc',mly'))

-- | accumulateCircuit except with a with no output, and only a returned accumulator.
accumulateCircuit0 :: acc -> Circuit (a,acc) acc -> Circuit a acc
{-# INLINE accumulateCircuit0 #-}
accumulateCircuit0 acc0 mly0 =
    accumulateFunction (acc0,mly0) f
    where f a (acc,Circuit cf) =
              let (acc',mly') = cf (a,acc)
               in (acc',(acc',mly'))

{-
mapCircuit :: [Circuit a b] -> Circuit a [b]
mapCircuit mlys = Circuit $ \a ->
    let fs = runCircuit <$> mlys
        (bs,mlys') = unzip [ f $! a | f <- fs ]
     in (bs,mapCircuit mlys')
     -}

--- Streaming ---


-- | Feeds a list of inputs into a Circuit automata and gathers a list of outputs.
stream :: [a] -> Circuit a b -> [b]
{-# INLINE stream #-}
stream as mly = snd $ mapAccumL runCircuit' mly as
    where runCircuit' (Circuit f) a = let (b,f') = f a in (f',b)

{-
-- | A convenience function which Feeds a list of inputs into a Circuit automata
-- and then applies a monadic operation before gathering the containered list of
-- outputs.
streamM :: Monad m => [a] -> Circuit a b -> (b -> m c) -> m [c]
streamM as mly fM = runT . supply as $ auto mly ~> autoM fM

-- | Same as streamM, but where we ignore the resulting containered list.
streamM_ :: Monad m => [a] -> Circuit a b -> (b -> m c) -> m ()
streamM_ as mly fM = runT_ . supply as $ auto mly ~> autoM fM
-}

--- Instances ---


instance C.Category Circuit where
    --id :: Circuit a a
    {-# INLINE id #-}
    id = Circuit $ \a -> (a,C.id)
    --(.) :: Circuit b c -> Circuit a b -> Circuit a c
    {-# INLINE (.) #-}
    (.) = dot
        where dot (Circuit crc1) (Circuit crc2) = Circuit $ \a ->
                  let (b, crcA2') = crc2 a
                      (c, crcA1') = crc1 b
                   in (c, crcA1' `dot` crcA2')

instance Arrow Circuit where
    --arr :: (a -> b) -> Circuit a b
    {-# INLINE arr #-}
    arr f = Circuit $ \a -> (f a, arr f)
    --first :: Circuit a b -> Circuit (a,c) (b,c)
    {-# INLINE first #-}
    first (Circuit crc) = Circuit $ \(a,c) ->
        let (b, crcA') = crc a
         in ((b,c), first crcA')

instance ArrowChoice Circuit where
    --left :: Circuit a b -> Circuit (Either a c) (Either b c)
    {-# INLINE left #-}
    left crcA@(Circuit crc) = Circuit $
        \case
          Left a ->
              let (b,crcA') = crc a
               in (Left b,left crcA')
          Right c -> (Right c,left crcA)

instance ArrowLoop Circuit where
    --loop :: Circuit (a,c) (b,c) -> Circuit a b
    {-# INLINE loop #-}
    loop (Circuit crc) = Circuit $ \a ->
        let ((b,c),crcA') = crc (a,c)
        in (b,loop crcA')
