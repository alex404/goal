{-# LANGUAGE DataKinds,ScopedTypeVariables #-}

import Goal.Core
import qualified Goal.Core.Vector.Boxed as B

import Goal.Probability

generatorM :: forall k r . KnownNat k => Proxy k -> Random r Double
generatorM _ = do
    let k = fromIntegral $ natValInt (Proxy :: Proxy k)
    x <- standardNormal
    return $ x + k

--generatorM' :: forall k r . (KnownNat k, KnownNat k', k' <= k) => Proxy k' -> Random r Double
--generatorM' _ = do
--    let k = fromIntegral $ natValInt (Proxy :: Proxy k)
--    x <- standardNormal
--    return $ x + k

main :: IO ()
main = do

    let xs :: B.Vector 100 Int
        xs = B.generateP natValInt

    (mxs :: B.Vector 100 Double) <- realize $ B.generatePM generatorM

    (xs' :: B.Vector 50 Int) <- realize $ subsampleVector xs

    print xs
    print mxs
    print xs'

