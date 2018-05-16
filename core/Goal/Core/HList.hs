{-# LANGUAGE UndecidableInstances #-}

module Goal.Core.HList where


--- Imports ---


import GHC.TypeLits
import qualified Goal.Core.Vector.Boxed as B


--- HLists ---


data HList :: [*] -> * where
    Null :: HList '[]
    (:+:) :: a -> HList as -> HList (a ': as)

infixr 6 :+:


--- Basic Type Families ---


type family ToTypeList as where
    ToTypeList (HList as) = as

type family Head as where
    Head (a ': as) = a

type family Tail as where
    Tail (a ': as) = as

type family Append as b where
    Append '[] b = '[b]
    Append (a ': as) b = a ': Append as b

type family ReverseAcc xs acc where
    ReverseAcc '[] acc = acc
    ReverseAcc (x ': xs) acc = ReverseAcc xs (x ': acc)

type family Reverse xs where
    Reverse xs = ReverseAcc xs '[]

type Last as = Head (Reverse as)

type Init as = Reverse (Tail (Reverse as))

--- Kind Function Type Families ---

type family Append3 (as :: [* -> * -> *]) (b :: * -> * -> *) where
    Append3 '[] b = '[b]
    Append3 (a ': as) b = a ': Append3 as b

type family Tail3 (as :: [* -> * -> *]) where
    Tail3 (a : as) = as

type family Head3 (as :: [* -> * -> *]) where
    Head3 (a ': as) = a

type family ReverseAcc3 (fs :: [* -> * -> *]) (acc :: [* -> * -> *]) where
    ReverseAcc3 '[] acc = acc
    ReverseAcc3 (x ': xs) acc = ReverseAcc3 xs (x ': acc)

type family Reverse3 (fs :: [* -> * -> *]) where
    Reverse3 fs = ReverseAcc3 fs '[]

type Last3 (fs :: [* -> * -> *]) = Head3 (Reverse3 fs)

type Init3 (fs :: [* -> * -> *]) = Reverse3 (Tail3 (Reverse3 fs))

--- Classes ---


--reverseAcc (x :+: xs) acc = reverseAcc xs (x :+: acc)

class Reversing (as :: [*]) where
    reverseAcc :: HList as -> HList bs -> HList (ReverseAcc as bs)

instance Reversing '[] where
    reverseAcc _ acc = acc

instance Reversing as => Reversing (a : as) where
    reverseAcc (x :+: xs) acc = reverseAcc xs (x :+: acc)

hReverse :: Reversing as => HList as -> HList (Reverse as)
hReverse as = reverseAcc as Null

class Appending as b where
    append :: HList as -> b -> HList (Append as b)

instance Appending '[] b where
    append _ b = b :+: Null

instance Appending as b => Appending (a ': as) b where
    append (a :+: as) b = a :+: append as b

class Homogeneous as a where
    homogenize :: HList (a ': as) -> [a]

instance Homogeneous '[] a where
    homogenize (a :+: Null) = [a]

instance Homogeneous as a => Homogeneous (a ': as) a where
    homogenize (a :+: as) = a : homogenize as

hZip :: KnownNat k => B.Vector k x -> B.Vector k y -> B.Vector k (HList [x,y])
hZip = B.zipWith (\x y -> x :+: y :+: Null)

hUnzip :: KnownNat k => B.Vector k (HList [x,y]) -> (B.Vector k x, B.Vector k y)
hUnzip xys = B.unzip $ B.map (\(x :+: y :+: Null) -> (x,y)) xys

hHead :: HList xs -> Head xs
hHead (x :+: _) = x

hLast :: Reversing xs => HList xs -> Last xs
hLast xs = hHead $ hReverse xs
