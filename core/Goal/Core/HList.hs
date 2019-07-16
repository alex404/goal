{-# LANGUAGE UndecidableInstances #-}
-- | Yet another implementation of HLists.
module Goal.Core.HList
    ( -- * Types
      HList ((:+:), Null)
      -- ** Families
    , Append
    , Reverse
    -- * Functions on type-lists
    , append
    , hReverse
    , hZip
    , hUnzip
    , hZip2
    , hUnzip2
    , hHead
    , hSingleton
    , hTail
    ) where


--- Imports ---

import Data.Kind

--- HLists ---


-- | The basic HList type.
data HList :: [Type] -> Type where
    Null :: HList '[]
    (:+:) :: a -> HList as -> HList (a ': as)

infixr 6 :+:


--- Basic Type Families ---

-- | A type-list with an element appended.
type family Append as b where
    Append '[] b = '[b]
    Append (a ': as) b = a ': Append as b

-- | An accumulator for reversing a type-list.
type family ReverseAcc xs acc where
    ReverseAcc '[] acc = acc
    ReverseAcc (x ': xs) acc = ReverseAcc xs (x ': acc)

-- | A reversed type-list.
type family Reverse xs where
    Reverse xs = ReverseAcc xs '[]

--- Kind Function Type Families ---

--- Classes ---


-- | Reversable HLists.
class Reversing as where
    -- | Accumulator for reversing an 'HList'.
    reverseAcc :: HList as -> HList bs -> HList (ReverseAcc as bs)

instance Reversing '[] where
    reverseAcc _ acc = acc

instance Reversing as => Reversing (a : as) where
    reverseAcc (x :+: xs) acc = reverseAcc xs (x :+: acc)

-- | Reverses an 'HList'.
hReverse :: Reversing as => HList as -> HList (Reverse as)
hReverse as = reverseAcc as Null

class Appending as b where
    -- | Appends an element to an 'HList'.
    append :: HList as -> b -> HList (Append as b)

instance Appending '[] b where
    append _ b = b :+: Null

instance Appending as b => Appending (a ': as) b where
    append (a :+: as) b = a :+: append as b

class Homogeneous as a where
    -- | Reduces a homogeneous 'HList' to a list.
    homogenize :: HList (a ': as) -> [a]

instance Homogeneous '[] a where
    homogenize (a :+: Null) = [a]

instance Homogeneous as a => Homogeneous (a ': as) a where
    homogenize (a :+: as) = a : homogenize as

-- | Zips a list of elements an 'HList's into a list of 'HList's.
hZip :: [x] -> [HList xs] -> [HList (x : xs)]
{-# INLINE hZip #-}
hZip = zipWith (:+:)

-- | Unzip a list of 'HList's into a list of head elements, and tail 'HList's.
hUnzip :: [HList (x : xs)] -> ([x], [HList xs])
{-# INLINE hUnzip #-}
hUnzip = unzip . map (\(x :+: xs) -> (x,xs))

-- | Zips two lists into a list of 'HList's over two types.
hZip2 :: [x] -> [y] -> [HList [x,y]]
{-# INLINE hZip2 #-}
hZip2 = zipWith (\x y -> x :+: y :+: Null)

-- | Unzip list of 'HList's over two types into a pair of lists.
hUnzip2 :: [HList [x,y]] -> ([x], [y])
{-# INLINE hUnzip2 #-}
hUnzip2 = unzip . map (\(x :+: y :+: Null) -> (x,y))

-- | Converts a value into a singleton 'HList'.
hSingleton :: x -> HList '[x]
{-# INLINE hSingleton #-}
hSingleton x = x :+: Null

-- | The first element of an 'HList'.
hHead :: HList (x : xs) -> x
{-# INLINE hHead #-}
hHead (x :+: _) = x

-- | The tail of an 'HList'.
hTail :: HList (x : xs) -> HList xs
{-# INLINE hTail #-}
hTail (_ :+: ys) = ys
