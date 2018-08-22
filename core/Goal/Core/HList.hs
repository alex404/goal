{-# LANGUAGE UndecidableInstances #-}

-- | Yet another implementation of HLists.
module Goal.Core.HList
    ( -- * Types
      HList ((:+:), Null)
      -- ** Families
    , ToTypeList
    , Head
    , Tail
    , Append
    , ReverseAcc
    , Reverse
    , Last
    , Init
    -- ** Function Kind Families
    , Head3
    , Tail3
    , Append3
    , ReverseAcc3
    , Reverse3
    , Last3
    , Init3
    -- * Functions on type-lists
    , append
    , Reversing
    , hReverse
    , hZip
    , hUnzip
    , hZip2
    , hUnzip2
    , hHead
    , hSingleton
    , hTail
    , hLast
    ) where


--- Imports ---


--- HLists ---


-- | The basic HList type.
data HList :: [*] -> * where
    Null :: HList '[]
    (:+:) :: a -> HList as -> HList (a ': as)

infixr 6 :+:


--- Basic Type Families ---

-- | Retrieve the type-list which defines an HList.
type family ToTypeList as where
    ToTypeList (HList as) = as

-- | The first type in the type-list.
type family Head as where
    Head (a ': as) = a

-- | The tail of the type-list.
type family Tail as where
    Tail (a ': as) = as

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

-- | The last element of a type-list.
type Last as = Head (Reverse as)

-- | All but the last element of a type-list.
type Init as = Reverse (Tail (Reverse as))

--- Kind Function Type Families ---

-- | A type-list (of bivariate type-functions) with an element appended.
type family Append3 (as :: [* -> * -> *]) (b :: * -> * -> *) where
    Append3 '[] b = '[b]
    Append3 (a ': as) b = a ': Append3 as b

-- | The tail of a type-list of bivariate type-functions.
type family Tail3 (as :: [* -> * -> *]) where
    Tail3 (a : as) = as

-- | The first type in a type-list of bivariate type-functions.
type family Head3 (as :: [* -> * -> *]) where
    Head3 (a ': as) = a

-- | An accumulator for reversing a type-list.
type family ReverseAcc3 (fs :: [* -> * -> *]) (acc :: [* -> * -> *]) where
    ReverseAcc3 '[] acc = acc
    ReverseAcc3 (x ': xs) acc = ReverseAcc3 xs (x ': acc)

-- | A reversed type-list of bivariate functions.
type family Reverse3 (fs :: [* -> * -> *]) where
    Reverse3 fs = ReverseAcc3 fs '[]

-- | The last element of a type-list of binvariate functions.
type Last3 (fs :: [* -> * -> *]) = Head3 (Reverse3 fs)

-- | All but the last element of a type-list of bivariate functions.
type Init3 (fs :: [* -> * -> *]) = Reverse3 (Tail3 (Reverse3 fs))

--- Classes ---


-- | Reversable HLists.
class Reversing (as :: [*]) where
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

-- | Zips two 'Vector's into a 'Vector' of length-2 'HList's.
hZip :: [x] -> [HList xs] -> [HList (x : xs)]
{-# INLINE hZip #-}
hZip = zipWith (:+:)

-- | Unzips a 'Vector' of length-2 'HList's into two 'Vector's.
hUnzip :: [HList (x : xs)] -> ([x], [HList xs])
{-# INLINE hUnzip #-}
hUnzip = unzip . map (\(x :+: xs) -> (x,xs))

-- | Zips two 'Vector's into a 'Vector' of length-2 'HList's.
hZip2 :: [x] -> [y] -> [HList [x,y]]
{-# INLINE hZip2 #-}
hZip2 = zipWith (\x y -> x :+: y :+: Null)

-- | Unzips a 'Vector' of length-2 'HList's into two 'Vector's.
hUnzip2 :: [HList [x,y]] -> ([x], [y])
{-# INLINE hUnzip2 #-}
hUnzip2 = unzip . map (\(x :+: y :+: Null) -> (x,y))

-- | The first element of an 'HList'.
hSingleton :: x -> HList '[x]
{-# INLINE hSingleton #-}
hSingleton x = x :+: Null

-- | The first element of an 'HList'.
hHead :: HList xs -> Head xs
{-# INLINE hHead #-}
hHead (x :+: _) = x
hHead _ = error "Invalid pattern match in hHead"

-- | The first element of an 'HList'.
hTail :: HList xs -> HList (Tail xs)
{-# INLINE hTail #-}
hTail (_ :+: ys) = ys
hTail _ = error "Invalid pattern match in hHead"

-- | The last element of an 'HList'.
hLast :: Reversing xs => HList xs -> Last xs
hLast xs = hHead $ hReverse xs
