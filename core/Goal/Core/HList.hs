module Goal.Core.HList where

data HList :: [*] -> * where
    Null :: HList '[]
    (:+:) :: a -> HList as -> HList (a ': as)

infixr 6 :+:

type family ToTypeList as where
    ToTypeList (HList as) = as

type family Head as where
    Head (a ': as) = a

type family Tail as where
    Tail (a ': as) = as

type family Last as where
    Last '[a] = a
    Last (a ': as) = Last as

type family Init as where
    Init '[a] = '[]
    Init (a ': as) = a ': Init as

class Homogeneous as a where
    homogenize :: HList (a ': as) -> [a]

instance Homogeneous '[] a where
    homogenize (a :+: Null) = [a]

instance Homogeneous as a => Homogeneous (a ': as) a where
    homogenize (a :+: as) = a : homogenize as
