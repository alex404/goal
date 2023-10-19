#! stack runghc

{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators,TypeFamilies #-}
--- Imports ---


-- Goal --

import Goal.Geometry


--- Globals ---

trng :: Cartesian # Symmetric (Euclidean 3) (Euclidean 3)
trng = fromTuple (1,2,3,4,5,6)

diag :: Cartesian # Diagonal (Euclidean 2) (Euclidean 2)
diag = fromTuple (1,2)

scl :: Cartesian # Scale (Euclidean 2) (Euclidean 2)
scl = singleton 2

tns :: Cartesian # Tensor (Euclidean 2) (Euclidean 3)
tns = fromTuple (1,2,2,3,2,1)

tns' :: Cartesian # Tensor (Euclidean 2) (Euclidean 2)
tns' = fromTuple (1,2,2,3)
--- Main ---


main :: IO ()
main = do
    let tns1 = toTensor trng
        tns2 = toTensor diag
        tns3 = toTensor scl
        shr = schurComplement tns1 tns tns2
    --print $ woodburyMatrix tns2 tns shr

    putStrLn "Scale tests:"
    print $ changeOfBasis tns3 tns'
    print $ changeOfBasis scl tns'

    putStrLn "Diag tests:"
    print $ changeOfBasis tns2 tns'
    print $ changeOfBasis diag tns'

    print $ schurComplement tns1 tns tns2
    print $ schurComplement tns1 tns diag

    print $ woodburyMatrix tns2 tns shr
    print $ woodburyMatrix diag tns shr
    print $ schurComplement tns2 (transpose tns) tns1
