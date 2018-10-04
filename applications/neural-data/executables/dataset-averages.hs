{-# LANGUAGE DataKinds #-}

import NeuralData

import Goal.Core
import qualified Goal.Core.Vector.Storable as S

import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV

import Data.Semigroup ((<>))

import qualified Data.List as L



--- Averaging ---

covsToVector :: CoefficientsOfVariation -> S.Vector 4 Double
covsToVector (CoefficientsOfVariation rmn1 rsd1 tcmn1 tcsd1) =
    S.fromTuple (rmn1,rsd1,tcmn1,tcsd1)

vectorToCOVs :: S.Vector 4 Double -> CoefficientsOfVariation
vectorToCOVs v =
    let [rmn1,rsd1,tcmn1,tcsd1] = S.toList v
     in CoefficientsOfVariation rmn1 rsd1 tcmn1 tcsd1

averageCoVs :: [[CoefficientsOfVariation]] -> [CoefficientsOfVariation]
averageCoVs cvsss =
    vectorToCOVs . average . map covsToVector <$> L.transpose cvsss


--- CLI ---


data AverageOpts = AverageOpts String String

cvOpts :: Parser AverageOpts
cvOpts = AverageOpts
    <$> strArgument
        ( help "Which data collection to average" )
    <*> strOption
        ( long "filter" <> short 'f' <> help "Filter out datasets which contain these strings" <> value "")

filterDatasets :: [String] -> [Dataset] -> [Dataset]
filterDatasets wrds dsts =
    let filterfun wrd (Dataset dststr) = not $ L.isSubsequenceOf wrd dststr
        foldfun wrd = filter (filterfun wrd)
     in foldr foldfun dsts wrds

runOpts :: AverageOpts -> IO ()
runOpts (AverageOpts clcstr@"patterson-2013" fltstr) = do

    let wrds = words fltstr
        csvpth = "projects/" ++ clcstr ++ "/datasets.csv"

    Right (_,dsts0) <- CSV.decodeByName <$> BS.readFile csvpth
    let dsts = filterDatasets wrds $ toList dsts0

    putStrLn "Averaging Datasets:"
    sequence_ $ print <$> dsts

    cvss <- forM dsts $ \(Dataset dststr) -> do

        let cvpth = "projects/" ++ clcstr ++ "/analysis/cv/" ++ dststr ++ ".csv"

        print cvpth

        Right (_,cvs) <- CSV.decodeByName <$> BS.readFile cvpth

        return (toList cvs)

    let cvs = averageCoVs cvss

    BS.writeFile ("projects/" ++ clcstr ++ "/analysis/cv/average.csv")
        $ CSV.encodeDefaultOrderedByName cvs

runOpts _ = error "Invalid project"
    --let filterFun (Dataset dststr',_) = not $ L.isSubsequenceOf "pooled" dststr'
    --    alcsvss' = fmap snd . filter filterFun $ zip dsts allcsvss

    --forM_ allcsvss $ \allcsvs -> do

    --    let (cvcvs,ficvs) = unzip allcsvs

    --    BS.writeFile (clcstr ++ "/analysis/cv/" ++ dststr ++ ".csv")
    --        $ CSV.encodeDefaultOrderedByName cvcvs


--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Compute average lines from existing analysis"

    runOpts =<< execParser opts


