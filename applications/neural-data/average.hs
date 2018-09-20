import NeuralData

import Goal.Core

import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV

import Data.Semigroup ((<>))

import qualified Data.List as L


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
        csvpth = clcstr ++ "/datasets.csv"

    Right (_,dsts0) <- CSV.decodeByName <$> BS.readFile csvpth
    let dsts = filterDatasets wrds $ toList dsts0

    (cvss,filss) <- fmap unzip . forM dsts $ \(Dataset dststr) -> do

        let cvpth = clcstr ++ "/analysis/cv/" ++ dststr ++ ".csv"
            fipth = clcstr ++ "/analysis/fi/" ++ dststr ++ ".csv"

        Right (_,cvs) <- CSV.decodeByName <$> BS.readFile cvpth
        Right (_,fis) <- CSV.decodeByName <$> BS.readFile fipth

        return (toList cvs,toList fis)

    let cvs = averageCoVs cvss
        fis = averageVMIs filss


    BS.writeFile (clcstr ++ "/analysis/cv/" ++ "average.csv")
        $ CSV.encodeDefaultOrderedByName cvs

    BS.writeFile (clcstr ++ "/analysis/fi/" ++ "average.csv")
        $ CSV.encodeDefaultOrderedByName fis

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

    print "foo"

    runOpts =<< execParser opts


