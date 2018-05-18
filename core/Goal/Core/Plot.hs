-- | This module provides some additional tools for working with the Chart library.
module Goal.Core.Plot
    ( -- * Module Exports
      module Graphics.Rendering.Chart
    , module Data.Colour
    , module Data.Colour.Names
    , module Data.Colour.SRGB.Linear
    , module Graphics.Rendering.Chart.Backend.Cairo
    , module Graphics.Rendering.Chart.Grid
    , module Graphics.Rendering.Chart.State
    , module Goal.Core.Plot.Contour
    -- * Plots
    -- ** PixMap
    , pixMapPlot
    -- ** Histograms
    , histogram
    , histogramPlot
    , histogramPlot0
    , logHistogramPlot
    , logHistogramPlot0
    -- * Layouts
    , goalLayout
    , goalLayoutLR
    , radiansAbscissa
    , radiansAbscissaLR
    -- ** PixMap
    , pixMapLayout
    -- ** Histogram
    , histogramLayout
    , histogramLayoutLR
    --, logHistogramLayout
    -- * Util
    , rgbaGradient
    , loopRadiansPlotData
    ) where


--- Imports ---

import Control.Monad
import Goal.Core.Util
import Data.List

import Control.Lens.Setter hiding (Identity)

import Numeric


-- Re-exports --

import Graphics.Rendering.Chart hiding (x0,y0,Point,Vector,Matrix,xy)
import Data.Colour
import Data.Colour.Names hiding (tan)
import Data.Colour.SRGB.Linear
import Graphics.Rendering.Chart.Backend.Cairo
import Graphics.Rendering.Chart.State
import Graphics.Rendering.Chart.Grid

-- Scientific --

import Goal.Core.Plot.Contour


--- Util ---

-- | Returns an ordered list of colours useful for plotting.
rgbaGradient
    :: (Double, Double, Double, Double) -- ^ Initial (R,G,B,A)
    -> (Double, Double, Double, Double) -- ^ Final (R,G,B,A)
    -> Int -- ^ Number of steps
    -> [AlphaColour Double] -- ^ List of colours
rgbaGradient (rmn,gmn,bmn,amn) (rmx,gmx,bmx,amx) n =
    zipWith (flip withOpacity) [amn,amn + astp .. amx]
    $ zipWith3 rgb [rmn,rmn + rstp .. rmx] [gmn,gmn + gstp .. gmx] [bmn,bmn + bstp .. bmx]
    where rstp = (rmx - rmn) / fromIntegral (n-1)
          gstp = (gmx - gmn) / fromIntegral (n-1)
          bstp = (bmx - bmn) / fromIntegral (n-1)
          astp = (amx - amn) / fromIntegral (n-1)

-- | Sets the x-axis to range from 0 to 2pi.
radiansAbscissa ::  EC (Layout Double y) ()
radiansAbscissa = do

    let axs = [0,pi/2,pi,3*pi/2,2*pi]
    layout_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layout_x_axis . laxis_override .=
        axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]

-- | Sets the x-axis to range from 0 to 2pi in an LR layout.
radiansAbscissaLR ::  EC (LayoutLR Double y1 y2) ()
radiansAbscissaLR = do

    let axs = [0,pi/2,pi,3*pi/2,2*pi]
    layoutlr_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layoutlr_x_axis . laxis_override .=
        axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]

-- | Given a list of points from 0 to less than 2pi, this function copies the
-- first point to the end of the list, to create the appearance of a wrapped
-- plot.
loopRadiansPlotData :: [(Double,x)] -> [(Double,x)]
loopRadiansPlotData ((x,y):xys) = ((x,y):xys) ++ [(x+2*pi,y)]
loopRadiansPlotData [] = []

-- | Some default plot settings (a matter of taste).
goalLayout :: EC (Layout x y) ()
goalLayout = do
    layout_left_axis_visibility .= AxisVisibility True False True
    layout_bottom_axis_visibility .= AxisVisibility True False True
    layout_y_axis . laxis_override .= axisGridHide
    layout_x_axis . laxis_override .= axisGridHide

-- | Some default plot settings (a matter of taste).
goalLayoutLR :: EC (LayoutLR x y z) ()
goalLayoutLR = do
    layoutlr_left_axis_visibility .= AxisVisibility True False True
    layoutlr_bottom_axis_visibility .= AxisVisibility True False True
    layoutlr_right_axis_visibility .= AxisVisibility True False True
    layoutlr_left_axis . laxis_override .= axisGridHide
    layoutlr_right_axis . laxis_override .= axisGridHide
    layoutlr_x_axis . laxis_override .= axisGridHide

--- Plots ---


-- PixMap --

-- | Returns a pixmap representation of a matrix style list of list of
-- 'AlphaColour' 'Double's. The first point of the first list will be in the
-- upper left corner of the plot. Be advised that the 'Layout' used to render
-- the pixMap should have a reversed y-axis.
pixMapPlot :: [[AlphaColour Double]] -> Plot Int Int
pixMapPlot clrss = foldr1 joinPlot $ do
    (r,clrs) <- zip [0..] clrss
    (c,clr) <- zip [0..] clrs
    return . toPlot . execEC $ do
        plot_fillbetween_style .= solidFillStyle clr
        plot_fillbetween_values .= [(c,(r,r+1)),(c+1,(r,r+1))]

-- | Sets the layout to be appropriate for a pixMap, which creates a box around
-- the pixmap (with the axis lines) and defines the y-axis as reversed for
-- matrix style coordinates.
pixMapLayout :: EC (Layout Int Int) ()
pixMapLayout = do
    layout_top_axis_visibility .= AxisVisibility True False False
    layout_left_axis_visibility .= AxisVisibility True False False
    layout_right_axis_visibility .= AxisVisibility True False False
    layout_bottom_axis_visibility .= AxisVisibility True False False
    layout_y_axis . laxis_reverse .= True

-- Histogram --

-- | Calculates a histogram, and returns the over and underflow. The x-centres
-- are the points in between each calculated bin.
histogram
    :: Int -- ^ Number of bins
    -> Double -- ^ Min range
    -> Double -- ^ Max range
    -> [[Double]] -- ^ Data set
    -> ([(Double,[Int])],[[Double]],[[Double]]) -- ^ ([(x-centre, y-counts)],underflow,overflow)
histogram n mn mx xss =
    let xss' = sort <$> xss
        (unds,xss'') = unzip $ span (< mn) <$> xss'
        (ovfs,xss''') = unzip $ span (> mx) . reverse <$> xss''
        xssin = reverse <$> xss'''
        bns = range mn mx (n+1)
        vls = transpose $ histogramRow (tail $ take n bns) <$> xssin
        stp = (head (tail bns) - head bns) / 2
        hsts = zip (subtract stp <$> tail bns) vls
     in (hsts,unds,ovfs)
    where histogramRow bns [] = replicate (length bns + 1) 0
          histogramRow [] xs = [length xs]
          histogramRow (bn:bns') xs =
              let (hds,xs') = span (< bn) xs
               in genericLength hds : histogramRow bns' xs'

-- | Creates a histogram out of a data set. The data set is a list of list of values, where
-- each sublist is a collection of data along an axis. Under and overflow is
-- returned as output. The bars are centered at the mid point between each pair
-- of bins.
histogramPlot
    :: Int -- ^ Number of bins
    -> Double -- ^ Min range
    -> Double -- ^ Max range
    -> [[Double]] -- ^ Data set
    -> EC (PlotBars Double Int) ([[Double]],[[Double]]) -- ^ New Plot
histogramPlot n mn mx xss = do
    let (hsts,unds,ovfs) = histogram n mn mx xss
    plot_bars_alignment .= BarsCentered
    plot_bars_values .= hsts
    return (unds,ovfs)

-- | Generates a histogram plot where the min and max bin value is taken from
-- the data set.
histogramPlot0 :: Int -> [[Double]] -> EC (PlotBars Double Int) ()
histogramPlot0 n xss = do
    let mx = maximum $ maximum <$> xss
        mn = minimum $ minimum <$> xss
    void $ histogramPlot n mn mx xss

-- | Good defaults for a histogram layout.
histogramLayout :: Int -> Double -> Double -> EC (Layout Double Int) ()
histogramLayout n mn mx =
    layout_x_axis . laxis_generate .= const (histogramAxis n mn mx)

-- | Good defaults for a histogram layoutLR.
histogramLayoutLR :: Int -> Double -> Double -> EC (LayoutLR Double Int b) ()
histogramLayoutLR n mn mx =
    layoutlr_x_axis . laxis_generate .= const (histogramAxis n mn mx)

histogramAxis :: (PlotValue a, RealFloat a) => Int -> a -> a -> AxisData a
histogramAxis n mn mx = do
    let bns = range mn mx (n+1)
        rng = abs $ maximum bns
    makeAxis (labelFun rng) (bns,[],bns)
        where labelFun rng = map (labelFunStep rng)
              labelFunStep rng x
                  | rng >= 1000 = showEFloat (Just 2) x ""
                  | rng <= 0.01 = showEFloat (Just 2) x ""
                  | otherwise = reverse . dropWhile (== '.') . dropWhile (== '0') . reverse $ showFFloat (Just 2) x ""


--- Logified Histograms ---

-- | A logified version of 'histogramPlot'. The y-axis represents the logarithm
-- of the original count + 1. This can help better visualize outliers.
logHistogramPlot
    :: Double -- ^ Log Base
    -> Int -- ^ Number of bins
    -> Double -- ^ Min range
    -> Double -- ^ Max range
    -> [[Double]] -- ^ Data set
    -> EC (PlotBars Double Double) ([[Double]],[[Double]]) -- ^ New Plot
logHistogramPlot bs n mn mx xss = do
    let (hsts,unds,ovfs) = histogram n mn mx xss
        hsts' = [(x, logBase bs . fromIntegral . (+1) <$> ys) | (x,ys) <- hsts]
    plot_bars_alignment .= BarsCentered
    plot_bars_values .= hsts'
    return (unds,ovfs)

-- | Generates a log-histogram plot where the min and max bin value is taken from
-- the data set.
logHistogramPlot0 :: Double -> Int -> [[Double]] -> EC (PlotBars Double Double) ()
logHistogramPlot0 bs n xss = do
    let mx = maximum $ maximum <$> xss
        mn = minimum $ minimum <$> xss
    void $ logHistogramPlot bs n mn mx xss

{-
-- | Base layout for a log-histogram.
logHistogramLayout :: PlotBars Double Double -> Layout Double Double -> Layout Double Double
logHistogramLayout pbrs lyt =
    let vls = pbrs ^. plot_bars_values
        bns = fst <$> vls
        stp = (head (tail bns) - head bns) / 2
        bns' = (head bns - stp) : map (+stp) bns
        cl = fromIntegral . ceiling . maximum . concat $ snd <$> vls
        rng = abs $ maximum bns'
        xLabelFun x
            | rng >= 1000 = showEFloat (Just 2) x ""
            | rng <= 0.01 = showEFloat (Just 2) x ""
            | otherwise = reverse . dropWhile (== '.') . dropWhile (== '0') . reverse $ showFFloat (Just 2) x ""
    in layout_plots %~ (plotBars pbrs:)
        $ layout_y_axis . laxis_generate .~ const (makeAxis yLabelFun ([0..cl],[],[0..cl]))
        $ layout_x_axis . laxis_generate .~ const (makeAxis xLabelFun (bns',[],bns'))
        $ lyt
    where lbs = 10
          yLabelFun 0 = show 0
          yLabelFun x = show (round lbs) ++ "e" ++ show (round x) ++ "-1"
          -}
