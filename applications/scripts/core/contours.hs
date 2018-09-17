import Goal.Core

--- Test ---

main :: IO ()
main = do
    let f x y = 0.2 * x + 0.2 * y + sin x + sin y
        sz = 7.5
        nstp = 1000
        rng = (-sz,sz,nstp)
        niso = 20
        cntrs = contours rng rng niso f
        clrs = rgbaGradient (0,0,1,1) (1,0,0,1) niso
        lyt = execEC $ do

            layout_title .= "Contours of a 2D Sin Curve"

            sequence $ do

                ((_,cntr),clr) <- zip cntrs clrs

                return . plot . liftEC $ do

                    plot_lines_style .= solidLine 3 clr
                    plot_lines_values .= cntr

    goalRenderableToSVG "core" "contours" 600 600 $ toRenderable lyt
