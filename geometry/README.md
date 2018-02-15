In Goal, nearly all numerical objects are described as points on manifolds.
A manifold may be a relatively sophisticated object, like a manifold of
probability distributions or neural networks, but may also indicate simpler
concepts such as vector spaces.
                                                                              
A manifold is a set of points which can be described locally by euclidean
spaces. In geometry, a point is an element of a manifold.  However,
arbitrary types of points (e.g. probability distributions) are difficult to
represent directly on a computer. Therefore points in Goal are always
represented by their coordinates (vectors of real numbers) in terms of a
given chart/coordinate system.
                                                                              
In Goal, charts are are represented by phantom types. Mathematically, charts
are maps between the manifold and the relevant cartesian coordinate system.
However, since we do not represent the points of a manifold explicitly,
we also cannot represent charts explicitly. As such, charts in Goal are used
simply to indicate how to interpret the coordinates of a given point.

Points from Euclidean spaces are the most simple kind of point in Goal, and
`fromList` is the most straightforward way of creating a point. For example,

    x :: Cartesian :#: Euclidean
    x = fromList (Euclidean 2) [1,2]

creates a 2-dimensional euclidean point in cartesian coordinates, and

    x :: Polar :#: Euclidean
    x = fromList (Euclidean 2) [2.23,1.11]

creates approximately the same point in polar coordinates. Coordinate
systems can in turn be automatically changed with instances of the
`transition` function.

This library also provides tools for working with differential geometry,
function spaces, convex analysis, and matrices. See the scripts provided in
this library and further goal-libraries for examples of their application.

*Scripts*

**coordinates**: Demonstrates the conversion between polar and cartesian coordinates.

**gradient-descent**: A simple example of gradient descent using the
differential geometry tools provided by goal-geometry.
