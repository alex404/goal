This library is built on the *machines* library, and uses Mealy automata in
order to describe processes of various kinds. We often work wish to work
with dynamical systems in Goal, and so processes are defined differently
depending on whether they are deterministic or stochastic. This is because
mathematically, differential equations are easy to define on manifolds,
whereas defining stochastic differential equations on manifolds is far less
straight forward.

Simple mechanical systems are described by their position and velocity at
any time. Formally, we define a mechanical system in terms of positions, and
then consider the tangent bundle of that system in order to consider
position and velocity together. A pendulum, for example, is a manifold of
angles parameterized by a mass and length

    m :: Double
    m = 1
    
    l :: Double
    l = 1
    
    pndl :: Pendulum
    pndl = Pendulum m l

and a particular position and velocity of the pendulum is defined in
generalized coordinates on the phase space as

    x0 :: Partials :#: PhaseSpace Pendulum
    x0 = fromList (Bundle pndl) [1,1]

where `PhaseSpace` is just a synonym for `Bundle Generalized` as defined in
*goal-geometry*, and the `Partials` coordinate system indicates that we're
working with velocities as opposed to momenta.

If we wish to simulate this physical system, we may apply the
`lagrangianFlow` function in order to create a `Flow` under a particular
`ForceField`

    f :: Gravity
    f = Gravity 9.81
    
    flw :: Flow (Partials :#: (PhaseSpace Pendulum))
    flw = lagrangianFlow f x0

We may then generate simulations of this flow by feeding it times through the `stream`
function

    ts :: [Double]
    ts = [0,0.05..]
    
    stps :: [Partials :#: PhaseSpace Pendulum]
    stps = stream ts flw

Alternatively, we may use the `langevinStep` function to define the
step-wise, stochastic dynamics, and use the `chain` function to create a
Markov chain, and thereby simulate a stochastic mechanical system.

Based on the Markov chain definitions, this library also provides Gibbs
sampling algorithms for the exponential family harmoniums defined in
*goal-probability*.

*Scripts*

**gibbs**: Demonstrates Gibbs sampling on a population code.

**ito-process**: Defines and animates a simple Ito process.

**markov-chain**: Defines and animates a simple Markov chain.

**pendulum-simulation**: Defines and animates a pendulum.

**pendulum-vector-field**: Displays the vector field of a pendulum.

**rk4**: Compares the Euler and Runga-Kutta algorithms for simulating
differential equations.
