/** @file       documentation.h
 *  @brief      Doxygen documentation for the main page.
 *  @mainpage   Bimetric 3+1 Toolkit
 *
 *  Bimetric 3+1 Toolkit implements the standard 3+1 evolution for
 *  spherically symmetric bimetric and GR spacetimes. <p>
 *  The project goal is the construction of a computational toolkit for the evolution of
 *  3+1 equations of **modified theories of gravity** in spherical symmetry, at the same
 *  time resembling the functionality of the EinsteinToolkit to ease transition to the
 *  full scale projects. The current implementation handles bimetric and GR spacetimes.
 *
 *  Note that the toolkit can handle *any* system of differential equations.
 *  For a pedagogical example how to use the toolkit, see
 *  <a href="http://qft.nu/_private_gowdy">Gowdy spacetimes solver</a>.
 *
 *  @version    0.28, inception 2018-04-23, last modified **2018-09-18 09:05**
 *  @author     Mikica Kocic
 *
 *  @par Source code:
 *  <a href="bimetric-ss.zip" target="_blank">bimetric-ss.zip</a><br>
 *  A3 overview: <a href="bimetric-ss-in-3+1-overview.pdf">Bimetric equations
 *  in standard 3+1 form in spherical symmetry</a><br>
 *
 *  @par Main file:
 *  bim-solver.cpp, see also @ref modList
 *
 *  @note The 3+1 decomposition is based on @cite Kocic:2018ddp.
 *  The numerical algorithms are based on
 *  @cite Press:2007nrec, @cite Baumgarte:2010nr, and @cite Shibata:2015nr.
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  ## Example simulation <a href="example.html">(dust collapse testbed)</a> ##
 *
 *  @htmlonly
 *  <div width="100%" style="background:white;"><iframe width="512" height="514"
 *  style="margin-left:100px;border: 1px solid gray; padding: 2px;
 *  border-radius: 5px; box-shadow: 1px 2px 5px 2px #888888;"
 *  src="https://www.youtube.com/embed/fenVNgHxnK4?rel=0&cc_load_policy=1"
 *  frameborder="0" allowfullscreen></iframe></span> </div>
 *  @endhtmlonly
 *
 *  The simulation is compliant with the results of @cite Nakamura:1980.
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  - - -
 *  # Main Features #
 *  @copydoc features
 *
 *  - - -
 *  # Implementation details #
 *  @copydoc details
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  - - -
 *  # Compiling the code #
 *
 *  The code can use OpenMP for almost all parallel loops on the spatial grid.
 *  These loops are indicated by `OMP_parallel_for`, which expands to:
 *  @code
 *     _Pragma("omp parallel for")   i.e., into: #pragma omp parallel for
 *  @endcode
 *
 *  To compile with enabled parallelism, use the `-fopenmp` switch:
 *  @code
 *     g++ -Wall -m64 -std=c++14 -O3 -fopenmp bim-solver.cpp -o bim-solver
 *  @endcode
 *
 *  The code can also be compiled to be executed on a cluster with the MPI support.
 *  To compile, use the `mpic++` compiler with the `-D_USEMPI` switch:
 *  @code
 *     mpic++  -D_USEMPI  -Wall -m64 ...
 *  @endcode
 *
 *  Then use `mpiexec` to run the resulting executable on a cluster.
 *
 *  @warning The maximal slicing is not compliant with MPI since the boundary value
 *  problem for the slicing differential equation requires access to the whole grid!
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  - - -
 *  # References #
 *  @copydoc citelist
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  @copyright  GNU General Public License (GPLv3).
 *
 *  @page example Testbed: GR dust collapse
 *
 *  The reference simulation is a dust collapse in GR. The initial conditions come from
 *  @cite Nakamura:1980. The initial rest-mass density profile is given in Figure 1 below.
 *  The slicing represented by the proper time of the Eulerian observer is given in
 *  Figure 2. The time variation of the radius of each Lagrange shell with the emergence
 *  of the apparent horizon is given in Figure 3 (cf. Fig. 4 in @cite Nakamura:1980).
 *
 *  - - -
 *  @image html fig-1.png "Figure 1."
 *         The initial rest-mass density profile of the dust.
 *  - - -
 *  @image html fig-2.png "Figure 2."
 *         The slicing represented by the proper time of the Eulerian observer.
 *         The collapse of the lapse before reaching the singularity is clearly visible.
 *  - - -
 *  @image html fig-3.png "Figure 3."
 *         The time variation of the radius of each Lagrange shell (solid lines).
 *         The dashed line shows the apparent horizon.
 *  - - -
 *  The animation of Figure 3 is shown below:
 *
 *  @htmlonly
 *  <div width="100%" style="background:white" align="center"><iframe width="512"
 *  height="514" style="margin-left:0;border: 1px solid gray; padding: 2px;
 *  border-radius: 5px; box-shadow: 1px 2px 5px 2px #888888;"
 *  src="https://www.youtube.com/embed/fenVNgHxnK4?rel=0&cc_load_policy=1"
 *  frameborder="0" allowfullscreen></iframe></span> </div><br>&nbsp;
 *  @endhtmlonly
 *
 *  @page example2 Bimetric dust collapse
 *
 *  Adding the second metric to the GR dust solution introduces a system of coupled ODE
 *  governing two conformal factors of two metrics. In such a case, the ODE for the
 *  initial data becomes:
 *  - - -
 *  @image html eq-1.png
 *  - - -
 *  Here we demanded that the matter distribution profile is the same as in GR case.
 *  This is a generalized
 *  <a href="https://en.wikipedia.org/wiki/Lane-Emden_equation">Lane-Emden</a>
 *  equation. In GR, the Lane-Emden equation occurs in the case of a polytropic fluid
 *  (which is a special case in the above equation for a fixed beta model).
 *  Hence, adding a second metric makes a pressureless fluid in one sector to appear
 *  as the influence of a nontrivial polytropic fluid on both metrics.
 *
 *  An example of the initial conditions is given in Figure 4.
 *  - - -
 *  @image html fig-4.png "Figure 4."
 *  - - -
 *  In Figure 4, the testbed GR initial conditions are shown in magenta for
 *  comparison. The bimetric initial conditions are depicted in blue/red and appear
 *  as a wavy departure from GR. The evolution is under development :)
 *  - - -
 *  @image html fig-5.png "Figure 5."
 *  - - -
 *
 *  @page features  Main Features
 *
 *  ## Formalism ##
 *
 *  - Bimetric equations in standard 3+1 form
 *    (with the evolution of p and the explicit bimetric lapse ratio).
 *
 *  - Matter equations for the perfect fluid in conservative form.
 *
 *  - Equations regularized for spherical symmetry.
 *
 *  ## Gauge setup ##
 *
 *  - Lapse:
 *    * Maximal slicing (boundary value problem at each slice, see maximalSlice.h).
 *    * Bona-Massó slicing condition and the K-driver (evolution).
 *    * Algebraic slicing (normal coordinates, and (1+erf)^{-2}).
 *
 *  - Shift:
 *    * Planned: Minimal distortion and Gamma-driver.
 *
 *  ## Boundary conditions ##
 *
 *  - Imposed parity conditions for local flattness at `r = 0` (on the inner boundary).
 *
 *  - Extrapolation (linear 2nd, or 4th order) or
 *    Sommerfeld outgoing wave radiative condition on the outer boundary.
 *
 *  ## Spatial discretization ##
 *
 *  - 2nd, 4th, or 6th order centered differences on a staggered grid
 *    (see finiteDifferences.h).
 *
 *  - Planned: Upwind differences on shift advection terms.
 *
 *  ## Temporal discretization ##
 *
 *  - Method of Lines:
 *    * Runge-Kutta: 2, 3, and 4 steps
 *    * Iterated Crank-Nicolson (ICN): 2 and 3 steps
 *    * Averaged ICN: 2 and 3 steps
 *    * Generic N-step algorithm with arbitrary coefficients (see methodOfLines.h)
 *
 *  - Kreiss-Oliger dissipation (2nd and 4th order) in all of the evolution equations.
 *
 *  - Courant-Friedrichs-Lewy factor (CFL) 0.5 as default.
 *
 *  ## Grid-driver code ##
 *
 *  - Uniform spatial grid (see gridDriver.h).
 *
 *  ## Numerics ##
 *
 *  - Implemented classes: Matrix, Vector, BandLUDecomposition (band-diagonal matrix
 *    LU decomposition with Crout factorization), CubicSpline (normal cubic spline
 *    interpolation), arbirary FD extrapolation, and arbitrary FD derivatives.
 *
 *  - 64-bit and 128-bit floating point
 *
 *  - C++ code with OpenMP and MPI support adapted for high-performance computing.
 *
 *  - Planned: Transition to Cactus
 *
 *  ## Horizons ##
 *
 *  - Apparent horizon finder.
 *
 *  ## Initial Data ##
 *
 *  - Minkowski GR (opt. with a gauge wave)
 *
 *  - Bimetric Minkowski (opt. with a gauge wave)
 *
 *  - GR collisionless matter (dust) with Gaussian profile (as the main testbed).
 *
 *  - Bimetric Gaussian dust with a "polytropic" conformally flat g and f.
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  @page details Implementation Details
 *
 *  ## Main Class ##
 *
 *   @li BimetricEvolve: Encapsulates a 3+1 evolution solver for bimetric spacetimes.
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  ## Modules across the source files ##
 *  @copydoc modList
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  ## Data Flow ##
 *  @copydoc dataFlow
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  - - -
 *  # Finite Difference Approximation #
 *  @copydoc fdm
 *
 *  @page modList  Module and File Overview
 *
 *  <pre style="margin-left:30px;">
 *     +-------------------------------------------------------+    +-----------------+
 *     |                  Bimetric Evolution                   |    |  Bimetric Model |
 *     +-------------------------------------------------------+    +-----------------+
 *     | bim-solver.cpp   eomEvolution.h    eomGauge.h         |    | bimetricModel.h |
 *     | maximalSlice.h   eomMatter.h       eomLapseRatios.h   |    +-----------------+
 *     |                  eomSources.h      eomMiscEquations.h |
 *     |                  eomConstraints.h                     |
 *     +-------------------------------------------------------+
 *  &nbsp;
 *  +-----------------+   +-------------------+   +----------------+   +---------------+
 *  |   Grid Driver   |   |   Initial Data    |   | MoL Integrator |   |  Output Data  |
 *  +-----------------+   +-------------------+   +----------------+   +---------------+
 *  | gridDriver.h    |   | gridInitialData.h |   |  integrator.h  |   | gridOutput.h  |
 *  | gridFunctions.h |   +-------------------+   +----------------+   +---------------+
 *  +-----------------+
 *         +---------------------+      +--------------+      +-----------------+
 *         |  Numerical Methods  |      |  Data Types  |      |    Utilities    |
 *         +---------------------+      +--------------+      +-----------------+
 *         | finiteDifferences.h |      | matrix.h     |      | trackUsedTime.h |
 *         | methodOfLines.h     |      | dataTypes.h  |      | mpiWorld.h      |
 *         | cubicSpline.h       |      +--------------+      | mpiDummyWorld.h |
 *         | bandSol.h           |                            | paramsHolder.h  |
 *         +---------------------+                            | slog.h          |
 *                                                            +-----------------+ </pre>
 *
 *  @page dataFlow  Data Flow
 *
 *  <pre style="margin-left:30px;">
 *  +---------------------------------+              +----------------------------+
 *  | Mathematica Notebook            |------------>>| Initial Data & Parameters  |
 *  | bimetric-ss-in-3+1 (cpp solver) |  Solve the   |  input.dat  (initial data) |
 *  +---------------------------------+  constraint  |  config.ini (parameters)   |
 *            |   Export the             equations   +----------------------------+
 *            |   C++ code                                         |
 *            V                                                    V
 *  +---------------------------------+  eom-std.h   +----------------------------+
 *  | eomEvolution.h eomMatter.h      |------------>>|       bim-solver.cpp       |
 *  | eomSources.h   eomConstraints.h |              +----------------------------+
 *  | eomGauge.h     eomLapseRatios.h |                            |
 *  | eomMiscEquations.h              |                            V
 *  +---------------------------------+              +----------------------------+
 *                                                   |         output.dat         |
 *     See also:  maximalSlice.h                     +----------------------------+
 *     &nbsp;
 *     Numerics:  methodOfLines.h, finiteDifferences.h,
 *                dataTypes.h, matrix.h, bandSol.h, cubicSpline.h  </pre>
 *
 *  <br> <!-- ----------------------------------------------------------------------- -->
 *  ## Regularization of spherically symmetric evolution codes ##
 *
 *  - The equations are regularized according to:
 *    @cite Alcubierre:2004gn, @cite Alcubierre:2010is,
 *    @cite Baumgarte:2013,  @cite CorderoCarrion:2012ic,
 *    @cite Mewes:2018szi, @cite Ruchlin:2017com, and @cite Ruiz:2007rs.
 *
 *  - A nice overview of the parity conditions at `r = 0` can be found in
 *    @cite Ruiz:2007rs.
 *
 *  - Also, parity conditions for components of vectors and tensors can be found in
 *    Table 1 in @cite Baumgarte:2013.
 *
 */
