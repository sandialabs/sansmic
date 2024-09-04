// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file odesolver.cpp
 * @brief The ODESolver class implementation
 * @details The function of this code is completely explained and documented in
 * the text, *Computer solution of ordinary differential equations, the initial
 * value problem* by L. F. Shampine and M. K. Gordon. Further details on use of
 * this code are available in `Solving ordinary differential equations with
 * ODE, STEP, and INTRP`, by L. F. Shampine and M. K. Gordon, SLA-73-1060.
 */

#include <iostream>
#include <vector>

#include "libsansmic.hpp"

/**
 * @details Create a new ODESolver object with `neqn` = 0 and no `f` object.
 */
sansmic::ODESolver::ODESolver(void) { neqn = 0; }

/**
 * @brief Create a new solver object with a specific ODE in mind.
 * @param n the number of equations
 * @param f the derivative function object
 */
sansmic::ODESolver::ODESolver(int n, Derivable *f) {
  neqn = n;
  dy_dx = f;
  yy = std::vector<double>(neqn + 1, 0.0);
  wt = std::vector<double>(neqn + 1, 0.0);
  p = std::vector<double>(neqn + 1, 0.0);
  yp = std::vector<double>(neqn + 1, 0.0);
  ypout = std::vector<double>(neqn + 1, 0.0);
  for (int i = 0; i < neqn + 1; i++) {
    std::vector<double> tmp = std::vector<double>(17, 0.0);
    phi.push_back(tmp);
  }
  alpha = std::vector<double>(13, 0.0);
  beta = std::vector<double>(13, 0.0);
  sig = std::vector<double>(14, 0.0);
  v = std::vector<double>(13, 0.0);
  w = std::vector<double>(13, 0.0);
  g = std::vector<double>(14, 0.0);
  psi = std::vector<double>(13, 0.0);
  iphase = 1;
  phase1 = true;
  x = 0;
  h = 0;
  hold = 0;
  istart = 1;
  start = true;
  told = 0;
  delsgn = 0;
  iwork1 = 0;
  iwork2 = 0;
  nornd = false;
  iwork3 = 0;
  iwork4 = 0;
  iwork5 = 0;
}

/**
 * @brief Get the number of equations this solver was configured for.
 * @returns number of equations
 */
int sansmic::ODESolver::num_eqn(void) { return neqn; }

/**
 * @brief Integrate a system of neqn first order ODEs
 * @param y solution vector at t
 * @param t independent variable
 * @param tout point at which solution is desired
 * @param relerr relative error tolerance for local error test
 * @param abserr absolute error tolerance for local error test
 * @param iflag status of integration
 * @details @verbatim embed:rst

.. rubric:: Abstract

Subroutine ode integrates a system of :math:`n_{eqn}`
first order ordinary differential equations of the form

.. math::
    \frac{d}{dt}y_i = f(t, y_1, y_2, \dots, y_{n})

with :math:`y_i` given at :math:`t`.

The subroutine integrates from :math:`t` to :math:`t_{out}`.  On return the
parameters in the call list are set for continuing the integration.
the user has only to define a new value :math:`t_{out}` and call
:cpp:func:`~sansmic::ode` again.

The differential equations are actually solved by a suite of codes
:cpp:func:`~sansmic::de`, :cpp:func:`~sansmic::step1`, and
:cpp:func:`~sansmic::intrp`. :cpp:func:`~sansmic::ode` allocates virtual
storage [as class variables] and calls :cpp:func:`~sansmic::de`
:cpp:func:`~sansmic::de` is a supervisor which
directs the solution.
It calls on the routines :cpp:func:`~sansmic::step1` and
:cpp:func:`~sansmic::intrp` to advance the integration and to interpolate at
output points. :cpp:func:`~sansmic::step1` uses a modified divided difference
form of the Adams PECE formulas and local extrapolation.  It adjusts the order
and step size to control the local error per unit step in a generalized sense.
Normally each call to :cpp:func:`~sansmic::step1` advances the solution one
step in the direction of :math:`t_{out}`.  For reasons of efficiency
:cpp:func:`~sansmic::de` integrates beyond :math:`t_{out}` internally, though
never beyond :math:`t+10(t_{out}-t)`, and calls :cpp:func:`~sansmic::intrp` to
interpolate the solution at :math:`t_{out}`.  An option is provided to stop the
integration at :math:`t_{out}` but it should be used only if it is impossible
to continue the integration beyond :math:`t_{out}`.

This code is completely explained and documented in the text,
`Computer solution of ordinary differential equations, the initial value
problem` by L. F. Shampine and M. K. Gordon. for usage details, also see
SLA-73-1060.

The parameters represent:

- ``f`` -- subroutine :math:`f(t,y,\dot{y})` to evaluate derivatives
  :math:`\dot{y}(i)=dy(i)/dt`
- ``neqn`` -- number of equations to be integrated, :math:`n_{eqn}`
- ``y[]`` -- solution vector at :math:`t`
- ``t`` -- independent variable :math:`t`
- ``tout`` -- point at which solution is desired :math:`t_{out}`
- ``relerr``, ``abserr`` -- relative and absolute error tolerances for local
  error test. At each step the code requires :math:`|\epsilon_{local}| \le
  |y|\times \epsilon_{rel} + \epsilon_{abs}` for each component of the local
  error and solution vectors
- ``iflag`` -- indicates status of integration
- ``work`` -- [object] to hold information internal to code which is necessary
  for subsequent calls


First call to ode:
    The user must provide storage in his calling program for the arrays in the
    call list, y(neqn), work(100+21*neqn), iwork(5),declare `f` in an external
    statement, supply the subroutine :math:`f(t,y,\dot{y})` to evaluate
    :math:`dy(i)/dt = \dot{y}(i) = f(t,y(1),y(2),...,y(n_{eqn}))` and
    initialize the parameters:

    - ``neqn`` -- number of equations to be integrated
    - ``y[]`` -- vector of initial conditions
    - ``t`` -- starting point of integration
    - ``tout`` -- point at which solution is desired
    - ``relerr``, ``abserr`` -- relative and absolute local error tolerances
    - ``iflag`` -- +1, -1.  indicator to initialize the code. Normal input is
      +1. The user should set ``iflag==-1`` only if it is impossible to
      continue the integration beyond :math:`t_{out}`. All parameters except
      ``f``, ``neqn`` and
      ``tout`` may be altered by the code on output so must be variables in the
      calling program.


Output from ode:
    - ``neqn`` -- unchanged
    - ``y[]`` -- solution at :math:`t`
    - ``t`` -- last point reached in integration; normal return has ``t ==
      tout``.
    - ``tout`` -- unchanged
    - ``relerr``, ``abserr`` -- normal return has tolerances unchanged.
      ``iflag == 3`` signals tolerances increased
    - ``iflag``:
        - 2 -- normal return.  integration reached :math:`t_{out}`
        - 3 -- integration did not reach :math:`t_{out}` because error
          tolerances too small. ``relerr``, ``abserr`` increased appropriately
          for continuing
        - 4 -- integration did not reach :math:`t_{out}` because more than
          500 steps needed
        - 5 -- integration did not reach :math:`t_{out}` because equations
          appear to be stiff
        - 6 -- integration did not reach :math:`t_{out}` because solution
          vanished making pure relative error impossible. Must use non-zero
          abserr  to continue.
        - 7 -- invalid input parameters (fatal error) the value of ``iflag``
          is returned negative when the input value is negative and the
          integration does not reach ``tout``, i.e., -3, -4, -5, -6.

Subsequent calls to ode:
    Subroutine returns with all information needed
    to continue the integration. If the integration reached :math:`t_{out}`,
    the user need only define a new :math:`t_{out}` and call again. If the
    integration did not reach  :math:`t_{out}` and the user wants to continue
    he just calls again. In the case  ``iflag==6`` , the user must also alter
    the error criterion. The output value of ``iflag`` is the appropriate
    input value for subsequent calls. The only situation in which it should be
    altered is to stop the integration internally at the new :math:`t_{out}`,
    i.e., change output
    ``iflag=2`` to input ``iflag=-2``. Error tolerances may be changed by the
    user before continuing. All other parameters must remain unchanged.

@endverbatim
 */
void sansmic::ODESolver::ode(std::vector<double> &y, double &t, double tout,
                             double &relerr, double &abserr, int &iflag) {
  bool start = false, phase1 = false, nornd = false;

  if (abs(iflag) < 2 || abs(iflag) > 6) {
    // go directly to DE
  } else {
    start = istart > 0.0;
    phase1 = iphase > 0.0;
    nornd = iwork2 != -1;
  }

  de(y, t, tout, relerr, abserr, iflag, phase1, start, iwork1, nornd, iwork3,
     iwork4, iwork5);

  istart = -1;
  if (start) istart = 1;

  iphase = -1;
  if (phase1) iphase = 1;

  iwork2 = -1;
  if (nornd) iwork2 = 1;

  return;
}

/**
 * @brief Differential equation solver
 * @param y solution vector at T
 * @param t independent variable
 * @param tout point at which solution is desired
 * @param relerr relative error tolerance
 * @param abserr absolute error tolerance
 * @param iflag status of integration
 * @param phase1 work variable
 * @param start work variable
 * @param ns work variable
 * @param nornd work variable
 * @param k work variable
 * @param kold work variable
 * @param isnold work variable
 * @details @verbatim embed:rst

.. rubric:: Abstract

ode merely allocates storage for de to relieve the user of the inconvenience of
a long call list and provides the public interface to the DE function.

This code is completely explained and documented in the text,
*Computer solution of ordinary differential equations, the initial
value problem* by L. F. Shampine and M. K. Gordon.

@endverbatim
 */
void sansmic::ODESolver::de(vector<double> &y, double &t, double tout,
                            double &relerr, double &abserr, int &iflag,
                            bool &phase1, bool &start, int &ns, bool &nornd,
                            int &k, int &kold, int &isnold) {
  int maxnum = 500;
  double u = 1.0e-30;
  double fouru = 4.0 * u;
  double eps;
  int isn;
  double del;
  double absdel;
  double tend;
  int nostep;
  int kle4;
  bool stiff = false;
  bool crash = false;
  double releps;
  double abseps;
  int errcode = 0;
  int noop = 0;

  eps = max(relerr, abserr);
  if ((neqn < 1) || (t == tout) || (relerr < 0.0 || abserr < 0.0) ||
      (eps < 0.0) || (iflag == 0)) {
    iflag = 7;
    errcode = sansmic::ODE_IFLAG_SEVEN;
    return;
  }
  isn = sign(1, iflag);
  iflag = abs(iflag);
  if (iflag == 1) {
    noop = 1;
  } else if (t != told) {
    iflag = 7;
    errcode = sansmic::ODE_IFLAG_SEVEN;
    return;
  } else if ((iflag >= 2 && iflag <= 5) || (iflag == 6 && abserr > 0)) {
    noop = 2;
  } else if (iflag != 6 && iflag != 7) {
    iflag = 7;
    errcode = sansmic::ODE_IFLAG_SEVEN;
    return;
  } else if (iflag == 6) {
    errcode = sansmic::ODE_IFLAG_SIX;
    return;
  } else if (iflag == 7) {
    errcode = sansmic::ODE_IFLAG_SEVEN;
    return;
  } else {
    return;
  }

  // ON EACH CALL SET INTERVAL OF INTEGRATION AND COUNTER FOR NUMBER OF
  // STEPS.  ADJUST INPUT ERROR TOLERANCES TO DEFINE WEIGHT VECTOR FOR
  // ROUTINE STEP1
  del = tout - t;
  absdel = dabs(del);
  tend = t + 10.0 * del;
  if (isn < 0) {
    tend = tout;
  }
  nostep = 0;
  kle4 = 0;
  stiff = false;
  releps = relerr / eps;
  abseps = abserr / eps;
  if ((iflag == 1) || (isnold < 0) || !(delsgn * del > 0.0)) {
    // ON START AND RESTART ALSO SET WORK VARIABLES X AND YY(*), STORE THE
    // DIRECTION OF INTEGRATION AND INITIALIZE THE STEP SIZE
    start = true;
    x = t;
    for (int L = 1; L <= neqn; L++) {
      yy[L] = y[L];
    }
    delsgn = sign(1.0, del);
    h = sign(max(fouru * dabs(x), dabs(tout - x)), tout - x);
  }

  while (true) {
    // IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN
    if (!(dabs(x - t) < absdel)) {
      intrp(tout, y, kold);
      iflag = 2;
      t = tout;
      told = t;
      isnold = isn;
      return;
    }

    // IF CANNOT GO PAST OUTPUT POINT AND SUFFICIENTLY CLOSE,
    // EXTRAPOLATE AND RETURN
    if (!(isn > 0 || dabs(tout - x) >= fouru * dabs(x))) {
      h = tout - x;
      dy_dx->func(x, yy, yp);
      for (int L = 1; L <= neqn; L++) {
        y[L] = yy[L] + h * yp[L];
      }

      iflag = 2;
      t = tout;
      told = t;
      isnold = isn;
      return;
    }

    // TEST FOR TOO MANY STEPS
    if (!(nostep < maxnum)) {
      iflag = isn * 4;
      if (stiff) {
        iflag = isn * 5;
      }
      for (int L = 1; L <= neqn; L++) {
        y[L] = yy[L];
      }
      t = x;
      told = t;
      isnold = 1;
      return;
    }

    h = sign(min(dabs(h), dabs(tend - x)), h);
    for (int L = 1; L <= neqn; L++) {
      wt[L] = releps * dabs(yy[L]) + abseps;
      if (wt[L] <= 0.0) {
        // RELATIVE ERROR CRITERION INAPPROPRIATE
        iflag = isn * 6;
        for (int L = 1; L <= neqn; L++) {
          y[L] = yy[L];
        }
        t = x;
        told = t;
        isnold = 1;

        errcode = sansmic::ODE_BAD_ABSERR;
        return;
      }
    }

    step1(eps, start, k, kold, crash, phase1, ns, nornd);

    // TEST FOR TOLERANCES TOO SMALL
    if (crash) {
      iflag = isn * 3;
      relerr = eps * releps;
      abserr = eps * abseps;
      for (int L = 1; L <= neqn; L++) {
        y[L] = yy[L];
      }
      t = x;
      told = t;
      isnold = 1;
      return;
    }

    // AUGMENT COUNTER ON NUMBER OF STEPS AND TEST FOR STIFFNESS
    nostep++;
    kle4++;
    if (kold > 4) {
      kle4 = 0;
    }
    if (kle4 >= 50) {
      stiff = true;
    }
  }
  return;
}

/**
 * @brief Take a single solver step, used indirectly through ODE
 * @param eps local error tolerance
 * @param start logical variable set true for first step
 * @param k appropriate order for next step
 * @param kold order used for last successful step
 * @param crash logical variable set true when no step can be taken
 * @param phase1 elimainate local retention of variables
 * @param ns elimainate local retention of variables
 * @param nornd elimainate local retention of variables
 * @details @verbatim embed:rst

.. rubric:: Abstract

Subroutine :cpp:func:`~sansmic::step1` is normally used indirectly through
subroutine :cpp:func:`~sansmic::ode`. Because :cpp:func:`~sansmic::ode`
suffices for most problems and is much easier to use, using it should be
considered before using :cpp:func:`~sansmic::step1` alone.

Subroutine :cpp:func:`~sansmic::step1` integrates a system of :math:`n_{eqn}`
first order ordinary differential equations one step, normally from :math:`x`
to :math:`x+h`, using a modified divided difference form of the Adams PECE
formulas. Local extrapolation is used to improve absolute stability and
accuracy. The code adjusts its order and step size to control the local error
per unit step in a generalized sense.  special devices are included to control
roundoff error and to detect when the user is requesting too much accuracy.

This code is completely explained and documented in the text, *Computer
solution of ordinary differential equations, the initial value problem* by L.
F. Shampine and M. K. Gordon. Further details on use of this code are available
in *Solving ordinary differential equations with ODE, STEP, and INTRP*, by L.
F. Shampine and M. K. Gordon, SLA-73-1060.

The parameters represent:
- ``f`` -- subroutine to evaluate derivatives
- ``neqn`` -- number of equations to be integrated
- ``y[]`` -- solution vector at :math:`x`
- ``x`` -- independent variable
- ``h`` -- appropriate step size for next step.  normally determined by code
- ``eps`` -- local error tolerance
- ``wt[]`` -- vector of weights for error criterion
- ``start`` -- logical variable set true for first step, false otherwise
- ``hold`` -- step size used for last successful step
- ``k`` -- appropriate order for next step (determined by code)
- ``kold`` -- order used for last successful step
- ``crash`` -- logical variable set true when no step can be taken, else false
- ``yp[]`` -- derivative of solution vector at :math:`x` after successful step

The arrays ``phi``, ``psi`` are required for the interpolation subroutine
:cpp:func:`~sansmic::intrp`. The array ``p`` is internal to the code. The
remaining nine variables and arrays are included in the call list only to
eliminate local retention of variables between calls.

Input to :cpp:func:`~sansmic::step1`
    First call
        The user must provide storage in his calling program for all arrays
        in the call list, namely

        .. code:: fortran

            dimension
            y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12), 1
            alpha(12),beta(12),sig(13),v(12),w(12),g(13)


        The user must also declare `start`, `crash`, `phase1` and `nornd`
        logical variables and `f` an external subroutine, supply the subroutine
        ``f(x,y,yp)``  to evaluate

        .. math::

            \frac{d}{dx} y_i = \dot{y}_i = f(x, y_1, y_2, \dots, y_n)


        and initialize only the following parameters.
        (*Note: this is provided by the PlumeRise class for the specific
        SANSMIC case.*)

        - ``neqn`` -- number of equations to be integrated
        - ``y[]`` -- vector of initial values of dependent variables
        - ``x`` -- initial value of the independent variable
        - ``h`` -- nominal step size indicating direction of integration and
          maximum size of step.  must be variable
        - ``eps`` -- local error tolerance per step.  must be variable
        - ``wt[]`` -- vector of non-zero weights for error criterion
        - ``start`` -- ``true``

        step1  requires that the l2 norm of the vector with components
        local error(l)/wt(l)  be less than  eps  for a successful step.  the
        array  wt  allows the user to specify an error test appropriate
        for their problem.  For example,

        wt(l) = 1.0
            specifies absolute error,

        wt(l) = abs(y(l))
            error relative to the most recent value of the l-th component of
            the solution, wt(l) = abs(yp(l)) error relative to the most recent
            value of the l-th component of the derivative, wt(l) =
            amax1(wt(l),abs(y(l))) error relative to the largest magnitude of
            l-th component obtained so far, wt(l) abs(y(l))*relerr/eps +
            abserr/eps specifies a mixed relative-absolute test where  relerr
            is relative error, abserr  is absolute error and  eps =
            amax1(relerr,abserr) .


    Subsequent calls
        Subroutine :cpp:func:`~sansmic::step1` is designed so that all
        information needed to continue the integration, including the step size
        h  and the order k , is returned with each step.  with the exception of
        the step size, the error tolerance, and the weights, none of the
        parameters should be altered. the array  wt  must be updated after each
        step to maintain relative error tests like those above.  normally the
        integration is continued just beyond the desired endpoint and the
        solution interpolated there with subroutine intrp . if it is impossible
        to integrate beyond the endpoint, the step size may be reduced to hit
        the endpoint since the code will not take a step larger than the  h
        input. changing the direction of integration, i.e., the sign of  h ,
        requires the user set  start = .true. before calling  step1 again. this
        is the only situation in which  start should be altered.


Output from STEP1
    Successful step
        the subroutine returns after each successful step with  start  and
        crash  set .false. .  x  represents the independent variable
        advanced one step of length  hold  from its value on input and  y
        the solution vector at the new value of  x .  all other parameters
        represent information corresponding to the new  x  needed to
        continue the integration.

    Unsuccessful step
        when the error tolerance is too small for the machine precision,
        the subroutine returns without taking a step and  crash = .true. .
        an appropriate step size and error tolerance for continuing are
        estimated and all other information is restored as upon input
        before returning.  to continue with the larger tolerance, the user
        just calls the code again.  a restart is neither required nor
        desirable.


@endverbatim
 */
void sansmic::ODESolver::step1(double &eps, bool &start, int &k, int &kold,
                               bool &crash, bool &phase1, int &ns,
                               bool &nornd) {
  int ifail, kp1, kp2, km1, km2, nsp1, im1, nsm2, nsp2, limit1, limit2, knew,
      ip1;
  double u, twou, fouru, round, sum, absh, p5eps, realns, reali, temp1, temp2,
      temp3, temp4, temp5, temp6, tau, errk, errkm1, errkm2, err, xold, rho,
      errkp1, hnew, r;

  double two[] = {2.0,   4.0,   8.0,    16.0,   32.0,   64.0,  128.0,
                  256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};
  double gstr[] = {0.500,   0.0833,  0.0417,  0.0264,  0.0188,  0.0143, 0.0114,
                   0.00936, 0.00789, 0.00679, 0.00592, 0.00524, 0.00468};

  bool raise_order = false, lower_order = false, calc_coeff = true;

  u = 1.0e-30;
  twou = 2.0 * u;
  fouru = 4.0 * u;

  /* *******************************
   *    *** BEGIN BLOCK 0 ***
   * CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
   * PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
   * STARTING STEP SIZE.
   ******************************* */
  // IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
  crash = true;
  if (dabs(h) < fouru * dabs(x)) {
    h = sign(fouru * dabs(x), h);
    return;
  }

  p5eps = 0.5 * eps;

  // IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE
  round = 0.0;
  for (int L = 1; L <= neqn; L++) {  // DO 10
    round = round + (yy[L] / wt[L]) * (yy[L] / wt[L]);
  }  // 10

  round = twou * sqrt(round);
  if (p5eps < round) {
    eps = 2.0 * round * (1.0 + fouru);
    return;
  }

  crash = false;
  g[1] = 1.0;
  g[2] = 0.5;
  sig[1] = 1.0;
  if (start != 0) {
    // INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP
    dy_dx->func(x, yy, yp);
    sum = 0.0;
    for (int L = 1; L <= neqn; L++) {
      phi[L][1] = yp[L];
      phi[L][2] = 0.0;
      sum = sum + (yp[L] / wt[L]) * (yp[L] / wt[L]);
    }
    sum = sqrt(sum);
    absh = dabs(h);
    if (eps < 16.0 * sum * h * h) {
      absh = 0.25 * sqrt(eps / sum);
    }
    h = sign(max(absh, fouru * dabs(x)), h);
    hold = 0.0;
    k = 1;
    kold = 0;
    start = false;
    phase1 = true;
    nornd = true;
    if (!(p5eps > 100.0 * round)) {
      nornd = false;
      for (int L = 1; L <= neqn; L++) {
        phi[L][15] = 0.0;
      }
    }
  }

  ifail = 0;
  /* *******************************
   * END BLOCK 0
   ******************************* */

  /* *******************************
   *    *** BEGIN BLOCK 1 ***
   * COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
   * THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED.
   ******************************* */
  do {
    kp1 = k + 1;
    kp2 = k + 2;
    km1 = k - 1;
    km2 = k - 2;

    // NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
    // ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE
    if (h != hold) {
      ns = 0;
    }
    if (ns <= kold) {
      ns++;
    }
    nsp1 = ns + 1;
    if (k >= ns)

    {  // COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
      // ARE
      // CHANGED
      beta[ns] = 1.0;
      realns = (double)ns;
      alpha[ns] = 1.0 / realns;
      temp1 = h * realns;
      sig[nsp1] = 1.0;
      if (k >= nsp1) {
        for (int i = nsp1; i <= k; i++) {
          im1 = i - 1;
          temp2 = psi[im1];
          psi[im1] = temp1;
          beta[i] = beta[im1] * psi[im1] / temp2;
          temp1 = temp2 + h;
          alpha[i] = h / temp1;
          reali = (double)i;
          sig[i + 1] = reali * alpha[i] * sig[i];
        }
      }

      psi[k] = temp1;
      // COMPUTE COEFFICIENTS G(*)
      // INITIALIZE V(*) AND SET W(*).
      bool no_upd_v = false;
      if (ns > 1)
        no_upd_v = false;
      else {
        for (int iq = 1; iq <= k; iq++) {
          temp3 = (iq) * (iq + 1);
          v[iq] = 1.0 / temp3;
          w[iq] = v[iq];
        }
        no_upd_v = true;
      }
      if (!no_upd_v) {
        // IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*)
        if (!(k <= kold))

        {
          temp4 = k * kp1;
          v[k] = 1.0 / temp4;
          nsm2 = ns - 2;
          if (!(nsm2 < 1)) {
            for (int j = 1; j <= nsm2; j++) {
              int i = k - j;
              v[i] = v[i] - alpha[j + 1] * v[i + 1];
            }
          }
        }

        // UPDATE V(*) AND SET W(*)
        limit1 = kp1 - ns;
        temp5 = alpha[ns];
        for (int iq = 1; iq <= limit1; iq++) {
          v[iq] = v[iq] - temp5 * v[iq + 1];
          w[iq] = v[iq];
        }
        g[nsp1] = w[1];
      }

      // COMPUTE THE G(*) IN THE WORK VECTOR W(*)
      nsp2 = ns + 2;
      if (kp1 >= nsp2) {
        for (int i = nsp2; i <= kp1; i++) {
          limit2 = kp2 - i;
          temp6 = alpha[i - 1];
          for (int iq = 1; iq <= limit2; iq++) {
            w[iq] = w[iq] - temp6 * w[iq + 1];
          }
          g[i] = w[1];
        }
      }
    }

    u = u + 0.0;  // CONTINUE
    /* *******************************
     * *** END BLOCK 1
     ******************************* */

    /* *******************************
     *    *** START BLOCK 2 ***
     * PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED
     * SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K,
     * K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
     ******************************* */
    // CHANGE PHI TO PHI STAR
    if (k >= nsp1) {
      for (int i = nsp1; i <= k; i++) {
        temp1 = beta[i];
        for (int L = 1; L <= neqn; L++) {
          phi[L][i] = temp1 * phi[L][i];
        }
      }
    }

    // PREDICT SOLUTION AND DIFFERENCES
    for (int L = 1; L <= neqn; L++) {
      phi[L][kp2] = phi[L][kp1];
      phi[L][kp1] = 0.0;
      p[L] = 0.0;
    }

    for (int j = 1; j <= k; j++) {  // do 230
      int i = kp1 - j;
      ip1 = i + 1;
      temp2 = g[i];
      for (int L = 1; L <= neqn; L++) {
        p[L] = p[L] + temp2 * phi[L][i];
        phi[L][i] = phi[L][i] + phi[L][ip1];
      }
    }
    if (nornd) {
      for (int L = 1; L <= neqn; L++) {
        p[L] = yy[L] + h * p[L];
      }
    } else {
      for (int L = 1; L <= neqn; L++) {
        tau = h * p[L] - phi[L][15];
        p[L] = yy[L] + tau;
        phi[L][16] = (p[L] - yy[L]) - tau;
      }
    }

    xold = x;
    x = x + h;
    absh = dabs(h);
    dy_dx->func(x, p, yp);
    // estimate errors at orders k, k-1, k-2
    errkm2 = 0.0;
    errkm1 = 0.0;
    errk = 0.0;
    for (int L = 1; L <= neqn; L++) {
      temp3 = 1.0 / wt[L];
      temp4 = yp[L] - phi[L][1];

      if (km2 > 0) {
        errkm2 = errkm2 + ((phi[L][km1] + temp4) * temp3) *
                              ((phi[L][km1] + temp4) * temp3);
      }
      if (km2 >= 0) {
        errkm1 = errkm1 +
                 ((phi[L][k] + temp4) * temp3) * ((phi[L][k] + temp4) * temp3);
      }
      errk = errk + (temp4 * temp3) * (temp4 * temp3);
    }

    if (km2 > 0) {
      errkm2 = absh * sig[km1] * gstr[km2] * sqrt(errkm2);
    }
    if (km2 >= 0) {
      errkm1 = absh * sig[k] * gstr[km1] * sqrt(errkm1);
    }
    temp5 = absh * sqrt(errk);

    err = temp5 * (g[k] - g[kp1]);
    errk = temp5 * sig[kp1] * gstr[k];
    knew = k;

    // test if order should be lowered
    if (km2 > 0) {
      if (max(errkm1, errkm2) <= errk) {
        knew = km1;
      }
    }
    if (km2 == 0) {
      if (errkm1 <= 0.5 * errk) {
        knew = km1;
      }
    }

    if (err <= eps)
      calc_coeff = false;
    else {
      /* *******************************
       * END BLOCK 2
       ******************************* */

      /* *******************************
       * BEGIN BLOCK 3
       * THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) .
       * IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS
       *MORE THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
       * TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR
       *MACHINE PRECISION
       ******************************* */
      // RESTORE X, PHI(*,*) AND PSI(*)
      phase1 = false;
      x = xold;
      for (int i = 1; i <= k; i++) {
        temp1 = 1.0 / beta[i];
        ip1 = i + 1;
        for (int L = 1; L <= neqn; L++) {
          phi[L][i] = temp1 * (phi[L][i] - phi[L][ip1]);
        }
      }
      if (!(k < 2)) {
        for (int i = 2; i <= k; i++) {
          psi[i - 1] = psi[i] - h;
        }
      }
      // ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP
      // SIZE
      ifail++;
      temp2 = 0.5;

      if ((ifail - 3 > 0) && (p5eps < 0.25 * errk)) {
        temp2 = sqrt(p5eps / errk);
      }
      if (ifail - 3 >= 0) knew = 1;

      h = temp2 * h;
      k = knew;
      if (!(dabs(h) >= fouru * dabs(x))) {
        crash = true;
        h = sign(fouru * dabs(x), h);
        eps = eps + eps;
        return;
      }
    }
  }
  /* *******************************
   * END BLOCK 3
   ******************************* */
  while (calc_coeff);

  /* *******************************
   *    *** BEGIN BLOCK 4 ***
   * THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE
   * THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE
   * DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP.
   ******************************* */
  kold = k;
  hold = h;
  // CORRECT AND EVALUATE
  temp1 = h * g[kp1];
  if (nornd != 0) {
    for (int L = 1; L <= neqn; L++) {
      yy[L] = p[L] + temp1 * (yp[L] - phi[L][1]);
    }
  } else {
    for (int L = 1; L <= neqn; L++) {
      rho = temp1 * (yp[L] - phi[L][1]) - phi[L][16];
      yy[L] = p[L] + rho;
      phi[L][15] = (yy[L] - p[L]) - rho;
    }
  }

  dy_dx->func(x, yy, yp);

  // UPDATE DIFFERENCES FOR NEXT STEP
  for (int L = 1; L <= neqn; L++) {
    phi[L][kp1] = yp[L] - phi[L][1];
    phi[L][kp2] = phi[L][kp1] - phi[L][kp2];
  }
  for (int i = 1; i <= k; i++) {
    for (int L = 1; L <= neqn; L++) {
      phi[L][i] = phi[L][i] + phi[L][kp1];
    }
  }

  // ESTIMATE ERROR AT ORDER K+1 UNLESS:
  // IN FIRST PHASE WHEN ALWAYS RAISE ORDER, ALREADY DECIDED TO LOWER ORDER,
  // STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE
  errkp1 = 0.0;
  if (knew == km1 || k == 12) {
    phase1 = 0;
  }

  if (phase1) {
    raise_order = true;
    lower_order = false;
  } else if (knew == km1) {
    raise_order = false;
    lower_order = true;
  } else if (kp1 > ns) {
    raise_order = false;
    lower_order = false;
  } else {
    for (int L = 1; L <= neqn; L++) {
      errkp1 = errkp1 + (phi[L][kp2] / wt[L]) * (phi[L][kp2] / wt[L]);
    }
    errkp1 = absh * gstr[kp1] * sqrt(errkp1);

    // USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER FOR
    // NEXT STEP
    if (k <= 1) {
      if (errkp1 >= 0.5 * errk) {
        raise_order = false;
        lower_order = false;
      } else {
        raise_order = true;
        lower_order = false;
      }
    } else {
      // k > 1 so else do next check
      if (errkm1 <= min(errk, errkp1)) {
        raise_order = false;
        lower_order = true;
      } else if (errkp1 >= errk || k == 12) {
        raise_order = false;
        lower_order = false;
      } else {
        // HERE ERKP1 .LT. ERK .LT. AMAX1(ERKM1,ERKM2) ELSE ORDER WOULD
        // HAVE BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
        raise_order = true;
        lower_order = false;
      }
    }
  }
  if (raise_order) {
    // RAISE ORDER
    k = kp1;
    errk = errkp1;
  } else if (lower_order) {
    // LOWER ORDER
    k = km1;
    errk = errkm1;
  }

  // WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
  hnew = h + h;

  if (!phase1 && (p5eps < errk * two[k + 1])) {
    hnew = h;

    if (p5eps < errk) {
      temp2 = k + 1;
      r = pow(p5eps / errk, 1.0 / temp2);
      hnew = absh * max(0.5, min(0.9, r));
      hnew = sign(max(hnew, fouru * dabs(x)), h);
    }
  }

  h = hnew;
  return;
  /* *******************************
   * END BLOCK 4
   ******************************* */
}

/**
 * @brief Approximate the solution near x by polymonial.
 * @param xout points at which solution is known
 * @param yout solution at xout
 * @param kold pass from sansmic::step1 - unchanged
 * @details @verbatim embed:rst

.. rubric:: Abstract

The methods in subroutine :cpp:func:`~sansmic::step1` approximate the solution
near :math:`x` by a polynomial. Subroutine :cpp:func:`~sansmic::intrp`
approximates the solution at :math:`x_{out}` by evaluating the polynomial
there. Information defining this polynomial is passed from
:cpp:func:`~sansmic::step1` so :cpp:func:`~sansmic::intrp` cannot be used
alone.

This code is completely explained and documented in the text, *Computer
solution of ordinary differential equations, the initial value problem* by L.
F. Shampine and M. K. Gordon. Further details on use of this code are available
in *Solving ordinary differential equations with ODE, STEP, and INTRP*, by L.
F. Shampine and M. K. Gordon, SLA-73-1060.

Input to :cpp:func:`~sansmic::intrp`
    The user provides storage in the calling program for the arrays in the call
    list
    - ``dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)``

    and defines
    - ``xout`` -- point at which solution is desired.

    The remaining parameters are defined by a previous call to
    :cpp:func:`~sansmic::step1` and should be passed to
    :cpp:func:`~sansmic::intrp` unchanged.


Output from :cpp:func:`~sansmic::intrp`
    - ``yout`` -- solution at :math:`x_{out}`
    - ``ypout`` -- derivative of solution at :math:`x_{out}`

    the remaining parameters are returned unaltered from their input values.
    Integration with :cpp:func:`~sansmic::step1` may be continued.

@endverbatim
 */
void sansmic::ODESolver::intrp(double xout, vector<double> &yout, int kold) {
  int ki, kip1, jm1;
  double hi, temp1, term, psijm1, gamma, eta, limit1, erm, temp2, temp3;
  vector<double> g = vector<double>(14, 0.0);
  vector<double> w = vector<double>(14, 0.0);
  vector<double> rho = vector<double>(14, 0.0);
  g[1] = 1.0;
  rho[1] = 1.0;

  hi = xout - x;
  ki = kold + 1;
  kip1 = ki + 1;

  // INITIALIZE W[*] FOR COMPUTING G[*]
  for (int i = 1; i <= ki; i++) {
    temp1 = (double)i;
    w[i] = 1.0 / temp1;
  }
  term = 0.0;

  // COMPUTE G[*]
  for (int j = 2; j <= ki; j++) {
    jm1 = j - 1;
    psijm1 = psi[jm1];
    gamma = (hi + term) / psijm1;
    eta = hi / psijm1;
    limit1 = kip1 - j;
    for (int i = 1; i <= limit1; i++) {
      w[i] = gamma * w[i] + eta * w[i + 1];
    }
    g[j] = w[1];
    rho[j] = gamma * rho[jm1];
    erm = psijm1;
  }

  // INTERPOLATE
  for (int L = 1; L <= neqn; L++) {
    ypout[L] = 0.0;
    yout[L] = 0.0;
  }

  for (int j = 1; j <= ki; j++) {
    int i = kip1 - j;
    temp2 = g[i];
    temp3 = rho[i];
    for (int L = 1; L <= neqn; L++) {
      yout[L] = yout[L] + temp2 * phi[L][i];
      ypout[L] = ypout[L] + temp3 * phi[L][i];
    }
  }

  for (int L = 1; L <= neqn; L++) {
    yout[L] = yy[L] + hi * yout[L];
  }

  return;
}
