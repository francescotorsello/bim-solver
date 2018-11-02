/**
 *  @file      lowess.h
 *  @brief     Locally weighted scatterplot smoothing (LOWESS).
 *  @authors   W.S. Cleveland, Hannes Roest
 *  @copyright Copyright (c) 2015. All rights reserved.
 *             This software is released under a three-clause BSD license.
 */

#ifndef _LOWESS_H_INCLUDED
#define _LOWESS_H_INCLUDED

#include <cstdlib>
#include <cmath>
#include <algorithm>    // std::min, std::max

/** Templated class for LOWESS (locally weighted scatterplot smoothing).
 *
 *  The lowess code is translated from `RATFOR` lowess code written by W. S.
 *  Cleveland as obtained from `NETLIB`.
 *  It is based on two functions written in ratfor: `lowest`
 *  and `lowess`. The code has since been refactored and commented further.
 */
template <typename Container, typename T>
class Lowess
{
    inline T pow2(T x) { return x * x;  }
    inline T pow3(T x) { return x * x * x;  }

    /// Return the median of a sequence of numbers defined by the random
    /// access iterators begin and end.  The sequence must not be empty
    /// (median is undefined for an empty set).
    ///
    /// The numbers must be convertible to double.
    ///
    template <class RandAccessIter>
    T median(RandAccessIter begin, RandAccessIter end)
    {
      std::size_t size = end - begin;
      std::size_t middleIdx = size / 2;
      RandAccessIter target = begin + middleIdx;
      std::nth_element(begin, target, end);

      if (size % 2 != 0)
      {
        //Odd number of elements
        return *target;
      }
      else
      {
        //Even number of elements
        double a = *target;
        RandAccessIter targetNeighbor = target - 1;
        targetNeighbor = std::max_element(begin, target);
        return (a + *targetNeighbor) / 2.0;
      }
    }

    /// Calculate weights for weighted regression.
    bool calculate_weights(const Container& x,
                           const size_t n,
                           const T current_x,
                           const bool use_resid_weights,
                           const size_t nleft,
                           const Container& resid_weights,
                           Container& weights,
                           size_t& nrt,
                           const T h)
    {
      T r;
      size_t j;

      T h9 = .999 * h;
      T h1 = .001 * h;
      T a = 0.0; // sum of weights

      // compute weights (pick up all ties on right)
      for (j = nleft; j < n; j++)
      {
        // Compute the distance measure, then apply the tricube
        // function on the distance to get the weight.
        // use_resid_weights will be False on the first iteration, then True
        // on the subsequent ones, after some residuals have been calculated.
        weights[j] = 0.0;
        r = std::abs(x[j] - current_x);
        if (r <= h9)
        {
          if (r > h1)
          {
            // small enough for non-zero weight
            // compute tricube function: ( 1 - (r/h)^3 )^3
            weights[j] = pow3(1.0 - pow3(r / h));
          }
          else
          {
            weights[j] = 1.0;
          }

          if (use_resid_weights)
          {
            weights[j] = resid_weights[j] * weights[j];
          }

          a += weights[j];
        }
        else if (x[j] > current_x)
        {
          // get out at first zero wt on right
          break;
        }
      }

      // rightmost pt (may be greater than nright because of ties)
      nrt = j - 1;
      if (a <= 0.0)
      {
        return false;
      }
      else
      {
        // normalize weights (make sum of w[j] == 1)
        for (j = nleft; j <= nrt; j++)
        {
          weights[j] = weights[j] / a;
        }

        return true;
      }
    }

    /// Calculate smoothed/fitted y-value by weighted regression.
    void calculate_y_fit(const Container& x,
                         const Container& y,
                         const T current_x,
                         const size_t n,
                         const size_t nleft,
                         const size_t nrt,
                         const T h,
                         T& ys,
                         Container& weights)
    {
      T range = x[n - 1] - x[0];

      if (h > 0.0)
      {
        // use linear fit

        // No regression function (e.g. lstsq) is called. Instead a "projection
        // vector" p_i_j is calculated, and y_fit[i] = sum(p_i_j * y[j]) = y_fit[i]
        // for j s.t. x[j] is in the neighborhood of x[i]. p_i_j is a function of
        // the weights, x[i], and its neighbors.
        // To save space, p_i_j is computed in place using the weight vector.

        // find weighted center of x values
        T sum_weighted_x = 0.0; // originally variable a
        for (size_t j = nleft; j <= nrt; j++)
        {
          sum_weighted_x += weights[j] * x[j];
        }

        T b = current_x - sum_weighted_x; // originally variable b
        T weighted_sqdev = 0.0; // originally variable c
        for (size_t j = nleft; j <= nrt; j++)
        {
          weighted_sqdev += weights[j] *
                            (x[j] - sum_weighted_x) * (x[j] - sum_weighted_x);
        }

        if (sqrt(weighted_sqdev) > .001 * range)
        {
          // points are spread out enough to compute slope
          b = b / weighted_sqdev;
          for (size_t j = nleft; j <= nrt; j++)
          {
            // Compute p_i_j in place
            weights[j] = weights[j] * (1.0 + b * (x[j] - sum_weighted_x));
          }
        }
      }

      ys = 0.0;
      for (size_t j = nleft; j <= nrt; j++)
      {
        ys += weights[j] * y[j];
      }
    }

    bool lowest(const Container& x,
                const Container& y,
                size_t n,
                T current_x, //xs
                T& ys,
                size_t nleft,
                size_t nright,
                Container& weights, // vector w
                bool use_resid_weights,  // userw
                const Container& resid_weights)
    {
      T h;
      size_t nrt; // rightmost pt (may be greater than nright because of ties)

      h = std::max(current_x - x[nleft], x[nright] - current_x);

      // Calculate the weights for the regression in this neighborhood.
      // Determine if at least some weights are positive, so a regression
      // is ok.
      bool fit_ok = calculate_weights(x, n, current_x, use_resid_weights,
                                      nleft, resid_weights,
                                      weights, nrt, h);
      if (!fit_ok)
      {
        return fit_ok;
      }

      // If it is ok to fit, run the weighted least squares regression
      calculate_y_fit(x, y, current_x, n, nleft, nrt, h, ys, weights);

      return fit_ok;
    }

    /// Find the indices bounding the k-nearest-neighbors of the current point.
    void update_neighborhood(const Container& x,
                             const size_t n,
                             const size_t i,
                             size_t& nleft,
                             size_t& nright)
    {
      T d1, d2;
      // A subtle loop. Start from the current neighborhood range:
      // [nleft, nright). Shift both ends rightwards by one
      // (so that the neighborhood still contains ns points), until
      // the current point is in the center (or just to the left of
      // the center) of the neighborhood. This neighborhood will
      // contain the ns-nearest neighbors of x[i].
      //
      // Once the right end hits the end of the data, hold the
      // neighborhood the same for the remaining x[i]s.
      while (nright < n - 1)
      {
        // move nleft, nright to right if radius decreases
        d1 = x[i] - x[nleft];
        d2 = x[nright + 1] - x[i];
        // if d1 <= d2 with x[nright+1] == x[nright], lowest fixes
        if (d1 <= d2) break;
        // radius will not decrease by move right
        nleft++;
        nright++;
      }
    }

    /// Update the counters of the local regression.
    void update_indices(const Container& x,
                        const size_t n,
                        const T delta,
                        size_t& i,
                        size_t& last,
                        Container& ys)
    {
      // For most points within delta of the current point, we skip the
      // weighted linear regression (which save much computation of
      // weights and fitted points). Instead, we'll jump to the last
      // point within delta, fit the weighted regression at that point,
      // and linearly interpolate in between.

      // the last point actually estimated
      last = i;

      // This loop increments until we fall just outside of delta distance,
      // copying the results for any repeated x's along the way.
      T cut = x[last] + delta;
      for (i = last + 1; i < n; i++)
      {
        // find close points
        if (x[i] > cut) break;

        // i one beyond last pt within cut
        if (x[i] == x[last])
        {
          // exact match in x
          // if tied with previous x-value, just use the already
          // fitted y, and update the last-fit counter.
          ys[i] = ys[last];
          last = i;
        }
      }


      // the next point to fit the regression at is either one prior to i (since
      // i should be the first point outside of delta) or it is "last + 1" in the
      // case that i never got incremented. This insures we always step forward.
      // -> back 1 point so interpolation within delta, but always go forward
      i = std::max(last + 1, i - 1);
    }

    /// Calculate smoothed/fitted y by linear interpolation between the current
    /// and previous y fitted by weighted regression.
    void interpolate_skipped_fits(const Container& x,
                                  const size_t i,
                                  const size_t last,
                                  Container& ys)
    {
      // skipped points -- interpolate
      T alpha;
      T denom = x[i] - x[last]; // non-zero - proof?
      for (size_t j = last + 1; j < i; j = j + 1)
      {
        alpha = (x[j] - x[last]) / denom;
        ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last];
      }
    }

    /// Calculate residual weights for the next `robustifying` iteration.
    void calculate_residual_weights(const size_t n,
                                    const Container& weights,
                                    Container& resid_weights)
    {
      T r;

      for (size_t i = 0; i < n; i++)
      {
        resid_weights[i] = std::abs(weights[i]);
      }

      // ***********************************
      // Compute pseudo-median (take average even if we have an odd number of
      // elements), following the original implementation. We could also use a
      // true median calculation here:
      // T cmad = 6.0 * median(resid_weights.begin(), resid_weights.end());
      // ***********************************

      size_t m1 = n / 2; // FORTRAN starts with one, CPP with zero
      // size_t m1 = 1 + n / 2; // original FORTRAN code
      // size_t m2 = n - m1 + 1; // see below, we don't explicitly sort but use max_element

      // Use nth element to find element m1, which produces a partially sorted
      // vector. This means we can get element m2 by looking for the maximum in the
      // remainder.
      typename Container::iterator it_m1 = resid_weights.begin() + m1;
      std::nth_element(resid_weights.begin(), it_m1, resid_weights.end());
      typename Container::iterator it_m2 = std::max_element(
        resid_weights.begin(), it_m1);
      T cmad = 3.0 * (*it_m1 + *it_m2);
      T c9 = .999 * cmad;
      T c1 = .001 * cmad;

      for (size_t i = 0; i < n; i++)
      {
        r = std::abs(weights[i]);
        if (r <= c1)
        {
          // near 0, avoid underflow
          resid_weights[i] = 1.0;
        }
        else if (r > c9)
        {
          // near 1, avoid underflow
          resid_weights[i] = 0.0;
        }
        else
        {
          resid_weights[i] = pow2(1.0 - pow2(r / cmad));
        }
      }
    }

public:
    int lowess(const Container& x,
               const Container& y,
               double frac,    // parameter f
               int nsteps,
               T delta,
               Container& ys,
               Container& resid_weights,   // vector rw
               Container& weights   // vector res
               )
    {
      bool fit_ok;
      size_t i, last, nleft, nright, ns;

      size_t n = x.size();
      if (n < 2)
      {
        ys[0] = y[0];
        return 1;
      }

      // how many points around estimation point should be used for regression:
      // at least two, at most n points
      size_t tmp = (size_t)(frac * (double)n);
      ns = std::max(std::min(tmp, n), (size_t)2);

      // robustness iterations
      for (int iter = 1; iter <= nsteps + 1; iter++)
      {
        // start of array in C++ at 0 / in FORTRAN at 1
        nleft = 0;
        nright = ns - 1;
        last = -1;          // index of prev estimated point
        i = 0;              // index of current point

        // Fit all data points y[i] until the end of the array
        do
        {
          // Identify the neighborhood around the current x[i]
          // -> get the nearest ns points
          update_neighborhood(x, n, i, nleft, nright);

          // Calculate weights and apply fit (original lowest function)
          fit_ok = lowest(x, y, n, x[i], ys[i], nleft, nright,
                          weights, (iter > 1), resid_weights);

          // if something went wrong during the fit, use y[i] as the
          // fitted value at x[i]
          if (!fit_ok) ys[i] = y[i];

          // If we skipped some points (because of how delta was set), go back
          // and fit them by linear interpolation.
          if (last < i - 1)
          {
            interpolate_skipped_fits(x, i, last, ys);
          }

          // Update the last fit counter to indicate we've now fit this point.
          // Find the next i for which we'll run a regression.
          update_indices(x, n, delta, i, last, ys);
        }
        while (last < n - 1);

        // compute current residuals
        for (i = 0; i < n; i++)
        {
          weights[i] = y[i] - ys[i];
        }

        // compute robustness weights except last time
        if (iter > nsteps) break;

        calculate_residual_weights(n, weights, resid_weights);
      }
      return 0;
    }
};

#endif // _LOWESS_H_INCLUDED
