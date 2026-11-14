// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_COOLING_GAUSSIANS_HPP
#define VOLUME_COOLING_GAUSSIANS_HPP

//#define VOLESTI_DEBUG

#include <iterator>
#include <vector>
#include <list>
#include <math.h>
#include <chrono>

#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/gaussian_helpers.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/math_helpers.hpp"


/////////////////// Helpers for random walks

template <typename WalkType>
struct update_delta
{
    template <typename NT>
    static void apply(WalkType, NT) {}
};

template <typename Polytope, typename RandomNumberGenerator>
struct update_delta<GaussianBallWalk::Walk<Polytope, RandomNumberGenerator>>
{
    template <typename NT>
    static void apply(GaussianBallWalk::Walk<Polytope, RandomNumberGenerator> walk,
                      NT delta)
    {
        walk.update_delta(delta);
    }
};


////////////////////////////// Algorithms

// Gaussian Anealling

// Implementation is based on algorithm from paper "A practical volume algorithm",
// Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society 2015
// Ben Cousins, Santosh Vempala

// Compute the first variance a_0 for the starting gaussian
template <typename Polytope, typename NT>
void get_first_gaussian(Polytope& P,
                        NT const& frac,
                        NT const& chebychev_radius,
                        NT const& error,
                        std::vector<NT> & a_vals)
{
    // if tol is smaller than 1e-6 no convergence can be obtained when float is used
    NT tol = std::is_same<float, NT>::value ? 0.001 : 0.0000001;

    std::vector <NT> dists = P.get_dists(chebychev_radius);
    NT lower = 0.0;
    NT upper = 1.0;

    // Compute an upper bound for a_0
    unsigned int i;
    const unsigned int maxiter = 10000;
    for (i= 1; i <= maxiter; ++i) {
        NT sum = 0.0;
        for (auto it = dists.begin(); it != dists.end(); ++it)
        {
            sum += std::exp(-upper * std::pow(*it, 2.0))
                    / (2.0 * (*it) * std::sqrt(M_PI * upper));
        }

        if (sum > frac * error)
        {
            upper = upper * 10;
        } else {
            break;
        }
    }

    if (i == maxiter) {
#ifdef VOLESTI_DEBUG
        std::cout << "Cannot obtain sharp enough starting Gaussian" << std::endl;
#endif
        return;
    }

    //get a_0 with binary search.  The tolerance is relative to `upper`, as in
    //Cousins-Vempala's reference code (ConvexBody.m, getStartingGaussian):
    //an absolute 1e-7 never terminates once a_0 exceeds ~1e9 (the spacing of
    //doubles near a_0 is then larger than 1e-7), i.e. for tiny bodies.
    while (upper - lower > tol * upper)
    {
        NT mid = (upper + lower) / 2.0;
        NT sum = 0.0;
        for (auto it = dists.begin(); it != dists.end(); ++it) {
            sum += std::exp(-mid * std::pow(*it, 2.0))
                    / (2.0 * (*it) * std::sqrt(M_PI * mid));
        }

        if (sum < frac * error) {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    a_vals.push_back((upper + lower) / NT(2.0));
}


// Compute a_{i+1} when a_i is given
template
<
    typename RandomPointGenerator,
    typename Polytope,
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
NT get_next_gaussian(Polytope& P,
                     Point &p,
                     NT const& a,
                     const unsigned int &N,
                     const NT &ratio,
                     const NT &C,
                     const unsigned int& walk_length,
                     RandomNumberGenerator& rng)
{

    NT last_a = a;
    NT last_ratio = 0.1;
    // Worst-case cap (valid for every convex body): a^n Z(a) is log-concave
    // (Kalai-Vempala), so for a' = a(1-alpha)
    //   E[Y^2]/E[Y]^2 = Z(2a'-a)Z(a)/Z(a')^2 <= ((1-alpha)^2/(1-2alpha))^n.
    // Keep this <= 1+C; the sample variance below cannot see far tails.
    const NT q = std::pow(NT(1) + C, NT(1) / NT(P.dimension()));
    const NT alpha_max = std::sqrt((q - NT(1)) * q) - (q - NT(1));
    const NT min_shrink = NT(1) - alpha_max;
    //k is needed for the computation of the next variance a_{i+1} = a_i * (1-1/d)^k
    NT k = 1.0;
    const NT tol = 0.00001;
    bool done=false;
    std::vector<NT> fn(N,NT(0.0));
    std::list<Point> randPoints;
    typedef typename std::vector<NT>::iterator viterator;

    //sample N points
    PushBackWalkPolicy push_back_policy;
    RandomPointGenerator::apply(P, p, last_a, N, walk_length, randPoints,
                                push_back_policy, rng);

    while (!done)
    {
        NT new_a = last_a * std::max(std::pow(ratio,k), min_shrink);

        auto fnit = fn.begin();
        for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit, fnit++)
        {
            *fnit = eval_exp(*pit, new_a)/eval_exp(*pit, last_a);
        }
        std::pair<NT, NT> mv = get_mean_variance(fn);

        // Compute a_{i+1}
        if (mv.second/(mv.first * mv.first)>=C || mv.first/last_ratio<1.0+tol
            || std::pow(ratio,k) <= min_shrink)
        {
            if (k != 1.0)
            {
                k = k / 2;
            }
            done = true;
        } else {
            k = 2 * k;
        }
        last_ratio = mv.first;
    }
    return last_a * std::max(std::pow(ratio, k), min_shrink);
}

// Compute the sequence of spherical gaussians
template
<
    typename WalkType,
    typename RandomPointGenerator,
    typename Polytope,
    typename NT,
    typename RandomNumberGenerator
>
void compute_annealing_schedule(Polytope& P,
                                NT const& ratio,
                                NT const& C,
                                NT const& frac,
                                unsigned int const& N,
                                unsigned int const& walk_length,
                                NT const& chebychev_radius,
                                NT const& error,
                                std::vector<NT>& a_vals,
                                RandomNumberGenerator& rng)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::VT VT;

    // Compute the first gaussian
    get_first_gaussian(P, frac, chebychev_radius, error, a_vals);

#ifdef VOLESTI_DEBUG
    std::cout<<"first gaussian computed\n"<<std::endl;
#endif

    NT a_stop = 0.0;
    const NT tol = 0.001;
    unsigned int it = 0;
    unsigned int n = P.dimension();
    const unsigned int totalSteps = ((int)150/((1.0 - frac) * error))+1;

    if (a_vals[0]<a_stop) a_vals[0] = a_stop;

#ifdef VOLESTI_DEBUG
    std::cout<<"Computing the sequence of gaussians..\n"<<std::endl;
#endif

    Point p(n);

    while (true)
    {
        // Compute the next gaussian
        NT next_a = get_next_gaussian<RandomPointGenerator>
                      (P, p, a_vals[it], N, ratio, C, walk_length, rng);

        NT curr_fn = 0;
        NT curr_its = 0;
        auto steps = totalSteps;

        WalkType walk(P, p, a_vals[it], rng);
        //TODO: test update delta here?

        update_delta<WalkType>
                ::apply(walk, 4.0 * chebychev_radius
                        / std::sqrt(std::max(NT(1.0), a_vals[it]) * NT(n)));

        // Compute some ratios to decide if this is the last gaussian
        for (unsigned  int j = 0; j < steps; j++)
        {
            walk.apply(P, p, a_vals[it], walk_length, rng);
            curr_its += 1.0;
            curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
        }

        // Remove the last gaussian.
        // Set the last a_i equal to zero
        if (next_a>0 && curr_fn/curr_its>(1.0+tol))
        {
            a_vals.push_back(next_a);
            it++;
        } else if (next_a <= 0)
        {
            a_vals.push_back(a_stop);
            it++;
            break;
        } else {
            a_vals[it] = a_stop;
            break;
        }
    }
}

// Relative standard error of the mean of a (possibly autocorrelated) stream,
// by non-overlapping batch means with ~sqrt(N) batches.
template <typename NT>
NT batch_means_relative_se(std::vector<NT> const& ys, NT const& mean)
{
    const std::size_t N = ys.size();
    const std::size_t nb = std::max<std::size_t>(2, (std::size_t) std::sqrt((double) N));
    const std::size_t bs = N / nb;
    if (bs == 0 || mean <= NT(0)) return std::numeric_limits<NT>::max();
    NT bm_mean = 0;
    std::vector<NT> bm(nb, NT(0));
    for (std::size_t b = 0; b < nb; ++b) {
        for (std::size_t i = 0; i < bs; ++i) bm[b] += ys[b * bs + i];
        bm[b] /= NT(bs);
        bm_mean += bm[b] / NT(nb);
    }
    NT var_bm = 0;
    for (std::size_t b = 0; b < nb; ++b) var_bm += (bm[b] - bm_mean) * (bm[b] - bm_mean);
    var_bm /= NT(nb - 1);
    return std::sqrt(var_bm / NT(nb)) / mean;
}

template <typename NT>
struct gaussian_annealing_parameters
{
    gaussian_annealing_parameters(unsigned int d)
        :   frac(0.1)
        ,   ratio(NT(1)-NT(1)/(NT(d)))
        ,   C(NT(2))
        ,   N(500 * ((int) C) + ((int) (d * d / 2)))
        ,   W(6*d*d+800)
    {}

    NT frac;
    NT ratio;
    NT C;
    unsigned int N;
    unsigned int W;
};

template
<
    typename WalkTypePolicy,
    typename Polytope,
    typename RandomNumberGenerator

>
double volume_cooling_gaussians(Polytope& Pin,
                                RandomNumberGenerator& rng,
                                double const& error = 0.1,
                                unsigned int const& walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::VT 	VT;
    typedef typename WalkTypePolicy::template Walk
                                              <
                                                    Polytope,
                                                    RandomNumberGenerator
                                              > WalkType;
    typedef GaussianRandomPointGenerator<WalkType> RandomPointGenerator;

    //const NT maxNT = std::numeric_limits<NT>::max();//1.79769e+308;
    //const NT minNT = std::numeric_limits<NT>::min();//-1.79769e+308;

    auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    unsigned int m = P.num_of_hyperplanes();
    gaussian_annealing_parameters<NT> parameters(P.dimension());
    //RandomNumberGenerator rng(n);

    // Consider Chebychev center as an internal point
    auto InnerBall = P.ComputeInnerBall();
    if (InnerBall.second < 0.0) return -1.0;

    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(c.getCoefficients());

    // The schedule and walk use absolute constants (a_0 search starts at
    // upper = 1 with absolute tol 1e-7; walk step 4r/sqrt(max(1,a) n)), so
    // they are only scale-correct for an O(1) inner radius.  Rescale the
    // (row-normalized, centred) body to unit inner radius, and undo the
    // scaling on the result: vol(P) = r^n vol(P / r).
    const NT scale = radius;
    P.normalize();
    P.set_vec(P.get_vec() / scale);
    radius = NT(1);
    P.set_InnerBall(std::pair<Point, NT>(Point(n), radius));

    // Computing the sequence of gaussians
#ifdef VOLESTI_DEBUG
    std::cout<<"\n\nComputing annealing...\n"<<std::endl;
    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
#endif

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = parameters.ratio;
    NT C = parameters.C;
    unsigned int N = parameters.N;

    compute_annealing_schedule
    <
        WalkType,
        RandomPointGenerator
    >(P, ratio, C, parameters.frac, N, walk_length, radius, error, a_vals, rng);

#ifdef VOLESTI_DEBUG
    std::cout<<"All the variances of schedule_annealing computed in = "
            << (double)clock()/(double)CLOCKS_PER_SEC-tstart2<<" sec"<<std::endl;
    auto j=0;
    for (auto avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++){
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;
#endif

    // Initialization for the approximation of the ratios
    unsigned int W = parameters.W;
    unsigned int mm = a_vals.size()-1;
    std::vector<NT> last_W2(W,0);
    std::vector<NT> fn(mm,0);
    std::vector<NT> its(mm,0);
    VT lamdas;
    lamdas.setZero(m);
    NT vol = std::pow(M_PI/a_vals[0], (NT(n))/2.0);
    Point p(n); // The origin is the Chebychev center of the Polytope
    unsigned int i=0;

    typedef typename std::vector<NT>::iterator viterator;
    viterator itsIt = its.begin();
    viterator avalsIt = a_vals.begin();
    viterator minmaxIt;

#ifdef VOLESTI_DEBUG
    std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    std::cout<<"computing ratios..\n"<<std::endl;
#endif

    //iterate over the number of ratios
    for (viterator fnIt = fn.begin();
         fnIt != fn.end();
         fnIt++, itsIt++, avalsIt++, i++)
    {
        //initialize convergence test
        bool done = false;
        NT curr_eps = error/std::sqrt((NT(mm)));
        NT min_val = std::numeric_limits<NT>::min();
        NT max_val = std::numeric_limits<NT>::max();
        unsigned int min_index = W-1;
        unsigned int max_index = W-1;
        unsigned int index = 0;
        unsigned int min_steps = 0;
        std::vector<NT> last_W = last_W2;
        std::vector<NT> ys;          // the ratio stream of this phase
        NT next_check = 0;

        // Set the radius for the ball walk
        WalkType walk(P, p, *avalsIt, rng);

        update_delta<WalkType>
                ::apply(walk, 4.0 * radius
                         / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n)));

        while (!done || (*itsIt)<min_steps)
        {
            walk.apply(P, p, *avalsIt, walk_length, rng);

            *itsIt = *itsIt + 1.0;
            NT y = eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
            ys.push_back(y);
            *fnIt = *fnIt + y;
            NT val = (*fnIt) / (*itsIt);

            last_W[index] = val;
            if (val <= min_val)
            {
                min_val = val;
                min_index = index;
            } else if (min_index == index)
            {
                minmaxIt = std::min_element(last_W.begin(), last_W.end());
                min_val = *minmaxIt;
                min_index = std::distance(last_W.begin(), minmaxIt);
            }

            if (val >= max_val)
            {
                max_val = val;
                max_index = index;
            } else if (max_index == index)
            {
                minmaxIt = std::max_element(last_W.begin(), last_W.end());
                max_val = *minmaxIt;
                max_index = std::distance(last_W.begin(), minmaxIt);
            }

            // The window test alone only bounds the drift of the running mean
            // over the last W steps; its delivered error scales like
            // sqrt(curr_eps), not curr_eps.  Additionally require the
            // batch-means relative standard error of this ratio to be
            // <= curr_eps/2, so the product of mm ratios has sd <= error/2.
            if ( (max_val-min_val)/max_val <= curr_eps/2.0 && *itsIt >= next_check )
            {
                if (batch_means_relative_se(ys, (*fnIt) / (*itsIt)) <= curr_eps/2.0)
                    done=true;
                else
                    next_check = NT(1.2) * (*itsIt);
            }

            index = index%W + 1;
            if (index == W) index = 0;
        }
#ifdef VOLESTI_DEBUG
        std::cout << "ratio " << i << " = " << (*fnIt) / (*itsIt)
                  << " N_" << i << " = " << *itsIt << std::endl;
#endif
        vol *= ((*fnIt) / (*itsIt));
    }

#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        auto steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
#endif

    return vol * std::pow(scale, NT(n));
}


template
<
    typename WalkTypePolicy = GaussianCDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b, double>,
    typename Polytope
>
double volume_cooling_gaussians(Polytope &Pin,
                                 double const& error = 0.1,
                                 unsigned int const& walk_length = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_cooling_gaussians<WalkTypePolicy>(Pin, rng, error, walk_length);
}


template
<
    typename WalkTypePolicy = GaussianCDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b, double>,
    typename Polytope
>
double volume_cooling_gaussians(Polytope &Pin,
                                Cartesian<double>::Point const& interior_point,
                                unsigned int const& walk_length = 1,
                                double const& error = 0.1)
{
    RandomNumberGenerator rng(Pin.dimension());
    Pin.set_interior_point(interior_point);

    return volume_cooling_gaussians<WalkTypePolicy>(Pin, rng, error, walk_length);
}

#endif // VOLUME_COOLING_GAUSSIANS_HPP
