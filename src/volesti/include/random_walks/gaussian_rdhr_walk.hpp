// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_RDHR_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_RDHR_WALK_HPP

#include "random_walks/gaussian_helpers.hpp"
#include "generators/boost_random_number_generator.hpp"


// Pick a point from the distribution exp(-a_i||x||^2) on the chord
template
<
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
void chord_random_point_generator_exp(Point &lower,
                                      Point & upper,
                                      const NT &a_i,
                                      Point &p,
                                      RandomNumberGenerator& rng)
{
    NT r, r_val, fn;
    Point bef = upper - lower;
    // pick from 1-dimensional gaussian if enough weight is inside polytope P
    if (a_i > EXP_CHORD_TOLERENCE && bef.squared_length() > NT(0))
    {
        // Exact truncated normal for every chord (ttc): the rectangle branch's
        // bound underflows to 0 on far-tail chords and then never accepts.
        Point a = -1.0 * lower;
        Point b = (1.0 / std::sqrt(bef.squared_length())) * bef;
        Point z = (a.dot(b) * b) + lower;
        // Signed positions of the chord ends along the unit direction b from
        // z (was (lower[0]-z[0])/b[0]: undefined when b[0] ~ 0), then an exact
        // truncated normal (the old rejection loop hung on far-tail chords).
        NT low_bd = (lower - z).dot(b);
        NT up_bd = (upper - z).dot(b);
        if (low_bd > up_bd) std::swap(low_bd, up_bd);
        const NT sd = NT(1) / std::sqrt(2.0 * a_i);
        r = sd * sample_truncated_std_normal(low_bd / sd, up_bd / sd, rng);
        p = (r * b) + z;

    // select using rejection sampling from a bounding rectangle
    } else {
        NT M = get_max(lower, upper, a_i);
        while (true) {
            r = rng.sample_urdist();//urdist(rng2);
            Point pef = r * upper;
            p = ((1.0 - r) * lower) + pef;
            r_val = M * rng.sample_urdist();//urdist(var.rng);
            fn = eval_exp(p, a_i);
            if (r_val < fn) {
                break;
            }
        }
    }
}

// Random directions hit-and-run walk with spherical Gaussian target distribution

struct GaussianRDHRWalk
{

    struct parameters {};
    parameters param;

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    Walk(Polytope&, Point const&, NT const&, RandomNumberGenerator&)
    {}

    Walk(Polytope&, Point const&, NT const&, RandomNumberGenerator&,
         parameters&)
    {}

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j = 0u; j < walk_length; ++j)
        {
            Point v = GetDirection<Point>::apply(p.dimension(), rng);
            std::pair <NT, NT> dbpair = P.line_intersect(p, v);

            NT min_plus = dbpair.first;
            NT max_minus = dbpair.second;
            Point upper = (min_plus * v) + p;
            Point lower = (max_minus * v) + p;

            chord_random_point_generator_exp(lower, upper, a_i, p, rng);
        }
    }
};

};


#endif // RANDOM_WALKS_GAUSSIAN_RDHR_WALK_HPP
