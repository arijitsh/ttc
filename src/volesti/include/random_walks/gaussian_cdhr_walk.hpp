// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_CDHR_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_CDHR_WALK_HPP

#include "sampling/sphere.hpp"
#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"

// Pick a point from the distribution exp(-a_i||x||^2) on the coordinate chord
template
<
    typename NT,
    typename RandomNumberGenerator
>
NT chord_random_point_generator_exp_coord(const NT &l,
                                          const NT &u,
                                          const NT &a_i,
                                          RandomNumberGenerator& rng)
{
    NT r, r_val, fn, dis;
    // pick from 1-dimensional gaussian if enough weight is inside polytope P
    if (a_i > EXP_CHORD_TOLERENCE)
    {
        // Exact truncated normal for every chord (ttc).  Both old branches
        // could loop forever on a chord in the far tail: the N(0,1)-until-
        // inside loop (acceptance ~ tail mass), and the rectangle rejection,
        // whose bound M = exp(-a l^2) underflows to 0 once a l^2 > 745 so that
        // `r_val < fn` (0 < 0) never holds.  See sample_truncated_std_normal.
        const NT sd = NT(1) / std::sqrt(2.0 * a_i);
        dis = sd * sample_truncated_std_normal(l / sd, u / sd, rng);

    // select using rejection sampling from a bounding rectangle
    }
    else
    {
        NT M = get_max_coord(l, u, a_i);
        while (true)
        {
            r = rng.sample_urdist();
            dis = (1.0 - r) * l + r * u;
            r_val = M * rng.sample_urdist();
            fn = std::exp(-a_i * dis * dis);
            if (r_val < fn)
            {
                break;
            }
        }
    }
    return dis;
}


// Coordinate directions hit-and-run walk with spherical Gaussian target distribution

struct GaussianCDHRWalk
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

    Walk(Polytope& P,
         Point const& p,
         NT const& a_i,
         RandomNumberGenerator &rng)
    {
        initialize(P, p, a_i, rng);
    }

    Walk(Polytope& P,
         Point const& p,
         NT const& a_i,
         RandomNumberGenerator& rng,
         parameters&)
    {
        initialize(P, p, a_i, rng);
    }

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
            auto rand_coord_prev = _rand_coord;
            _rand_coord = rng.sample_uidist();
            // ttc: recompute the slacks b - A p from scratch every 256
            // steps; the incremental update drifts over the long chains of
            // the cooling schedule.
            std::pair <NT, NT> bpair = (++_steps % 256 == 0)
                    ? P.line_intersect_coord(_p, _rand_coord, _lamdas)
                    : P.line_intersect_coord(_p, _p_prev, _rand_coord,
                                             rand_coord_prev, _lamdas);
            NT dis = chord_random_point_generator_exp_coord
                        (_p[_rand_coord] + bpair.second,
                         _p[_rand_coord] + bpair.first,
                         a_i,
                         rng);
            _p_prev = _p;
            _p.set_coord(_rand_coord, dis);
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           Point const& p,
                           NT const& a_i,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _rand_coord = rng.sample_uidist();
        _p = p;
        std::pair <NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord, _lamdas);
        NT dis = chord_random_point_generator_exp_coord
                    (_p[_rand_coord] + bpair.second,
                     _p[_rand_coord] + bpair.first,
                     a_i, rng);
        _p_prev = p;
        _p.set_coord(_rand_coord, dis);
    }

    unsigned int _rand_coord;
    unsigned long _steps = 0;
    Point _p;
    Point _p_prev;
    typename Point::Coeff _lamdas;
};

};

#endif // RANDOM_WALKS_GAUSSIAN_CDHR_WALK_HPP
