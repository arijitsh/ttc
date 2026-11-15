#ifndef GAUSSIAN_HELPERS_HPP
#define GAUSSIAN_HELPERS_HPP

#define EXP_CHORD_TOLERENCE 0.00000001

// evaluate the pdf of point p
template <typename Point, typename NT>
NT eval_exp(Point const& p, NT const& a)
{
    return std::exp(-a * p.squared_length());
}
template <typename Point, typename NT, typename MT>
NT eval_exp(Point const& p, MT const& inv_covariance_matrix, NT const& a_next, NT const& a_curr)
{
    Eigen::Matrix<NT, Eigen::Dynamic, 1> dist_vector = p.getCoefficients();
    NT mahalanobis_dist = dist_vector.transpose() * inv_covariance_matrix * dist_vector;
    NT log_ratio = (a_curr - a_next) * mahalanobis_dist;
    return std::exp(log_ratio);
}

template <typename Point, typename NT>
NT get_max(Point const& l, Point const& u, NT const& a_i)
{
    NT res;
    Point a = -1.0 * l;
    Point bef = u - l;
    Point b = (1.0 / std::sqrt((bef).squared_length())) * bef;
    Point z = (a.dot(b) * b) + l;
    NT low_bd = (l[0] - z[0]) / b[0], up_bd = (u[0] - z[0]) / b[0];
    if (low_bd * up_bd > 0)
    {
        //if(std::signbit(low_bd)==std::signbit(up_bd)){
        res = std::max(eval_exp(u, a_i), eval_exp(l, a_i));
    }
    else
    {
        res = eval_exp(z, a_i);
    }

    return res;
}


template <typename NT>
NT get_max_coord(NT const& l, NT const& u, NT const& a_i)
{
    const NT zero = 0;
    if (l < zero && u > zero)
    {
        return NT(1);
    }
    return std::max(std::exp(-a_i * l * l), std::exp(-a_i * u * u));
}

// Exact sample from the standard normal truncated to [alpha, beta] (ttc).
// The Gaussian chord samplers used to draw N(0,1) until it landed in the
// chord whenever the chord was >= 2 sd long, without checking that the chord
// carries any mass: a chord in the far tail (e.g. [8, 100] sd) has acceptance
// ~e^-32 and the walk hung (ttc: 2-D simplex, --seed 10).  Every branch below
// has acceptance bounded away from 0:
//   * 0 in [alpha, beta], width >= 2:  N(0,1) rejection (acceptance >= 0.47);
//   * narrow interval:                 uniform rejection under the maximum of
//                                      the density on the interval (>= e^-1);
//   * one-sided tail (alpha >= 0):     Robert's (1995) exponential-proposal
//                                      rejection, optimal rate
//                                      lambda = (alpha + sqrt(alpha^2 + 4)) / 2.
template <typename NT, typename RandomNumberGenerator>
NT sample_truncated_std_normal(NT alpha, NT beta, RandomNumberGenerator& rng)
{
    if (beta <= NT(0))
    {
        return -sample_truncated_std_normal(-beta, -alpha, rng);
    }
    if (alpha < NT(0))
    {
        if (beta - alpha >= NT(2))
        {
            while (true)
            {
                const NT r = rng.sample_ndist();
                if (r >= alpha && r <= beta) return r;
            }
        }
        while (true)  // |x| <= 2 on the interval: acceptance >= e^-2
        {
            const NT x = alpha + (beta - alpha) * rng.sample_urdist();
            if (rng.sample_urdist() <= std::exp(-x * x / NT(2))) return x;
        }
    }
    // 0 <= alpha < beta
    if ((beta - alpha) * (beta + alpha) <= NT(2))
    {
        while (true)  // density ratio on the interval >= e^-1
        {
            const NT x = alpha + (beta - alpha) * rng.sample_urdist();
            if (rng.sample_urdist() <= std::exp(-(x * x - alpha * alpha) / NT(2))) return x;
        }
    }
    const NT lambda = (alpha + std::sqrt(alpha * alpha + NT(4))) / NT(2);
    while (true)
    {
        const NT z = alpha - std::log(NT(1) - rng.sample_urdist()) / lambda;
        if (z > beta) continue;
        if (rng.sample_urdist() <= std::exp(-(z - lambda) * (z - lambda) / NT(2))) return z;
    }
}

#endif // GAUSSIAN_HELPERS_HPP
