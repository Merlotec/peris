//
// Created by ncbmk on 10/17/24.
//

#ifndef PERIS_H
#define PERIS_H

#include <limits>
#include <string>

#ifndef __ssize_t_defined
typedef long int ssize_t;
#endif

/// Pareto efficient relative investment solver (PERIS).
namespace peris {
    template<typename A>
    concept AgentConcept = requires(A agent, double price, double quality)
    {
        { agent.item_id() } -> std::same_as<int>;
        { agent.income() } -> std::same_as<double>;
        { agent.utility(price, quality) } -> std::same_as<double>;
        { agent.debug_info() } -> std::same_as<std::string>;
    };

    template<typename I>
    concept ItemConcept = requires(I item)
    {
        { item.quality() } -> std::same_as<double>;
    };

    /// Represents a single agent in the model, who wishes to maximize utility.
    template<typename A, typename I>
        requires AgentConcept<A>&& ItemConcept<I>
    struct Allocation {
        /// The item being allocated.
        I item;

        /// The agent being allocated to this item. Note - this can change, hence is not const.
        A agent;

        /// The current allocation price.
        double price;

        /// The current allocation utility for the agent.
        double utility;

        /// Indicates a household that cannot be at this index (given a previous invalidation).
        /// If it allocated here we know we would experience an infinite loop.
        bool doublecross = false;

        /// Out of phase agents are not indifferent to the agent below, and hence cannot be shifted down allocation one while maintaining the same utility.
        bool out_of_phase = false;

        double favourite = 0;

        double quality() const {
            return item.quality();
        }

        void set_price(double price) {
            this->price = price;
            recalculate_utility();
        }

        void recalculate_utility() {
            utility = agent.utility(price, quality());
        }
    };

    template<typename A>
        requires AgentConcept<A>
    double indifferent_quality(A& agent, double price, double u_0, double y_min, double y_max, double epsilon = 1e-4,
        int max_iter = 100) {
        double lower = y_min;
        double upper = y_max;
        double mid = 0.0f;
        int iter = 0;

        while (iter < max_iter) {
            mid = (lower + upper) / 2.0f;
            double u_mid = agent.utility(price, mid);
            double diff = u_mid - u_0;

            if (std::abs(diff) < epsilon)
                return mid;

            if (diff > 0)
                upper = mid;
            else
                lower = mid;

            iter++;
        }

        // Return NaN (not a number) if solution was not found within the tolerance of epsilon
        return std::numeric_limits<double>::quiet_NaN();
    }

    template<typename A>
        requires AgentConcept<A>
    double indifferent_price(A& agent, double quality, double u_0, double x_min, double x_max, double epsilon = 1e-4,
        int max_iter = 100) {
        double lower = x_min;
        double upper = x_max;
        double mid = 0.0f;
        int iter = 0;
        double u_mid;
        while (iter < max_iter) {
            mid = (lower + upper) / 2.0f;
            u_mid = agent.utility(mid, quality);
            double diff = u_mid - u_0;

            if (std::abs(diff) < epsilon)
                return mid;

            // Because the utility function is increasing in quality we swap this from the solver for quality.
            if (diff > 0)
                lower = mid;
            else
                upper = mid;

            iter++;
        }

        // Return NaN if solution was not found within tolerance
        return std::numeric_limits<double>::quiet_NaN();
    }
}

#endif // PERIS_H