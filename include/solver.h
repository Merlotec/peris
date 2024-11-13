//
// Created by ncbmk on 10/22/24.
//
// Updated documentation and comments for clarity.
//

#ifndef SOLVER_H
#define SOLVER_H

#include "peris.h"
#include "render.h"
#include <cassert>
#include <optional>
#include <algorithm>

#ifndef __ssize_t_defined
typedef long int ssize_t;
#endif

#define BUFFER_MARGIN_FACTOR 10.f

namespace peris {

    /// Enumeration of possible results from the solver functions
    enum SolutionResult : int {
        err_cycle = -4,
        err_budget_constraint = -3,
        err_nan = -2,
        err_unknown = -1,
        terminated = 0,
        success = 1,
        repeat = 2,
    };

    struct SolutionTick {
        ssize_t agent_to_promote = -1;
        ssize_t agent_to_displace = -1;

        size_t i;

        SolutionResult result;

        SolutionTick(SolutionResult result, size_t i) : result(result), i(i) {}

        static SolutionTick promote(size_t i, ssize_t a) {
            auto tick = SolutionTick(SolutionResult::repeat, i);
            tick.agent_to_promote = a;
            return tick;
        }

        static SolutionTick displace(size_t i, ssize_t a) {
            auto tick = SolutionTick(SolutionResult::repeat, i);
            tick.agent_to_displace = a;
            return tick;
        }

        static SolutionTick doublecross(size_t i, size_t displace, size_t promote) {
            auto tick = SolutionTick(SolutionResult::repeat, i);
            tick.agent_to_displace = displace;
            tick.agent_to_promote = promote;
            return tick;
        }
    };

    /// @brief The Solver class calculates Pareto efficient allocations of items to agents.
    ///        Each agent receives exactly one item. The goal is to determine prices such that
    ///        no agent would prefer any other agent's allocation at the given prices.
    ///
    /// @tparam A The agent type, must satisfy AgentConcept (e.g., have income(), utility() methods)
    /// @tparam I The item type, must satisfy ItemConcept (e.g., have quality() method)
    ///
    template<typename A, typename I>
        requires AgentConcept<A>&& ItemConcept<I>
    class Solver {
        /// Vector containing the current allocations of items to agents
        /// Each Allocation contains an agent, item, price, and utility
        std::vector<Allocation<A, I> > allocations;

        /// @brief Swaps the agents between allocations at indices a and b,
        ///        and recalculates the utility for the affected allocations.
        ///
        /// @param a Index of the first allocation
        /// @param b Index of the second allocation
        ///
        static void swap_agents(std::vector<Allocation<A, I>>& allocations, size_t a, size_t b) {
            auto agent_a = allocations[a].agent;
            auto agent_b = allocations[b].agent;

            // Swap the agent objects.
            allocations[a].agent = agent_b;
            // Calculate utility for new agent.
            allocations[a].recalculate_utility();

            // Repeat for other agent
            allocations[b].agent = agent_a;
            allocations[b].recalculate_utility();
        }

        /// @brief Moves the agent at index 'a' to index 'b' in the allocations vector,
        ///        shifting agents between positions 'b' and 'a-1' up by one position.
        ///        This effectively inserts the agent at position 'a' into position 'b',
        ///        pushing other agents forward in the vector.
        ///
        /// @param a Source index of the agent to move (must be greater than b)
        /// @param b Destination index where the agent will be placed
        ///
        static void displace_up(std::vector<Allocation<A, I>>& allocations, size_t a, size_t b) {
            assert(b < a); // Ensure that the source index 'a' is greater than the destination index 'b'.

            // Temporarily store the agent at position 'b' as it will be overridden.
            auto free_agent = allocations[b].agent;

            // Move the agent from position 'a' to position 'b'.
            allocations[b].agent = allocations[a].agent;
            allocations[b].recalculate_utility(); // Update utility after changing the agent.

            // Shift agents from position 'b+1' to 'a' up by one position.
            // This loop moves each agent into the position of the previous agent.
            for (size_t i = b + 1; i <= a; ++i) {
                // Store the current agent to be moved in the next iteration.
                auto agent_buffer = allocations[i].agent;

                // Move the 'free_agent' into the current position.
                allocations[i].agent = free_agent;
                allocations[i].recalculate_utility(); // Update utility after changing the agent.

                // Update 'free_agent' for the next iteration.
                free_agent = agent_buffer;
            }
            // After the loop, 'free_agent' holds the agent originally at position 'a',
            // which has already been moved to position 'b', so it can be discarded.
        }

        /// @brief Moves the agent at index 'a' to index 'b' in the allocations vector,
        ///        shifting agents between positions 'a+1' and 'b' down by one position.
        ///        This effectively inserts the agent at position 'a' into position 'b',
        ///        pulling other agents backward in the vector.
        ///
        /// @param a Source index of the agent to move (must be less than b)
        /// @param b Destination index where the agent will be placed
        ///
        static void displace_down(std::vector<Allocation<A, I>>& allocations, size_t a, size_t b) {
            assert(b > a); // Ensure that the destination index 'b' is greater than the source index 'a'.

            // Temporarily store the agent at position 'b' as it will be overridden.
            auto free_agent = allocations[b].agent;

            // Move the agent from position 'a' to position 'b'.
            allocations[b].agent = allocations[a].agent;
            allocations[b].recalculate_utility(); // Update utility after changing the agent.

            // Shift agents from position 'b-1' down to 'a' by one position.
            // This loop moves each agent into the position of the next agent.
            for (ssize_t i = b - 1; i >= (ssize_t)a; --i) {
                // Store the current agent to be moved in the next iteration.
                auto agent_buffer = allocations[i].agent;

                // Move the 'free_agent' into the current position.
                allocations[i].agent = free_agent;
                allocations[i].recalculate_utility(); // Update utility after changing the agent.

                // Update 'free_agent' for the next iteration.
                free_agent = agent_buffer;
            }
            // After the loop, 'free_agent' holds the agent originally at position 'a',
            // which has already been moved to position 'b', so it can be discarded.
        }

        static double calculate_efficient_price(const Allocation<A, I>& l, const Allocation<A, I>& a, const double epsilon, const int max_iter) {

            double max_price = l.agent.income() - epsilon;
            double min_price = l.price == 0 ? epsilon : l.price - epsilon;

            double efficient_price = indifferent_price(l.agent, a.quality(), l.utility,
                min_price, max_price, epsilon, max_iter);

            if (std::isnan(efficient_price)) {
                // If unable to compute efficient price, use max price if it improves utility
                if (l.agent.utility(max_price, a.quality()) > l.utility) {
                    efficient_price = max_price;
                }
            }
            return efficient_price;
        }

        static double calculate_max_price(const Allocation<A, I>& l, const Allocation<A, I>& a, const double epsilon, const int max_iter) {

            double max_price = a.agent.income() - epsilon;
            double min_price = std::max(epsilon, l.price - epsilon);

            double u_lower = a.agent.utility(l.price, l.quality());

            double efficient_price = indifferent_price(a.agent, a.quality(), u_lower,
                min_price, max_price, epsilon, max_iter);

            if (std::isnan(efficient_price)) {
                // If unable to compute efficient price, use max price if it improves utility
                // TODO: confirm
                if (a.agent.utility(max_price, a.quality()) > u_lower) {
                    efficient_price = max_price;
                }
            }
            return efficient_price;
        }

        static double calculate_pulled_price(const Allocation<A, I>& l, const Allocation<A, I>& a, const double epsilon, const int max_iter) {

            double max_price = l.agent.income() - epsilon;
            double min_price = std::max(epsilon, l.price - epsilon);

            double efficient_price = indifferent_price(a.agent, l.quality(), a.utility,
                min_price, max_price, epsilon, max_iter);

            // if (isnan(efficient_price)) {
            //     // If unable to compute efficient price, use max price if it improves utility
            //     if (a.agent.utility(max_price, l.quality()) > a.utility) {
            //         efficient_price = max_price;
            //     }
            // }
            return efficient_price;
        }

    public:
        ///
        /// @brief Constructs the Solver with given agents and items.
        ///
        ///        - Asserts that the number of agents equals the number of items.
        ///        - Sorts agents in order of increasing income.
        ///        - Sorts items in order of increasing quality.
        ///        - Initializes Allocations vector pairing each agent with an item,
        ///          setting an initial guess price based on the agent's income.
        ///
        /// @param agents Vector of agents
        /// @param items  Vector of items
        /// @param guess_factor Initial guess factor for prices (used as price = guess_factor * agent.income())
        ///
        Solver(std::vector<A> agents, std::vector<I> items, double guess_factor) {
            // Ensure that there is one item per agent (numbers of each are the same).
            assert(agents.size() == items.size());

            // Sort agents by income (increasing order), so lower-income agents are first.
            std::sort(agents.begin(), agents.end(), [](A a, A b) { return a.income() < b.income(); });
            // Sort items by quality (increasing order), so lower-quality items are first.
            std::sort(items.begin(), items.end(), [](I a, I b) { return a.quality() < b.quality(); });

            allocations.reserve(items.size());
            // Combine the items and agents into allocations.
            for (size_t i = 0; i < items.size(); i++) {
                // Set the initial guess price to an arbitrary guess according to the function p_i = guess_factor * y_i
                const double guess_price = guess_factor * agents[i].income();
                Allocation<A, I> allocation = {
                    .item = items[i],
                    .agent = agents[i],
                    .price = guess_price,
                    .utility = agents[i].utility(guess_price, items[i].quality())
                };

                allocations.push_back(allocation);
            }
        }

        ///
        /// @brief Returns a copy of the current allocations vector.
        ///
        /// @return A vector of Allocations containing the current allocation state.
        ///
        std::vector<Allocation<A, I>> get_allocations() {
            return std::vector<Allocation<A, I>>(this->allocations);
        }

        ///
        /// @brief Solves the allocation model to achieve Pareto efficiency among agents.
        ///        The algorithm assigns items to agents in a way that no agent can be made better off
        ///        without making another agent worse off. It iteratively adjusts allocations and prices,
        ///        possibly swapping agents to improve overall efficiency.
        ///
        /// @param render_state Pointer to a RenderState object for visualization (can be nullptr)
        /// @param epsilon The tolerance for numerical approximations (default is 1e-5)
        /// @param max_iter The maximum number of iterations for convergence (default is 200)
        ///
        /// @return A SolutionResult indicating success or type of error.
        ///
        SolutionResult solve(RenderState<A, I>* render_state, double epsilon = 1e-5, int max_iter = 400) {
            // If there are no agents, return success.
            if (allocations.empty()) {
                return SolutionResult::success;
            }
            //return align(render_state, 0, epsilon, max_iter);
            // // Perform initial alignment.
            // SolutionResult res;
            // if ((res = align(render_state, 0, epsilon, max_iter)) < 0) {
            //     // Exit command.
            //     return res;
            // }
            // std::cout << "Allocated" << std::endl;
            // SolutionResult pres;
            // while ((pres = push(epsilon, max_iter)) == SolutionResult::repeat) {
            //     if (render_state != nullptr) {
            //         if (!render_state->draw_allocations(this->allocations, -1)) {
            //             return SolutionResult::terminated;
            //         }
            //     }
            // }
            //
            // if (pres != SolutionResult::success) {
            //     return pres;
            // }
            //
            // return SolutionResult::success;

            // If there is a double crossing agent that is being sorted, this should be its index.
            // This lets us know not to promote the agent when we experience the double cross.
            ssize_t doublecross_idx = -1;
            int doublecross_id = -1;

            size_t i = 0;
            size_t reserve = 0;
            size_t to_skip = 0;

            bool pause = false;
            bool tick_visual = false;
            bool just_displace = false;

            SolutionTick tick = SolutionTick(repeat, 0);

            std::vector<Allocation<A, I>> savestate{};
            size_t try_offset = 0;

            bool try_swap = false;

            bool enable_doublecross = false;

            while (true) {
                if (render_state != nullptr) {
                    if (to_skip > 0) {
                        --to_skip;
                    }
                    else {
                        auto rc = render_state->draw_allocations(this->allocations, tick.i);
                        if (pause) {
                            tick_visual = false;
                            if (rc == RenderCommand::pause) {
                                pause = false;
                            }
                            else if (rc == RenderCommand::tick) {
                                tick_visual = true;
                            }
                        }
                        else {
                            if (rc == RenderCommand::pause) {
                                pause = true;
                            }
                        }

                        if (rc == RenderCommand::enable_doublecross) {
                            enable_doublecross = !enable_doublecross;
                        }

                        if (rc == RenderCommand::skip) {
                            to_skip = 100;
                        }

                        if (pause && !tick_visual) {
                            continue;
                        }
                    }
                }

                if (i >= allocations.size()) {
                    return SolutionResult::success;
                    if (just_displace) {
                        return SolutionResult::success;
                    } else {
                        i = 0;
                        just_displace = true;
                        continue;
                    }
                }

                // Align
                tick = try_align(allocations, i, allocations.size() - reserve, 0.f, epsilon, max_iter);

                if (tick.result == SolutionResult::success) {
                    if (reserve == 0) {
                        return SolutionResult::success;
                        if (just_displace) {
                            return SolutionResult::success;
                        } else {
                            i = 0;
                            just_displace = true;
                            continue;
                        }
                    }
                    else {
                        // Perform reserve allocation.
                        // Move reserve to correct place.
                        const size_t r = allocations.size() - reserve;
                        // Find optimal location - this assumes everything below r is perfectly allocated.
                        const size_t pref = most_preferred(allocations, r, false, epsilon);
                        // TODO: maybe check end of non-reserve to see if we should just slot it in?

                        assert(pref < r);

                        displace_up(allocations, r, pref);
                        doublecross_idx = pref; // So that we know what to use as the current double-cross target.
                        doublecross_id = allocations[pref].agent.item_id();
                        allocations[pref].doublecross = true;
                        i = pref + 1;
                        savestate = std::vector<Allocation<A, I>>(allocations); // Clone existing state.
                        try_offset = 0; // Reset this because we are allocating a new reserve agent.
                        --reserve;
                        continue;
                    }
                }
                else if (tick.result != SolutionResult::repeat) {
                    return tick.result;
                }

                bool increase_offset = false;

                if (tick.agent_to_displace != -1) {
                    if (tick.agent_to_promote != -1 && !just_displace) {
                        if (tick.agent_to_promote <= doublecross_idx) {
                            //assert(allocations[tick.agent_to_promote].agent.item_id() == doublecross_id);
                            // Push back because we know that this is optimal position.
                            // auto res = pull_back(allocations, tick.i, epsilon, max_iter);
                            // if (res <= 0) {
                            //     return res;
                            // }
                            // i = tick.i + 1;

                            // Reset to previous state and reallocate.
                            // Restore savestate.
                            //increase_offset = true;
                            if (enable_doublecross) {
                                pause = true;
                            }
                            i = tick.i + 1;
                        }
                        else {
                            // try_swap = !try_swap;
                            // if (try_swap) {
                            //     displace_up(allocations, tick.i, tick.agent_to_displace);
                            //     i = tick.agent_to_displace;
                            //     continue;
                            // } else {
                            //     displace_down(allocations, tick.agent_to_promote, allocations.size() - 1);p
                            //     i = tick.agent_to_promote;
                            //     continue;
                            // }
                            if (tick.agent_to_promote < doublecross_idx) {
                                --doublecross_idx;
                            }
                            displace_down(allocations, tick.agent_to_promote, allocations.size() - 1);
                            i = tick.agent_to_promote;

                            // displace_up(allocations, tick.i, tick.agent_to_displace);
                            // i = tick.agent_to_displace;
                            // if (tick.agent_to_displace > tick.agent_to_promote) {
                            //     displace_up(allocations, tick.i, tick.agent_to_displace);
                            //     displace_down(allocations, tick.agent_to_promote, allocations.size() - 1);
                            //     i = tick.agent_to_promote;
                            // } else if (tick.agent_to_displace < tick.agent_to_promote) {
                            //     displace_up(allocations, tick.i, tick.agent_to_displace);
                            //     displace_down(allocations, tick.agent_to_promote + 1, allocations.size() - 1);
                            //     i = tick.agent_to_displace;
                            // } else {
                            //     displace_down(allocations, tick.agent_to_promote, allocations.size() - 1);
                            //     i = tick.agent_to_promote;
                            // }

                            ++reserve;
                        }
                    }
                    else {
                        assert(tick.agent_to_displace < tick.i); // The agent to displace should be at a lower index.
                        // Displace the current agent 'a' to position 'agent_to_displace', shifting other agents accordingly.
                        // if (doublecross_idx != -1 && tick.agent_to_displace < doublecross_idx + try_offset) {
                        //     // We are displacing an offset agent so this will ruin our structure and potentially cause a loop.
                        //
                        //     if (tick.i == doublecross_idx) {
                        //         i = tick.i + 1;
                        //     } else {
                        //         increase_offset = true;
                        //     }
                        // } else {
                        //     displace_up(allocations, tick.i, tick.agent_to_displace);
                        //     i = tick.agent_to_displace;
                        // }
                        displace_up(allocations, tick.i, tick.agent_to_displace);
                        i = tick.agent_to_displace;
                    }

                    if (increase_offset) {
                        // if (tick.agent_to_promote != -1 && tick.i + 1 < doublecross_idx + try_offset) {
                        //    // return SolutionResult::err_unknown;
                        //     if (tick.agent_to_promote == doublecross_idx) {
                        //         // Ignore - within epsilon (must be).
                        //         i = tick.i + 1;
                        //     } else {
                        //         displace_up(allocations, tick.i, tick.agent_to_displace);
                        //         i = tick.agent_to_displace;
                        //         continue;
                        //     }
                        // }

                        // if (false && doublecross_idx + try_offset + 1 >= allocations.size() - reserve - 1) { // We are on the 'frontier'.
                        //     i = allocations.size() - reserve; // Allocate next reserve agent.
                        //     continue;
                        // } else {
                        //     allocations = std::vector<Allocation<A, I>>(savestate);
                        //     ++try_offset; // Try the next offset.
                        //     if (try_offset > 0) {
                        //         auto res = pull_to(allocations, doublecross_idx + 1, try_offset, epsilon, max_iter);
                        //         if (res < 1) {
                        //             return res;
                        //         }
                        //     }
                        //
                        //     i = doublecross_idx + try_offset;
                        //     if (i >= allocations.size()) {
                        //         return SolutionResult::success;
                        //     }
                        // }
                    }
                }
                else if (tick.agent_to_promote != -1) {
                    // if (tick.agent_to_promote > doublecross_idx) {
                    //     displace_down(allocations, tick.agent_to_promote, allocations.size() - 1);
                    //     i = tick.agent_to_promote;
                    //     ++reserve;
                    // } else {
                    //     i = tick.i + 1;
                    // }
                    i = tick.i + 1;
                }

            }

        }

    private:
        ///
        /// @brief Finds the index of the most preferred allocation for agent at index 'i' among allocations.
        ///
        ///        It goes through other allocations and computes the utility the agent at index 'i' would get
        ///        from those allocations (after adjusting the price slightly to avoid division by zero or
        ///        price equality issues). It returns the index of the allocation that provides the maximum utility.
        ///
        /// @param allocations
        /// @param i Index of the agent/allocation to consider.
        /// @param search_above If true, search only for allocations with higher indices than 'i'
        /// @param epsilon Tolerance used for price adjustments
        ///
        /// @return The index of the most preferred allocation for agent at index 'i'
        ///
        static size_t most_preferred(std::vector<Allocation<A, I>>& allocations, size_t i, bool search_above, double epsilon) {
            assert(i < allocations.size());
            double u_max;
            if (allocations[i].agent.income() < allocations[i].price) {
                u_max = allocations[i].agent.utility(allocations[i].price, allocations[i].item.quality());
            }
            else {
                u_max = NAN;
            }
            size_t i_max = i;

            size_t limit = search_above ? allocations.size() : i;

            // Iterate over allocations to find the most preferred one for agent i
            for (size_t j = 0; j < limit; ++j) {
                if (i != j) {
                    // Compute the utility of agent i if they were allocated allocation j
                    const double u_alt = allocations[i].agent.utility(allocations[j].price, allocations[j].item.quality());
                    if (u_alt > u_max + epsilon || (std::isnan(u_max) && !std::isnan(u_alt))) {
                        u_max = u_alt;
                        i_max = j;
                    }
                }
            }

            return i_max;
        }

        ///
        /// @brief Aligns the allocations starting from index 'i' to achieve Pareto efficiency.
        ///
        ///        The align function is the core of the algorithm. It assigns allocations and adjusts prices
        ///        such that each agent prefers their own allocation over any other allocation,
        ///        i.e., no agent would prefer another allocation at the given prices.
        ///
        ///        It may involve adjusting prices, displacing agents (moving them up or down in the allocation list),
        ///        to ensure that allocations are Pareto efficient.
        ///
        /// @param render_state Pointer to a RenderState for visualization (can be nullptr)
        /// @param i Starting index for alignment
        /// @param epsilon Numerical tolerance for calculations
        /// @param max_iter Maximum number of iterations for numerical methods
        ///
        /// @return A SolutionResult indicating success or type of error.
        ///
        SolutionResult align(RenderState<A, I>* render_state, size_t i, double epsilon, int max_iter = 100) {
            // Initialize the first allocation if starting from index 0
            if (i == 0) {
                // Set the price of the first allocation to zero.
                // This essentially 'anchors' the algorithm so that all other allocations can be distributed around this value.
                allocations[0].set_price(0.0f);

                // We have allocated agent 0 so we can start at agent 1.
                i = 1;
            }

            // Keeps track of the highest index that has so far been reached.
            // This is used for rendering so that we only re-render when the next allocation has been reached.
            size_t head = i;

            std::vector<int> reserve_set{};

            size_t reserve_count = 0;

            size_t repeat = 0;
            size_t repeat_n = 0;

            bool break_loop = false;

            bool pause = false;
            bool tick = false;
            bool swap_always = false;

            size_t to_skip = 0;

            std::optional<int> immune_id{};

            // Iterate through each agent starting from index i
            while (i < allocations.size()) {
                // Visualization code if render_state is provided
                if (render_state != nullptr) {
                    if (to_skip > 0) {
                        --to_skip;
                    }
                    else {
                        auto rc = render_state->draw_allocations(this->allocations, i);
                        if (pause) {
                            tick = false;
                            if (rc == RenderCommand::pause) {
                                pause = false;
                            }
                            else if (rc == RenderCommand::tick) {
                                tick = true;
                            }
                        }
                        else {
                            if (rc == RenderCommand::pause) {
                                pause = true;
                            }
                        }

                        if (rc == skip) {
                            // Skip n elements
                            to_skip = head * allocations.size();
                        }
                        else if (rc == RenderCommand::enable_swap_always) {
                            swap_always = !swap_always;
                        }
                    }
                }

                if (pause && !tick) {
                    continue;
                }
                head = std::max(i, head);

                // Initialize variables to keep track of agents to displace or promote
                ssize_t agent_to_displace = -1;
                ssize_t agent_to_promote = -1;

                Allocation<A, I>& a = allocations[i];     // Current allocation
                Allocation<A, I>& l = allocations[i - 1]; // Previous allocation

                // if (i >= allocations.size() - reserve_count) {
                //     std::vector<int> new_reserve_set{};
                //     for (size_t j = i; j < allocations.size(); ++j) {
                //         new_reserve_set.push_back(allocations[j].agent.item_id());
                //     }
                //     std::sort(new_reserve_set.begin(), new_reserve_set.end());
                //
                //     // If we are in a cycle, allocate the lowest item in the cycle and do not allow any switches for agents pushed out by this double-crossing agent.
                //     // In other words just keep the bad agent where it is.
                //     // To do this we mark the id of the agent to be 'immune' from promotion.
                //     if (reserve_set == new_reserve_set) {
                //         //immune_id = a.agent.item_id();
                //         std::cout << "n = " << new_reserve_set.size() << std::endl;
                //         if (reserve_set.size() == repeat_n) {
                //             ++repeat;
                //         } else {
                //             repeat = 0;
                //         }
                //         repeat_n = reserve_set.size();
                //     } else {
                //         immune_id = std::optional<int>();
                //         repeat = 0;
                //     }
                //
                //     if (repeat > 3) {
                //         immune_id = a.agent.item_id();
                //     }
                //
                //     reserve_set = std::vector<int>(new_reserve_set);
                //     --reserve_count;
                // }

                double efficient_price = calculate_efficient_price(l, a, epsilon, max_iter);

                // Check if the efficient price exceeds the current agent's income
                if (efficient_price + epsilon > a.agent.income()) {
                    // Find the most preferred allocation for the current agent among those below
                    size_t new_i = most_preferred(i, false, epsilon); // Do not search above because not yet allocated.
                    agent_to_displace = new_i;
                }

                if (agent_to_displace == -1) {

                    for (ssize_t j = i - 1; j >= 0; --j) {
                        const Allocation<A, I>& prev = allocations[j];
                        // Check if the previous agent prefers the current allocation at efficient price
                        if (prev.agent.income() > efficient_price) {
                            if (j < i - 1 && prev.agent.utility(efficient_price, a.quality()) > prev.utility + epsilon) {
                                // Mark this agent to promote
                                agent_to_promote = j;
                                // displace this agent down there.
                                //agent_to_displace = j;
                                // Calculate the new efficient price by making this 'prev' agent indifferent.
                                double new_price = calculate_efficient_price(prev, a, epsilon, max_iter);

                                if (new_price > efficient_price) {
                                    efficient_price = new_price;
                                }
                            }
                        }
                    }

                    // Check if the efficient price exceeds the current agent's income
                    if (efficient_price + epsilon > a.agent.income()) {
                        // Find the most preferred allocation for the current agent among those below
                        size_t new_i = most_preferred(i, false, epsilon); // Do not search above because not yet allocated.
                        agent_to_displace = new_i;
                    }

                    if (agent_to_displace == -1) {
                        // Check if this agent prefers any previous allocations
                        double u_max = a.agent.utility(efficient_price, a.quality());

                        if (std::isnan(u_max))
                            return SolutionResult::err_nan;
                        for (ssize_t j = i - 1; j >= 0; --j) {
                            const Allocation<A, I>& prev = allocations[j];
                            // Check if the agent can afford the previous allocation
                            if (a.agent.income() > prev.price + epsilon) {
                                double u_prev = a.agent.utility(prev.price, prev.quality());
                                if (std::isnan(u_prev))
                                    return SolutionResult::err_nan;
                                if (u_prev > u_max + epsilon) {
                                    // The current agent 'a' prefers 'prev''s allocation; mark 'prev' as the agent to displace.
                                    u_max = u_prev;
                                    agent_to_displace = j;
                                }
                            }
                        }

                        // If no displacement is needed, update the current allocation's price and utility.
                        if (agent_to_displace == -1) {
                            if (efficient_price + epsilon > a.agent.income()) {
                                // Find best place to move agent to.
                                size_t new_i = most_preferred(i, false, epsilon); // Do not search above because not yet allocated.
                                agent_to_displace = new_i;
                            }
                            else {
                                // Update the allocation with the efficient price and utility
                                double efficient_utility = a.agent.utility(efficient_price, a.quality());
                                if (std::isnan(efficient_utility))
                                    return SolutionResult::err_nan;

                                a.price = efficient_price;
                                a.utility = efficient_utility;
                            }
                        }
                    }
                }

                if (agent_to_displace >= 0) {
                    // Handle displacement or promotion of agents
                    if (agent_to_promote >= 0 && !swap_always) {
                        // Displace the agent to the reserve area
                        // Check if this agent is immune from being promoted.
                        if (allocations[agent_to_promote].agent.item_id() != immune_id) {
                            displace_down(agent_to_promote, allocations.size() - 1);
                            i = std::max(agent_to_promote, 1L);
                            ++reserve_count;
                        }
                        else {
                            ++i;
                        }
                    }
                    else {
                        assert(agent_to_displace < i); // The agent to displace should be at a lower index.
                        // Displace the current agent 'a' to position 'agent_to_displace', shifting other agents accordingly.
                        displace_up(i, agent_to_displace);
                        i = std::max(agent_to_displace, 1L);
                    }
                }
                else {
                    // The current allocation is successful, so we can move on to the next one.
                    ++i;
                }
            }
            return SolutionResult::success;
        }

        static SolutionTick try_align(std::vector<Allocation<A, I>>& allocations, size_t i, size_t end, double p0, double epsilon, int max_iter = 100) {
            // Initialize the first allocation if starting from index 0
            if (i == 0) {
                // Set the price of the first allocation to zero.
                // This essentially 'anchors' the algorithm so that all other allocations can be distributed around this value.
                allocations[0].set_price(p0);

                // We have allocated agent 0 so we can start at agent 1.
                i = 1;
            }

            // Keeps track of the highest index that has so far been reached.
            // This is used for rendering so that we only re-render when the next allocation has been reached.
            size_t head = i;

            // std::vector<int> reserve_set{};
            //
            // size_t reserve_count = 0;
            //
            // bool break_loop = false;

            // Iterate through each agent starting from index i
            while (i < end) {
                head = std::max(i, head);

                // if (i >= allocations.size() - reserve_count) {
                //     std::vector<int> new_reserve_set{};
                //     for (size_t j = i; j < allocations.size(); ++j) {
                //         new_reserve_set.push_back(allocations[j].agent.item_id());
                //     }
                //     std::sort(new_reserve_set.begin(), new_reserve_set.end());
                //
                //     break_loop = reserve_set == new_reserve_set;
                //
                //     if (break_loop) {
                //         std::cout << "breakl" << std::endl;
                //     }
                //
                //     reserve_set = new_reserve_set;
                //     --reserve_count;
                // }

                // Initialize variables to keep track of agents to displace or promote
                ssize_t agent_to_displace = -1;
                ssize_t agent_to_promote = -1;

                Allocation<A, I>& a = allocations[i];     // Current allocation
                Allocation<A, I>& l = allocations[i - 1]; // Previous allocation

                double efficient_price = calculate_efficient_price(l, a, epsilon, max_iter);

                // Check if the efficient price exceeds the current agent's income
                if (efficient_price + epsilon > a.agent.income()) {
                    // Find the most preferred allocation for the current agent among those below
                    size_t new_i = most_preferred(allocations, i, false, epsilon); // Do not search above because not yet allocated.
                    agent_to_displace = new_i;
                }

                if (agent_to_displace == -1) {

                    for (ssize_t j = i - 1; j >= 0; --j) {
                        const Allocation<A, I>& prev = allocations[j];
                        // Check if the previous agent prefers the current allocation at efficient price
                        if (prev.agent.income() > efficient_price) {
                            if (j < i - 1 && prev.agent.utility(efficient_price, a.quality()) > prev.utility + epsilon) {
                                // Mark this agent to promote
                                agent_to_promote = j;

                                // Calculate the new efficient price by making this 'prev' agent indifferent.
                                double new_price = calculate_efficient_price(prev, a, epsilon, max_iter);

                                if (new_price > efficient_price) {
                                    efficient_price = new_price;
                                }
                            }
                        }
                    }

                    // Check if the efficient price exceeds the current agent's income
                    if (efficient_price + epsilon > a.agent.income()) {
                        // Find the most preferred allocation for the current agent among those below
                        size_t new_i = most_preferred(allocations, i, false, epsilon); // Do not search above because not yet allocated.
                        agent_to_displace = new_i;
                    }

                    if (agent_to_displace == -1) {
                        // Check if this agent prefers any previous allocations
                        double u_max = a.agent.utility(efficient_price, a.quality()) + epsilon;

                        if (std::isnan(u_max))
                            return SolutionTick(SolutionResult::err_nan, i);
                        for (ssize_t j = i - 1; j >= 0; --j) {
                            const Allocation<A, I>& prev = allocations[j];
                            // Check if the agent can afford the previous allocation
                            if (a.agent.income() > prev.price + epsilon) {
                                double u_prev = a.agent.utility(prev.price, prev.quality());
                                if (std::isnan(u_prev))
                                    return SolutionTick(SolutionResult::err_nan, i);
                                if (u_prev > u_max) {
                                    // The current agent 'a' prefers 'prev''s allocation; mark 'prev' as the agent to displace.
                                    u_max = u_prev;
                                    agent_to_displace = j;
                                }
                            }
                        }
                        double efficient_utility = a.agent.utility(efficient_price, a.quality());
                        // If no displacement is needed, update the current allocation's price and utility.
                        if (agent_to_displace == -1) {
                            if (efficient_price + epsilon > a.agent.income()) {
                                // Find best place to move agent to.
                                size_t new_i = most_preferred(allocations, i, false, epsilon); // Do not search above because not yet allocated.
                                agent_to_displace = new_i;
                                return SolutionTick::displace(new_i, i);
                            }
                            else {
                                // Update the allocation with the efficient price and utility
                                if (std::isnan(efficient_utility))
                                    return SolutionTick(SolutionResult::err_nan, i);
                            }
                        }

                        if (!std::isnan(efficient_utility)) {
                            // Set price and util anyway.
                            a.price = efficient_price;
                            a.utility = efficient_utility;
                        }
                    }
                }

                if (agent_to_displace >= 0) {
                    // Handle displacement or promotion of agents
                    if (agent_to_promote >= 0) {
                        // Displace the agent to the reserve area
                        return SolutionTick::doublecross(i, agent_to_displace, agent_to_promote);
                    }
                    else {
                        assert(agent_to_displace < i); // The agent to displace should be at a lower index.
                        // Displace the current agent 'a' to position 'agent_to_displace', shifting other agents accordingly.
                        return SolutionTick::displace(i, agent_to_displace);
                    }
                }
                else {
                    if (agent_to_promote >= 0) {
                        return SolutionTick::promote(i, agent_to_promote);
                    }
                    // The current allocation is successful, so we can move on to the next one.
                    ++i;
                }
            }
            return SolutionTick(SolutionResult::success, i);
        }

        static SolutionResult pull_back(std::vector<Allocation<A, I>>& allocations, const size_t i, double epsilon, int max_iter = 100) {
            for (size_t j = i;; --j) {
                // Start with i-1.
                Allocation<A, I>& a = allocations[j];     // Current allocation
                Allocation<A, I>& l = allocations[j - 1]; // Previous allocation

                // Only need to pull back if we prefer last.
                if (a.agent.utility(l.price, l.quality()) > a.utility + epsilon) {
                    // We want to swap!!!
                    // Shift l up so that we dont want to switch to it.
                    const double pulled_price = calculate_pulled_price(l, a, epsilon, max_iter);

                    // if (isnan(pulled_price))
                    //     return SolutionResult::err_nan;
                    // if (pulled_price > l.agent.income() - epsilon) {
                    //     return SolutionResult::err_budget_constraint;
                    // }
                    // l.set_price(pulled_price);

                    if (std::isnan(pulled_price)) {
                        if (j >= 2) {
                            displace_up(allocations, j - 1, j - 2);
                            j += 1;
                        }
                        else {
                            return SolutionResult::err_nan;
                        }
                    }
                    else {
                        l.set_price(pulled_price);
                    }
                }
                else {
                    return SolutionResult::success;
                }
            }
        }

        static SolutionResult pull_to(std::vector<Allocation<A, I>>& allocations, const size_t i, const size_t n, double epsilon, int max_iter = 100) {
            assert(i >= 1);
            for (size_t j = i; j < i + n; ++j) {
                // Start with i-1.
                Allocation<A, I>& a = allocations[j];     // Current allocation
                Allocation<A, I>& l = allocations[j - 1]; // Previous allocation

                // We want to swap!!!
                // Shift a up so that we dont want to switch to it.
                const double max_price = calculate_max_price(l, a, epsilon, max_iter);

                // if (isnan(pulled_price))
                //     return SolutionResult::err_nan;
                // if (pulled_price > l.agent.income() - epsilon) {
                //     return SolutionResult::err_budget_constraint;
                // }
                // l.set_price(pulled_price);
                if (max_price > a.agent.income()) {
                    return SolutionResult::err_budget_constraint;
                }

                if (std::isnan(max_price)) {
                    return SolutionResult::err_nan;
                }
                a.set_price(max_price);
            }
        }


        ///
        /// @brief Adjusts prices to ensure that no agent prefers another allocation.
        ///
        ///        The push method iteratively increases prices where necessary to prevent agents
        ///        from preferring allocations assigned to other agents. It ensures that the prices
        ///        are set such that each agent's utility from their own allocation is at least as
        ///        high as any other allocation they could afford.
        ///
        /// @param epsilon Numerical tolerance for calculations
        /// @param max_iter Maximum number of iterations for numerical methods
        ///
        /// @return A SolutionResult indicating whether any prices were updated or an error occurred.
        ///
        SolutionResult push(double epsilon, int max_iter) {
            size_t updated = 0;
            // Start from the top allocation and move backwards
            for (ssize_t i = allocations.size() - 1; i >= 0; --i) {
                Allocation<A, I>& a = allocations[i];     // Current allocation
                double efficient_price = a.price;

                // Check if any earlier agents prefer the current allocation at 'efficient_price'
                for (ssize_t j = allocations.size() - 1; j >= 0; --j) {
                    if (i == j) {
                        continue;
                    }
                    Allocation<A, I>& other = allocations[j];

                    // Ensure that 'other' can afford the current allocation
                    if (efficient_price + epsilon < other.agent.income()) {
                        // If agent 'other' prefers the current allocation at 'efficient_price'
                        if (other.agent.utility(efficient_price, a.quality()) > other.utility + epsilon) {
                            // Update the efficient price to prevent 'other' from preferring it
                            double max_price = other.agent.income() - epsilon;
                            double min_price = other.price == epsilon;

                            double new_price = indifferent_price(other.agent, a.quality(), other.utility, min_price, max_price, epsilon, max_iter);
                            if (std::isnan(new_price)) {
                                if (other.agent.utility(max_price, a.quality()) > other.utility) {
                                    efficient_price = max_price;
                                }
                                else {
                                    return SolutionResult::err_nan;
                                }
                            }
                            else if (new_price > efficient_price) {
                                efficient_price = new_price;
                            }
                        }
                    }
                }

                // Update to reflect new efficient price
                if (efficient_price > a.price + epsilon) {
                    // Update efficient price and utility
                    if (efficient_price + epsilon > a.agent.income()) {
                        return SolutionResult::err_budget_constraint;
                    }
                    a.set_price(efficient_price);
                    ++updated;
                }
            }

            if (updated > 0) {
                return SolutionResult::repeat;
            }
            else {
                return SolutionResult::success;
            }
        }



    public:
        size_t equilibriate_demand(double delta, double epsilon = 1e-5) {

            for (size_t i = 0; i < allocations.size(); i++) {
                allocations[i].favourite = 0.0;
            }

            struct FavReg {
                size_t i;
                double u;
            };

            for (size_t i = 0; i < allocations.size(); i++) {
                auto a = allocations[i];
                double u_max = 0.0;
                std::vector<FavReg> i_max{};
                for (size_t j = 0; j < allocations.size(); j++) {
                    // Determine whether
                    Allocation<A, I>& other = allocations[j];
                    if (a.agent.income() > other.price) {
                        double u = a.agent.utility(other.price, other.quality());
                        if (u > u_max) {
                            u_max = u;

                            i_max.erase(std::remove_if(i_max.begin(), i_max.end(),
                                [u_max, epsilon](FavReg f) {
                                    return f.u + epsilon < u_max;
                                }),
                                i_max.end());

                            i_max.push_back({ .i = j, .u = u });
                        }
                        else if (u + epsilon > u_max) {
                            i_max.push_back({ .i = j, .u = u });
                        }
                    }
                }

                for (auto x : i_max) {
                    allocations[x.i].favourite += 1.0 / ((double) i_max.size());
                }
            }

            size_t movements = 0;

            for (size_t i = 1; i < allocations.size(); i++) {
                if (allocations[i].favourite >= 2) {
                    // Increase price
                    allocations[i].price += delta * allocations[i].favourite; // Total movement right needs to match total movement left.
                    allocations[i].set_price(allocations[i].price);
                    ++movements;
                }
                else if (allocations[i].favourite <= 0) {
                    allocations[i].price -= delta;
                    allocations[i].set_price(allocations[i].price);
                    ++movements;
                }
            }

            return movements;

        }

        double max_difference(std::vector<Allocation<A, I>>& other_allocations) {
            assert(other_allocations.size() == allocations.size());

            double max = 0.0;

            for (size_t i = 0; i < allocations.size(); ++i) {
                double diff = std::abs(allocations[i].price - other_allocations[i].price);
                if (diff > max) {
                    diff = max;
                }
            }
            return max;
        }

        // Reassigns in order of income from lowest to highest.
        bool reassign_ordered() {
            std::vector<A> ord = get_agents();
            std::sort(ord.begin(), ord.end(), [](A a, A b) { return a.income() < b.income(); });

            std::vector<size_t> avail(allocations.size());

            // Populate the vector so each element equals its index
            for (size_t i = 0; i < avail.size(); ++i) {
                avail[i] = i;
            }

            for (auto a : ord) {
                // Allocate in order
                // Find fav agent that is still unassigned.
                double u_max = 0.0;
                ssize_t j_max = -1;
                for (size_t j = 0; j < avail.size(); ++j) {
                    size_t i = avail[j];
                    if (a.income() > allocations[i].price) {
                        double u = a.utility(allocations[i].price, allocations[i].quality());
                        // Check if util is max:
                        if (u > u_max) {
                            j_max = j;
                            u_max = u;
                        }
                    }
                }

                if (j_max >= 0) {
                    allocations[avail[j_max]].agent = a;
                    allocations[avail[j_max]].recalculate_utility();
                    avail.erase(avail.begin() + j_max);
                }
                else {
                    return false;
                }
            }
            return true;
        }

        std::vector<I> get_items() {
            std::vector<I> s;
            s.reserve(allocations.size());

            for (auto a : allocations) {
                s.push_back(a.item);
            }

            return s;
        }

        std::vector<A> get_agents() {
            std::vector<A> s;
            s.reserve(allocations.size());

            for (auto a : allocations) {
                s.push_back(a.agent);
            }

            return s;
        }

        // Will fail if we are not in equilibrium!
        bool reassign(double epsilon = 1e-5) {
            auto new_allocations = std::vector<Allocation<A, I>>(allocations);

            std::vector<std::vector<size_t>> favourites(allocations.size());

            struct FavReg {
                size_t i;
                double u;
            };

            for (size_t i = 0; i < allocations.size(); i++) {
                auto a = allocations[i];
                double u_max = 0.f;
                std::vector<FavReg> i_max{};
                for (size_t j = 0; j < allocations.size(); j++) {
                    // Determine whether
                    Allocation<A, I>& other = allocations[j];
                    if (a.agent.income() > other.price) {
                        double u = a.agent.utility(other.price, other.quality());
                        if (u > u_max) {
                            if (u > u_max + epsilon) {
                                i_max.clear();
                            }
                            u_max = u;
                            i_max.erase(std::remove_if(i_max.begin(), i_max.end(),
                                [u_max, epsilon](FavReg f) {
                                    return f.u + epsilon < u_max;
                                }),
                                i_max.end());
                            i_max.push_back({ .i = j, .u = u });
                        }
                        else if (u + epsilon > u_max) {
                            i_max.push_back({ .i = j, .u = u });
                        }
                    }
                }

                for (auto x : i_max) {
                    favourites[x.i].push_back(i);
                }
            }

            for (auto& v : favourites) {
                if (v.empty())
                    return false;
            }

            std::vector<std::vector<size_t>> favclone(favourites);

            std::vector<size_t> allocs(favourites.size());

            size_t allocated = 0;
            size_t min_alloc = 1;
            while (allocated < allocations.size()) {
                size_t new_allocated = 0;
                for (size_t j = 0; j < favourites.size(); ++j) {
                    std::vector<size_t>& favs = favourites[j];
                    if (!favs.empty() && (favs.size() <= 1 || min_alloc > 1)) {
                        min_alloc = 1;
                        // Agent j prefers allocation i.
                        size_t i = favs[0];
                        new_allocations[j].agent = allocations[i].agent;
                        new_allocations[j].recalculate_utility();
                        allocs[j] = i;
                        allocated += 1;
                        new_allocated += 1;
                        
                        for (size_t w = 0; w < favourites.size(); ++w) {
                            if (w != j) {
                                std::vector<size_t>& other = favourites[w];
                                size_t s = other.size();
                                auto remove = other.erase(std::remove(other.begin(), other.end(), i), other.end());
                                if (s != other.size() && other.size() == 0) {
                                    // Another agent had this as their most preferred, so the allocation must be invalid!!
                                    return false;
                                }
                            }
                        }
                        favs.clear();
                    }
                }

                // Now check for uniqueness.
                
                for (size_t j = 0; j < favourites.size(); ++j) {
                    std::vector<size_t>& favs = favourites[j];
                    for (size_t f : favs) {
                        bool unique = true;
                        for (size_t w = 0; w < favourites.size(); ++w) {
                            if (w != j) {
                                std::vector<size_t>& other = favourites[w];
                                if (std::find(other.begin(), other.end(), f) != other.end()) {
                                    unique = false;
                                    break;
                                }
                            }
                        }

                        if (unique) {
                            new_allocations[j].agent = allocations[f].agent;
                            new_allocations[j].recalculate_utility();
                            favs.clear();
                            ++new_allocated;
                        }
                    }
                }

                if (new_allocated == 0)
                    ++min_alloc;
            }
            

            allocations = new_allocations;
            return true;
        }

        ///
        /// @brief Verifies that the current allocation is a valid solution.
        ///
        ///        Checks that no agent would prefer any other agent's allocation at the given prices.
        ///        Ensures that each agent's utility is correct and that no agent can improve their utility
        ///        by switching to another allocation they can afford.
        ///
        /// @param epsilon Numerical tolerance for comparisons (default is 1e-6)
        ///
        /// @return True if the current allocation is valid; false otherwise.
        ///
        bool verify_solution(const double epsilon = 1e-5) const {
            bool valid = true;
            for (size_t i = 0; i < allocations.size(); ++i) {
                double u = allocations[i].agent.utility(allocations[i].price, allocations[i].item.quality());
                if (u != allocations[i].utility) {
                    std::cout << "Agent " << i << " has utility mismatch!" << std::endl;
                    return false;
                }

                for (size_t j = 0; j < allocations.size(); ++j) {
                    if (i != j) {
                        if (allocations[j].agent.item_id() == allocations[i].agent.item_id()) {
                            std::cout << "Agent " << i << " has the same item_id as " << j << "; item_id= " << allocations[j].agent.item_id() << std::endl;
                            valid = false;
                        }
                        // Compute the utility agent i would get from allocation j
                        const double u_alt = allocations[i].agent.utility(allocations[j].price, allocations[j].item.quality());
                        if (u_alt > u + epsilon) {
                            std::cout << "Agent " << i << " prefers allocation " << j << "; " << u_alt << ">" << u << std::endl;
                            valid = false;
                        }
                    }
                }
            }
            return valid;
        }

        ///
        /// @brief Draws the current allocations using the provided RenderState.
        ///
        /// @param render_state Pointer to a RenderState object for visualization.
        ///
        /// @return True if drawing was successful; false if terminated.
        ///
        RenderCommand draw(RenderState<A, I>* render_state) {
            return render_state->draw_allocations(this->allocations, -1);
        }

        ///
        /// @brief Performs a linear regression of price on quality.
        ///
        ///        Calculates the best-fit line of the form price = a + b * quality
        ///        using least squares regression, and outputs the regression coefficients
        ///        and the coefficient of determination (R^2).
        ///
        void regress_price_on_quality() {
            // Ensure there are enough data points
            if (allocations.size() < 2) {
                std::cerr << "Not enough data points to perform regression." << std::endl;
                return;
            }

            size_t n = allocations.size();
            double sum_x = 0.0f;   // Sum of qualities
            double sum_y = 0.0f;   // Sum of prices
            double sum_xx = 0.0f;  // Sum of qualities squared
            double sum_xy = 0.0f;  // Sum of quality * price

            // Calculate sums required for regression coefficients
            for (const auto& alloc : allocations) {
                double x = alloc.quality();
                double y = alloc.price;
                sum_x += x;
                sum_y += y;
                sum_xx += x * x;
                sum_xy += x * y;
            }

            double x_bar = sum_x / n;
            double y_bar = sum_y / n;

            double Sxy = sum_xy - n * x_bar * y_bar;
            double Sxx = sum_xx - n * x_bar * x_bar;

            if (Sxx == 0.0f) {
                std::cerr << "Cannot compute regression coefficients; division by zero." << std::endl;
                return;
            }

            double b = Sxy / Sxx;          // Slope
            double a = y_bar - b * x_bar;  // Intercept

            std::cout << "Regression result: price = " << a << " + " << b << " * quality" << std::endl;

            // Calculate the coefficient of determination (R^2)
            double ss_tot = 0.0f;
            double ss_res = 0.0f;
            for (const auto& alloc : allocations) {
                double x = alloc.quality();
                double y = alloc.price;
                double y_pred = a + b * x;
                ss_tot += (y - y_bar) * (y - y_bar);
                ss_res += (y - y_pred) * (y - y_pred);
            }

            double r_squared = 1 - (ss_res / ss_tot);
            std::cout << "Coefficient of determination (R^2): " << r_squared << std::endl;
        }
    };
}

#endif //SOLVER_H