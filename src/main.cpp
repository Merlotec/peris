#include <iostream>
#include "world.h"
#include "peris.h"
#include <random>
#include <algorithm>
#include <chrono>
#include <assert.h>

World create_world(const size_t school_count, const size_t house_count) {

    // Number of school places must be exact, so school_count divide into house_count perfectly.
    assert(house_count % school_count == 0);

    const ssize_t school_capacity = house_count / school_count;

    // Initialize with capacities
    std::vector<School> schools;
    schools.reserve(school_count);

    std::vector<House> houses;
    houses.reserve(house_count);

    std::vector<Household> households;
    households.reserve(house_count);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<double> school_quality_distribution(0.8, 0.28);
    std::normal_distribution<double> ability_distribution(0, 1.0);
    std::normal_distribution<double> aspiration_distribution(0.59, 0.12);

    constexpr double mean_household_inc = 100.0;

    double cv = 0.4; // Adjusts variance and skewness of distribution
    double variance = (cv * mean_household_inc) * (cv * mean_household_inc);
    double standard_deviation = std::sqrt(variance);

    // Calculate the parameters m and s for the underlying normal distribution
    double sigma_squared = std::log((variance / (mean_household_inc * mean_household_inc)) + 1.0);
    double sigma = std::sqrt(sigma_squared);
    double mu = std::log(mean_household_inc) - (sigma_squared / 2.0);

    // Construct the lognormal distribution with parameters mu and sigma
    std::lognormal_distribution<double> household_income_distribution(mu, sigma);

    std::uniform_real_distribution<double> location_axis_distribution(-1.0, 1.0);

    std::cout << "Created distributions" << std::endl;
    for (int i = 0; i < school_count; ++i) {
        // Sample
        const double quality = school_quality_distribution(gen);
        const double x = location_axis_distribution(gen);
        const double y = location_axis_distribution(gen);

        const School school = { .capacity = school_capacity, .x = x, .y = y, .quality = quality, .attainment = -1.f, .num_pupils = 0 };

        schools.push_back(school);
    }

    std::sort(schools.begin(), schools.end(), [](const School& a, const School& b) {
        return a.quality < b.quality;
        });

    std::cout << "Created schools" << std::endl;

    for (int i = 0; i < house_count; ++i) {
        // Sample
        const double x = location_axis_distribution(gen);
        const double y = location_axis_distribution(gen);

        const House house = { .x = x, .y = y, .school = -1 };
        houses.push_back(house);
    }

    std::cout << "Created houses" << std::endl;

    std::vector<House> allocated_houses;
    allocated_houses.reserve(houses.size());

    for (size_t sc = 0; sc < schools.size(); ++sc) {
        School& school = schools[sc];
        std::sort(houses.begin(), houses.end(), [school](const House& a, const House& b) {
            const double a_dis = (a.x - school.x) * (a.x - school.x) + (a.y - school.y) * (a.y - school.y);
            const double b_dis = (b.x - school.x) * (b.x - school.x) + (b.y - school.y) * (b.y - school.y);

            return a_dis > b_dis;
            });

        size_t n = std::min((size_t)school.capacity, houses.size());
        for (size_t i = 0; i < n; ++i) {
            House h = houses.back();
            h.set_school(sc, schools[sc].quality);
            allocated_houses.push_back(h);
            houses.pop_back();
        }
    }
    std::cout << "Allocated houses: " << allocated_houses.size() << std::endl;
    assert(allocated_houses.size() == house_count);

    for (int i = 0; i < house_count; ++i) {
        const double income = household_income_distribution(gen);
        const double ability = ability_distribution(gen);
        const double aspiration = aspiration_distribution(gen);// > 0.7 ? 0.7 : 0.8;

        const Household household = { .id = i, .inc = income, .ability = ability, .aspiration = aspiration, .school = -1, .house = -1 };
        households.push_back(household);
    }

    std::cout << "Created households: " << households.size() << std::endl;

    return World(households, schools, allocated_houses);
}

int main() {
    constexpr size_t school_count = 200;
    constexpr size_t house_count = 200;

    std::cout << "Creating with " << school_count << " schools" << " and " << house_count << " houses" << std::endl;


    World world = create_world(school_count, house_count);
    while (!world.validate()) {
        world = create_world(school_count, house_count);

        std::cout << "Created world" << std::endl;
    }




    // Setup render window to draw visuals.
    peris::RenderState<Household, House> render_state{};

    //auto solver = world.solver();

    std::vector<peris::Allocation<Household, House>> allocations{};
    peris::Solver<Household, House>::align_right(world.households, world.houses, allocations, 0, &render_state);
    auto solver = peris::Solver(allocations);

    solver.draw(&render_state);

    //auto pres = solver.solve(&render_state);

    //std::cout << "Solver finished with code " << pres << std::endl;

    if (peris::Solver<Household, House>::verify_solution(allocations)) {
        std::cout << "Verification successful" << std::endl;
    }
    else {
        std::cout << "Verification failed!" << std::endl;
    }

    // Run regression to calculate gsp:

    solver.regress_price_on_quality();

    double delta = 0.05;

    peris::RenderCommand cmd;
    while ((cmd = solver.draw(&render_state)) >= 0) {
        if (cmd == peris::tick) {
            bool valid = false;
            size_t i;
            size_t mv;
            std::vector<peris::Allocation<Household, House>> prev = solver.get_allocations();

            for (i = 0; i < 100; ++i) {
                mv = solver.equilibriate_demand(delta);
                if (mv == 0) {
                    valid = true;
                    break;
                }
                delta *= 0.95;
            }
            delta = std::max(0.05, solver.max_difference(prev) * 0.1);

            
            if (valid) {
                std::cout << "Valid allocation (" << i << ") ticks" << std::endl;
            }
            else {
                std::cout << "Allocation not yet valid (" << mv << ") final movements..." << std::endl;

            }
            if (peris::Solver<Household, House>::verify_solution(prev)) {
                std::cout << "Solution verified!" << std::endl;
            }
            // if (solver.reassign()) {
            //     std::cout << "Successfully reassigned agents!" << std::endl;
            //     if (solver.verify_solution()) {
            //         std::cout << "Solution verified!" << std::endl;
            //     }
            // } else {
            //     std::cout << "Failed to reassign agents!" << std::endl;
            // }
        }
    }
    return 0;
}