/*
   Taken from:
   The Computer Language Benchmarks Game
   https://salsa.debian.org/benchmarksgame-team/benchmarksgame/

   An implementation pretty much from scratch, with inspiration from the Rust
   version, which used the idea of saving some of the ingredients of the
   compution in an array instead of recomputing them.
   
   contributed by cvergu
   slightly modified by bmmeijers
*/

#define _USE_MATH_DEFINES // https://docs.microsoft.com/en-us/cpp/c-runtime-library/math-constants?view=msvc-160

#include <cmath>
#include <iostream>
#include <fstream>
#include <tuple>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>


// these values are constant and not allowed to be changed
const double SOLAR_MASS = 4 * M_PI * M_PI;
const double DAYS_PER_YEAR = 365.24;
const unsigned int BODIES_COUNT = 5;
const char *FILE_NAME = "locations.csv";
const unsigned int ITERATIONS_BEFORE_SAVING_TO_FILE = 1e6;


class vector3d {
public:
    double x, y, z;

    double norm() const noexcept {
        return x * x + y * y + z * z;
    }

    double magnitude(double dt) const noexcept {
        double sum = norm();
        return dt / (sum * sqrt(sum));
    }
};

vector3d operator+(vector3d v1, vector3d v2) {
    return vector3d{
            v1.x + v2.x, v1.y + v2.y, v1.z + v2.z
    };
}

vector3d operator-(vector3d v1, vector3d v2) {
    return vector3d{
            v1.x - v2.x, v1.y - v2.y, v1.z - v2.z
    };
}

vector3d &operator+=(vector3d &v1, vector3d v2) {
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;

    return v1;
}

vector3d &operator-=(vector3d &v1, vector3d v2) {
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;

    return v1;
}

vector3d &operator*=(vector3d &v, double mag) {
    v.x *= mag;
    v.y *= mag;
    v.z *= mag;

    return v;
}

vector3d operator*(vector3d v, double mag) {
    return vector3d{
            v.x * mag, v.y * mag, v.z * mag
    };
}

vector3d operator/(vector3d v, double mag) {
    return vector3d{
            v.x / mag, v.y / mag, v.z / mag
    };
}


class body {
public:
    std::string name;
    vector3d position;
    vector3d velocity;
    double mass;
};


void advance(body state[BODIES_COUNT], double dt) {
    /*
     * We precompute the quantity (r_i - r_j)
     */
    // 2D array (to hold: BODIES_COUNT x BODIES_COUNT elements)
    vector3d rij[BODIES_COUNT][BODIES_COUNT];

    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            rij[i][j] = state[i].position - state[j].position;
        }
    }

    double magnitudes[BODIES_COUNT][BODIES_COUNT];

    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            magnitudes[i][j] = rij[i][j].magnitude(dt);
        }
    }

    /*
     * Compute the new speed using v_i = a_i dt, where
     * a_i = \sum_{j \neq i} m_j (r_i - r_j)/|r_i - r_j|^3
     */
    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            vector3d dist = rij[i][j];
            double mag = magnitudes[i][j];
            state[i].velocity -= dist * (state[j].mass * mag);
            state[j].velocity += dist * (state[i].mass * mag);
        }
    }

    /*
     * Compute the new positions
     */
    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        state[i].position += state[i].velocity * dt;
    }
}

void offset_momentum(body state[BODIES_COUNT]) {
    vector3d &sun_velocity = state[0].velocity;

    for (unsigned int i = 1; i < BODIES_COUNT; ++i) {
        sun_velocity -= state[1].velocity * state[1].mass / SOLAR_MASS;
    }
}

double energy(const body state[BODIES_COUNT]) {
    double energy = 0;

    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        const body &body1 = state[i];
        energy += 0.5 * body1.mass * body1.velocity.norm();
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            const body &body2 = state[j];
            vector3d r12 = body1.position - body2.position;
            energy -= body1.mass * body2.mass / sqrt(r12.norm());
        }
    }

    return energy;
}

void create_csv(std::ofstream &csv) {
    csv << "# name of body; position x; position y; position z\n";
}

void append_csv(std::ofstream &csv, std::vector<std::tuple<std::string, vector3d>> &locations) {
    std::ostringstream csv_steam;
    for (auto location: locations) {
        auto position = get<1>(location);
        csv_steam << get<0>(location) << ";" << position.x << ";" << position.y << ";" << position.z << "\n";
    }
    std::string var = csv_steam.str();
    csv << var;
}

void append_csv(std::ofstream &csv, body state[BODIES_COUNT]) {
    std::ostringstream csv_steam;
    for (int i = 0; i < BODIES_COUNT; ++i) {
        csv_steam << state[i].name << ";" << state[i].position.x << ";" << state[i].position.y << ";"
                  << state[i].position.z << "\n";
    }
    std::string var = csv_steam.str();
    csv << var;
}


body state[] = {
        // Sun
        {
                .name = "sun",
                .position = {
                        0,
                        0,
                        0
                },
                .velocity = {
                        0,
                        0,
                        0
                },
                .mass = SOLAR_MASS
        },
        // Jupiter
        {
                .name = "jupiter",
                .position = {
                        4.84143144246472090e+00,
                        -1.16032004402742839e+00,
                        -1.03622044471123109e-01
                },
                .velocity = {
                        1.66007664274403694e-03 * DAYS_PER_YEAR,
                        7.69901118419740425e-03 * DAYS_PER_YEAR,
                        -6.90460016972063023e-05 * DAYS_PER_YEAR
                },
                .mass = 9.54791938424326609e-04 * SOLAR_MASS
        },
        // Saturn
        {
                .name = "saturn",
                .position = {
                        8.34336671824457987e+00,
                        4.12479856412430479e+00,
                        -4.03523417114321381e-01
                },
                .velocity = {
                        -2.76742510726862411e-03 * DAYS_PER_YEAR,
                        4.99852801234917238e-03 * DAYS_PER_YEAR,
                        2.30417297573763929e-05 * DAYS_PER_YEAR
                },
                .mass = 2.85885980666130812e-04 * SOLAR_MASS
        },
        // Uranus
        {
                .name = "uranus",
                .position = {
                        1.28943695621391310e+01,
                        -1.51111514016986312e+01,
                        -2.23307578892655734e-01
                },
                .velocity = {
                        2.96460137564761618e-03 * DAYS_PER_YEAR,
                        2.37847173959480950e-03 * DAYS_PER_YEAR,
                        -2.96589568540237556e-05 * DAYS_PER_YEAR
                },
                .mass = 4.36624404335156298e-05 * SOLAR_MASS
        },
        // Neptune
        {
                .name = "neptune",
                .position = {
                        1.53796971148509165e+01,
                        -2.59193146099879641e+01,
                        1.79258772950371181e-01
                },
                .velocity = {
                        2.68067772490389322e-03 * DAYS_PER_YEAR,
                        1.62824170038242295e-03 * DAYS_PER_YEAR,
                        -9.51592254519715870e-05 * DAYS_PER_YEAR
                },
                .mass = 5.15138902046611451e-05 * SOLAR_MASS
        }
};

void simulate_bodies(const unsigned int n, const char *filename) {
    offset_momentum(state);
    std::cout << energy(state) << std::endl;
    std::ofstream csv;
    csv.open(filename, std::ofstream::trunc);
    create_csv(csv);
    append_csv(csv, state);
    std::vector<std::tuple<std::string, vector3d>> locations;
    for (int i = 0; i < n; ++i) {
        if (i % ITERATIONS_BEFORE_SAVING_TO_FILE == 0) {
            append_csv(csv, locations);
            locations.clear();
        }
        advance(state, 0.01);
        for (auto &body: state) {
            auto body_location = std::make_tuple(body.name, body.position);
            locations.push_back(body_location);
        }
    }
    append_csv(csv, locations);
    csv.close();
    std::cout << energy(state) << std::endl;
}

void time_function() {
    int iterations[] = {5000};
    for (auto iteration: iterations) {
        std::string filename = "locations-" + std::to_string(iteration) + ".csv";
        auto t1 = std::chrono::high_resolution_clock::now();
        simulate_bodies(iteration, filename.c_str());
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> ms_double = t2 - t1;
        std::cout << iteration << " iterations, " << ms_double.count() / 1000 << "s\n";
    }
}

int main(int argc, char **argv) {
    time_function();
    return EXIT_SUCCESS;
}

//int main(int argc, char **argv) {
//    if (argc != 2) {
//        std::cout << "This is " << argv[0] << std::endl;
//        std::cout << "Call this program with an integer as program argument" << std::endl;
//        std::cout << "(to set the number of iterations for the n-body simulation)." << std::endl;
//        return EXIT_FAILURE;
//    } else {
//        const unsigned int n = atoi(argv[1]);
//        simulate_bodies(n, FILE_NAME);
//        return EXIT_SUCCESS;
//    }
//}
