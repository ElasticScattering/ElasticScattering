#pragma once

typedef struct ScatteringParameters
{
    int mode;
    int dim;
    int impurity_count;
    int integrand_steps;

    int is_clockwise;
    int is_incoherent;
    int is_diag_regions;

    unsigned int impurity_seed;

    double region_size;
    double region_extends;
    double impurity_density;
    int impurity_grid_dim;
    
    double particle_speed;
    double impurity_radius;
    double tau;
    double temperature;

    double alpha;
    double phi;
    double magnetic_field;
    double angular_speed;

    friend bool operator==(const ScatteringParameters& lhs, const ScatteringParameters& rhs) {
        return (lhs.mode == rhs.mode && lhs.impurity_seed == rhs.impurity_seed &&
            lhs.region_size == rhs.region_size && lhs.dim == rhs.dim &&
            lhs.particle_speed == rhs.particle_speed &&
            lhs.impurity_count == rhs.impurity_count && lhs.impurity_radius == rhs.impurity_radius &&
            lhs.alpha == rhs.alpha && lhs.phi == rhs.phi && lhs.temperature == rhs.temperature &&
            lhs.magnetic_field == rhs.magnetic_field && lhs.tau == rhs.tau &&
            lhs.integrand_steps == rhs.integrand_steps && lhs.is_clockwise == rhs.is_clockwise &&
            lhs.region_extends == rhs.region_extends && lhs.is_clockwise == rhs.is_clockwise &&
            lhs.is_diag_regions == rhs.is_diag_regions && lhs.is_incoherent == rhs.is_incoherent);
    }

    friend bool operator!=(const ScatteringParameters& lhs, const ScatteringParameters& rhs) {
        return !(lhs == rhs);
    }
} ScatteringParameters;
