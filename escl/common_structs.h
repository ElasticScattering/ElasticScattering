typedef struct
{
	int particle_row_count;
	int particle_count;
	int impurity_count;

	double region_size;
	double particle_max_lifetime;
	double particle_speed;      // v
	double particle_mass;       // m
	double impurity_radius;     // r
	double impurity_radius_sq;  // r^2
	double tau;

	double alpha;
	double phi;
	double magnetic_field;      // B
	double angular_speed;       // w
} SimulationParameters;

typedef struct
{
	double2 position;
	double radius;
} Circle;