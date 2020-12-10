#ifndef REGION_ORBIT_DETAILS
    double smod(const double a, const double b)
    {
        return a - b * floor(a / b);
    }

    double GetBoundTime(const double phi, const double alpha, const double w, const bool is_incoherent, const bool is_diag_region, const bool is_electron, const bool is_future)
    {
        if (!is_incoherent) return DBL_MAX;

        const double v = is_diag_region ? (phi + alpha - PI / 4.0) : (phi + alpha);
        const double remaining = smod(v, PI * 0.5);

        // Oud
        //double dphi = ((!is_electron && is_future) || (is_electron && !is_future)) ? remaining : (2.0 * alpha - remaining);
        double dphi = (is_electron != is_future) ? remaining : (2.0 * alpha - remaining);
        return dphi / w;
    }

    double GetBoundAngle(const double phi, const double alpha, const bool clockwise)
    {
        const int multiple = (int)(phi + alpha / (PI / 2.0));
        const double bound1 = multiple * PI / 2.0 - alpha;
        const double bound2 = multiple * PI / 2.0 + alpha;

        double dangle1 = GetCrossAngle(phi, bound1, clockwise);
        double dangle2 = GetCrossAngle(phi, bound2, clockwise);

        return (dangle1 < dangle2) ? bound1 : bound2;
    }

    double2 GetCyclotronOrbitCenter(const double2 p, const double2 velocity, const double radius, const double vf, const bool is_electron)
    {
        double2 shift = { velocity.y, -velocity.x };
        shift = shift * radius / vf;

        return is_electron ? (p - shift) : (p + shift);
    }

    bool CirclesCross(const Orbit *orbit, const double2 p2, const double r2)
    {
        const double2 q = orbit->center - p2;
        const double dist_squared = q.x * q.x + q.y * q.y;

        const double r_add = orbit->radius + r2;
        const double r_min = orbit->radius - r2;

        //Oud
        //return (dist_squared >= r_add * r_add || dist_squared <= r_min * r_min) ? false : true;
        return !(dist_squared >= r_add * r_add || dist_squared <= r_min * r_min);
    }

    double4 GetCrossPoints(const Orbit *orbit, const double2 p2, const double r2)
    {
        const double2 q = orbit->center - p2;

        const double dist_squared = dot(q, q);
        const double dist = sqrt(dist_squared);
        const double xs = (dist_squared + orbit->radius_squared - r2 * r2) / (2.0 * dist);
        const double ys = sqrt(orbit->radius_squared - xs * xs);

        const double2 u = (p2 - orbit->center) / dist;

        double4 points = {
            orbit->center.x + u.x * xs + u.y * ys,
            orbit->center.y + u.y * xs + -u.x * ys,

            orbit->center.x + u.x * xs + -u.y * ys,
            orbit->center.y + u.y * xs + u.x * ys
        };

        return points;
    }

    double GetPhi(const double2 pos, const Orbit *orbit)
    {
        double p = (pos.x - orbit->center.x) / orbit->radius;

        p = (p > 1.0) ? 1.0 : p;
        p = (p < -1.0) ? -1.0 : p;

        double phi = acos(p);
        phi = (pos.y < orbit->center.y) ? PI2 - phi : phi;

        return phi;
    }

    double GetCrossAngle(const double p, const double q, const bool clockwise)
    {
        const double g = clockwise ? (p - q) : (q - p);
        return smod(g, PI2);
    }

    double GetFirstCrossTime(const double2 pos, const Orbit *orbit, const double2 ip, const double ir, const double w)
    {
        const double4 cross_points = GetCrossPoints(orbit, ip, ir);

        const double2 p1 = { cross_points.x, cross_points.y };
        const double2 p2 = { cross_points.z, cross_points.w };

        const double phi0 = GetPhi(pos, orbit);
        const double phi1 = GetPhi(p1, orbit);
        const double phi2 = GetPhi(p2, orbit);

        const double t1 = GetCrossAngle(phi0, phi1, orbit->clockwise) / w;
        const double t2 = GetCrossAngle(phi0, phi2, orbit->clockwise) / w;

        return min(t1, t2);
    }
#endif

#ifndef REGION_IMPURITY_GRID_DETAILS

	double GetAngle(double2 pos, const Orbit *orbit) {
		if (abs(pos.y - orbit->center.y) > EPSILON) {
			double sign = (pos.x < orbit->center.x) ? -1.0 : 1.0;
			double angle = sign * asin((pos.y - orbit->center.y) / orbit->radius);
			return fmod(angle, PI2);
		}
		else {
			return (pos.x > orbit->center.x) ? 0 : PI;
		}
	}

	bool PointInSegment(double point, double l0, double l1)
	{
		return (point > l0 && point < l1) || (point < l0 && point > l1);
	}

	bool GetFirstBoundaryIntersect(const double2 p1, const double2 p2, const Orbit *orbit, const double L, const double start_phi, Intersection *intersection) {
		double2 u = (p2 - p1) / L;
		double projection_distance = u.x * (orbit->center.x - p1.x) + u.y * (orbit->center.y - p1.y);
		double2 proj = p1 + (u * projection_distance);
		double proj_circle_distance_sq = pow(proj.x - orbit->center.x, 2) + pow(proj.y - orbit->center.y, 2);
	
		// Test if line enters circle.
		if (proj_circle_distance_sq >= orbit->radius_squared) {
			return false; 
		}
	
		// Calculate segment to the edge of the circle.
		double2 to_edge = u * sqrt(orbit->radius_squared - proj_circle_distance_sq);

		double2 i1 = proj + to_edge;
		double2 i2 = proj - to_edge;

		// Determine if the intersection happen on the line segment.
		bool horizontal_line = abs(u.x) > abs(u.y);

		Intersection in1;
		in1.position = i1;
		in1.incident_angle = GetAngle(i1, orbit);
		in1.dphi = GetCrossAngle(start_phi, in1.incident_angle, orbit->clockwise);
	
		Intersection in2;
		in2.position = i2;
		in2.incident_angle = GetAngle(i2, orbit);
		in2.dphi = GetCrossAngle(start_phi, in2.incident_angle, orbit->clockwise);

		bool i1_valid = horizontal_line ? PointInSegment(i1.x, p1.x, p2.x) : PointInSegment(i1.y, p1.y, p2.y);
		bool i2_valid = horizontal_line ? PointInSegment(i2.x, p1.x, p2.x) : PointInSegment(i2.y, p1.y, p2.y);
	
		//@Refactor, kan dit simpeler omdat clockwise al in GetCrossAngle zit?
		if (i1_valid && i2_valid)
		{
			bool phi1_lower = in1.dphi < in2.dphi;
			if (orbit->clockwise) intersection = (phi1_lower) ? &in1 : &in2;
			else             intersection = (phi1_lower) ? &in2 : &in1;
		}
		else if (i1_valid) intersection = &in1;
		else               intersection = &in2;

		return true;
	}

	bool GetNextCell(const Orbit *orbit,
		const int current_cell,
		const double2 current_cell_lowleft, 
		const Intersection last_intersection,
		const double L,
		const int cells_per_row,
		int *next_cell, 
		Intersection* next_intersection)
	{
		double2 low_left = current_cell_lowleft;
		double2 low_right = low_left + double2(L, 0);
		double2 top_right = low_left + double2(L, L);
		double2 top_left = low_left + double2(0, L);

		Intersection left, right, up, down;

		double last_phi = last_intersection.dphi;
		bool up_hit = GetFirstBoundaryIntersect(top_left, top_right, orbit, L, last_phi, &up);
		bool down_hit = GetFirstBoundaryIntersect(low_left, low_right, orbit, L, last_phi, &down);
		bool right_hit = GetFirstBoundaryIntersect(low_right, top_right, orbit, L, last_phi, &right);
		bool left_hit = GetFirstBoundaryIntersect(low_left, top_left, orbit, L, last_phi, &left);

		if (!up_hit && !down_hit && !right_hit && !left_hit)
			return false;

		// Select [one] from all valid intersects.
		//	double dphi = GetCrossAngle(start_intersect.incident_angle, i.incident_angle, orbit->clockwise);
		int closest_cell_index = current_cell; // ???
		Intersection closest_intersection = last_intersection; //
	
		if (up_hit) {
			// Test if this intersection is closer than what we have
			double dphi = GetCrossAngle(last_intersection.incident_angle, up.incident_angle, orbit->clockwise);
			int next_cell_candidate = current_cell - cells_per_row;
			if (dphi < closest_intersection.dphi && next_cell_candidate >= 0) {
				closest_intersection.dphi = dphi;
				closest_cell_index = current_cell - cells_per_row;
			}
		}

		if (down_hit) {
			// Test if this intersection is closer than what we have
			double dphi = GetCrossAngle(last_intersection.incident_angle, up.incident_angle, orbit->clockwise);
			int next_cell_candidate = current_cell + cells_per_row; // @Todo, dit controleren..
			if (dphi < closest_intersection.dphi && next_cell_candidate >= 0) {
				closest_intersection.dphi = dphi;
				closest_cell_index = current_cell - cells_per_row;
			}
		}
		//etc ...

		if (closest_cell_index == current_cell) {
			return false;
		}

		// The intersection can not put us out of bounds.
		//@Todo: doe dit met cell_x, cell_y, tijdens de 4-test, of achteraf?

		return true;
	}

	////////////////////////////

	int get_cell_index(const v2 pos, const v2 range, const int cells_per_row)
	{
		return to_index(to_grid(pos.x, pos.y, range, cells_per_row), cells_per_row);
	}

	CellRange get_cell_id(const v2 pos, const v2 range, const int cells_per_row)
	{
		CellRange cell_range;

		cell_range.start = to_index(to_grid(pos.x, pos.y, range, cells_per_row), cells_per_row);
		cell_range.end = cell_range.start + 1;

		return cell_range;
	}

	////////////////////////


	int to_grid(const double x, const double2 range, const int cells_per_row)
	{
		return (int)((x - range.x) / (range.y - range.x) * (cells_per_row));
	}

	v2i to_grid(const double x, const double y, const double2 range, const int cells_per_row)
	{
		return {
			(int)((x - range.x) / (range.y - range.x) * (cells_per_row)),
			(int)((y - range.x) / (range.y - range.x) * (cells_per_row))
		};
	}

	v2 to_world(const int cell_index, const int cells_per_row, const double2 spawn_range)
	{

		double x = cell_index % cells_per_row;
		double y = floor(cell_index / cells_per_row);
		double2 low_left = (double2)(x, y) * (spawn_range.y - spawn_range.x) + spawn_range.x;

		return low_left;
	}

	bool within_bounds(v2i p, const int cells_per_row) {
		return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
	}

	int to_index(const int x, const int y, const int cells_per_row) {
		return y * cells_per_row + x;
	}

	int to_index(const v2i p, const int cells_per_row) {
		return to_index(p.x, p.y, cells_per_row);
	}

	double remap(double x, double s0, double s1, double t0, double t1)
	{
		return t0 + (x - s0) / (s1 - s0) * (t1 - t0);
	}
#endif

#ifndef REGION_LIFETIME_TRACING
	double SingleLifetime(const Particle* p, const Orbit *orbit, const double2 impurity, const double impurity_radius, const double angular_speed)
	{
		double2 d = p->starting_position - impurity;
		if ((impurity_radius * impurity_radius) > dot(d, d))
			return 0;

		if (CirclesCross(orbit, impurity, impurity_radius))
		{
			return GetFirstCrossTime(p->starting_position, orbit, impurity, impurity_radius, angular_speed);
		}

		return DBL_MAX;
	}

	double TraceOrbit(Particle* p, const Orbit *orbit, BUFFER_ARGS)
	{
		Intersection entry_point;
		entry_point.position = p->starting_position;
		entry_point.dphi = 0; // ????

		// Early exit als laatste intersectie een grotere hoek heeft dan bound_angle.

		double lifetime = DBL_MAX;
		bool hit = false;
		while (!hit) {
			int impurity_start = cell_indices[p->cell_index];
			int impurity_end = cell_indices[p->cell_index + 1]; // null?

			for (int i = impurity_start; i < impurity_end; i++) {
				double t = SingleLifetime(p, orbit, impurities[i], sp->impurity_radius, sp->angular_speed);

				lifetime = (t < lifetime) ? t : lifetime;
				hit = (lifetime < DBL_MAX);
			}

			if (!hit) {
				double2 cell_pos = to_world(p->cell_index, sp->cells_per_row, sp->impurity_spawn_range);

				int next_cell;
				Intersection next_intersection;
				bool success = GetNextCell(orbit, p->cell_index, cell_pos, entry_point, sp->cell_size, sp->cells_per_row, &next_cell, &next_intersection);

				if (!success) {
					// The End?
					break;
				}

				p->cell_index = next_cell;
				entry_point = next_intersection;
			}
		}

		return lifetime;
	}

	Orbit MakeOrbit(const double2 pos, const double phi, ScatteringParameters* sp)
	{
		const bool clockwise = sp->is_clockwise == 1;
		const bool incoherent = sp->is_incoherent == 1;
		const bool diag_regions = sp->is_diag_regions == 1;

		const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, incoherent, diag_regions, clockwise, false);
		const double bound_angle = GetBoundAngle(phi, sp->alpha, clockwise);
		const double bound_phi = GetCrossAngle(phi, bound_angle, clockwise);

		const v2 vel = (double2)(cos(phi), sin(phi)) * sp->particle_speed;
		const double orbit_radius = sp->particle_speed / sp->angular_speed;
		const double2 center = GetCyclotronOrbitCenter(pos, vel, orbit_radius, sp->particle_speed, clockwise);

		Orbit orbit(center, orbit_radius, clockwise, bound_time, bound_phi);

		return orbit;
	}
#endif
