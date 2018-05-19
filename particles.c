#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "particles.h"
#include "consts.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
      _a > _b ? _a : _b; })


void print_particle(particle p) {
    printf("Particle index %d:\n", p.index);
    printf("Coordinates %lf %lf %lf:\n", p.x_c, p.y_c, p.z_c);
    printf("Velocities %lf %lf %lf:\n", p.x_v, p.y_v, p.z_v);

    printf("Accelerations %lf %lf %lf:\n", p.acc_x, p.acc_y, p.acc_z);
    printf(
        "New Accelerations %lf %lf %lf:\n\n",
        p.new_acc_x, p.new_acc_y, p.new_acc_z
    );
}


void particle_to_line(particle p) {
    printf(
        "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n",
        p.x_c, p.y_c, p.z_c, p.x_v, p.y_v, p.z_v
    );
}


void print_p_set(particle_set p_set) {
    printf("Particle set with number %d\n\n", p_set.number);
    for(int i = 0; i < p_set.number; i ++) {
        print_particle(p_set.particles[i]);
    }
}

void print_to_lines_p_set(particle_set p_set) {
    for(int i = 0; i < p_set.number; i ++) {
        particle_to_line(p_set.particles[i]);
    }
}


void free_particle_set(particle_set p_set) {
    free(p_set.particles);
}


particle_set new_particle_set(int number) {
    particle_set p_set;
    p_set.number = number;
    p_set.particles = malloc(sizeof(particle) * number);
    return p_set;
}

double safe_value(double value) {
    if (fabs(value) < MIN_NUMB)
        return MIN_NUMB;
    return value;
}


double power(double value, int times) {
    if (times == 0) {
        return 1;
    }
    return power(value, times - 1) * value;
}


double norm_distance(particle *particle_a, particle *particle_b) {
    double x = fabs(particle_a->x_c - particle_b->x_c);
    double y = fabs(particle_a->y_c - particle_b->y_c);
    double z = fabs(particle_a->z_c - particle_b->z_c);
    return max(MIN_NUMB, sqrt(x * x + y * y + z * z));
}


double calc_potential(particle *p1, particle *p2, particle *p3) {
    double rij = norm_distance(p1, p2);
    double rik = norm_distance(p1, p3);
    double rkj = norm_distance(p2, p3);

    double srij = rij * rij;
    double srik = rik * rik;
    double srkj = rkj * rkj;

    double first_part = (1 / power((rij * rik * rkj), 3));
    double second_part_l = (
        3.0 * (-srij + srik + srkj)
        * (srij - srik + srkj) * (srij + srik - srkj)
    );
    double second_part_m = 8.0 * power((rij * rik * rkj), 5);
    return first_part + second_part_l / second_part_m;
}


void update_acceleration(particle *p1, particle *p2, particle *p3) {
    // calculating x
    double old_x = p1->x_c;
    //printf("MIN_NUMB %E\n", MIN_NUMB);
    //printf("h_bb_x %E\n", SQRT_E_CONST * p1->x_c);
    double h_x = safe_value(SQRT_E_CONST * p1->x_c);

    double up_x = old_x + h_x;
    double down_x = old_x - h_x;
    //printf("h_x %E\n", h_x);
    //printf("up_x %E\n", up_x);
    //printf("down_x %E\n", down_x);

    p1->x_c = up_x;
    double potential_x1 = calc_potential(p1, p2, p3);

    p1->x_c = down_x;
    double potential_x2 = calc_potential(p1, p2, p3);

    p1->x_c = old_x;

    //printf("Before update %E\n", p1->new_acc_x);
    //printf("potential_x1 %E\n", potential_x1);
    //printf("potential_x2 %E\n", potential_x2);
    //printf("Diff %E\n", potential_x1 - potential_x2);
    p1->new_acc_x -= (potential_x1 - potential_x2) / (up_x - down_x);
    //printf("After update %E\n", p1->new_acc_x);


    // calculating y
    double old_y = p1->y_c;
    double h_y = safe_value(SQRT_E_CONST * p1->y_c);

    double up_y = old_y + h_y;
    double down_y = old_y - h_y;

    p1->y_c = up_y;
    double potential_y1 = calc_potential(p1, p2, p3);

    p1->y_c = down_y;
    double potential_y2 = calc_potential(p1, p2, p3);

    p1->y_c = old_y;
    p1->new_acc_y -= (potential_y1 - potential_y2) / (up_y - down_y);


    // calculating z
    double old_z = p1->z_c;
    double h_z = safe_value(SQRT_E_CONST * p1->z_c);

    double up_z = old_z + h_z;
    double down_z = old_z - h_z;

    p1->z_c = up_z;
    double potential_z1 = calc_potential(p1, p2, p3);

    p1->z_c = down_z;
    double potential_z2 = calc_potential(p1, p2, p3);

    p1->z_c = old_z;
    p1->new_acc_z -= (potential_z1 - potential_z2) / (up_z - down_z);
}

void update_distance(particle *p, double timedelta) {
    p->x_c += (
        p->x_v * timedelta
        // + (1.0 / 2) * p->acc_x * timedelta * timedelta
        + (1.0 / 2) * p->acc_x * timedelta
    );

    p->y_c += (
        p->y_v * timedelta
        // + ( 1.0 / 2) * p->acc_y * timedelta * timedelta
        + ( 1.0 / 2) * p->acc_y * timedelta
    );

    p->z_c += (
        p->z_v * timedelta
        // + (1.0 / 2) * p->acc_z * timedelta * timedelta
        + (1.0 / 2) * p->acc_z * timedelta
    );
}

void reset_acceleration(particle *self) {
    self->acc_x = self->new_acc_x;
    self->acc_y = self->new_acc_y;
    self->acc_z = self->new_acc_z;

    self->new_acc_x = 0.0;
    self->new_acc_y = 0.0;
    self->new_acc_z = 0.0;
}

void velocity_update(particle *p, double timedelta) {
    p->x_v += (1.0 / 2) * (p->acc_x + p->new_acc_x) * timedelta;
    p->y_v += (1.0 / 2) * (p->acc_y + p->new_acc_y) * timedelta;
    p->z_v += (1.0 / 2) * (p->acc_z + p->new_acc_z) * timedelta;
    reset_acceleration(p);
}

void calc_accelerations(particle_set p_set) {
    for (int i = 0; i < p_set.number; i ++) {
        for (int j = i + 1; j < p_set.number; j ++) {
            for (int k = j + 1; k < p_set.number; k ++) {
                particle *part_a = &(p_set.particles[i]);
                particle *part_b = &(p_set.particles[j]);
                particle *part_c = &(p_set.particles[k]);

                update_acceleration(part_a, part_b, part_c);
                update_acceleration(part_b, part_a, part_c);
                update_acceleration(part_c, part_a, part_b);
            }
        }
    }
}

void update_velocity(particle_set p_set, double timedelta) {
    for (int i = 0; i < p_set.number; i ++) {
        velocity_update(&(p_set.particles[i]), timedelta);
    }
}


void update_distances(particle_set p_set, double timedelta) {
    for (int i = 0; i < p_set.number; i ++) {
        update_distance(&(p_set.particles[i]), timedelta);
    }
}

void reset_accelerations(particle_set p_set) {
    for (int i = 0; i < p_set.number; i ++) {
        reset_acceleration(&(p_set.particles[i]));
    }
}
