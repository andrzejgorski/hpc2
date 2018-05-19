#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "particles.h"


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


void print_p_set(particle_set p_set) {
    printf("Particle set with number %d\n\n", p_set.number);
    for(int i = 0; i < p_set.number; i ++) {
        print_particle(p_set.particles[i]);
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
