#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particles.h"
#include "files.h"


particle_set read_particle(char* filename) {
    FILE* file;

    int lines = count_lines_in_file(filename);
    int index = 0;

    particle_set p_set = new_particle_set(lines);

    file = fopen(filename, "r");
    if (file) {
        while (true) {
            if (fscanf(file, "%lf", &(p_set.particles[index].x_c)) == EOF)
                break;
            if (fscanf(file, "%lf", &(p_set.particles[index].y_c)) == EOF)
                break;
            if (fscanf(file, "%lf", &(p_set.particles[index].z_c)) == EOF)
                break;

            if (fscanf(file, "%lf", &(p_set.particles[index].x_v)) == EOF)
                break;
            if (fscanf(file, "%lf", &(p_set.particles[index].y_v)) == EOF)
                break;
            if (fscanf(file, "%lf", &(p_set.particles[index].z_v)) == EOF)
                break;

            p_set.particles[index].index = index++;
        }
        fclose(file);
    }
    return p_set;
}


int main (int argc, char** argv) {
    if (argc <= 5) {
        fprintf(stderr, "Error: Program require more arguments\n");
    }

    char *input_filename = argv[1];
    char *output_filename = argv[2];
    int stepcount = atoi(argv[3]);
    int deltatime = atoi(argv[4]);

    // if present puts the result after each step i
    // (counted from 1) in a file particles_out_i.txt
    bool debug_mode;

    if (argc == 6 && strcmp(argv[5], "-v") == 0) {
        debug_mode = true;
    } else {
        debug_mode = false;
    }

    particle_set p_set = read_particle(input_filename);
    // print_p_set(p_set);

    update_distances(p_set, deltatime);
    calc_accelerations(p_set);
    update_velocity(p_set, deltatime);

    print_to_lines_p_set(p_set);

    free_particle_set(p_set);

    // printf("Input filename %s\n", input_filename);
    // printf("Output filename %s\n", output_filename);
    // printf("Stepcount %d\n", stepcount);
    // printf("Deltatime %d\n", deltatime);
    // printf("argc %d\n", argc);
    // printf("Debug mode %s\n", argv[5]);
    // printf("Debug mode %d\n", debug_mode);
}
