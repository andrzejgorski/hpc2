
typedef int bool;
#define true 1
#define false 0

struct __particle {
    unsigned int index;

    double x_c;
    double y_c;
    double z_c;

    double x_v;
    double y_v;
    double z_v;

    double acc_x;
    double acc_y;
    double acc_z;

    double new_acc_x;
    double new_acc_y;
    double new_acc_z;

};

typedef struct __particle particle;

struct __particle_set {
    int number;
    particle *particles;
};
typedef struct __particle_set particle_set;

void print_particle(particle p);
void print_p_set(particle_set p_set);
void print_to_lines_p_set(particle_set p_set);
void free_particle_set(particle_set p_set);
particle_set new_particle_set(int number);

void update_distances(particle_set p_set, double timedelta);
void update_velocity(particle_set p_set, double timedelta);
void calc_accelerations(particle_set p_set);
void reset_accelerations(particle_set p_set);
