import itertools
import math
from functools import partial
from copy import deepcopy


E0 = 1
sqrt_e_const = 4.69041575982343e-08
MIN_NUMB = 0.1 ** 10


def logger(string, *args):
    pass
    # print(string % args)


def safe_value(value):
    if abs(value) < MIN_NUMB:
        return MIN_NUMB
    return value


class Particle(object):
    count = 0
    def __init__(self, x_c, y_c, z_c, x_v=None, y_v=None, z_v=None):
        self.x_c = float(x_c)
        self.y_c = float(y_c)
        self.z_c = float(z_c)

        self.x_v = float(x_v)
        self.y_v = float(y_v)
        self.z_v = float(z_v)

        self.acc_x = 0.0
        self.acc_y = 0.0
        self.acc_z = 0.0

        self.new_acc_x = 0.0
        self.new_acc_y = 0.0
        self.new_acc_z = 0.0

        self.index = Particle.count
        Particle.count += 1

    def __repr__(self):
        # return str(self.index)
        return (
            ("%.16f" % self.x_c) + ' ' +
            ("%.16f" % self.y_c) + ' ' +
            ("%.16f" % self.z_c) + ' ' +
            ("%.16f" % self.x_v) + ' ' +
            ("%.16f" % self.y_v) + ' ' +
            ("%.16f" % self.z_v)
        )

    def shift_x(self, value):
        new_particle = deepcopy(self)
        new_particle.x_c += value
        return new_particle

    def shift_y(self, value):
        new_particle = deepcopy(self)
        new_particle.y_c += value
        return new_particle

    def shift_z(self, value):
        new_particle = deepcopy(self)
        new_particle.z_c += value
        return new_particle

    def with_new_coordinates(self, x, y, z):
        new_particle = deepcopy(self)
        new_particle.x_c = x
        new_particle.y_c = y
        new_particle.z_c = z
        return new_particle

    def update_acceleration(self, particle_2, particle_3):
        logger(
            'In %s: updating particles %s and %s',
            self.index, particle_2.index, particle_3.index
        )
        logger(
            'Before update the coordinates are x: %s, y: %s, z: %s',
            self.new_acc_x, self.new_acc_y, self.new_acc_z
        )
        potential_func = partial(calc_potential, particle_2, particle_3)
        h_x = safe_value(sqrt_e_const * self.x_c)
        self.new_acc_x -= (
            2 * (
                potential_func(self.shift_x(h_x))
                - potential_func(self.shift_x(-h_x))
            ) / ((self.x_c + h_x) - (self.x_c - h_x))
        )

        h_y = safe_value(sqrt_e_const * self.y_c)
        self.new_acc_y -= (
            2 * (
                potential_func(self.shift_y(h_y))
                - potential_func(self.shift_y(-h_y))
            ) / ((self.y_c + h_y) - (self.y_c - h_y))
        )

        h_z = safe_value(sqrt_e_const * self.z_c)
        self.new_acc_z -= (
            2 * (
                potential_func(self.shift_z(h_z))
                - potential_func(self.shift_z(-h_z))
            ) / ((self.z_c + h_z) - (self.z_c - h_z))
        )
        logger(
            'After update the coordinates are x: %s, y: %s, z: %s',
            self.new_acc_x, self.new_acc_y, self.new_acc_z
        )

    def reset_acceleration(self):
        self.acc_x = self.new_acc_x
        self.acc_y = self.new_acc_y
        self.acc_z = self.new_acc_z

        self.new_acc_x = 0.0
        self.new_acc_y = 0.0
        self.new_acc_z = 0.0

    def update_distance(self, timedelta):
        self.x_c += (
            self.x_v * timedelta
            + 1.0 / 2 * self.acc_x * timedelta ** 2
            # + 1.0 / 2 * self.acc_x * timedelta
        )

        self.y_c += (
            self.y_v * timedelta
            + 1.0 / 2 * self.acc_y * timedelta ** 2
            # + 1.0 / 2 * self.acc_y * timedelta
        )

        self.z_c += (
            self.z_v * timedelta
            + 1.0 / 2 * self.acc_z * timedelta ** 2
            # + 1.0 / 2 * self.acc_z * timedelta
        )

    def velocity_update(self, timedelta):
        # After deleting 1/2 the results are more or less identical.
        self.x_v += (1.0 / 2) * (self.acc_x + self.new_acc_x) * timedelta
        self.y_v += (1.0 / 2) * (self.acc_y + self.new_acc_y) * timedelta
        self.z_v += (1.0 / 2) * (self.acc_z + self.new_acc_z) * timedelta
        # self.x_v += (self.acc_x + self.new_acc_x) * timedelta
        # self.y_v += (self.acc_y + self.new_acc_y) * timedelta
        # self.z_v += (self.acc_z + self.new_acc_z) * timedelta
        self.reset_acceleration()


class ParticleSet(object):

    def __init__(self, particles):
        self.particles = particles
        self.calc_accelerations()
        self.reset_acceleration()

    def calc_accelerations(self):
        for i in xrange(len(self.particles)):
            for j in xrange(i + 1, len(self.particles)):
                for k in xrange(j + 1, len(self.particles)):
                    part_a = self.particles[i]
                    part_b = self.particles[j]
                    part_c = self.particles[k]

                    part_a.update_acceleration(part_b, part_c)
                    part_b.update_acceleration(part_a, part_c)
                    part_c.update_acceleration(part_a, part_b)

    def update_distances(self, timedelta):
        for particle in self.particles:
            particle.update_distance(timedelta)

    def __repr__(self):
        return '\n'.join(str(p) for p in self.particles)

    def velocity_update(self, timedelta):
        for particle in self.particles:
            particle.velocity_update(timedelta)

    def reset_acceleration(self):
        for particle in self.particles:
            particle.reset_acceleration()


def norm_distance(particle_a, particle_b):
    x = abs(particle_a.x_c - particle_b.x_c)
    y = abs(particle_a.y_c - particle_b.y_c)
    z = abs(particle_a.z_c - particle_b.z_c)
    return max(MIN_NUMB, math.sqrt(x ** 2 + y ** 2 + z ** 2))


def calc_potential(particle_i, particle_j, particle_k):
    rij = norm_distance(particle_i, particle_j)
    rik = norm_distance(particle_i, particle_k)
    rkj = norm_distance(particle_k, particle_j)

    srij = rij ** 2
    srik = rik ** 2
    srkj = rkj ** 2
    first_part = (1 / (rij * rik * rkj) ** 3)
    second_part_l = (
        3.0 * (-srij + srik + srkj)
        * (srij - srik + srkj) * (srij + srik - srkj)
    )
    second_part_m = 8.0 * (rij * rik * rkj) ** 5
    return first_part + second_part_l / second_part_m


def read_particle(input_line):
    values = input_line.split(' ')
    particle = Particle(
        values[0], values[1], values[2], values[3], values[4], values[5]
    )
    return particle


def read_input(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        return [read_particle(x) for x in lines]


def run_algorithm():
    timedelta = 0.5
    particles_set = ParticleSet(read_input('tests/part4.txt'))

    particles_set.update_distances(timedelta)
    particles_set.calc_accelerations()
    particles_set.velocity_update(timedelta)
    print(particles_set)


if __name__ == "__main__":
    run_algorithm()
