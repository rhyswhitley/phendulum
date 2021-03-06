#!/usr/bin/env python2

import matplotlib.pyplot as plt
import datetime, time

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

class spring(object):
    """
    This class simulates the dynamics of a simple pendulum whose momentum is
    driven by some arbitrary, external force through a viscous fluid, such that
    it experiences a drag force that holds it at a new equilibrium point. The idea
    is that the pendulum represents plant phenology and the external force may
    be water, temperature, etc. driving such phenology.
    """

    def __init__(self, time, force_env=1):
        self.mass = 1
        self.force_env = force_env
        self.dt = 1e-2
        self.k_hooke = 1
        self.k_drag = 2
        self.time = time

    def hookes_law(self, x):
        # Hooke's Law for a spring
        return -self.k_hooke*x
    def drag(self, v):
        # drag on pendulum
        return -self.k_drag*v

    def dynamic(self):
        # set zeros (could have empty list) -- ugly either way
        displ       = [0]*self.time
        accel       = [0]*self.time
        veloc       = [0]*self.time
        force_drag  = [0]*self.time
        force_resist= [0]*self.time
        t_force_env = [0]*self.time
        # instantaneous vectors
        for t in range(self.time-1):
            # point past which the environmental force is removed
            if t < (self.time/2):
                t_force_env[t] = self.force_env
            else:
                t_force_env[t] = 0
            # forces on the pendulum
            force_drag[t] = self.drag( veloc[t] )
            force_resist[t] = self.hookes_law( displ[t] )
            force_total = t_force_env[t] + force_resist[t] + force_drag[t]
            # vectors of the pendulum
            accel[t+1] = force_total/self.mass
            veloc[t+1] = veloc[t] + accel[t]*self.dt
            displ[t+1] = displ[t] + veloc[t]*self.dt
        # return wave as tuple
        return {'a':accel, 'v':veloc, 'x':displ,
                'Fe':t_force_env, 'Fr':force_resist, 'Fd':force_drag}

    def plot_spring_dynamics(self):
        motion = self.dynamic()

        #fig = plt.figure( figsize=(4,8) )
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        # leaf phenology
        ax1.plot( motion['x'], color='green', label='$x(t)$', linestyle='-', lw=2 )
        # phenology vectors
        ax2.plot( motion['a'], color='red', label='$a(t)$', linestyle='-', lw=2 )
        ax2.plot( motion['v'], color='blue', label='$v(t)$', linestyle='-', lw=2 )
        # forces
        ax3.plot( motion['Fe'], color='black', label='$F_{e}(t)$', linestyle='--', lw=2 )
        ax3.plot( motion['Fr'], color='lightgray', label='$F_{r}(t)$', linestyle='-', lw=2 )
        ax3.plot( motion['Fd'], color='gray', label='$F_{d}(t)$', linestyle='-', lw=2 )
        # legends
        ax1.legend( loc=1 )
        ax2.legend( loc=1 )
        ax3.legend( loc=1 )
        # axis labels
        ax1.set_ylabel( 'phenology', fontsize=14 )
        ax2.set_ylabel( 'vector', fontsize=14 )
        ax3.set_ylabel( 'force', fontsize=14 )
        ax3.set_xlabel( 'time (t)', fontsize=14 )
        # adjust spacing
        plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1,
                            hspace=0.1, wspace=0.4)
        # print file
        plt.savefig("pendulum_sig.pdf")


tseries = 2000
x = spring( tseries )
x.plot_spring_dynamics()


