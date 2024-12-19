# import figure_properties as fp

import numpy as np
import matplotlib.pyplot as plt
from sleep_cycles import compute_fourier_sleep, ttx

def_colors = {'ret': '#4d9221',
              'fet': '#c51b7d'}


def plot_cycles(ax, age=25):
    s_vals, f_vals, rem_thres = compute_fourier_sleep(age, debug=False)
    angs = np.arange(0, 2*np.pi, 2*np.pi/(24/np.diff(ttx)[0]))
    night = np.zeros_like(angs)
    night[np.where(s_vals < 0)] = 1
    pot_rem = np.zeros_like(angs)
    pot_rem[np.where(f_vals >= rem_thres)] = 1
    r_rem = night*pot_rem
    ax.plot(angs, s_vals)  # awake state
    ax.plot(angs, f_vals)  # underlying rythm
    ax.plot(angs, r_rem)  # rem sleep
    # ax.set_rlim(-2, 2)
    return ax


def plot_background(ax, age=25, show_cutoffs=False, fill_in=False):
    s_vals, f_vals, rem_thres = compute_fourier_sleep(age, debug=False)
    angs = np.arange(0, 2*np.pi, 2*np.pi/(24/np.diff(ttx)[0]))
    night = np.zeros_like(angs)
    night[np.where(s_vals < 0)] = 1
    pot_rem = np.zeros_like(angs)
    pot_rem[np.where(f_vals >= rem_thres)] = 1
    r_rem = night*pot_rem

    slp_on = angs[np.where(s_vals < 0)[0][1]]
    ax.plot([0, 0], [0, 2],
            color="black", linewidth=2)
    ax.plot([slp_on, slp_on], [0, 2],
            color="black", linewidth=2)

    r_circle = np.ones_like(angs)
    r_circle0 = np.zeros_like(angs)

    r_nrem_zone = np.maximum(r_circle0, night)
    r_rem_zone = np.minimum(r_circle, r_rem)

    ax.fill_between(angs, 0, r_nrem_zone, alpha=0.5, color='k')
    ax.fill_between(angs, 0, r_rem_zone, alpha=0.5, color='white')
    ax.set_rlim(0, 1.6)
    ax.text(slp_on, 1.2, "Sleep\non")
    ax.text(0, 1.2, "Awake")
    ax.set_title('Age: ' + str(age), loc='left')
    if show_cutoffs:
        r_circ_ret = np.ones_like(angs)*0.4
        ax.plot(angs, r_circ_ret, def_colors['ret'], ls='--', lw=1.5)
        r_circ_fet = np.ones_like(angs)*1.6
        ax.plot(angs, r_circ_fet, def_colors['fet'], ls='--', lw=1.5)
        if fill_in:
            ax.fill_between(angs, 0, r_circ_ret, alpha=0.5, color=def_colors['ret'])
            ax.fill_between(angs, r_circ_fet, 1, alpha=0.5, color=def_colors['fet'])
    return ax


def clean_ax(ax):
    ax.set_rticks([])
    ax.set_thetagrids([])
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2],
                  ['0 hrs', '6 hrs', '12 hrs', '18 hrs'])
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi/2.0)
    return ax


if __name__ == '__main__':

    ax = plt.subplot(111, projection='polar')
    # for age in [2, 25, 60, 90]:
    #     ax = plot_background(age=age)
    #     plt.savefig('age_clock_'+str(age)+'.png', dpi=150)
    #     plt.clf()

    age = 50
    # ax = plot_background(ax, age=age)
    ax = plot_cycles(ax, age=age)

    # ax = plot_background(age=age, show_cutoffs=True)
    # plt.savefig('age_clock_coff_'+str(age)+'.png', dpi=150)

    # age = 25
    # ax = plot_background(age=age, show_cutoffs=True, fill_in=True)
    # plt.savefig('age_clock_coff_ret'+str(age)+'.png', dpi=150)

    # age = 25
    # ax = plot_background(ax, age=age, show_cutoffs=True, fill_in=True)
    # plt.savefig('age_clock_coff_ret+fet'+str(age)+'.png', dpi=150)
    clean_ax(ax)
    plt.show()

