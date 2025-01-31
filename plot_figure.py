import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.transforms import ScaledTranslation

from pwm_fourier import cmpt_fourier, cmpt_ss_durs, cmpt_aging_sleep, cmpt_tot_prop_sleep
# from pwm_fourier import cmpt_prop_sleep


def neat_axs(ax_list, r=True, t=True, b=True):
    for ax in ax_list:
        # ax.get_xaxis().set_visible(False)
        ll = []
        if r:
            ll.append('right')
        if t:
            ll.append('top')
        if b:
            ll.append('bottom')
        ax.spines[ll].set_visible(False)
    return ax_list


def draw_labels(ax, xx, tt, rem=True):
    # labels for rem cycles
    if rem:
        l_strs = ['R1', 'R2', 'R3', 'R4']
        y_offset = 0.05
        lab_diff = -0.02
    else:
        l_strs = ['Q1', 'Q2', 'Q3', 'Q4']
        y_offset = -0.39
        lab_diff = 0.03
    whens = []
    for ii in xx:  # places small markers
        ax.plot([tt[ii], tt[ii]],
                [y_offset, y_offset+0.01], c='k', lw=0.5)
        whens.append(tt[ii])
    whens.append(1)
    if rem:
        ax.plot([tt[-1], tt[-1]],   # places small markers
                [y_offset, y_offset+0.01], c='k', lw=0.5)
        ax.plot([tt[xx[0]], tt[-1]],
                [y_offset+0.005, y_offset+0.005], c='k', lw=0.5)
    else:
        ax.plot([tt[xx[0]], tt[xx[-1]]],
                [y_offset+0.005, y_offset+0.005], c='k', lw=0.5)
    for iii in range(4):
        ax.text(x=(whens[iii] + whens[iii+1])/2,
                y=y_offset+lab_diff, s=l_strs[iii],
                va='center', ha='center')
    return ax
        

def plot_left(ax1, ax2, ax3, ax4):
    tt, circ_rhy, clck_24, sleep_pressure, n3_cutoff, yy = cmpt_fourier(total_width,
                                                                        duty)
    sleepstate, s_n1n2, s_n2n3, s_n3n2, s_nr, s_rn1 = cmpt_ss_durs(yy, n3_cutoff,
                                                                   circ_rhy, sleep_pressure)
    # circ_rhy and clck_24 plots
    l1, = ax1.plot(tt, circ_rhy, label='C(A,t)', color='C2')
    l2, = ax1.plot(tt, clck_24, label='K(t)', color='C1')
    ax1.text(0.2, 0.55, s='C(A,t)', color=l1.get_color())
    ax1.text(0.2, -0.4, s='K(t)', color=l2.get_color())
    ax1.set_xticks([0, 0.7, 1])
    ax1.set_xticklabels([0, 'D', 24])
    ax1.set_xlabel('Awake since (hrs)')
    ax1.hlines(0, xmin=0, xmax=1, ls='--', color='gray', lw=0.5)
    ax1.set_ylim(-1.2, 1.2)
    # sleep pressure and its fourier
    rect = patches.Rectangle((0.68, -0.4), 0.32, 0.55, linewidth=0.1,
                             edgecolor='white', facecolor='lightgray')
    ax2.hlines(0, xmin=0, xmax=1, ls='--', color='gray', lw=0.5)
    ax2.add_patch(rect)
    l3, = ax2.plot(tt, circ_rhy+clck_24, label='S(A,t)', c='red')
    l4, = ax2.plot(tt, yy, c='k', lw=1., label='Z(A,t)')

    # ax2.legend(frameon=False, ncols=2)
    ax2.text(0.1, 0.55, s='S(A,t)', color=l3.get_color())
    ax2.text(0.1, -0.4, s='Z(A,t)', color=l4.get_color())
    ax2.set_xticks([0, 0.7, 1])
    ax2.set_xticklabels([0, 'D', 24])
    ax2.set_xlabel('Awake since (hrs)')
    ax2.set_ylim(-1.2, 1.2)
    ax2.set_yticks([-1, 0, 1])
    ax2.set_yticklabels([])
    ax2.text(-.15, -0.05, 'Sleep', transform=ax2.transAxes,
             rotation=90, va='bottom', ha='center')
    ax2.text(-.15, 0.55, 'Awake', transform=ax2.transAxes,
             rotation=90, va='bottom', ha='center')
    # zoom in of ax2 sleep part
    ax3.hlines(0, xmin=0, xmax=1, ls='--', color='gray', lw=0.5)
    ax3.plot(tt, sleep_pressure, c='r')
    ax3.plot(tt, yy, c='k')
    draw_labels(ax3, s_nr[0], tt)  # REM label line
    draw_labels(ax3, s_n1n2[0], tt, rem=False)  # non-REM label line
    ax3.fill_between(tt, sleep_pressure, yy, where=(sleepstate == 0),
                     color='C0', alpha=0.8)
    ax3.fill_between(tt, n3_cutoff, yy, where=(sleepstate == -3),
                     color='k', alpha=1)
    ax3.text(-.05, 0.5, 'Sleep', transform=ax3.transAxes,
             rotation=90, va='bottom', ha='center')
    ax3.text(-.05, 0.75, 'Awake', transform=ax3.transAxes,
             rotation=90, va='bottom', ha='center')
    ax3.set_yticks([0])
    ax3.set_yticklabels([])
    ax3.set_ylim(-0.4, 0.15)
    # hypnogram
    ax4.plot(tt, sleepstate)
    ax4.tick_params(which='both', length=7)
    ax4.set_yticks([-3, -2, -1, 0, 1])
    ax4.set_yticklabels(['N3', 'N2    ', 'N1', 'R     ', 'A'])
    # setting limits
    ax3.set_xlim(0.68, 1)
    ax4.set_xlim(0.68, 1)
    neat_axs([ax1, ax2], b=False)
    neat_axs([ax3, ax4])
    ax3.set_xticks([])
    ax4.set_xticks([])
    return ax1, ax2, ax3, ax4


# def plot_right(ax1, ax2, ax3):
#     ages, sleepmins, remmins, nremmins = cmpt_aging_sleep(dt_age=1)
#     # ax1.fill_between(ages, np.ones_like(sleepmins)*24, np.zeros_like(sleepmins),
#     #                  color='white', alpha=1, label='Awake')
#     ax1.fill_between(ages, np.array(sleepmins)/60, np.zeros_like(sleepmins),
#                      color='C0', alpha=0.8, label='REM')
#     ax1.fill_between(ages, np.array(nremmins)/60, np.zeros_like(nremmins),
#                      color='lightgray', alpha=1, label='non-REM')
#     ax1.text(x=55, y=17, s='Awake', va='center', ha='center')
#     ax1.legend(frameon=False, ncols=2)
#     ax1.set_ylim(0, 24)
#     ax1.set_xlim(0, 100)
#     ax1.set_yticks([4, 8, 12, 16, 20, 24])
#     ax1.set_ylabel('Sleep per day (hrs)')
#     ax1.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
#     ax1.set_xlabel('Age (yrs)')
#     test_ages = [20, 50, 70]
#     # dt_mins = 0.001*24*60
#     rtimes, ntimes = cmpt_prop_sleep(test_ages)
#     print(rtimes, ntimes)
#     for ii, xx in enumerate(test_ages):
#         ax2.plot(np.arange(len(rtimes[ii])-1), np.diff(rtimes[ii])*24*60)
#     for ii, xx in enumerate(test_ages):
#         ax3.plot(np.arange(len(ntimes[ii])-1), np.diff(ntimes[ii])*24*60)
#     neat_axs([ax2, ax3])
#     return ax1, ax2, ax3

def plot_right(ax1, ax2):
    ages, sleepmins, remmins, nremmins = cmpt_aging_sleep(dt_age=1)
    # ax1.fill_between(ages, np.ones_like(sleepmins)*24, np.zeros_like(sleepmins),
    #                  color='white', alpha=1, label='Awake')
    ax1.fill_between(ages, np.array(sleepmins)/60, np.zeros_like(sleepmins),
                     color='C0', alpha=0.8, label='REM')
    ax1.fill_between(ages, np.array(nremmins)/60, np.zeros_like(nremmins),
                     color='lightgray', alpha=1, label='non-REM')
    # ax1.text(x=55, y=17, s='Awake', va='center', ha='center')
    # ax1.legend(frameon=False, ncols=2)
    colors = {'REM': 'C0', 'non-REM': 'lightgray', 'Awake': 'white'}
    alphas = {'REM': 0.8, 'non-REM': 1, 'Awake': 1}
    edgeC = {'REM': 'w', 'non-REM': 'w', 'Awake': 'k'}
    labels = list(colors.keys())
    handles = [plt.Rectangle((0, 0), 0.2, 0.2, facecolor=colors[label],
                             alpha=alphas[label], edgecolor=edgeC[label]) for label in labels]
    ax1.legend(handles, labels, frameon=False, ncols=3, borderpad=0.1, columnspacing=1, handletextpad=0.2)
    ax1.set_ylim(0, 24)
    ax1.set_xlim(0, 100)
    ax1.set_yticks([4, 8, 12, 16, 20, 24])
    ax1.set_ylabel('Sleep per day (hrs)')
    ax1.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax1.set_xlabel('Age (yrs)')
    test_ages = [20, 40, 60, 80]
    r_, n1_, n2_, n3_ = cmpt_tot_prop_sleep(test_ages)
    for ii, ag in enumerate(test_ages):
        ax2 = plot_total_bar(ax2, r_[ii], n1_[ii],
                             n2_[ii], n3_[ii], ii)
    colors = {'REM': 'C0', 'N1': 'darkgray', 'N2': 'gray', 'N3': 'k'}
    alphas = {'REM': 0.8, 'N1': 1, 'N2': 1, 'N3': 1}
    labels = list(colors.keys())
    handles = [plt.Rectangle((0, 0), 0.2, 0.2, color=colors[label],
                             alpha=alphas[label]) for label in labels]
    ax2.legend(handles, labels, frameon=False, ncols=4, borderpad=0.1, columnspacing=1, handletextpad=0.2)
    ax2.set_ylabel('Sleep per day (hrs)')
    ax2.set_xticks(range(len(test_ages)))
    ax2.set_xticklabels(test_ages)
    ax2.set_xlabel('Age (yrs)')
    ax2.set_ylim(0, 10)
    # ax2.set_xlim(-1, len(test_ages)+2)
    neat_axs([ax2])
    return ax1, ax2


def plot_total_bar(ax, r_, n1_, n2_, n3_, x=0):
    total = n3_ + n2_ + n1_ + r_
    width = 0.6
    ax.bar(x, n3_, color='k', width=width)
    ax.bar(x, n2_, bottom=n3_, color='gray', width=width)
    ax.bar(x, n1_, bottom=n3_+n2_, color='darkgray', width=width)
    ax.bar(x, r_, bottom=n3_+n2_+n1_, color='C0', width=width, alpha=0.8)
    ax.text(x, n3_ / 2, f'{n3_/total*100:.0f}%', ha='center',
            va='center', color='white')
    ax.text(x, n3_ + (n2_/2), f'{n2_/total*100:.0f}%', ha='center',
            va='center', color='white')
    ax.text(x, n3_ + n2_ + (n1_/2), f'{n1_/total*100:.0f}%', ha='center',
            va='center', color='white')
    ax.text(x, n3_ + n2_ + n1_ + (r_/2), f'{r_/total*100:.0f}%',
            ha='center', va='center', color='white')
    return ax


def plot_subplotlabels(fig, axs):
    for label, ax in axs.items():
        ll = label.lower()
        ax.text(
            0.0, 1.0, ll, transform=(ax.transAxes + ScaledTranslation(-30/72, 3/72, fig.dpi_scale_trans)),
            fontsize='large', va='bottom', weight='bold')
    return axs

if __name__ == '__main__':
    left = [['A', 'B'],
            ['A', 'B'],
            ['C', 'C'],
            ['C', 'C'],
            ['C', 'C'],
            ['C', 'C'],
            ['C', 'C'],
            ['D', 'D']]
    right = [['E',  'E'],
             ['E',  'E'],
             ['E',  'E'],
             ['E',  'E'],
             ['F',  'F'],
             ['F',  'F'],
             ['F',  'F'],
             ['F',  'F']]
    
    # right = [['E',  'E'],
    #          ['E',  'E'],
    #          ['E',  'E'],
    #          ['E',  'E'],
    #          ['F',  'G'],
    #          ['F',  'G'],
    #          ['H',  'G'],
    #          ['H',  'G']]
    # plot_right(ax_dict['E'], ax_dict['F'], ax_dict['H'])
    grid_spec = [[left, right]]
    fig, ax_dict = plt.subplot_mosaic(grid_spec, figsize=(7, 4.5), layout='constrained')
    total_width, duty = 1, 0.7
    plot_left(ax_dict['A'], ax_dict['B'], ax_dict['C'], ax_dict['D'])
    plot_right(ax_dict['E'], ax_dict['F'])
    plot_subplotlabels(fig, ax_dict)
    # plt.tight_layout(pad=0.5)
    #plt.show()
    # print('Done')
    plt.savefig('hypnogram.png', dpi=300)
    #plt.savefig('hypnogram.svg')
