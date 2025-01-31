import numpy as np
import matplotlib.pyplot as plt


def fetch_period(age):
    lls = [1, 2, 3, 4, 6, 8, 12, 24]  # First 2 yrs of age
    # number of hrs in a sleep-awake cycle
    if age < 2:
        return lls[int(age/0.25)]  # cycle per 24 hrs
    else:
        return 24


def fetch_duty(age):
    # correction for the amount of overall sleep per day over age
    # 0yr -> 16 hrs sleep per 24hrs  D(A) = 0.33
    # 2yr -> 12 hrs sleep per 24hrs  D(A) = 0.5
    # 16yr -> 8 hrs sleep per 24hrs  D(A) = 0.66
    # 80yr -> 6 hrs sleep per 24hrs  D(A) = 0.75
    if age <= 2:
        return 0.33333 + (age)*(0.1666667/2)
    elif age <= 16:
        return 0.5 + (age-2)*(0.16666/14)
    else:
        return 0.66667 + (age-16)*(0.083333/64)

    
def pwm(tt, v, duty, total_width):
    ff = np.ones_like(tt)*v
    ff[np.where(tt > (duty*total_width))] = 0
    return ff


def ramp(tt, v, total_width):
    ff = -v + (tt*v/total_width)
    return ff


def ramp_n3(tt, v, total_width, offset=0.025):
    ff = -(v + offset) + (tt*v/total_width)
    return ff


def comp_coeffs_sqr(n, v, duty):
    an = v*np.sin(2*np.pi*n*duty)/(n*np.pi)
    bn = 2*v*np.sin(np.pi*n*duty)**2/(n*np.pi)
    return an, bn


def comp_coeffs_ramp(n, v):
    an = 0
    bn = -v/(n*np.pi)
    return an, bn


def sqr_parts(v, duty, max_n, tt, total_width):
    a0_s = v*duty
    val = np.zeros_like(tt)
    for ii in range(max_n):
        cos_comp = np.cos(2*np.pi*(ii+1)*tt/total_width)
        sin_comp = np.sin(2*np.pi*(ii+1)*tt/total_width)
        an_s, bn_s = comp_coeffs_sqr(ii+1, v, duty)
        val += an_s*cos_comp + bn_s*sin_comp
    return val + a0_s


def ramp_parts(v, max_n, tt, total_width):
    a0_r = v/2 - 1
    val = np.zeros_like(tt)
    for ii in range(max_n):
        cos_comp = np.cos(2*np.pi*(ii+1)*tt/total_width)
        sin_comp = np.sin(2*np.pi*(ii+1)*tt/total_width)
        an_r, bn_r = comp_coeffs_ramp(ii+1, v)
        val += an_r*cos_comp + bn_r*sin_comp
    return val + a0_r


def cmpt_ss_durs(yy, n3_cutoff, circ_rhy, sleep_pressure):
    slp = np.diff(yy)  # slope
    drctn = np.sign(slp)  # direction of slope
    dir_change = np.diff(drctn)  # acceleration / deceleration
    # N3 state compute
    bb = np.zeros_like(yy)
    bb[np.argwhere(yy < n3_cutoff)] = 1
    bb[np.where(circ_rhy > 0)] = 0  # reset those in awake state
    in_n3 = np.where(bb > 0)
    s_n2n3 = [in_n3[0][0]]  # N2 -> N3 switch
    s_n3n2 = []  # N3 -> N2 switch
    for ij in np.where(np.diff(in_n3) > 1)[1]:
        s_n3n2.append(in_n3[0][ij])
        s_n2n3.append(in_n3[0][ij+1])
    s_n3n2.append(in_n3[0][-1])
    # REM sleep and N1 state
    cc = np.zeros_like(yy)
    cc[np.argwhere(dir_change == -2)+1] = 1  # end of REM
    cc[np.where(circ_rhy > 0)] = 0  # reset those in awake state
    s_rn1 = np.where(cc > 0)  # REM -> N1
    sleepstate = np.ones_like(yy)
    sleepstate[np.where(yy < 0)] = -1  # N1   (set everything as N1 to begin)
    sleepstate[(np.hstack((0, drctn)) >= 0) & (yy < 0)] = 0  # rem
    sleepstate[(yy < sleep_pressure) & (yy < 0)] = -2  # N2
    sleepstate[in_n3] = -3  # N3
    # N2 state
    aa = np.zeros_like(yy)
    aa[np.where(yy < sleep_pressure)] = 1
    aa[np.where(circ_rhy > 0)] = 0
    d = np.diff(aa)
    idxx, = d.nonzero()
    s_n1n2 = np.where(d > 0)  # N1 -> N2
    s_nr = np.where(d < 0)  # N2/N3 -> REM
    return sleepstate, s_n1n2, s_n2n3, s_n3n2, s_nr, s_rn1


def cmpt_circ_clck(total_width, duty):
    tt = np.arange(0, total_width, 0.001)
    v = 1  # height of modulation - dont change this!
    # duty = 0.7  # fraction of total_width the input is on
    circ_rhy = pwm(tt, v, duty, total_width)
    clck_24 = ramp(tt, v, total_width)
    return tt, v, circ_rhy, clck_24


def cmpt_sleep_press(v, duty, tt, total_width, N=15):
    yy1 = sqr_parts(v, duty, N, tt, total_width)
    yy2 = ramp_parts(v, N, tt, total_width)
    return yy1 + yy2


def cmpt_fourier(total_width=1, duty=0.8):
    # total_width = 1  # 2*np.pi
    tt, v, circ_rhy, clck_24 = cmpt_circ_clck(total_width, duty)
    sleep_pressure = circ_rhy + clck_24
    n3_cutoff = ramp_n3(tt, v, total_width)
    yy = cmpt_sleep_press(v, duty, tt, total_width, 15)
    return tt, circ_rhy, clck_24, sleep_pressure, n3_cutoff, yy


def cmpt_time_in(tt, s_n1n2, s_nr, s_n2n3, s_n3n2, s_rn1):
    n1times = []  # in minutes
    n2times = []
    n3starttimes = []
    n3endtimes = []
    rtimes = []  # in minutes
    for ii in s_n1n2:
        n2times.append(tt[ii]*24*60)
        # plt.vlines(tt[ii], ymax=1, ymin=-1, color='r')
    for ii in s_nr:
        rtimes.append(tt[ii]*24*60)
        # plt.vlines(tt[ii], ymax=1, ymin=-1, color='g')
    for ii in s_n2n3:
        n3starttimes.append(tt[ii]*24*60)
        # plt.vlines(tt[ii], ymax=1, ymin=-1, color='k')
    for ii in s_n3n2:
        n3endtimes.append(tt[ii]*24*60)
        # plt.vlines(tt[ii], ymax=1, ymin=-1, color='purple')
    for ii in s_rn1:
        n1times.append(tt[ii]*24*60)
        # plt.vlines(tt[ii], ymax=1, ymin=-1, color='gold')
    return n1times, n2times, n3starttimes, n3endtimes, rtimes


def cmpt_aging_sleep(dt_age=1):
    dt_mins = 0.001*24*60
    ages = np.arange(0.1, 100, dt_age)
    sleepmins = []
    remmins = []
    nremmins = []
    for ii in ages:
        P = fetch_period(ii)
        num_cycles = int(24/P)  # cycles per night
        d = fetch_duty(ii)
        # print(ii, P, num_cycles, d)
        tt, circ_rhy, clck_24, sleep_P, n3_co, yy = cmpt_fourier(total_width=P/24,
                                                                 duty=d)
        sleepstate, s_n1n2, s_n2n3, s_n3n2, s_nr, s_rn1 = cmpt_ss_durs(yy,
                                                                       n3_co,
                                                                       circ_rhy,
                                                                       sleep_P)
        sleepmins.append(len(np.where(sleepstate != 1)[0])*dt_mins*num_cycles)
        remmins.append(len(np.where(sleepstate == 0)[0])*dt_mins*num_cycles)
        nremmins.append(len(np.where(sleepstate < 0)[0])*dt_mins*num_cycles)
    return ages, sleepmins, remmins, nremmins


def cmpt_tot_prop_sleep(ages):
    dt_hrs = 0.001*24
    rmins = []
    n1mins = []
    n2mins = []
    n3mins = []
    for ii in ages:
        P = fetch_period(ii)
        num_cycles = int(24/P)  # cycles per night
        d = fetch_duty(ii)
        # print(ii, P, num_cycles, d)
        tt, circ_rhy, clck_24, sleep_P, n3_co, yy = cmpt_fourier(total_width=P/24,
                                                                 duty=d)
        sleepstate, s_n1n2, s_n2n3, s_n3n2, s_nr, s_rn1 = cmpt_ss_durs(yy,
                                                                       n3_co,
                                                                       circ_rhy,
                                                                       sleep_P)
        rmins.append(len(np.where(sleepstate == 0)[0])*dt_hrs*num_cycles)
        n1mins.append(len(np.where(sleepstate == -1)[0])*dt_hrs*num_cycles)
        n2mins.append(len(np.where(sleepstate == -2)[0])*dt_hrs*num_cycles)
        n3mins.append(len(np.where(sleepstate == -3)[0])*dt_hrs*num_cycles)
    return rmins, n1mins, n2mins, n3mins

        
def cmpt_prop_sleep(ages):
    dt_mins = 0.001*24*60
    rtimes = []
    n2times = []
    for ii in ages:
        P = fetch_period(ii)
        # num_cycles = int(24/P)  # cycles per night
        d = fetch_duty(ii)
        # print(ii, P, num_cycles, d)
        tt, circ_rhy, clck_24, sleep_P, n3_co, yy = cmpt_fourier(total_width=P/24,
                                                                 duty=d)
        sleepstate, s_n1n2, s_n2n3, s_n3n2, s_nr, s_rn1 = cmpt_ss_durs(yy,
                                                                       n3_co,
                                                                       circ_rhy,
                                                                       sleep_P)
        # rem cycle lengths
        when_rem = []
        for ii in s_nr[0]:
            when_rem.append(tt[ii])
        #when_rem.append(1)
        # non-rem cycle lengths
        when_nrem = []
        for ii in s_n1n2[0]:
            when_nrem.append(tt[ii])
        #when_nrem.append(1)
        rtimes.append(np.array(when_rem))
        n2times.append(np.array(when_nrem))
    return rtimes, n2times




if __name__ == '__main__':
    # cmpt_aging_sleep()
    print('Hello')
    # tt, circ_rhy, clck_24, sleep_P, n3_cutoff, yy = cmpt_fourier(total_width=0.25,
    #                                                              duty=0.6)
    # sleepstate, s_n1n2, s_n2n3, s_n3n2, s_nr, s_rn1 = cmpt_ss_durs(yy,
    #                                                                n3_cutoff,
    #                                                                circ_rhy,
    #                                                                sleep_P)
    # n1times, n2times, n3starttimes, n3endtimes, rtimes = cmpt_time_in(tt,
    #                                                                   s_n1n2,
    #                                                                   s_nr,
    #                                                                   s_n2n3,
    #                                                                   s_n3n2,
    #                                                                   s_rn1)
    # plt.plot(tt, sleep_P)
    # plt.plot(tt, sleepstate)
    # plt.show()

    
    # cmpt_pwr()
    # plot_total_bar(sleepstate)

    # plt.plot(tt, pwm(tt, v, duty))
    # plt.plot(tt, ramp(tt, v))
    # plt.plot(tt, sleep_pressure)
    # plt.plot(tt, ramp_n3(tt, v))

    # print(np.diff(n2times))
    # print(np.diff(rtimes))
    # print(np.diff(n3times))
