
import numpy as np
import matplotlib.pyplot as plt

lls = [1, 2, 3, 4, 6, 8, 12, 24]  # First 2 yrs of age
ttx = np.arange(0, 24, 0.01)  # time in hours


def angle_to_hrs(ang):
    # 0 -> 0hrs; 2*pi -> 24hrs
    return 24*ang / (2*np.pi)


def fetch_freq(year):
    # number of hrs in a sleep-awake cycle
    if year < 2:
        return lls[int(year/0.25)]  # full cycle per 24 hrs
    else:
        return 24


def fetch_duty(year):
    # correction for the amount of overall sleep per day over age
    # 0yr -> 16 hrs sleep per 24hrs  D(A) = 0.33
    # 2yr -> 12 hrs sleep per 24hrs  D(A) = 0.5
    # 16yr -> 8 hrs sleep per 24hrs  D(A) = 0.66
    # 80yr -> 6 hrs sleep per 24hrs  D(A) = 0.75
    if year <= 2:
        return 0.33333 + (year)*(0.1666667/2)
    elif year <= 16:
        return 0.5 + (year-2)*(0.16666/14)
    else:
        return 0.6666 + (year-16)*(0.083333/64)

    
def fetch_bias(year):
    # correction for the amount of overall sleep per day over age
    # 0yr -> 16 hrs per 24hrs
    # 2yr -> 12 hrs per 24hrs
    # 16yr -> 8 hrs per 24hrs
    # 80yr -> 6 hrs per 24hrs
    if year <= 2:
        return -0.5 + 0.5*(year/2)
    elif year <= 16:
        return (year-2)*(0.427)/14
    else:
        return 0.427 + ((year-16)*(0.28)/(80-16))


def kth_term(freq, k, ttx=ttx):
    # kth term of the fourier mode for a sqare wave
    # https://en.wikipedia.org/wiki/Square_wave
    ff = 2*k - 1
    return np.sin(2*np.pi*ff*(ttx)/freq) / ff


def find_awake_sleep_frac(f_vals, debug=True):
    awake_time = len(np.where(f_vals > 0)[0])  # awake
    sleep_time = len(f_vals) - awake_time
    frac_awake = awake_time/(awake_time + sleep_time)
    if debug:
        print('Frac of awake/sleep:', frac_awake)
        print('Sleep hours: ', (1-frac_awake)*24)
    return frac_awake


def find_rem_sleep_frac(f_vals, rem_thres, s_vals, debug=False):
    sleep_idx = np.where(s_vals < 0)  # when in sleep cycle
    rem_idx = np.where(f_vals[sleep_idx] >= rem_thres[sleep_idx])  # if in rem
    sleep_time = len(sleep_idx[0])
    rem_time = len(rem_idx[0])
    frac_rem = rem_time/(sleep_time + rem_time)
    if debug:
        print('Frac of sleep/REM:', frac_rem)
    return frac_rem


def start_with_awake(f_mode, year, bias):
    # move it so that it begins with awake mode
    zero_crossings = np.where(np.diff(np.sign(f_mode + bias)))[0]
    if len(zero_crossings) == 1:  # may happen due to tt least count
        correction_shift = 0
    elif year <= 2:
        correction_shift = zero_crossings[0] + 1
    else:
        correction_shift = zero_crossings[-1]
    vals = np.roll(f_mode + bias, -1*correction_shift)
    return vals


def compute_sine_sleep(year, ttx=ttx, debug=False):
    bias = fetch_bias(year)   # hours per sleep cycle
    freq = fetch_freq(year)  # fixes amount of sleep per day
    offset = np.arcsin(bias)
    s_vals = np.sin((2*ttx*np.pi/freq) - offset) + bias
    if debug:
        print('Age (yrs): ', year)
        plt.plot(ttx, s_vals)  # starts with awake
    return s_vals


def compute_fourier_sleep(year, ttx=ttx, debug=False):
    ''' this is just fourier transform with the 5th term amplified
    for better version consider FT of pwm,  (TODO?)
    https://liraeletronica.weebly.com/uploads/4/9/3/5/4935509/spectral_analysis_of_a_pwm_signal.pdf
    '''
    freq = fetch_freq(year)  # fixes amount of sleep bouts per day
    bias = fetch_bias(year)   # hours per sleep cycle
    s_vals = compute_sine_sleep(year, ttx, debug=debug)
    # # REM mode of sleep - the first 5 fourier modes of a square wave
    f_mode = (4/np.pi)*(kth_term(freq, 1, ttx=ttx) + kth_term(freq, 2, ttx=ttx) +
                        kth_term(freq, 3, ttx=ttx) + kth_term(freq, 4, ttx=ttx) +
                        3*kth_term(freq, 5, ttx=ttx))
    f_vals = start_with_awake(f_mode, year, bias)
    # a threshold for when it is in rem mode
    rem_thres = np.zeros_like(ttx) - 0.75 + bias + ((bias < 0)*1.5*bias)
    if debug:
        print('Age (yrs): ', year)
        plt.plot(ttx, f_vals, 'r.')
        plt.plot(ttx, [rem_thres]*len(ttx), 'k')
    return s_vals, f_vals, rem_thres


if __name__ == '__main__':
    # year = 0.5
    # compute_sine_sleep(year)

    year = 12

    print(fetch_bias(year))
    
    # s_vals, f_vals, rem_thres = compute_fourier_sleep(year, debug=True)
    # frac_awake = find_awake_sleep_frac(s_vals, debug=False)
    # frac_rem = find_rem_sleep_frac(f_vals, rem_thres, s_vals, debug=False)
    # sl_hrs = 24 - (frac_awake*24)
    # print('Sleep hrs: ', sl_hrs)
    # print('REM hrs: ', sl_hrs*frac_rem)


    # ages = np.arange(0, 100, 0.1)
    # sleephrs = []
    # remhrs = []
    # nremhrs = []
    # biass = []
    # rem_thress = []
    # for ii in ages:
    #     s_vals, f_vals, rem_thres = compute_fourier_sleep(ii, debug=False)
    #     frac_awake = find_awake_sleep_frac(s_vals, debug=False)
    #     frac_rem = find_rem_sleep_frac(f_vals, rem_thres, s_vals, debug=False)
    #     sleephrs.append(24 - (frac_awake*24))
    #     remhrs.append(sleephrs[-1]*frac_rem)
    #     nremhrs.append(sleephrs[-1]-remhrs[-1])
    #     biass.append(fetch_bias(ii))
    #     rem_thress.append(rem_thres)
    # plt.plot(np.array(ages), nremhrs)
    # plt.plot(np.array(ages), sleephrs)
    # # biass = np.array(biass)
    # # plt.plot(np.array(ages), biass)
    # # plt.plot(np.array(ages), np.array(rem_thress))
    # plt.xlim(0, 100)
    # plt.ylim(0, 24)             
    # plt.ylabel('Sleep hours')
    # plt.xlabel('Age (yrs')

    plt.show()
