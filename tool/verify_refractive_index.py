#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

from ErwinJr2.material import rIdx


def testGaAs():
    # Mid IR refractive index for GaAs as a function of wavelength in um _gaas
    index = np.array([
        [2.066, 3.3378], [2.254, 3.3306], [2.480, 3.3240], [2.755, 3.3180],
        [3.100, 3.3125], [3.542, 3.3075], [4.133, 3.3027], [4.275, 3.3017],
        [4.428, 3.3008], [4.592, 3.2998], [4.768, 3.2988], [4.959, 3.2978],
        [5.166, 3.2968], [5.391, 3.2954], [5.636, 3.2946], [5.904, 3.2934],
        [6.199, 3.2921], [6.526, 3.2907], [6.888, 3.2891], [7.293, 3.2874],
        [7.749, 3.2854], [8.266, 3.2831], [8.856, 3.2803], [9.537, 3.2769],
        [10.33, 3.2727], [11.27, 3.2671], [12.40, 3.2597], [13.78, 3.2493],
        [15.50, 3.2336], [17.71, 3.2081]
    ])
    x = np.linspace(0.97, 19, 100)
    plt.plot(x, rIdx['GaAs'](x))
    plt.plot(index[:, 0], index[:, 1], '*')


def testInAs():
    index = np.array([
        [3.100, 3.539], [3.757, 3.497], [3.875, 3.491], [4.000, 3.485],
        [4.133, 3.479], [4.275, 3.473], [4.428, 3.467], [4.592, 3.461],
        [4.769, 3.457], [4.959, 3.451], [5.166, 3.445], [5.391, 3.441],
        [5.636, 3.436], [5.904, 3.432], [6.199, 3.427], [6.526, 3.423],
        [6.888, 3.420], [7.293, 3.416], [7.749, 3.412], [8.266, 3.408],
        [8.856, 3.404], [9.537, 3.400], [10.33, 3.394], [11.27, 3.388],
        [12.40, 3.381], [13.78, 3.372], [15.50, 3.358], [17.71, 3.338]
    ])
    x = np.linspace(3, 20, 100)
    plt.plot(x, rIdx['InAs'](x))
    plt.plot(index[:, 0], index[:, 1], '*')


def testInP():
    index = np.array([
        [2.066, 3.129], [2.254, 3.121], [2.479, 3.114], [2.754, 3.107],
        [3.099, 3.101], [3.541, 3.095], [4.131, 3.089], [4.958, 3.083],
        [6.197, 3.074], [6.886, 3.069], [7.746, 3.062], [8.853, 3.053],
        [10.33, 3.038], [12.39, 3.012], [13.77, 2.990], [14.58, 2.975]
    ])
    x = np.linspace(2, 15, 100)
    plt.plot(x, rIdx['InP'](x))
    plt.plot(index[:, 0], index[:, 1], '*')


def testAlAs():
    x = np.linspace(0.56, 10, 100)
    plt.plot(x, rIdx['AlAs'](x))
    # plt.plot(x, rIdx['AlGaAs'](x, 1))


def testAu():
    realindex = np.array([
        [3, 0.704], [4, 1.25], [5, 1.95], [6, 2.79], [7, 3.79], [8, 4.93],
        [10, 7.62], [12, 10.8], [14, 14.5], [16, 18.7], [18, 23.3],
    ])
    imagindex = np.array([
        [3, 21.8], [4, 29], [5, 36.2], [6, 43.4], [7, 50.5], [8, 57.6],
        [10, 71.5], [12, 85.2], [14, 98.6], [16, 112], [18, 124],
    ])
    x = np.linspace(2.5, 20, 100)
    index = rIdx['Au'](x)
    plt.plot(x, index.real, color='C0', label="n")
    plt.plot(realindex[:, 0], realindex[:, 1], '*', color='C0')
    plt.plot(x, index.imag, color='C1', label="k")
    plt.plot(imagindex[:, 0], imagindex[:, 1], '*', color='C1')
    plt.legend()


def testSi():
    x = np.linspace(2, 20, 100)
    n = rIdx['SiNx'](x)
    plt.plot(x, n.real, label="n")
    plt.plot(x, n.imag, label="k")
    plt.legend()
    plt.ylim(-1, 5)


if __name__ == '__main__':
    plt.figure().gca().set_title("GaAs")
    testGaAs()
    plt.figure().gca().set_title("InAs")
    testInAs()
    plt.figure().gca().set_title("InP")
    testInP()
    plt.figure().gca().set_title("AlAs")
    testAlAs()
    plt.figure().gca().set_title("Au")
    testAu()
    plt.figure().gca().set_title("SiNx")
    testSi()
    plt.show()
