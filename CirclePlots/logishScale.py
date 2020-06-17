import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def getThDisp(th):
    thVal = th
    negFac = 1
    if (thVal<90) and (thVal>0):
        thDispVal = 90*np.log(thVal+negFac)/np.log(90+negFac)
    elif (thVal>90) and (thVal<180):
        thDispVal = 180-90*np.log(180-thVal+negFac)/np.log(90+negFac)
    elif (thVal>180) and (thVal<270):
        thDispVal = 180+90*np.log(270-thVal+negFac)/np.log(90+negFac)
    elif (thVal>270) and (thVal<360):
        thDispVal = 360-90*np.log(360-thVal+negFac)/np.log(90+negFac)
    else:
        thDispVal = thVal
    return thDispVal

def fillData(part,thMin,thMax,rVal,r,thDisp):
    rVal = np.log(rVal)/np.log(10)
    degToRad = np.pi/180
    n = 1
    n2 = 10
    thMinIn = thMin
    thMaxIn = thMax
    if part=='e':
        thMin = thMinIn
        thMax = thMaxIn
    elif part =="p":
        # thMax = 360-thMin
        # thMin = 360-thMax
        thMin = int(thMinIn/6*n2)
        thMax = int(thMaxIn/6*n2)
        # thMin = 360*n2-int(thMaxIn/6*n2)
        # thMax = 360*n2-int(thMinIn/6*n2)
        # print (thMin,360*n2-thMinIn,thMax,360*n2-thMaxIn)
    # else:
    #     print("Error!")
    #     return 0
    print(thMin,thMax)
    for i in range(thMin*n,thMax*n):
        # thVal = (i/n)*90
        thVal = i/n
        if part=='p':
            thVal = i/(n2*n)
            # thVal = i/(n*6)
            print(rVal,thVal)
        thDispVal = getThDisp(thVal)
        # th.append(thVal*degToRad)
        thDisp.append(thDispVal*degToRad)
        r.append(rVal)
    return 1



thp = []
thpDisp = []
rp = []

the = []
theDisp = []
re = []

rVal = 18
part = 'e'
thMin = 170
thMax = 178
fillData(part,thMin,thMax,rVal,re,theDisp)

rVal = 10
part = 'e'
thMin = 164
thMax = 176
fillData(part,thMin,thMax,rVal,re,theDisp)


rVal = 5
part = 'e'
thMin = 109
thMax = 170
fillData(part,thMin,thMax,rVal,re,theDisp)


rVal = 275
part = 'p'
thMin = 0 # 0
thMax = 4 # 0.67
fillData(part,thMin,thMax,rVal,rp,thpDisp)

rVal = 100
part = 'p'
thMin = 2 # 0.33
thMax = 10 # 1.67
fillData(part,thMin,thMax,rVal,rp,thpDisp)

rVal = 41
part = 'p'
thMin = 7 # 1.17
thMax = 17 # 2.83
fillData(part,thMin,thMax,rVal,rp,thpDisp)

rVal = 40
part = 'p'
thMin = 15 # 2.5
thMax = 65 # 10.83
fillData(part,thMin,thMax,rVal,rp,thpDisp)

rVal = 39
part = 'p'
thMin = 70
thMax = 77
fillData(part,thMin,thMax,rVal,rp,thpDisp)

# # fig = plt.figure()
# # ax = fig.add_subplot(111, projection='polar')
# # c = ax.scatter(theDisp, re, 'b.', label='electron')
# # c = ax.scatter(thpDisp, rp, 'r.', label='electron')
# plt.axes(polar='true')
# plt.plot(theDisp, re, 'b.',label='electron')
# plt.plot(thpDisp, rp, 'r.',label='proton')
# # plt.xticks(np.pi/180. * np.linspace(360,  0, 4*4, endpoint=False),['$90^{-1}$\N{DEGREE SIGN}', '$90^{-1/2}$\N{DEGREE SIGN}', '$90^{0}$\N{DEGREE SIGN}', '$90^{1/2}$\N{DEGREE SIGN}', '$90^{1}$\N{DEGREE SIGN}', '180-$90^{1/2}$\N{DEGREE SIGN}', '180-$90^{0}$\N{DEGREE SIGN}', '180-$90^{-1/2}$\N{DEGREE SIGN}', '180-$90^{-1}$\N{DEGREE SIGN}', '180-$90^{-1/2}$\N{DEGREE SIGN}', '180-$90^{0}$\N{DEGREE SIGN}', '180-$90^{1/2}$\N{DEGREE SIGN}', '$90^{1}$\N{DEGREE SIGN}', '$90^{1/2}$\N{DEGREE SIGN}', '$90^{0}$\N{DEGREE SIGN}', '$90^{-1/2}$\N{DEGREE SIGN}'])
# plt.xticks(np.pi/180. * np.linspace(360,  0, 4*4, endpoint=False),['0.011\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '$90$\N{DEGREE SIGN}', '90.11\N{DEGREE SIGN}', '91\N{DEGREE SIGN}', '99.49\N{DEGREE SIGN}', '180\N{DEGREE SIGN}', '99.49\N{DEGREE SIGN}', '91\N{DEGREE SIGN}', '90.11\N{DEGREE SIGN}', '90\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}'])
# plt.yticks(np.linspace(0,5,6),['$10^0$','$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])
# ax.set_rlabel_position(135)
# plt.set_rlabel_position(0)

ax = plt.subplot(111, projection='polar')
ax.plot(theDisp, re, 'b.',label='electron')
ax.plot(thpDisp, rp, 'r.',label='proton')
ax.set_thetamin(0)
ax.set_thetamax(180)
# plt.xticks(np.pi/180. * np.linspace(360,  0, 4*4, endpoint=False),['0.011\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '$90$\N{DEGREE SIGN}', '170.51\N{DEGREE SIGN}', '179\N{DEGREE SIGN}', '179.895\N{DEGREE SIGN}', '179.989\N{DEGREE SIGN}', '179.895\N{DEGREE SIGN}', '179\N{DEGREE SIGN}', '170.51\N{DEGREE SIGN}', '90\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}','0.011\N{DEGREE SIGN}'])
plt.xticks(np.pi/180. * np.linspace(180,  0, 2*4+1, endpoint=True),['179.989\N{DEGREE SIGN}', '179.895\N{DEGREE SIGN}', '179\N{DEGREE SIGN}', '170.51\N{DEGREE SIGN}', '90\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}','0.011\N{DEGREE SIGN}'])
# ax.set_rmax(2)
# ax.set_thetalim(-np.pi, np.pi)
ax.set_rticks(np.linspace(0,2.5,6, endpoint=True))  # Less radial ticks
ax.set_rlabel_position(0)  # Move radial labels away from plotted line
ax.text(-0.4, 2.5, 'log$_{10}$E\n(log(GeV))', fontsize=12)
# ax.grid(True)

ax.set_title("Elastic e-p Scattering", va='bottom')
plt.legend()
plt.show()


# plt.xticks(np.arange(0, 360, 45))
# plt.xticklabels(['$90^{-1}~^{\deg}$', '$90^{-1/2}~^{\deg}$', '$90^{0}~^{\deg}$', '$90^{1/2}~^{\deg}$', '$90^{1}~^{\deg}$', '', 'E', ''])
# plt.ylim(1E-1,500)
# plt.yscale('log')
# plt.polar(th, r, 'b.')
