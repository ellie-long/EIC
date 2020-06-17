import numpy as np
import scipy as sp
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit

# momLog = True
momLog = False

evtLog = True
# evtLog = False

electronOnly = True
# electronOnly = False

# protonOnly = True
protonOnly = False

dataFilename = './data/elastic_theta_momentum.txt'
# dataFilename = './data/5by41-electron.dat'

# [0]          | [1]               | [2]                | [3]                | [4]
# Angle [Rad]  |  Bin Width [Rad]  |  Momentum [GeV/c]  |  Bin Width[GeV/c]  |  Events/10fb^-1

# Note: All datafile lines that aren't data must have a '!' in them for the code to work properly
# For electron and proton tagging to work, the line before the electron data must include "Electron"
#                                          and the line before the proton data must include "Proton"


def getThDisp(th):
    thVal = th
    negFac = 1
    if (thVal<90) and (thVal>0):
        thDispVal = 90*np.log(thVal+negFac)/np.log(90+negFac)
    elif (thVal>90) and (thVal<180):
        thDispVal = 180-90*np.log(180-thVal+negFac)/np.log(90+negFac)
    elif (thVal>180) and (thVal<270):
        thDispVal = 360-90*np.log(180-thVal+negFac)/np.log(90+negFac)
    elif (thVal>270) and (thVal<360):
        thDispVal = 360-90*np.log(360-thVal+negFac)/np.log(90+negFac)
    else:
        thDispVal = thVal
    thDispVal = 360-thDispVal
    return thDispVal

def fillDataLog(part,thMin,thMax,rVal,r,thDisp):
    rVal = np.log(rVal)/np.log(10)
    degToRad = np.pi/180
    n = 10
    n2 = 10
    thMinIn = thMin
    thMaxIn = thMax
    if part=='e':
        thMin = thMinIn*n2
        thMax = thMaxIn*n2
    elif part =="p":
        thMin = int(thMinIn/6*n2)
        thMax = int(thMaxIn/6*n2)
    for i in range(thMin*n,thMax*n):
        thVal = i/(n2*n)
        thDispVal = getThDisp(thVal)
        thDispVal = 360-thDispVal
        thDisp.append(thDispVal*degToRad)
        r.append(rVal)
    return 1

def fillData(part,thMin,thMax,rVal,r,thDisp):
    rVal = np.log(rVal)/np.log(10)
    degToRad = np.pi/180
    n = 1
    n2 = 1
    thMinIn = thMin
    thMaxIn = thMax
    if part=='e':
        thMin = thMinIn
        thMax = thMaxIn
    elif part =="p":
        thMin = int(thMinIn/6*n2)
        thMax = int(thMaxIn/6*n2)
    for i in range(thMin*n,thMax*n):
        thVal = i/n
        if part=='p':
            thVal = i/(n2*n)
        thDispVal = thVal
        thDisp.append(thDispVal*degToRad)
        r.append(rVal)
    return 1

def fillIncData(part,thMin,thMax,rVal,r,thDisp):
    rValIn = rVal
    degToRad = np.pi/180
    n = 1
    n2 = 1
    thMinIn = thMin
    thMaxIn = thMax
    if part=='e':
        thMin = thMinIn
        thMax = thMaxIn
    elif part =="p":
        thMin = int(thMinIn/6*n2)
        thMax = int(thMaxIn/6*n2)
    for i in range(thMin*n,thMax*n):
        thVal = i/n
        if part=='p':
            thVal = i/(n2*n)
        thDispVal = thVal
        thDisp.append(thDispVal*degToRad)
        if thVal<=90:
            rVal = np.log(thVal*rValIn/90)/np.log(10)
        elif thVal>90:
            rVal = np.log((180-thVal)*rValIn/90)/np.log(10)
        r.append(rVal)
    return 1

def fillIncDataLog(part,thMin,thMax,rVal,r,thDisp):
    rValIn = rVal
    # rVal = np.log(rVal)/np.log(10)
    degToRad = np.pi/180
    n = 10
    n2 = 10
    thMinIn = thMin
    thMaxIn = thMax
    if part=='e':
        thMin = thMinIn*n2
        thMax = thMaxIn*n2
    elif part =="p":
        thMin = int(thMinIn/6*n2)
        thMax = int(thMaxIn/6*n2)
    for i in range(thMin*n,thMax*n):
        thVal = i/(n2*n)
        thDispVal = getThDisp(thVal)
        thDispVal = 360-thDispVal
        thDisp.append(thDispVal*degToRad)
        if thVal<=90:
            rVal = np.log(thVal**2*rValIn/90**2)/np.log(10)
        elif thVal>90:
            rVal = np.log((180-thVal)**2*rValIn/90**2)/np.log(10)
        r.append(rVal)
    return 1


degToRad = np.pi/180
radToDeg = 180/np.pi

data = open(dataFilename,'r')
lines = data.readlines()
data.close()
angle = []
angWidth = []
angMin = []
angMax = []
angleDisp = []
angDispWidth = []
angDispMin = []
angDispMax = []
momentum = []
momWidth = []
momMin = []
momMax = []
momentumLog = []
momWidthLog = []
momMinLog = []
momMaxLog = []
numEvts = []
part = []
partVal = 'u'
maxEvt = 0
minEvt = 1E9
for line in lines:
    p = line.split()
    if "Electron" in line:
        partVal = 'e'
    if "Proton" in line:
        partVal = 'p'
    if ("!" not in line):
        angleVal = float(p[0])*radToDeg
        angWidthVal = float(p[1])*radToDeg
        angMinVal = angleVal - angWidthVal
        angMaxVal = angleVal + angWidthVal
        momentumVal = float(p[2])
        momWidthVal = float(p[3])
        momMinVal = momentumVal - momWidthVal
        momMaxVal = momentumVal + momWidthVal
        momentumLogVal = np.log(momentumVal)/np.log(10)
        momMinLogVal = np.log(momMinVal)/np.log(10)
        momMaxLogVal = np.log(momMaxVal)/np.log(10)
        numEvtsVal = float(p[4])
        if numEvtsVal > maxEvt:
            maxEvt = numEvtsVal
        if numEvtsVal < minEvt:
            minEvt = numEvtsVal
        angDispMinVal = getThDisp(angMaxVal)*degToRad
        angDispMaxVal = getThDisp(angMinVal)*degToRad
        angDispWidthVal = angDispMaxVal - angDispMinVal
        # angleDispVal = getThDisp(angleVal)*degToRad
        angleDispVal = angDispMinVal + (angDispWidthVal/2)
        # print(angWidthVal,angDispWidthVal*radToDeg)
        if (electronOnly and protonOnly):
            print("Error! Can't have both electronOnly and protonOnly set to True!")
            exit()
        if (electronOnly and (partVal == 'e')):
            angle.append(angleVal*degToRad)
            angleDisp.append(angleDispVal)
            angWidth.append(angWidthVal*degToRad)
            angDispWidth.append(angDispWidthVal)
            angMin.append(angMinVal)
            angDispMin.append(getThDisp(angMinVal)*degToRad)
            angMax.append(angMaxVal)
            angDispMax.append(getThDisp(angMaxVal)*degToRad)
            momentum.append(momentumVal)
            momWidth.append(momWidthVal)
            momMin.append(momMinVal)
            momMax.append(momMaxVal)
            momentumLog.append(momentumLogVal)
            momMinLog.append(momMinLogVal)
            momMaxLog.append(momMaxLogVal)
            momWidthLog.append(momMaxLogVal-momMinLogVal)
            numEvts.append(numEvtsVal)
            part.append(partVal)
        if (protonOnly and (partVal == 'p')):
            angle.append(angleVal*degToRad)
            angleDisp.append(angleDispVal)
            angWidth.append(angWidthVal*degToRad)
            angDispWidth.append(angDispWidthVal)
            angMin.append(angMinVal)
            angDispMin.append(getThDisp(angMinVal)*degToRad)
            angMax.append(angMaxVal)
            angDispMax.append(getThDisp(angMaxVal)*degToRad)
            momentum.append(momentumVal)
            momWidth.append(momWidthVal)
            momMin.append(momMinVal)
            momMax.append(momMaxVal)
            momentumLog.append(momentumLogVal)
            momMinLog.append(momMinLogVal)
            momMaxLog.append(momMaxLogVal)
            momWidthLog.append(momMaxLogVal-momMinLogVal)
            numEvts.append(numEvtsVal)
            part.append(partVal)
        if ((not electronOnly) and (not protonOnly)):
            angle.append(angleVal*degToRad)
            angleDisp.append(angleDispVal)
            angWidth.append(angWidthVal*degToRad)
            angDispWidth.append(angDispWidthVal)
            angMin.append(angMinVal)
            angDispMin.append(getThDisp(angMinVal)*degToRad)
            angMax.append(angMaxVal)
            angDispMax.append(getThDisp(angMaxVal)*degToRad)
            momentum.append(momentumVal)
            momWidth.append(momWidthVal)
            momMin.append(momMinVal)
            momMax.append(momMaxVal)
            momentumLog.append(momentumLogVal)
            momMinLog.append(momMinLogVal)
            momMaxLog.append(momMaxLogVal)
            momWidthLog.append(momMaxLogVal-momMinLogVal)
            numEvts.append(numEvtsVal)
            part.append(partVal)
alphas = []
alphasLog = []
colorMap = [[],[]]
cFactor = maxEvt - minEvt
minEvtLog = np.log(minEvt)/np.log(10)
maxEvtLog = np.log(maxEvt)/np.log(10)
cFactorLog = np.log(maxEvt)/np.log(10)-np.log(minEvt)/np.log(10)
colorMap = np.empty((2,len(numEvts)))
colorMapLog = np.empty((2,len(numEvts)))
for i in range(0,len(numEvts)):
    alphaVal = (numEvts[i]-minEvt)/cFactor
    alphaLogVal = (np.log(numEvts[i])/np.log(10)-np.log(minEvt)/np.log(10))/cFactorLog
    alphas.append(alphaVal)
    alphasLog.append(alphaLogVal)
    colorMap[0][i] = numEvts[i]-minEvt
    colorMap[1][i] = alphaVal
    colorMapLog[0][i] = np.log(numEvts[i])/np.log(10)-np.log(minEvt)/np.log(10)
    colorMapLog[1][i] = alphaLogVal


cmapDisp = cm.get_cmap('viridis')
fig = plt.figure()
ax = plt.subplot(111, projection='polar')

if momLog:
    ax.set_rmax(math.ceil(max(momMaxLog)))
    ax.set_rticks(np.linspace(0,math.ceil(max(momMaxLog)),math.ceil(max(momMaxLog)), endpoint=False))  # Less radial ticks
    if evtLog:
        ax.bar(angle,momWidthLog,width=angWidth,bottom=momMinLog,color=cmapDisp(alphasLog))
        ax.bar(angleDisp,momWidthLog,width=angDispWidth,bottom=momMinLog,color=cmapDisp(alphasLog))
    else:
        ax.bar(angle,momWidthLog,width=angWidth,bottom=momMinLog,color=cmapDisp(alphas))
        ax.bar(angleDisp,momWidthLog,width=angDispWidth,bottom=momMinLog,color=cmapDisp(alphas))
    ax.set_rmin(0)
    ax.set_rmax(math.ceil(max(momMaxLog)))
    if protonOnly:
        ax.set_rticks(np.linspace(0,math.ceil(max(momMaxLog)),2*math.ceil(max(momMaxLog)), endpoint=False))  # Less radial ticks
        ax.set_rmin(min(momMaxLog)-0.2)
        ax.set_rmax(max(momMaxLog)+0.2)
        ax.text(0.08, (max(momMaxLog)+0.2)*(1.01), 'log$_{10}$E\'\n(log(GeV))', fontsize=12)
        ax.text(2.5, (max(momMaxLog)+0.2)*(1.15), 'Linear\nAngle', fontsize=12)
        ax.text(3.75, (max(momMaxLog)+0.2)*(1.15), 'Log\nAngle', fontsize=12)
    else:
        ax.text(0.08, math.ceil(max(momMaxLog))*(1.01), 'log$_{10}$E\'\n(log(GeV))', fontsize=12)
        ax.text(2.5, math.ceil(max(momMaxLog))+1, 'Linear\nAngle', fontsize=12)
        ax.text(3.75, math.ceil(max(momMaxLog))+1, 'Log\nAngle', fontsize=12)


else:
    ax.set_rmax(math.ceil(max(momMax)))
    # ax.set_rticks(np.linspace(0,math.ceil(max(momMax)),math.ceil(max(momMax))/50, endpoint=False))  # Less radial ticks
    if evtLog:
        ax.bar(angle,momWidth,width=angWidth,bottom=momMin,color=cmapDisp(alphasLog))
        ax.bar(angleDisp,momWidth,width=angDispWidth,bottom=momMin,color=cmapDisp(alphasLog))
    else:
        ax.bar(angle,momWidth,width=angWidth,bottom=momMin,color=cmapDisp(alphas))
        ax.bar(angleDisp,momWidth,width=angDispWidth,bottom=momMin,color=cmapDisp(alphas))
    ax.set_rmin(0)
    ax.set_rmax(math.ceil(max(momMax)))
    ax.text(2.5, math.ceil(max(momMax))*1.3, 'Linear\nAngle', fontsize=12)
    ax.text(3.75, math.ceil(max(momMax))*1.3, 'Log\nAngle', fontsize=12)
    ax.text(0.08, math.ceil(max(momMax))*(1.01), 'E\' (GeV)', fontsize=12)


ax.set_thetamin(0)
ax.set_thetamax(180)
# plt.xticks(np.pi/180. * np.linspace(360,  0, 4*4, endpoint=False),['0.011\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '$90$\N{DEGREE SIGN}', '170.51\N{DEGREE SIGN}', '179\N{DEGREE SIGN}', '179.895\N{DEGREE SIGN}', '179.989\N{DEGREE SIGN}', '179.895\N{DEGREE SIGN}', '179\N{DEGREE SIGN}', '170.51\N{DEGREE SIGN}', '90\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}','0.011\N{DEGREE SIGN}'])
plt.xticks(np.pi/180. * np.linspace(360,  0, 4*4, endpoint=False),['0.011\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '$90$\N{DEGREE SIGN}', '170.51\N{DEGREE SIGN}', '179\N{DEGREE SIGN}', '179.895\N{DEGREE SIGN}', '179.989\N{DEGREE SIGN}', '157.5\N{DEGREE SIGN}', '135\N{DEGREE SIGN}', '112.5\N{DEGREE SIGN}', '90\N{DEGREE SIGN}', '67.5\N{DEGREE SIGN}', '45\N{DEGREE SIGN}', '22.5\N{DEGREE SIGN}','0\N{DEGREE SIGN}'])
# plt.xticks(np.pi/180. * np.linspace(180,  0, 2*4+1, endpoint=True),['179.989\N{DEGREE SIGN}', '179.895\N{DEGREE SIGN}', '179\N{DEGREE SIGN}', '170.51\N{DEGREE SIGN}', '90\N{DEGREE SIGN}', '9.49\N{DEGREE SIGN}', '1\N{DEGREE SIGN}', '0.105\N{DEGREE SIGN}','0.011\N{DEGREE SIGN}'])
# ax.set_thetalim(-np.pi, np.pi)
ax.set_rlabel_position(0)  # Move radial labels away from plotted line

ax.set_title("Elastic e-p Scattering", va='bottom')

ax3 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
if evtLog:
    ax3.text(-10, maxEvtLog*(1+0.02), '10^n Events', fontsize=12)
    cNorm = mpl.colors.Normalize(vmin=minEvtLog, vmax=maxEvtLog)
else:
    ax3.text(-10, maxEvt*(1+0.05), '# Events', fontsize=12)
    cNorm = mpl.colors.Normalize(vmin=minEvt, vmax=maxEvt)
cb1 = mpl.colorbar.ColorbarBase(ax3, norm=cNorm)
plt.show()



# plt.xticks(np.arange(0, 360, 45))
# plt.xticklabels(['$90^{-1}~^{\deg}$', '$90^{-1/2}~^{\deg}$', '$90^{0}~^{\deg}$', '$90^{1/2}~^{\deg}$', '$90^{1}~^{\deg}$', '', 'E', ''])
# plt.ylim(1E-1,500)
# plt.yscale('log')
# plt.polar(th, r, 'b.')
