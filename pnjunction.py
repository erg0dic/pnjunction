# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 07:02:38 2019

@author: Irtaza

Code for generating plots from collected data
"""
import numpy as np
from scipy.stats import linregress
#from scipy.stats import polyfit
import matplotlib.pyplot as plt

def idf(slope, T=293):
    return 1.6e-19/(slope*T*1.38e-23)




def id_compute_pseudo(data, low, high, voltl, volth, error, T, plot=True):
    "computes ideality given data array, slice paramenters and temperature"
    vals = linregress(data[:,0][low:high], np.log(data[:,1][low:high]))
    if plot==True:
#        plt.figure(low)
#        plt.semilogy(data[:,0], data[:,1], '.')
#        plt.semilogy(data[:,0][low:high], data[:,1][low:high], '.')
        plt.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), '-', 
                     label = r'$T= {0} K, n = {1} \pm {2}$'.format(T, round(idf(vals[0], T), 3), round(error*idf(vals[0], T), 3)))
        plt.xlabel('V / V')
        plt.ylabel(r'$ \log{I}$ / A')
        plt.legend(loc=0, fontsize=14)
    return idf(vals[0], T)


def optimizer(filename, temperature, resolution, nid_limitl=0, nid_limitu=100, inter_points=5):
    "automatically computes a range of values for n_id profile"
    data = np.loadtxt(filename)
    orig_len = len(data)
    
    data = data[[i for i, item in enumerate(data[:,1]) if item > 0]] #  get rid of negative values
    
    offset = orig_len - len(data)
    
    nids = []
    low = 0
    high = inter_points
    while high < len(data):
        vals = linregress(data[:,0][low:high], np.log(data[:,1][low:high]))
        
        slope = vals[0]
        intercept = vals[1]
        
        nid = idf(slope, T=temperature)
        sat_current = np.exp(intercept)
        
        
        if nid < nid_limitu and nid > nid_limitl:
            nids.append([nid, sat_current, low+offset, high+offset])
        
        low += resolution
        high += resolution
        
    return nids
    
ger = ["g_2681.txt", "g_2707.txt", "g_2732.txt", "g_2748.txt", "g_2831.txt", "g_3012.txt", "g_3153.txt"  ]    
    
ger_temps = [268.1, 270.7, 273.2, 274.8, 283.1, 301.2, 315.3]    
    
 
nids_ger = []   
for i in range(len(ger)):
    nids_ger.append(optimizer(ger[i], ger_temps[i], 1, 1, 1.5, 5))
    
    
    
points = [[1.4070139160940418, 2.9647672790925197e-07, 1, 6],
          [1.1392530345771246, 2.1006694397411531e-07, 9, 14],
          [1.1933457684217084, 7.632029572686902e-13, 104, 109],
          [1.2765428944479726, 3.434695151460126e-07, 20, 25],
          [1.0950976851280503, 1.2289620336482655e-06, 16, 21],
          [1.131751362398514, 1.2127407207597237e-06, 28, 33]]

   
egs = []
   
led4 = ["led4_264.3.txt", "led4_268.0.txt", "led4_274.1.txt", "led4_278.3.txt", "led4_282.5.txt", 
        "led4_288.3.txt", "led4_338.txt", "led4_2976.txt", "led4_3084.txt", "led4_3184.txt", "led4_3283.txt", 
        "led4_3478.txt", "led4_3573.txt", "led4_3706.txt" ]

led4_temps = [264.3, 268, 274.1, 278.3, 282.5, 288.3, 338, 297.6, 308.4, 318.4, 328.3, 347.8, 357.3, 370.6 ]


nids_led4 = []   
for i in range(len(led4)):
    nids_led4.append(optimizer(led4[i], led4_temps[i], 1, 1, 1.5, 10))

points = np.array([nids_led4[i][0][1] for i in range(len(nids_led4))])
nids = np.array([nids_led4[i][0][0] for i in range(len(nids_led4))])
val = linregress(1/(nids*np.array(led4_temps)), np.log(points))
egs.append(val)

fig, ax = plt.subplots()
ax.semilogy(1/(nids*np.array(led4_temps)), points, '.')
ax.semilogy(np.linspace(0.0026, 0.0038, 1000),np.exp(np.linspace(0.0026, 0.0038, 1000)*val[0]+val[1]), 
            label = r'Best fit band gap $E_g = {0} \pm {1}$ eV'.format(round(-val[0]*1.38e-23/1.6e-19, 3), 
                                               0.1*round(-val[0]*1.38e-23/1.6e-19, 0) ))
ax.set_xlabel(r'$1/\eta T$ / K$^{-1}$')
ax.set_ylabel(r'$ \log{I_0}$ / A')
ax.grid(True,which="both",ls="-", color='0.75')
ax.legend(loc=0, framealpha=0.8, prop={'size':10})
ax.errorbar(1/(nids*np.array(led4_temps)), points, xerr = 0.01*1/(((nids**2)*np.array(led4_temps))), fmt='.')
ax.set_title(r"Log Saturation current $I_0$ vs 1/$\eta$T for LED 4")

led2 = ["led2_264.5.txt", "led2_267.4.txt", "led2_272.5.txt", "led2_279.3.txt", "led2_282.2.txt", 
        "led2_287.9.txt", "led2_338.txt", "led2_2976.txt", "led2_3084.txt", "led2_3185.txt", "led2_3283.txt",
        "led2_3478.txt", "led2_3573.txt", "led2_3706.txt"]

led2_temps = [264.5, 267.4, 272.5, 279.3, 282.2, 287.9, 338, 297.6, 308.4, 318.5, 328.3, 347.8, 357.3, 370.6]
    
nids_led2 = [] 
for i in range(len(led2)):
    nids_led2.append(optimizer(led2[i], led2_temps[i], 1, 1, 1.5, 10))

points = np.array([nids_led2[i][0][1] for i in range(len(nids_led2))])
nids = np.array([nids_led2[i][0][0] for i in range(len(nids_led2))])
val = linregress(1/(nids*np.array(led2_temps)), np.log(points))

fig, ax = plt.subplots()
ax.semilogy(1/(nids*np.array(led2_temps)), points, '.')
ax.semilogy(np.linspace(0.0026, 0.0038, 1000),np.exp(np.linspace(0.0026, 0.0038, 1000)*val[0]+val[1]), 
            label = r'Best fit band gap $E_g = {0} \pm {1}$ eV'.format(round(-val[0]*1.38e-23/1.6e-19, 3), 
                                               0.12*round(-val[0]*1.38e-23/1.6e-19, 0) ))
ax.set_xlabel(r'$1/\eta T$ / K$^{-1}$')
ax.set_ylabel(r'$ \log{I_0}$ / A')
ax.grid(True,which="both",ls="-", color='0.75')
ax.legend(loc=0, framealpha=0.8, prop={'size':10})
ax.errorbar(1/(nids*np.array(led2_temps)), points, xerr = 0.01*1/(((nids**2)*np.array(led2_temps))), fmt='.')
ax.set_title(r"Log Saturation current $I_0$ vs 1/$\eta$T for LED 2")



led6 = ["led6_267.5.txt", "led6_274.1.txt", "led6_279.7.txt", "led6_282.5.txt", "led6_287.5.txt", "led6_338.txt"
        , "led6_2643.txt", "led6_2976.txt", "led6_3085.txt", "led6_3185.txt", "led6_3283.txt",
        "led6_3477.txt", "led6_3572.txt", "led6_3703.txt"]

led6_temps = [267.5, 274.1, 279.7, 282.5, 287.5, 338, 264.3, 297.6, 308.5, 318.5, 328.3, 347.7, 357.2, 370.3]
    

nids_led6 = [] 
for i in range(len(led6)):
    nids_led6.append(optimizer(led6[i], led6_temps[i], 1, 1, 1.5, 10))
    
points = np.array([nids_led6[i][0][1] for i in range(len(nids_led6))])
nids = np.array([nids_led6[i][0][0] for i in range(len(nids_led6))])
val = linregress(1/(nids*np.array(led6_temps)), np.log(points))

fig, ax = plt.subplots()
ax.semilogy(1/(nids*np.array(led6_temps)), points, '.')
ax.semilogy(np.linspace(0.0026, 0.0038, 1000),np.exp(np.linspace(0.0026, 0.0038, 1000)*val[0]+val[1]), 
            label = r'Best fit band gap $E_g = {0} \pm {1}$ eV'.format(round(-val[0]*1.38e-23/1.6e-19, 3), 
                                               0.12*round(-val[0]*1.38e-23/1.6e-19, 0) ))
ax.set_xlabel(r'$1/\eta T$ / K$^{-1}$')
ax.set_ylabel(r'$ \log{I_0}$ / A')
ax.grid(True,which="both",ls="-", color='0.75')
ax.legend(loc=0, framealpha=0.8, prop={'size':10})
ax.errorbar(1/(nids*np.array(led6_temps)), points, xerr = 0.01*1/(((nids**2)*np.array(led6_temps))), fmt='.')
ax.set_title(r"Log Saturation current $I_0$ vs 1/$\eta$T for LED 6")

for vals in egs:
    print(vals[0]*1.38e-23/1.6e-19)


from mpl_toolkits.axes_grid1.inset_locator import inset_axes    
    
def id_compute(data, low, high, voltl, volth, error, T, plot=True):
    "computes ideality given data array, slice paramenters and temperature"
    vals = linregress(data[:,0][low:high], np.log(data[:,1][low:high]))
    if plot==True:
#        plt.figure(low)
#        plt.semilogy(data[:,0], data[:,1], '.')
#        plt.semilogy(data[:,0][low:high], data[:,1][low:high], '.')
        ax.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'k-', linewidth=2.5
                     )
#    axins.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'b-', linewidth=2.5,
 #                     )
    #axins.tick_params(labelleft=False, labelbottom=False)
    return [voltl, volth, vals[0], vals[1]]

### germanium ####

semis = ["g_803.txt",
"g_2707.txt",
"germanium_29444444.txt",
"g_3012.txt",
"g_3153.txt"]
#plt.figure(100)
fig, ax = plt.subplots()
opts = []
i0_pts = []
temps = [80.3, 270.7, 294.4, 301.2, 315.3]
axins = inset_axes(ax, width=2, height=1.3, loc=2)
axins2 = inset_axes(ax, width=1.3, height=0.9, loc=4)
for i in range(len(semis)):
    x = np.loadtxt(semis[i])
    points = optimizer(semis[i], temps[i], 1, 1, 2, 5)
    d = id_compute(x, points[0][2], points[0][3], x[:,0][points[0][2]]-0.04, x[:,0][points[0][2]]+0.08, 0.15, temps[i])
    ax.semilogy(x[:,0], x[:,1], '.-', label = r'$T= {0} K, n = {1} \pm {2}$'.format(temps[i], round(points[0][0], 3), round(0.15*points[0][0], 3)))
    axins.semilogy(x[:,0], x[:,1], '-')
    axins2.plot(x[:,0], x[:,1], '-', linewidth=1.5)
    axins.semilogy(np.linspace(d[0], d[1], 1000),  np.exp(d[2]*np.linspace(d[0], d[1], 1000)+d[3]), 'k-', linewidth=2.5,
                     )
    i0_pts.append([points[0][0], points[0][1]])
    
axins.tick_params(labelleft=False, labelbottom=False)
axins2.tick_params(labelleft=False, labelbottom=False)
ax.set_xlabel('V / V')
ax.set_ylabel(r' I / A')
ax.grid(True,which="both",ls="-", color='0.75')
ax.legend(loc=8, framealpha=0.8, prop={'size':10})
ax.set_title(r"Ge IV plot various T")   
##### silicon   #######

 
semis = ["si2_805.txt", "si2_2724.txt", "si2_285.txt", "si2_294.txt", 
       "si2_3002.txt", "si2_3097.txt", "si2_3149.txt"]
#plt.figure(100)
fig, ax = plt.subplots()
i0_pts = []
temps = [80.5, 272.4, 285, 294, 300.2, 309.7, 314.9]
axins = inset_axes(ax, width=2, height=1.3, loc=2)
axins2 = inset_axes(ax, width=1.3, height=0.9, loc=4)
for i in range(len(semis)):
    x = np.loadtxt(semis[i])
    points = optimizer(semis[i], temps[i], 1, 1, 2, 5)
    d = id_compute(x, points[0][2], points[0][3], x[:,0][points[0][2]]-0.04, x[:,0][points[0][2]]+0.08, 0.15, temps[i])
    ax.semilogy(x[:,0], x[:,1], '.-', label = r'$T= {0} K, n = {1} \pm {2}$'.format(temps[i], round(points[0][0], 3), round(0.15*points[0][0], 3)))
    axins.semilogy(x[:,0], x[:,1], '-')
    axins2.plot(x[:,0], x[:,1], '-', linewidth=1.5)
    axins.semilogy(np.linspace(d[0], d[1], 1000),  np.exp(d[2]*np.linspace(d[0], d[1], 1000)+d[3]), 'k-', linewidth=1.5,
                     )
    i0_pts.append([points[0][0], points[0][1]])
    
axins.tick_params(labelleft=False, labelbottom=False)
axins2.tick_params(labelleft=False, labelbottom=False)
ax.set_xlabel('V / V')
ax.set_ylabel(r'I / A')
ax.grid(True,which="both",ls="-", color='0.75')
ax.legend(loc=8, framealpha=0.8, prop={'size':10})
ax.set_title(r"Si IV plot various T") 

#################### some code to reverse engineer E_g values for Si, Ge, GaAs

#i0_pts = np.array(i0_pts)
#temps = np.array(temps)
#val = linregress(1/(i0_pts[:,0][1:]*temps[1:]), np.log(i0_pts[:,1][1:]))
#
#plt.figure(100)
#plt.semilogy(1/(i0_pts[:,0][1:]*temps[1:]), i0_pts[:,1][1:], '.')
#plt.semilogy(np.linspace(0.0024, 0.0038, 1000),np.exp(np.linspace(0.0024, 0.0038, 1000)*val[0]+val[1]))

#### gaas ####

semis = ["gaas_805.txt", "gaas_2736.txt", "gaas_2835.txt", "gaas_293.txt",  "gaas_3024.txt", 
      "gaas_3135.txt"]
#plt.figure(100)
fig, ax = plt.subplots()
opts = []
i0_pts = []
temps = [80.5, 273.6, 283.5, 293, 302.4, 313.5]
axins = inset_axes(ax, width=2, height=1.3, loc=2)
axins2 = inset_axes(ax, width=1.3, height=0.9, loc=4)
for i in range(len(semis)):
    x = np.loadtxt(semis[i])
    if i == 0:
        points = optimizer(semis[i], temps[i], 1, 1, 2, 10)
    elif i != 0:
        points = optimizer(semis[i], temps[i], 1, 1, 2, 5)
    minp = 2
    index = 0
    value = []
    for j in range(len(points)):
        if points[j][0] <= minp:
            minp = points[j][0]
            index = j
        if j == len(points)-1:
            value.append(points[index])
    print(value)
    d = id_compute(x, value[0][2], value[0][3], x[:,0][value[0][2]]-0.04, x[:,0][value[0][2]]+0.08, 0.15, temps[i])
    ax.semilogy(x[:,0], x[:,1], '.', label = r'$T= {0} K, n = {1} \pm {2}$'.format(temps[i], round(value[0][0], 3), round(0.15*value[0][0], 3)))
    axins.semilogy(x[:,0], x[:,1], '.')
    axins2.plot(x[:,0], x[:,1], '-', linewidth=1.5)
    axins.semilogy(np.linspace(d[0], d[1], 1000),  np.exp(d[2]*np.linspace(d[0], d[1], 1000)+d[3]), 'k-', linewidth=1.5,
                     )
    i0_pts.append([points[0][0], points[0][1]])
    
axins.tick_params(labelleft=False, labelbottom=False)
axins2.tick_params(labelleft=False, labelbottom=False)
ax.set_xlabel('V / V')
ax.set_ylabel(r'I / A')
ax.grid(True,which="both",ls="-", color='0.75')
ax.legend(loc=8, framealpha=0.8, prop={'size':10})
ax.set_title(r"GaAs IV plot various T")

##### series resistance effect

nits = ["g_803.txt", "si2_805.txt", "gaas_805.txt" ]

fig, ax = plt.subplots()
names = ["Germanium", "Silicon", "Gallium Arsenide"]
for i in range(len(nits)):
    x = np.loadtxt(nits[i])
    if i == 0:
        x = x[0:200]
    points = optimizer(nits[i], 80.5, 1, 1.8, 2.1, 5)
    def id_compute(data, low, high, voltl, volth, error, T, plot=True):
        "computes ideality given data array, slice paramenters and temperature"
        vals = linregress(data[:,0][low:high], np.log(data[:,1][low:high]))
        if plot==True:
#        plt.figure(low)
#        plt.semilogy(data[:,0], data[:,1], '.')
#        plt.semilogy(data[:,0][low:high], data[:,1][low:high], '.')
            plt.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'k-', linewidth=1.5
                     )
#    axins.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'b-', linewidth=2.5,
 #                     )
    #axins.tick_params(labelleft=False, labelbottom=False)
        return [voltl, volth, vals[0], vals[1]]
    plt.semilogy(x[:,0], x[:,1], '.', label= r'{0} at 80 K at n = {1} $\pm$ {2}'.format(names[i], round(points[0][0], 3), round(0.15*points[0][0], 3)))
    d = id_compute(x, points[0][2], points[0][3], x[:,0][points[0][2]]-0.04, x[:,0][points[0][2]]+0.15, 0.15, 80.5)
    def id_compute(data, low, high, voltl, volth, error, T, plot=True):
        "computes ideality given data array, slice paramenters and temperature"
        vals = linregress(data[:,0][low:high], np.log(data[:,1][low:high]))
        if plot==True:
#        plt.figure(low)
#        plt.semilogy(data[:,0], data[:,1], '.')
#        plt.semilogy(data[:,0][low:high], data[:,1][low:high], '.')
            plt.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'r--', linewidth=1
                     )
#    axins.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'b-', linewidth=2.5,
 #                     )
    #axins.tick_params(labelleft=False, labelbottom=False)
        return [voltl, volth, vals[0], vals[1]]    
    
    points = optimizer(nits[i], 80.5, 1, 30, 100, 5)
    d = id_compute(x, points[0][2], points[0][3], x[:,0][points[0][2]]-0.5, x[:,0][points[0][2]]+0.3, 0.15, 80.5)    
    plt.legend(loc=0, framealpha=0.8, prop={'size':10})
    def id_compute(data, low, high, voltl, volth, error, T, plot=True):
        "computes ideality given data array, slice paramenters and temperature"
        vals = linregress(data[:,0][low:high], np.log(data[:,1][low:high]))
        if plot==True:
#        plt.figure(low)
#        plt.semilogy(data[:,0], data[:,1], '.')
#        plt.semilogy(data[:,0][low:high], data[:,1][low:high], '.')
            plt.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'b--', linewidth=1
                     )
#    axins.semilogy(np.linspace(voltl, volth, 1000), np.exp(vals[0]*np.linspace(voltl, volth, 1000)+vals[1]), 'b-', linewidth=2.5,
 #                     )
    #axins.tick_params(labelleft=False, labelbottom=False)
        return [voltl, volth, vals[0], vals[1]]    

    points = optimizer(nits[i], 80.5, 1, 1, 1.6, 3)
    d = id_compute(x, points[0][2], points[0][3], x[:,0][points[0][2]]-0.1, x[:,0][points[0][2]]+0.3, 0.15, 80.5)
plt.xlabel('V / V')
plt.ylabel('I / A')
plt.grid(True,which="both",ls="-", color='0.75')


############# liquid nitrogen temperature IV

nits = ["led2_76.txt", "led4_76.txt", "led6_76.txt" ]
fig, ax = plt.subplots()
for i in range(len(nits)):
    x = np.loadtxt(nits[i])
    ax.semilogy(x[:,0], x[:,1], '.', label = 'Led {0} at T = 76 K'.format(2*(i+1)))
plt.xlabel('V / V')
plt.ylabel(r' I / A')
plt.legend(loc=0, framealpha=0.8, prop={'size':10})


####### reverse bias ###

nits = ["bzx_27_294.txt", "bzx_27_804.txt"]
temps = [294, 80.4, 294, 80.5]
vbs = [0.77, 1.00]
fig, ax = plt.subplots()
axins = inset_axes(ax, width=2, height=1.3, loc=5)
for i in range(len(nits)):
    x = np.loadtxt(nits[i])
    axins.semilogy(x[:,0], x[:,1], '-')
    ax.plot(x[:,0], x[:,1], '.', label = 'BZX {0} at T = {1} K'.format((int(nits[i][4])+0.1*int(nits[i][5])), temps[i]))
    val = linregress(x[:,0][-10:], x[:,1][-10:])
    ax.plot(np.linspace(x[:,0][-40], x[:,0][-1]+0.01, 100), (val[0]*np.linspace(x[:,0][-40], x[:,0][-1]+0.01, 100)+val[1]), 'k-', linewidth=1.5, label = r'$V_b$ = {0} $\pm$ {1} at T = {2} K'.format(vbs[i], round(0.03*vbs[i], 1),temps[i] ))
ax.set_xlabel('V / V')
ax.set_ylabel(r'I / A')
ax.legend(loc=0, framealpha=0.8, prop={'size':10})

##################### reverse bias characteristics

nits = ["bzx_91_294.txt", "bzx_805_91_zoomedin1.txt"]
temps = [294, 80.4, 294, 80.5]
vbs = [7.85, 8.74]
fig, ax = plt.subplots()
axins = inset_axes(ax, width=2, height=1.3, loc=10)
i=0
x = np.loadtxt(nits[i])
axins.semilogy(x[:,0], x[:,1], '-')
ax.plot(x[:,0], x[:,1], '.', label = 'BZX {0} at T = {1} K'.format((int(nits[i][4])+0.1*int(nits[i][5])), temps[i]))
val = linregress(x[:,0][700:720], x[:,1][700:720])
ax.plot(np.linspace(x[:,0][720]-0.07, x[:,0][-1]+0.01, 100), (val[0]*np.linspace(x[:,0][720]-0.07, x[:,0][-1]+0.01, 100)+val[1]), 'k-', linewidth=1.5, label = r'$V_b$ = {0} $\pm$ {1} at T = {2} K'.format(vbs[i], round(0.03*vbs[i], 2),temps[i] ))
i=1
x = np.loadtxt(nits[i])
axins.semilogy(x[:,0], x[:,1], '-')
ax.plot(x[:,0], x[:,1], '.', label = 'BZX {0} at T = {1} K'.format((int(nits[i][4])+0.1*int(nits[i][5])), temps[i]))
val = linregress(x[:,0][-20:], x[:,1][-20:])
ax.plot(np.linspace(x[:,0][380]-0.01, x[:,0][-1]+0.01, 100), (val[0]*np.linspace(x[:,0][380]-0.01, x[:,0][-1]+0.01, 100)+val[1]), 'k-', linewidth=1.5, label = r'$V_b$ = {0} $\pm$ {1} at T = {2} K'.format(vbs[i], round(0.03*vbs[i], 2),temps[i] ))
ax.set_xlabel('V / V')
ax.set_ylabel(r'I / A')
ax.legend(loc=0, framealpha=0.8, prop={'size':10})
