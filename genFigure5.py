
import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import scipy
from scipy.io import loadmat
from scipy.stats import gamma
import random
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plot_mode = 'Var'
plot_mode = 'CV'

# Create a figure
fig = plt.figure(figsize=(7.08, 6))

# Define grid specifications for the layout
gs = gridspec.GridSpec(3, 2, figure=fig)

# Create subplots based on the grid specifications
ax1 = fig.add_subplot(gs[0, :1])
ax2 = fig.add_subplot(gs[0, 1:])
ax3 = fig.add_subplot(gs[1, :1])
ax4 = fig.add_subplot(gs[1, 1:])
ax5 = fig.add_subplot(gs[2, :1])
ax6 = fig.add_subplot(gs[2, 1:])


# The current injection, tested to generate equivalent mean ISIs of 20, 25, 30, 35, 40 ms
BK_Iinj = {}

BK_Iinj[50] = []
BK_Iinj[250] = []


BK_Iinj[50] = [-0.1316, -0.1013, -0.047, 0.06]
BK_Iinj[250] = [-0.0906, -0.0545, 0.0078, 0.136]


# Do a similar Iinj range for one of the plots
BK_Iinj_range = {}

BK_Iinj_range[50] = []
BK_Iinj_range[250] = BK_Iinj_range[50]


BK_Iinj_range[50] = [-0.2, -0.1, 0.0, 0.1, 0.2]
BK_Iinj_range[250] = [-0.2, -0.1, 0.0, 0.1, 0.2]


meanISI = {}
meanISI[50] = []
meanISI[250] = []

meanISI_range = {}
meanISI_range[50] = []
meanISI_range[250] = []

medianISI = {}
medianISI[50] = []
medianISI[250] = []

medianISI_range = {}
medianISI_range[50] = []
medianISI_range[250] = []

stdISI = {}
stdISI[50] = []
stdISI[250] = []

varISI = {}
varISI[50] = []
varISI[250] = []

stdISI_range = {}
stdISI_range[50] = []
stdISI_range[250] = []

varISI_range = {}
varISI_range[50] = []
varISI_range[250] = []


ex_mean_idx = 0
ex_i_idx = 5
for BK_g in BK_Iinj.keys():
    idx=0
    for Iinj in  BK_Iinj[BK_g] + BK_Iinj_range[BK_g]:
        ISI = []
        file=loadmat(f'./NEURON/purkinje_pf_BKbuffer_{BK_g}pS_final_10nm/simdata/fig5_{BK_g}pS_final_10nm_{Iinj}inj/spikes.mat')
        sl=file['tspikelist'][0]
        # Sort into 500 individual traces
        spikes={}
        trial=-1
        splast=100000
        for spike in sl:
            if spike < splast:
                trial +=1
                spikes[trial] = []
            spikes[trial].append(spike)
            splast=spike

        for i in spikes.keys():
            s=np.array(spikes[i])
            ISI.extend(np.diff( s[(s>=100) & (s<600)] ))
        mean_str =  '%s' % float('%.3g' % np.mean(ISI))
        if idx == ex_mean_idx:
            if BK_g==50: ax=ax1
            else: ax=ax2
            keys = list(spikes.keys())
            random.shuffle(keys)
            new_i=0
            for i in keys:
                ax.plot(spikes[i], [new_i]*len(spikes[i]), 'ko', ms=0.5, markerfacecolor='black')
                new_i+=1
                s=np.array(spikes[i])
            ax.set_xlabel('Time (ms)')
            ax.set_ylabel('Trial number')
            ax.set_xlim(100, 600)
            ax.set_title(f'BK{BK_g}, I_inj: {Iinj}nA, Mean ISI: {mean_str}ms', fontweight = 'bold')

        fit_alpha, fit_loc, fit_scale= gamma.fit(ISI, np.mean(ISI)-5) # -5 sometimes needed to capture the correct mean

        fit_beta = 1.0/fit_scale
        # Gamma distriobution mean is a/b or a*scale plus we also have to add the shift 'location'
        gamma_mean = fit_alpha/fit_beta + fit_loc
        # Variance is a*scale^2 or a/b^2
        gamma_var = fit_alpha/(fit_beta*fit_beta)
        # STD is then sqrt of vat
        gamma_std = np.sqrt(gamma_var)
        
        mean_str =  '%s' % float('%.2g' % np.mean(ISI))
        std_str = '%s' % float('%.3g' % gamma_std)
        var_str = '%s' % float('%.3g' % np.var(ISI))
        cv_str = '%s' % float('%.2g' % (np.std(ISI)/np.mean(ISI)))

        if plot_mode == 'Var': lab = f'BK{BK_g} var: {var_str} ms$^2$'
        else: lab = f'BK{BK_g} CV: {cv_str}'
        
        if idx == ex_mean_idx:
            ax=ax3
            ax.hist(ISI, bins=range(0,300,4), label= lab, density=True, alpha=0.7)
            x = np.linspace (0, 200, 200)
            y1 = gamma.pdf(x, fit_alpha, scale=fit_scale, loc=fit_loc)
            ax.set_xlim(5, 75)
            ax.set_xlim(5, 125)
            ax.set_ylim(0, 0.035)
            #ax.set_ylabel('Probability density')
            ax.set_xlabel('ISI (ms)')
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax.text(0.30, 0.90, f'Equal mean = {mean_str}ms', transform=ax.transAxes, fontsize=7, verticalalignment='top', bbox=props)
            ax.legend(loc='center right')
        elif idx == ex_i_idx:
            ax=ax4
            ax.hist(ISI, bins=range(0,300,4), label= lab, density=True, alpha=0.7)
            x = np.linspace (0, 200, 200)
            y1 = gamma.pdf(x, fit_alpha, scale=fit_scale, loc=fit_loc)
            ax.set_xlim(5, 75)
            ax.set_xlim(5, 125)
            ax.set_ylim(0, 0.035)
            #ax.set_ylabel('Probability density')
            ax.set_xlabel('ISI (ms)')
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax.text(0.30, 0.90, f'Equal I_inj = {Iinj}nA', transform=ax.transAxes, fontsize=7, verticalalignment='top', bbox=props)
            ax.legend(loc='center right')

        if Iinj in BK_Iinj[BK_g]:
            meanISI[BK_g].append(float(round(np.mean(ISI),2)))
            medianISI[BK_g].append(float(round(np.median(ISI),2)))
            stdISI[BK_g].append(float(round(np.std(ISI),2)))
            varISI[BK_g].append(float(round(np.std(ISI)*np.std(ISI),2)))
        
        if Iinj in BK_Iinj_range[BK_g]:
            meanISI_range[BK_g].append(float(round(np.mean(ISI),2)))
            medianISI_range[BK_g].append(float(round(np.median(ISI),2)))
            stdISI_range[BK_g].append(float(round(np.std(ISI),2)))
            varISI_range[BK_g].append(float(round(np.std(ISI)*np.std(ISI),2)))

        idx+=1
    
for BK_g in BK_Iinj.keys():
    if plot_mode == 'Var':
        if BK_g == 50: ax5.bar(meanISI[BK_g], varISI[BK_g], label=f'BK{BK_g}', width=-1, align ='edge')
        else: ax5.bar(meanISI[BK_g], varISI[BK_g], label=f'BK{BK_g}', align ='edge')
        ax5.set_ylabel('ISI variance (ms$^2$)')
    else:
        if BK_g == 50: ax5.bar(meanISI[BK_g], np.array(stdISI[BK_g])/np.array(meanISI[BK_g]), label=f'BK{BK_g}', width=-1, align ='edge')
        else: ax5.bar(meanISI[BK_g], np.array(stdISI[BK_g])/np.array(meanISI[BK_g]), label=f'BK{BK_g}', align ='edge')
        ax5.set_ylabel('ISI CV')
    
    ax5.set_xlabel('Mean ISI (ms)')
    
ax5.legend()

shift=-0.002
shift=0.007
for BK_g in BK_Iinj.keys():
    if plot_mode == 'Var':
        if BK_g == 50: ax6.bar(np.array(BK_Iinj_range[BK_g]), varISI_range[BK_g], label=f'BK{BK_g}', width=-0.01, align ='edge')
        else: ax6.bar(np.array(BK_Iinj_range[BK_g])+shift, varISI_range[BK_g], label=f'BK{BK_g}', width=-0.01, align ='edge')
        ax6.set_ylabel('ISI variance (ms$^2$)')
    else:
        if BK_g == 50: ax6.bar(np.array(BK_Iinj_range[BK_g]), np.array(stdISI_range[BK_g])/np.array(meanISI_range[BK_g]), label=f'BK{BK_g}', width=-0.01, align ='edge')
        else: ax6.bar(np.array(BK_Iinj_range[BK_g])+shift, np.array(stdISI_range[BK_g])/np.array(meanISI_range[BK_g]), label=f'BK{BK_g}', width=-0.01, align ='edge')
        ax6.set_ylabel('ISI CV')
    
    ax6.set_xlabel('Iinj (nA)')
    shift+=0.003
    ax6.set_xticks(BK_Iinj_range[BK_g])
ax6.legend()
ax6.set_ylim(0,1.05)


ax1.text(-0.1,1.1,'a', fontweight = 'bold', transform=ax1.transAxes)
ax3.text(-0.1,1.1,'b', fontweight = 'bold', transform=ax3.transAxes)
ax4.text(-0.1,1.1,'c', fontweight = 'bold', transform=ax4.transAxes)
ax5.text(-0.1,1.1,'d', fontweight = 'bold', transform=ax5.transAxes)
ax6.text(-0.1,1.1,'e', fontweight = 'bold', transform=ax6.transAxes)


plt.subplots_adjust(wspace=0.3, hspace=0.5)
plt.savefig('Figures/Figure5.pdf', dpi=500)
plt.close()
