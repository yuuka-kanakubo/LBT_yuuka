import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm
from scipy.special import zeta
import itertools


C_2g = 3
C_2q = 4/3
c_g = 5.6
c_q = 5.8
alphas=0.3
N_C = 3
N_F = 3
C_A = 3
t_ti = 1.0/0.1973 #GEV^-1


def DebyeMass2(T):
	return (N_C + N_F/2)*4*np.pi*alphas*T*T/3

def qhat_g(E, T):
	"""
	Eq.(25) in arXiv: 2306.13742
        """ 
	return (C_2g/np.pi)*42.0*zeta(3)*alphas*alphas*T*T*T*np.log(c_g*E*T/(4.0*DebyeMass2(T)))

def qhat_q(E, T):
	""" 
	Eq.(25) in arXiv: 2306.13742
	""" 
	return (C_2q/np.pi)*42.0*zeta(3)*alphas*alphas*T*T*T*np.log(c_q*E*T/(4.0*DebyeMass2(T)))

def Gamma_inel_q(E, T):
	above_4 = (t_ti * (alphas * C_A * qhat_q(E, T)) / (2.0 * np.sqrt(DebyeMass2(T)))) * \
		(1.0 + (np.sqrt(DebyeMass2(T))/E) * np.log(np.sqrt(DebyeMass2(T))/E))
	below_4 = (t_ti * t_ti * (alphas * C_A * qhat_q(E, T)) / (4.0 * np.pi)) * \
		np.log(E/np.sqrt(DebyeMass2(T)))
	# Use np.where to choose
        ## This is in [GeV]
	return np.where(E > 4, above_4/0.1973, below_4/0.1973)

def Gamma_inel_g(E, T):
	above_4 = (t_ti * (alphas * C_A * qhat_g(E, T)) / (2.0 * np.sqrt(DebyeMass2(T)))) * \
		(1.0 + (np.sqrt(DebyeMass2(T))/E) * np.log(np.sqrt(DebyeMass2(T))/E))
	below_4 = (t_ti * t_ti * (alphas * C_A * qhat_g(E, T)) / (4.0 * np.pi)) * \
		np.log(E/np.sqrt(DebyeMass2(T)))
	# Use np.where to choose
        ## This is in [GeV]
        ##  ==>converting into [fm^-1]
	return np.where(E > 4, above_4/0.1973, below_4/0.1973)


def data_dummy(axis):
	return [0.,0.,0.,0]



def data_Gamma_inel_q(axis):
	data = [[5, 0.618712], [10, 1.10645], [15, 1.41], [20, 1.62896], [25, 
  1.79958], [30, 1.939], [35, 2.05667], [40, 2.15834], [45, 
  2.24776], [50, 2.3275], [55, 2.39942], [60, 2.46489], [65, 
  2.52495], [70, 2.58039], [75, 2.63188], [80, 2.67993], [85, 
  2.72496], [90, 2.76731], [95, 2.80729], [100, 2.84515], [105, 
  2.88108], [110, 2.91529], [115, 2.94791], [120, 2.97909], [125, 
  3.00895], [130, 3.0376], [135, 3.06512], [140, 3.09161], [145, 
  3.11713], [150, 3.14175], [155, 3.16553], [160, 3.18854], [165, 
  3.2108], [170, 3.23238], [175, 3.25331], [180, 3.27363], [185, 
  3.29337], [190, 3.31257], [195, 3.33125], [200, 3.34945], [205, 
  3.36717], [210, 3.38445], [215, 3.40131], [220, 3.41777], [225, 
  3.43385], [230, 3.44956], [235, 3.46492], [240, 3.47996], [245, 
  3.49466], [250, 3.50907], [255, 3.52317], [260, 3.537], [265, 
  3.55055], [270, 3.56384], [275, 3.57688], [280, 3.58967], [285, 
  3.60223], [290, 3.61457], [295, 3.62669], [300, 3.6386], [305, 
  3.6503], [310, 3.66181], [315, 3.67313], [320, 3.68426], [325, 
  3.69522], [330, 3.70601], [335, 3.71663], [340, 3.72708], [345, 
  3.73738], [350, 3.74752], [355, 3.75752], [360, 3.76737], [365, 
  3.77709], [370, 3.78667], [375, 3.79611], [380, 3.80543], [385, 
  3.81462], [390, 3.82369], [395, 3.83263], [400, 3.84147], [405, 
  3.85019], [410, 3.8588], [415, 3.86731], [420, 3.87571], [425, 
  3.884], [430, 3.8922], [435, 3.9003], [440, 3.90831], [445, 
  3.91622], [450, 3.92404], [455, 3.93177], [460, 3.93942], [465, 
  3.94698], [470, 3.95446], [475, 3.96185], [480, 3.96917], [485, 
  3.97641], [490, 3.98357], [495, 3.99066], [500, 3.99767], [505, 
  4.00462], [510, 4.01149], [515, 4.01829], [520, 4.02503], [525, 
  4.0317], [530, 4.03831], [535, 4.04485], [540, 4.05133], [545, 
  4.05775], [550, 4.06411], [555, 4.07041], [560, 4.07665], [565, 
  4.08283], [570, 4.08897], [575, 4.09504], [580, 4.10106], [585, 
  4.10703], [590, 4.11294], [595, 4.11881], [600, 4.12462], [605, 
  4.13039], [610, 4.1361], [615, 4.14177], [620, 4.14739], [625, 
  4.15297], [630, 4.1585], [635, 4.16398], [640, 4.16942], [645, 
  4.17482], [650, 4.18017], [655, 4.18549], [660, 4.19076], [665, 
  4.19599], [670, 4.20118], [675, 4.20633], [680, 4.21145], [685, 
  4.21652], [690, 4.22156], [695, 4.22656], [700, 4.23152], [705, 
  4.23645], [710, 4.24134], [715, 4.2462], [720, 4.25102], [725, 
  4.2558], [730, 4.26055], [735, 4.26527], [740, 4.26996], [745, 
  4.27462], [750, 4.27924], [755, 4.28383], [760, 4.28839], [765, 
  4.29292], [770, 4.29743], [775, 4.3019], [780, 4.30634], [785, 
  4.31075], [790, 4.31514], [795, 4.31949], [800, 4.32382], [805, 
  4.32812], [810, 4.33239], [815, 4.33664], [820, 4.34086], [825, 
  4.34505], [830, 4.34922], [835, 4.35336], [840, 4.35748], [845, 
  4.36156], [850, 4.36564], [855, 4.36968], [860, 4.3737], [865, 
  4.3777], [870, 4.38167], [875, 4.38562], [880, 4.38954], [885, 
  4.39344], [890, 4.39732], [895, 4.40119], [900, 4.40502], [905, 
  4.40883], [910, 4.41263], [915, 4.4164], [920, 4.42015], [925, 
  4.42388], [930, 4.42759], [935, 4.43128], [940, 4.43495], [945, 
  4.4386], [950, 4.44223], [955, 4.44584], [960, 4.44943], [965, 
  4.45301], [970, 4.45656], [975, 4.46009], [980, 4.46361], [985, 
  4.46711], [990, 4.47059], [995, 4.47405], [1000, 4.47749]]
	data = np.array(data)
        ## T = 0.3, quark inelastic Gamma [fm^-1].
        ## Calculated in mathematica.

	if axis == "x":
		return data[:, 0:1]
	elif axis == "y":
		return data[:, 1:2]




def data_Gamma_inel_g(axis):

	data = [[5, 1.57137], [10, 2.64947], [15, 3.30156], [20, 3.77026], [25, 
  4.13589], [30, 4.43536], [35, 4.68873], [40, 4.90815], [45, 
  5.10154], [50, 5.27437], [55, 5.4305], [60, 5.57287], [65, 
  5.70363], [70, 5.82454], [75, 5.93695], [80, 6.04196], [85, 
  6.14047], [90, 6.23322], [95, 6.32086], [100, 6.4039], [105, 
  6.48279], [110, 6.55792], [115, 6.62964], [120, 6.69823], [125, 
  6.76396], [130, 6.82705], [135, 6.88769], [140, 6.94608], [145, 
  7.00236], [150, 7.05669], [155, 7.1092], [160, 7.15999], [165, 
  7.20918], [170, 7.25688], [175, 7.30315], [180, 7.34808], [185, 
  7.39176], [190, 7.43425], [195, 7.47559], [200, 7.51588], [205, 
  7.55513], [210, 7.59342], [215, 7.63079], [220, 7.66728], [225, 
  7.70292], [230, 7.73776], [235, 7.77183], [240, 7.80518], [245, 
  7.83782], [250, 7.86981], [255, 7.90111], [260, 7.93179], [265, 
  7.96188], [270, 7.9914], [275, 8.02036], [280, 8.04879], [285, 
  8.07671], [290, 8.10412], [295, 8.13106], [300, 8.15753], [305, 
  8.18357], [310, 8.20917], [315, 8.23433], [320, 8.2591], [325, 
  8.28347], [330, 8.30747], [335, 8.3311], [340, 8.35438], [345, 
  8.3773], [350, 8.39987], [355, 8.42213], [360, 8.44407], [365, 
  8.4657], [370, 8.48703], [375, 8.50807], [380, 8.52882], [385, 
  8.54929], [390, 8.56949], [395, 8.58944], [400, 8.60913], [405, 
  8.62856], [410, 8.64775], [415, 8.66671], [420, 8.68543], [425, 
  8.70392], [430, 8.7222], [435, 8.74026], [440, 8.75812], [445, 
  8.77576], [450, 8.7932], [455, 8.81044], [460, 8.82749], [465, 
  8.84436], [470, 8.86104], [475, 8.87754], [480, 8.89387], [485, 
  8.91002], [490, 8.92601], [495, 8.94182], [500, 8.95747], [505, 
  8.97297], [510, 8.98832], [515, 9.0035], [520, 9.01854], [525, 
  9.03343], [530, 9.04818], [535, 9.06279], [540, 9.07726], [545, 
  9.09159], [550, 9.10579], [555, 9.11985], [560, 9.13379], [565, 
  9.14761], [570, 9.1613], [575, 9.17487], [580, 9.18831], [585, 
  9.20164], [590, 9.21486], [595, 9.22796], [600, 9.24095], [605, 
  9.25383], [610, 9.2666], [615, 9.27926], [620, 9.29182], [625, 
  9.30428], [630, 9.31664], [635, 9.3289], [640, 9.34106], [645, 
  9.35312], [650, 9.36509], [655, 9.37696], [660, 9.38875], [665, 
  9.40044], [670, 9.41204], [675, 9.42356], [680, 9.43499], [685, 
  9.44633], [690, 9.45759], [695, 9.46877], [700, 9.47987], [705, 
  9.49089], [710, 9.50184], [715, 9.51268], [720, 9.52346], [725, 
  9.53416], [730, 9.54479], [735, 9.55535], [740, 9.56583], [745, 
  9.57625], [750, 9.58659], [755, 9.59685], [760, 9.60705], [765, 
  9.61718], [770, 9.62725], [775, 9.63726], [780, 9.6472], [785, 
  9.65707], [790, 9.66689], [795, 9.67663], [800, 9.68631], [805, 
  9.69594], [810, 9.70551], [815, 9.715], [820, 9.72443], [825, 
  9.73381], [830, 9.74314], [835, 9.75241], [840, 9.76162], [845, 
  9.77076], [850, 9.77988], [855, 9.78893], [860, 9.79792], [865, 
  9.80687], [870, 9.81575], [875, 9.82459], [880, 9.83338], [885, 
  9.84211], [890, 9.8508], [895, 9.85945], [900, 9.86802], [905, 
  9.87656], [910, 9.88505], [915, 9.89349], [920, 9.90189], [925, 
  9.91025], [930, 9.91855], [935, 9.9268], [940, 9.93503], [945, 
  9.9432], [950, 9.95132], [955, 9.95941], [960, 9.96745], [965, 
  9.97546], [970, 9.98341], [975, 9.99132], [980, 9.9992], [985, 
  10.007], [990, 10.0148], [995, 10.0226], [1000, 10.0303]]

	data = np.array(data)
        ## T = 0.3, gluon inelastic Gamma [GeV].
        ## Calculated in mathematica.

	if axis == "x":
		return data[:, 0:1]
	elif axis == "y":
		return data[:, 1:2]




def Gamma_el_g(E, T):
	return 0.0
def Gamma_el_q(E, T):
	return 0.0

# List your data files and their plot settings
data_files = [
{"filename": "DATA_qhat_q.dat",       
	"xlabel": "$E_0$", 
	"ylabel": "$\hat{q} \ \mathrm{[GeV^{3}}]$",
	"title": "qhat of quarks",
	"data":data_dummy, 
	"func": qhat_q,
	"description": "(t-ti)=1 (fm), T=0.3(GeV)"},
{"filename": "DATA_qhat_g.dat",       
	"xlabel": "$E_0$", 
	"ylabel": "$\hat{q} \ \mathrm{[GeV^{3}]}$",      
	"title": "qhat of gluons",                  
	"data":data_dummy, 
	"func": qhat_g,
	"description": "(t-ti)=1 (fm), T=0.3(GeV)"},
{"filename": "DATA_Gamma_inel_q.dat", 
	"xlabel": "$E_0$", 
	"ylabel": "$\Gamma_{inel} \ \mathrm{[fm^{-1}]}$", 
	"title": "Inelastic scat. rate of quarks", 
	"data":data_Gamma_inel_q, 
	"func": Gamma_inel_q,
	"description": "(t-ti)=1 (fm), T=0.3(GeV)"},
{"filename": "DATA_Gamma_inel_g.dat", 
	"xlabel": "$E_0$", 
	"ylabel": "$\Gamma_{inel} \ \mathrm{[fm^{-1}]}$", 
	"title": "Inelastic scat. rate of gluons",  
	"data":data_Gamma_inel_g, 
	"func": Gamma_inel_g,
	"description": "(t-ti)=1 (fm), T=0.3(GeV)"},
{"filename": "DATA_Gamma_el_q.dat", 
	"xlabel": "$E_0$", 
	"ylabel": "$\Gamma_{el} \ \mathrm{[fm^{-1}]}$", 
	"title": "Elastic scat. rate of quarks",       
	"data":data_dummy, 
	"func": Gamma_el_q,
	"description": "T=0.3(GeV)"},
{"filename": "DATA_Gamma_el_g.dat", 
	"xlabel": "$E_0$", 
	"ylabel": "$\Gamma_{el} \ \mathrm{[fm^{-1}]}$", 
	"title": "Elastic scat. rate of gluons",       
	"data":data_dummy, 
	"func": Gamma_el_g,
	"description": "T=0.3(GeV)"},
{"filename": "DATA_channel_10GeVqJET_vsT.dat", 
	"xlabel": "$T (medium)$", 
	"ylabel": "$ \Gamma_{q+a->b+c}/\Gamma^{q}_{el} $", 
	"title": "Probability of each channel ",       
	"data":data_dummy, 
	"func": Gamma_el_g,
	"description": "E=10(GeV)"},
{"filename": "DATA_channel_100GeVqJET_vsT.dat", 
	"xlabel": "$T (medium)$", 
	"ylabel": "$ \Gamma_{q+a->b+c}/\Gamma^{q}_{el} $", 
	"title": "Probability of each channel ",       
	"data":data_dummy, 
	"func": Gamma_el_g,
	"description": "E=100(GeV)"},
{"filename": "DATA_channel_1000GeVqJET_vsT.dat", 
	"xlabel": "$T (medium)$", 
	"ylabel": "$ \Gamma_{q+a->b+c}/\Gamma^{q}_{el} $", 
	"title": "Probability of each channel ",       
	"data":data_dummy, 
	"func": Gamma_el_g,
	"description": "E=1000(GeV)"},
{"filename": "DATA_channel_10GeVgJET_vsT.dat", 
	"xlabel": "$T (medium)$", 
	"ylabel": "$ \Gamma_{g+a->b+c}/\Gamma^{g}_{el} $", 
	"title": "Probability of each channel ",       
	"data":data_dummy, 
	"func": Gamma_el_g,
	"description": "E=10(GeV)"},
{"filename": "DATA_channel_100GeVgJET_vsT.dat", 
	"xlabel": "$T (medium)$", 
	"ylabel": "$ \Gamma_{g+a->b+c}/\Gamma^{g}_{el} $", 
	"title": "Probability of each channel ",       
	"data":data_dummy, 
	"func": Gamma_el_g,
	"description": "E=100(GeV)"},
{"filename": "DATA_channel_1000GeVgJET_vsT.dat", 
	"xlabel": "$T (medium)$", 
	"ylabel": "$ \Gamma_{g+a->b+c}/\Gamma^{g}_{el} $", 
	"title": "Probability of each channel ",       
	"data":data_dummy, 
	"func": Gamma_el_g,
	"description": "E=1000(GeV)"},
]


n_files = len(data_files)

fig, axes = plt.subplots(n_files, 1, figsize=(6, 4*n_files))  # One column, n_files rows
T_fixed = 0.3

if n_files == 1:
    axes = [axes]  # make iterable

for ax, (i, data_info) in zip(axes, enumerate(data_files)):
    data = np.loadtxt(data_info["filename"])

    #Reading DATAfiles.
    if data_info["filename"]=="DATA_channel_10GeVqJET_vsT.dat" \
    or data_info["filename"]=="DATA_channel_100GeVqJET_vsT.dat" \
    or data_info["filename"]=="DATA_channel_1000GeVqJET_vsT.dat" \
    or data_info["filename"]=="DATA_channel_10GeVgJET_vsT.dat" \
    or data_info["filename"]=="DATA_channel_100GeVgJET_vsT.dat" \
    or data_info["filename"]=="DATA_channel_1000GeVgJET_vsT.dat":
         x = data[:,0]

         if "qJET" in data_info["filename"]:
	         CHANNEL=[
	                 "qg2gq",
	                 "qq2qq",
	                 "qqbar2qqbar",
	                 "qqbar2qqbar_flavor_exchange",
	                 "qqbar2qqbar_flavor_kept",
	                 "qqbar2gg"
                 ]
         else:
	         CHANNEL=[
	                 "gg2gg",
	                 "qg2qg",
	                 "gq2gq"
                 ]
         ls = itertools.cycle(('dotted', 'dashed', 'solid', 'dashdot')) 
         mk = itertools.cycle(('s', 'o', 'x', '^', '*')) 
         ax.set_ylim([1e-5, 1])
         for ii in range(len(CHANNEL)):
              y = data[:,ii+1]
              ax.plot(x, y, marker=next(mk), linestyle=next(ls), ms=1.5, color=cm.Accent(ii/len(CHANNEL)), alpha=0.75, label=CHANNEL[ii])
              ax.set_yscale("log")
    else:
         x = data[:,0]
         y = data[:,1]
         ax.plot(x, y, marker='o', linestyle='-', ms=1.5, color="goldenrod", alpha=0.5, label="LBT")

    xpos = np.mean(x)
    ypos = np.mean(y)

    if data_info["filename"]=="DATA_qhat_q.dat" \
	or data_info["filename"]=="DATA_qhat_g.dat":
	    ax.plot(x, [data_info["func"](x_, T_fixed) for x_ in x], linestyle='dotted', lw=1, color="black", alpha=0.75, label="analytical")

    elif data_info["filename"]=="DATA_Gamma_inel_q.dat" \
        or data_info["filename"]=="DATA_Gamma_inel_g.dat" \
        or data_info["filename"]=="DATA_Gamma_el_q.dat" \
        or data_info["filename"]=="DATA_Gamma_el_g.dat":
	    ax.plot(data_info["data"]("x"), data_info["data"]("y"), linestyle='--', lw=1, color="black", alpha=0.75, label="Mathematica")

    #ax.text(ax.get_xlim()[0], ax.get_ylim()[0], data_info["description"], ha='left', va='bottom',
    #		    bbox=dict(facecolor='white', edgecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))
    ax.set_xlabel(data_info["xlabel"])
    ax.set_ylabel(data_info["ylabel"])
    ax.set_title(data_info["title"] + "  " + data_info["description"])
    ax.grid(True)
    ax.legend(fontsize='small', loc='best')

plt.tight_layout()
plt.savefig('TestPlots.pdf')
plt.close()


for info in data_files:
    file_path = info["filename"]

    if os.path.exists(file_path):
	    os.remove(file_path)
	    print(f"File '{file_path}' deleted successfully.")
    else:
	    print(f"ERROR: File '{file_path}' not found.")
