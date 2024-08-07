from tkinter import *
from tkinter import messagebox
from tkinter import Toplevel
import numpy as np
from tkinter import messagebox
# for putting graph in tkinter window
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import interp1d #for interpolation of youngs modulus
import csv

'''-------------------interpolating files for et-val -----------------'''
file1 = "./carbon_steel_low_alloy.csv"
file2 = "./series_3xx_high_alloy.csv"
file3 = "./wrought_70_30_copper_nickel.csv"
file4 = "./nickel_chromium_molybdenum_iron_alloys.csv"
file5 = "./high_strength_bolting.csv"

''''--------------tkinter declarations----------------------'''
root = Tk()
root.title("Remaining Life Calculator")
root.geometry("1400x700")

'''---------------- variables used ------------------------------------'''
sa = StringVar()  #stress amp user input
sa2 = StringVar() #addiional stress amp user input

sigma_ys = StringVar() #sigma(ys) for copper-nickel option
sigma_uts = StringVar() #signma(uts) for carbon option
selectionvar = IntVar()  # radio selection
t = StringVar()  # for time
temp_val = StringVar() #temp val entered by user
temp_val2 = StringVar() # additional temp val entered by user

# for option high strength 2 more user inputs
s = StringVar()
maxm_stress = StringVar()

#var1 n var2 to check if user has checked options (multiple range and weld option respectively)
# for 1st checkbox -> additional temp n sa vals checked
var1 = IntVar()
# for 2nd checkbox -> weld checked
var2 = IntVar()

kf = StringVar() #var for kf var entered by user
clicked = StringVar() #var to check if user chooses options from below
options = ['B', 'C', 'D', 'E', 'F', 'F2', 'G', 'G2', 'W1', 'X', 'S1', 'S2', 'TJ']
clicked.set(options[0]) # intially option set at 'B'

# for final outputs->fatigue results
damage_per_cycle_var = StringVar()
cycles_to_failure_var = StringVar()
remaining_life_var = StringVar()

damage_per_cycle_var2 = StringVar()
cycles_to_failure_var2 = StringVar()
remaining_life_var2 = StringVar()

# for final outputs->BS results
damage_bs = StringVar()
n_bs = StringVar()
life_bs = StringVar()

bs_n_vals = [] #list to store n vals for plotting graph
n_bs_val = 0
damage_bs_val = 0
life_bs_val = 0

''' -------------------------------getting interpolated youngs modulus value --------------------------------------'''

def get_et_val(selection):
    file_data = 0
    if selection == 1:
        file_data = np.genfromtxt(file1, delimiter=',', skip_header=1)
    elif selection == 2:
        file_data = np.genfromtxt(file2, delimiter=',', skip_header=1)
    elif selection == 3:
        file_data = np.genfromtxt(file3, delimiter=',', skip_header=1)
    elif selection == 4:
        file_data = np.genfromtxt(file4, delimiter=',', skip_header=1)
    elif selection == 5:
        file_data = np.genfromtxt(file5, delimiter=',', skip_header=1)

    temperatures = file_data[:, 0]
    young_moduli = file_data[:, 1]

    # Handle NaN values in Young's moduli data
    nan_mask = np.isnan(young_moduli)
    if np.any(nan_mask):
        f = interp1d(temperatures[~nan_mask], young_moduli[~nan_mask], kind='linear', fill_value='extrapolate')
    else:
        f = interp1d(temperatures, young_moduli, kind='linear', fill_value='extrapolate')
    # Interpolate Young's modulus for user input temperature

    try:
        user_temp = float(temp_val.get())
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter valid numeric value for temperature")
        return

    interp_young_modulus = f(user_temp)
    return interp_young_modulus


''' ---------------------------------------------CONDITION A: carbon_steel_low_alloy ----------------------------------------'''

def cond_a(y, sigmays):
    temp = 10 ** y
    sq = y ** 2
    cube = y ** 3
    xval: float = 0

    if sigmays <= 552:
        if temp >= 20:
            xval = -4706.5245 + (1813.6228 * y) + (6785.5644 / y) - (368.12404 * sq) - (5133.7345 / sq) + (
                    30.708204 * cube) + (1596.1916 / cube)
        else:
            num = 38.1309 - (60.1705 * sq) + (25.0352 * sq * sq)
            denom = 1 + (1.80224 * sq) - (4.68904 * (sq ** 2)) + (2.26536 * (sq ** 3))
            xval = num / denom
    elif 793 <= sigmays <= 896:
        if temp >= 43:
            num = 5.37689 - (5.25401 * y) + (1.14427 * sq)
            denon = 1 - (0.960816 * y) + (0.291399 * sq) - (0.0562968 * cube)
            xval = num / denon
        else:
            num = -9.41749 + (14.7982 * y) - (5.94 * sq)
            denot = 1 - (3.46282 * y) + (3.63495 * sq) - (1.21849 * cube)
            xval = num / denot
    n = np.exp(xval * np.log(10))
    return n


def n_a(et_val):
    try:
        sa_val = float(sa.get())
        sigma_val = float(sigma_uts.get())
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter valid numeric value for SA and Sigma(uts)")
        return

    if kf.get():
        kf_val = float(kf.get())
        sa_val = sa_val * kf_val

    ay = np.log10(28.3e3 * (sa_val / et_val))
    an = cond_a(ay, sigma_val)
    return an


def a_plot(et_val):
    sigma_val = float(sigma_uts.get())

    sa_vals = np.linspace(40, 5000, 10000)
    logy = np.log10(28.3e3 * (sa_vals / et_val))
    an_vals = [cond_a(y, sigma_val) for y in logy]

    fig = Figure(figsize=(4.5, 4.5))
    plot = fig.add_subplot(1, 1, 1)
    plot.plot(an_vals, sa_vals)
    plot.set_xscale('log')
    plot.set_yscale('log')
    plot.set_xlim(1e1, 1e11)
    plot.set_ylim(10, 10000)
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    plot.yaxis.set_major_formatter(formatter)
    plot.set_xticks([10 ** i for i in range(1, 12)])
    plot.set_xticklabels([f'1e{i}' for i in range(1, 12)])  # Set x-axis tick labels
    plot.minorticks_on()
    plot.set_xlabel('N', fontsize=6, fontfamily='Arial')
    plot.set_ylabel('Sa', fontsize=6, fontfamily='Arial')
    # plot.legend(prop={'size': 6, 'family': 'Arial'})
    plot.grid(True, 'major', 'both', color='green')
    plot.grid(True, 'minor', 'both', color='grey')
    plot.set_title(
        'Fatigue Curve for Nickel-Chromium-Molybdenum-Iron, Alloys X, G, C-4, and C-276 for temperatures not exceeding 427°C',
        fontsize=5, fontfamily='Arial')
    # Add interactive cursor
    fig.suptitle('S-N Curve', fontsize=7, fontfamily='Times New Roman')
    for item in plot.get_xticklabels() + plot.get_yticklabels():
        item.set_fontsize(5)
        item.set_fontfamily('Arial')

    for item in [plot.xaxis.label, plot.yaxis.label]:
        item.set_fontsize(7)
        item.set_fontfamily('Arial')

    for widget in graph_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()


''' -----------------------------CONDITION B----------------------------------------'''


def cond_b(y):
    temp = 10 ** y
    sq = y ** 2
    logv = np.log(y)
    # xval: float = 0

    if temp >= 14.4:
        num = 17.0181 - 19.8713 * y + 4.21366 * sq
        deno = 1 - 0.1720606 * y - 0.633592 * sq
        if deno == 0:
            return ZeroDivisionError
        xval = num / deno
    else:
        deno = - 0.331096 + 4.3261 * logv / sq
        if deno == 0:
            return ZeroDivisionError
        xval = 1 / deno

    n = np.exp(xval * np.log(10))
    return n


def n_b(et_val):
    try:
        sa_val = float(sa.get())
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter numeric value for SA .")
        return

    if kf.get():
        kf_val = float(kf.get())
        sa_val = sa_val * kf_val

    by = np.log10(28.3e3 * (sa_val / et_val))
    bn = cond_b(by)
    return bn


def b_plot(et_val):
    sa_vals = np.linspace(70, 10000, 100000)
    logy = np.log10(28.3e3 * (sa_vals / et_val))
    bn_vals = [cond_b(y) for y in logy]

    fig = Figure(figsize=(4.5, 4.5))
    plot = fig.add_subplot(1, 1, 1)
    plot.plot(bn_vals, sa_vals)
    plot.set_xscale('log')
    plot.set_yscale('log')
    plot.set_xlim(1e1, 1e11)
    plot.set_ylim(10, 10000)
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    plot.yaxis.set_major_formatter(formatter)
    plot.set_xticks([10 ** i for i in range(1, 12)])
    plot.set_xticklabels([f'1e{i}' for i in range(1, 12)])  # Set x-axis tick labels
    plot.minorticks_on()
    plot.set_xlabel('N', fontsize=6, fontfamily='Arial')
    plot.set_ylabel('Sa', fontsize=6, fontfamily='Arial')
    plot.grid(True, 'major', 'both', color='green')
    plot.grid(True, 'minor', 'both', color='grey')
    plot.set_title(
        'Fatigue Curve for Series 3xx High Alloy Steels, Nickel-Chromium-Iron Alloy, Nickel-Iron-Chromium Alloy,\n '
        'and Nickel-Copper Alloy for Temperatures not Exceeding 800°F', fontsize=5, fontfamily='Arial')
    fig.suptitle('S-N Curve', fontsize=7, fontfamily='Times New Roman')
    for item in plot.get_xticklabels() + plot.get_yticklabels():
        item.set_fontsize(5)
        item.set_fontfamily('Arial')

    for item in [plot.xaxis.label, plot.yaxis.label]:
        item.set_fontsize(7)
        item.set_fontfamily('Arial')

    for widget in graph_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()


''' -----------------------------CONDITION C----------------------------------------'''


def cond_c(y, sigma_val):
    temp = 10 ** y
    sq = y ** 2
    cube = y ** 3
    logv = np.log(y)
    xval = 0

    if sigma_val <= 134:
        xval = -2.632 + (0.1186 / logv) + (15.12 * logv / sq) + (7.087 / sq)
    elif sigma_val == 207:
        if temp >= 24.5:
            xval = 8.580044 - (1.889784 * y) - (8.261383 * logv / y)
        else:
            xval = 5.89029 - (0.2280247 * y) - (6.649501 * logv / y)
    elif sigma_val == 310:
        if temp >= 46:
            xval = -884.9989 + (8936.214 / y) - (36034.88 / sq) + (72508.69 / cube) - 72703.36 / (
                        sq * sq) + 29053.66 / (
                           sq * cube)
        else:
            xval = -17.50197 + 109.168 / y - 236.7921 / sq + 257.9938 / cube - 137.1654 / (sq * sq) + 28.55546 / (
                    sq * cube)

    n = np.exp(xval * np.log(10))
    return n


def n_c(et_val):
    try:
        sa_val = float(sa.get())
        sigma_val = float(sigma_ys.get())
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter valid numeric values for SA and Sigma(ys)")
        return


    if kf.get():
        kf_val = float(kf.get())
        sa_val = sa_val * kf_val

    cy = np.log10(20.0e3 * (sa_val / et_val))
    cn = cond_c(cy, sigma_val)
    return cn


def c_plot(et_val):
    sigma_val = float(sigma_ys.get())

    sa_vals = np.linspace(60, 4000, 100000)
    logy = np.log10(20.0e3 * (sa_vals / et_val))
    cn_vals = [cond_c(y, sigma_val) for y in logy]

    fig = Figure(figsize=(4.5, 4.5))
    plot = fig.add_subplot(1, 1, 1)
    plot.plot(cn_vals, sa_vals)
    plot.set_xscale('log')
    plot.set_yscale('log')
    plot.set_xlim(1e1, 1e6)
    plot.set_ylim(10, 10000)
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    plot.yaxis.set_major_formatter(formatter)
    plot.set_xticks([10 ** i for i in range(1, 7)])
    plot.set_xticklabels([f'1e{i}' for i in range(1, 7)])  # Set x-axis tick labels
    plot.minorticks_on()
    plot.set_xlabel('N', fontsize=6, fontfamily='Arial')
    plot.set_ylabel('Sa', fontsize=6, fontfamily='Arial')
    # plot.legend(prop={'size': 6, 'family': 'Arial'})
    plot.grid(True, 'major', 'both', color='green')
    plot.grid(True, 'minor', 'both', color='grey')
    plot.set_title(
        '– Fatigue Curve for Wrought 70 Copper-Nickel For Temperatures Not Exceeding 371 °C – 134 (18 ) s ys = MPa ksi 134 s ys = MPa',
        fontsize=6, fontfamily='Arial')
    # Add interactive cursor
    fig.suptitle('S-N Curve', fontsize=8, fontfamily='Times New Roman')
    for item in plot.get_xticklabels() + plot.get_yticklabels():
        item.set_fontsize(5)
        item.set_fontfamily('Arial')

    for item in [plot.xaxis.label, plot.yaxis.label]:
        item.set_fontsize(7)
        item.set_fontfamily('Arial')

    for widget in graph_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()


''' -----------------------------CONDITION D----------------------------------------'''


def cond_d(y):
    temp = 10 ** y
    sq = y ** 2

    # xval = 0
    if temp >= 35.9:
        num = -42.08579 + (12.514054 * y)
        deno = 1 - (4.3290016 * y) + (0.60540862 * sq)
        if deno == 0:
            return ZeroDivisionError
        xval = num / deno
    else:
        num = 9.030556 - (8.1906623 * y)
        deno = 1 - (0.36077181 * y) - (0.47064984 * sq)
        if deno == 0:
            return ZeroDivisionError
        xval = num / deno
    n = 10 ** xval
    return n


def n_d(et_val):
    try:
        sa_val = float(sa.get())
    except ValueError:
        messagebox.showerror("Invalid Input","Please enter valid numeric values for SA")
        return

    if kf.get():
        kf_val = float(kf.get())
        sa_val = sa_val * kf_val

    dy = np.log10(28.3e3 * (sa_val / et_val))
    dn = cond_d(dy)
    return dn


def d_plot(et_val):
    sa_vals = np.linspace(100, 5000, 100000)
    logy = np.log10(28.3e3 * (sa_vals / et_val))
    dn_vals = [cond_d(y) for y in logy]

    fig = Figure(figsize=(4.5, 4.5))
    plot = fig.add_subplot(1, 1, 1)
    plot.plot(dn_vals, sa_vals)
    plot.set_xscale('log')
    plot.set_yscale('log')
    plot.set_xlim(1e1, 1e8)
    plot.set_ylim(10, 10000)
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    plot.yaxis.set_major_formatter(formatter)
    plot.set_xticks([10 ** i for i in range(1, 9)])
    plot.set_xticklabels([f'1e{i}' for i in range(1, 9)])  # Set x-axis tick labels
    plot.minorticks_on()
    plot.set_xlabel('N', fontsize=6, fontfamily='Arial')
    plot.set_ylabel('Sa', fontsize=6, fontfamily='Arial')
    plot.grid(True, 'major', 'both', color='green')
    plot.grid(True, 'minor', 'both', color='grey')
    plot.set_title(
        'Fatigue Curve for Nickel-Chromium-Molybdenum-Iron, Alloys X, G, C-4, and C-276 for temperatures not exceeding 427°C',
        fontsize=5, fontfamily='Arial')
    fig.suptitle('S-N Curve', fontsize=7, fontfamily='Times New Roman')
    for item in plot.get_xticklabels() + plot.get_yticklabels():
        item.set_fontsize(5)
        item.set_fontfamily('Arial')

    for item in [plot.xaxis.label, plot.yaxis.label]:
        item.set_fontsize(7)
        item.set_fontfamily('Arial')

    for widget in graph_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()


''' -----------------------------CONDITION E----------------------------------------'''


def cond_e(y, s_val, stress):
    temp = 2.7 * s_val
    sq = y ** 2
    cube = y ** 3

    # xval = 0

    if stress <= temp:
        xval = 3.75565644 - (75.58638 / y) + (403.70774 / sq) - (830.40346 / cube) + (772.53426 / (sq * sq)) - (
                267.75105 / (sq * cube))
    else:
        xval = -9.0006161 + (51.928295 / y) - (86.121576 / sq) + (73.1573 / cube) - (29.945507 / (sq * sq)) + (
                4.7332046 / (sq * cube))
    n = np.exp(xval * np.log(10))
    return n


def n_e(et_val):
    try:
        sa_val = float(sa.get())
        s_val = float(s.get())
        stress_val = float(maxm_stress.get())

    except ValueError:
        messagebox.showerror("Invalid Input","Please enter valid numeric values for SA, s and stress")
        return

    if kf.get():
        kf_val = float(kf.get())
        sa_val = sa_val * kf_val

    ey = np.log10(30e3 * (sa_val / et_val))
    en = cond_e(ey, s_val, stress_val)
    return en


def e_plot(et_val):
    s_val = float(s.get())
    stress_val = float(maxm_stress.get())

    sa_vals = np.linspace(30, 9000, 10000)
    logy = np.log10(30e3 * (sa_vals / et_val))
    en_vals = [cond_e(y, s_val, stress_val) for y in logy]

    fig = Figure(figsize=(4.5, 4.5))
    plot = fig.add_subplot(1, 1, 1)
    plot.plot(en_vals, sa_vals)
    plot.set_xscale('log')
    plot.set_yscale('log')
    plot.set_xlim(1e1, 1e6)
    plot.set_ylim(10, 10000)
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    plot.yaxis.set_major_formatter(formatter)
    plot.set_xticks([10 ** i for i in range(1, 7)])
    plot.set_xticklabels([f'1e{i}' for i in range(1, 7)])  # Set x-axis tick labels
    plot.minorticks_on()
    plot.set_xlabel('N', fontsize=6, fontfamily='Arial')
    plot.set_ylabel('Sa', fontsize=6, fontfamily='Arial')
    # plot.legend(prop={'size': 6, 'family': 'Arial'})
    plot.grid(True, 'major', 'both', color='green')
    plot.grid(True, 'minor', 'both', color='grey')
    plot.set_title('High strength bolting for temperatures not exceeding 371°C (700°F)', fontsize=5, fontfamily='Arial')
    fig.suptitle('S-N Curve', fontsize=7, fontfamily='Times New Roman')
    for item in plot.get_xticklabels() + plot.get_yticklabels():
        item.set_fontsize(5)
        item.set_fontfamily('Arial')

    for item in [plot.xaxis.label, plot.yaxis.label]:
        item.set_fontsize(7)
        item.set_fontfamily('Arial')

    for widget in graph_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()


''' -----------------------------CONDITIONS OVER----------------------------------------'''

'''-------------------------------BS PLOT------------------------------------'''
def bs_graph():
    sr_vals = np.linspace(10, 500, 10000)
    bs_n_vals_dict = {}

    for option in options:
        logc, m, sd = 0, 0, 0
        if option == 'B':
            logc, m, sd = 15.3697, 4, 0.1821
        elif option == 'C':
            logc, m, sd = 14.0344, 3.5, 0.2041
        elif option == 'D':
            logc, m, sd = 12.6008, 3, 0.2095
        elif option == 'E':
            logc, m, sd = 12.5171, 3, 0.2509
        elif option == 'F':
            logc, m, sd = 12.2371, 3, 0.2183
        elif option == 'F2':
            logc, m, sd = 12.0902, 3, 0.2279
        elif option == 'G':
            logc, m, sd = 11.7526, 3, 0.1793
        elif option == 'G2':
            logc, m, sd = 11.5918, 3, 0.1952
        elif option == 'W1':
            logc, m, sd = 11.3979, 3, 0.2140
        elif option == 'X':
            logc, m, sd = 11.9684, 3, 0.2134
        elif option == 'S1':
            logc, m, sd = 16.7710, 5, 0.2350
        elif option == 'S2':
            logc, m, sd = 16.5965, 5, 0.3900
        elif option == 'TJ':
            logc, m, sd = 12.9420, 3, 0.2330

        bs_n_vals_dict[option] = []
        for srt in sr_vals:
            logsr = np.log10(srt)
            logn = logc - (3 * sd) - m * logsr
            n_bs_val_dict = 10 ** logn
            if n_bs_val_dict >= 10 ** 7:
                n_bs_val_dict = 10 ** 20
            bs_n_vals_dict[option].append(n_bs_val_dict)

    fig = Figure(figsize=(9, 9))  # Adjusted figure size for better visualization
    plot = fig.add_subplot(1, 1, 1)

    for key in options:
        plot.plot(bs_n_vals_dict[key], sr_vals, label=key)

    plot.set_yscale('log')
    plot.set_xscale('log')
    plot.set_xlim(1e4, 1e8)
    plot.set_ylim(10, 500)
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    plot.yaxis.set_major_formatter(formatter)
    plot.set_xticks([10 ** i for i in range(4, 9)])
    plot.set_xticklabels([f'1e{i}' for i in range(4, 9)])
    plot.minorticks_on()
    plot.set_xlabel('N', fontsize=6, fontfamily='Arial')
    plot.set_ylabel('Sr', fontsize=6, fontfamily='Arial')
    plot.grid(True, 'major', 'both', color='green')
    plot.grid(True, 'minor', 'both', color='grey')
    plot.set_title('Mean Sr-N Curve', fontsize=8, fontfamily='Arial')
    plot.legend(fontsize='6')

    for item in plot.get_xticklabels() + plot.get_yticklabels():
        item.set_fontsize(5)
        item.set_fontfamily('Arial')

    for widget in graph_frame2.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=graph_frame2)
    canvas.draw()
    canvas.get_tk_widget().pack()


""" ------------------------------ calc for the BS class values----------------------------"""


# calc values corresponding to class n calls graph fnc to plot the graph
def bs_class_calc():
    global n_bs_val
    global damage_bs_val
    global life_bs_val

    try:
        sa_val = float(sa.get())
        time = int(t.get())
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter valid numeric values for SA and time")
        return

    selected_option = clicked.get()

    logc = 0
    m = 0
    sd = 0

    if selected_option == 'B':
        logc = 15.3697
        m = 4
        sd = 0.1821

    elif selected_option == 'C':
        logc = 14.0344
        m = 3.5
        sd = 0.2041

    elif selected_option == 'D':
        logc = 12.6008
        m = 3
        sd = 0.2095

    elif selected_option == 'E':
        logc = 12.5171
        m = 3
        sd = 0.2509

    elif selected_option == 'F':
        logc = 12.2371
        m = 3
        sd = 0.2183

    elif selected_option == 'F2':
        logc = 12.0902
        m = 3
        sd = 0.2279

    elif selected_option == 'G':
        logc = 11.7526
        m = 3
        sd = 0.1793

    elif selected_option == 'G2':
        logc = 11.5918
        m = 3
        sd = 0.1952

    elif selected_option == 'W1':
        logc = 11.3979
        m = 3
        sd = 0.2140

    elif selected_option == 'X':
        logc = 11.9684
        m = 3
        sd = 0.2134

    elif selected_option == 'S1':
        logc = 16.7710
        m = 5
        sd = 0.2350

    elif selected_option == 'S2':
        logc = 16.5965
        m = 5
        sd = 0.3900

    elif selected_option == 'TJ':
        logc = 12.9420
        m = 3
        sd = 0.2330

    sr = sa_val * 2
    logn = logc - (3 * sd) - m * np.log10(sr)
    n_bs_val = 10 ** logn
    damage_bs_val = float((1 / n_bs_val) * 100)
    life_bs_val = 1 / (damage_bs_val * 0.01) * (time / 365)

    damage_bs.set(f'{damage_bs_val:6f}')
    life_bs.set(f'{life_bs_val:6f}')
    n_bs.set(f'{n_bs_val:6f}')

    bs_graph()


"""------ changes to be made when reset button is clicked----------------------------"""
def on_reset():
    sa.set("")
    sa2.set("")
    temp_val2.set("")
    temp_val.set("")
    sigma_ys.set("")
    sigma_uts.set("")
    s.set("")
    t.set("")
    maxm_stress.set("")
    kf.set("")
    clicked.set(options[0])

    for widget in answer_frame.winfo_children():
        widget.destroy()
    for widget in graph_frame.winfo_children():
        widget.destroy()
    for widget in graph_frame2.winfo_children():
        widget.destroy()
    for widget in frame3.winfo_children():
        widget.destroy()

    answer_frame.config(highlightbackground="#f0f0f0",bg="#f0f0f0")
    graph_frame.config(highlightbackground="#f0f0f0",bg="#f0f0f0")
    graph_frame2.config(highlightbackground="#f0f0f0",bg="#f0f0f0")

    for widget in frame_weld.winfo_children():
        if not isinstance(widget, Checkbutton):
            widget.destroy()
    for widget in row2.winfo_children():
        widget.destroy()
    var1.set(None)
    var2.set(None)
    selectionvar.set(2)
    btn.config(text='SUBMIT', command=fatigue_calc)


'''-------------------------------Metal (Radio) selections----------------------------'''
def sigma_fnc(frame):
    if selectionvar.get() == 1:
        # Clear any previous entry fields in the subframe
        for widget in frame.winfo_children():
            widget.destroy()

        Label(frame, text="Enter Sigma(uts) value:", font=('calibre', 12, 'bold')).grid(row=0, column=0)
        Entry(frame, highlightbackground="black",highlightthickness='1',textvariable=sigma_uts, font=('calibre', 12, 'normal')).grid(row=0, column=1)

    elif selectionvar.get() == 3:
        for widget in frame.winfo_children():
            widget.destroy()

        Label(frame, text="Enter Sigma(ys) value:", font=('calibre', 12, 'bold')).grid(row=0, column=0)
        Entry(frame,highlightbackground="black",highlightthickness='1', textvariable=sigma_ys, font=('calibre', 12, 'normal')).grid(row=0, column=1)

    elif selectionvar.get() == 5:
        for widget in frame.winfo_children():
            widget.destroy()

        Label(frame, text="Enter S value:", font=('calibre', 12, 'bold')).grid(row=0, column=0, sticky='e')
        Entry(frame, highlightbackground="black",highlightthickness='1',textvariable=s, font=('calibre', 12, 'normal')).grid(row=0, column=1)

        Label(frame, text="Enter Maximum Nominal Stress value:", font=('calibre', 12, 'bold')).grid(row=3, column=0,
                                                                                                    pady='10')
        Entry(frame,highlightbackground="black",highlightthickness='1', textvariable=maxm_stress, font=('calibre', 12, 'normal')).grid(row=3, column=1, pady='10')

    else:
        for widget in frame.winfo_children():
            widget.destroy()


"""-------------------------displaying answers frame----------------------------------"""
#  fnc to call plot n calc functions above for option chosen by user
def fatigue_calc():
    for widget in answer_frame.winfo_children():
        widget.destroy()

    try:
        time = int(t.get())
    except ValueError:
        messagebox.showerror("Invalid Input","Please enter valid numeric value for time")
        return

    sel_val = int(selectionvar.get())

    output = None

    if sel_val == 1:
        output = n_a(get_et_val(sel_val))
        a_plot(get_et_val(sel_val))

    elif sel_val == 2:
        output = n_b(get_et_val(sel_val))
        b_plot(get_et_val(sel_val))

    elif sel_val == 3:
        output = n_c(get_et_val(sel_val))
        c_plot(get_et_val(sel_val))

    elif sel_val == 4:
        output = n_d(get_et_val(sel_val))
        d_plot(get_et_val(sel_val))

    elif sel_val == 5:
        output = n_e(get_et_val(sel_val))
        e_plot(get_et_val(sel_val))

    btn.config(text="RESUBMIT", command=on_reset)

    damage_per_cycle = float((1 / output) * 100)
    cycles_to_failure = output
    remaining_life = (1 / (damage_per_cycle * 0.01)) * (time / 365)

    damage_per_cycle_var.set(f'{damage_per_cycle:.6f}')
    cycles_to_failure_var.set(f'{cycles_to_failure:.6f}')
    remaining_life_var.set(f'{remaining_life:.6f}')

    if var1.get():
        cycles_to_failure2 = output
        damage_per_cycle2 = float((1 / output) * 100)
        damage_per_cycle_total = damage_per_cycle + damage_per_cycle2
        remaining_life_total = float((1 / (damage_per_cycle_total * 0.01)) * (time / 365))

        damage_per_cycle_var2.set(f'{damage_per_cycle_total:.6f}')
        cycles_to_failure_var2.set(f'{cycles_to_failure2:.6f}')
        remaining_life_var2.set(f'{remaining_life_total:.6f}')

    final_ans()

# setting vals to display in ans frame
def final_ans():
    ans_fatigue = Frame(answer_frame)

    Label(ans_fatigue, text='Fatigue Results: API-579/ASME-VIII-Div2', font=('calibre', 12, 'bold')).pack(side='top',padx=2,pady=6)

    Label(ans_fatigue, text='No. of Cycles to Failure:', font=('calibre', 10, 'bold')).pack(side='top', anchor='w', padx=2, pady=2)
    Label(ans_fatigue, textvariable=cycles_to_failure_var, font=('calibre', 10, 'normal')).pack(side='top',anchor='e',padx=2, pady=2)

    Label(ans_fatigue, text='Damage per Cycle (in %):', font=('calibre', 10, 'bold')).pack(side='top', anchor='w',padx=2, pady=2)
    Label(ans_fatigue, textvariable=damage_per_cycle_var, font=('calibre', 10, 'normal')).pack(side='top',anchor='e', padx=2, pady=2)

    Label(ans_fatigue, text='Remaining Life (in years):', font=('calibre', 10, 'bold')).pack(side='top', anchor='w',padx=2, pady=2)
    Label(ans_fatigue, textvariable=remaining_life_var, font=('calibre', 10, 'normal')).pack(side='top', anchor='e',padx=2, pady=2)


    if var1.get():
        Label(ans_fatigue, text='Total Damage (in %):', font=('calibre', 10, 'bold')).pack(side='top', anchor='w',padx=2, pady=2)
        Label(ans_fatigue, textvariable=damage_per_cycle_var2, font=('calibre', 10, 'normal')).pack(side='top',anchor='e',padx=2, pady=2)
        Label(ans_fatigue, text='Cumulative Remaining Life (in years):', font=('calibre', 10, 'bold')).pack(side='top', anchor='w',padx=2, pady=2)
        Label(ans_fatigue, textvariable=remaining_life_var2, font=('calibre', 10, 'normal')).pack(side='top',anchor='e',padx=2, pady=2)

    ans_fatigue.pack(side='left', pady=5)

    if var2.get():
        bs_class_calc()
        ans_bs = Frame(answer_frame,highlightbackground="black",highlightthickness='1')
        Label(ans_bs, text='BS7608', font=('calibre', 12, 'bold')).pack(side='top',padx=10,pady=10)

        Label(ans_bs, text='No. of Cycles to Failure:', font=('calibre', 10, 'bold')).pack(side='top', anchor='w',padx=2, pady=2)
        Label(ans_bs, textvariable=n_bs, font=('calibre', 10, 'normal')).pack(side='top', anchor='e',padx=2, pady=2)

        Label(ans_bs, text='Damage per Cycle (in %):', font=('calibre', 10, 'bold')).pack(side='top', anchor='w',padx=2, pady=2)
        Label(ans_bs, textvariable=damage_bs, font=('calibre', 10, 'normal')).pack(side='top', anchor='e', padx=2, pady=2)

        Label(ans_bs, text='Remaining Life (in years):', font=('calibre', 10, 'bold')).pack(side='top', anchor='w',padx=2, pady=2)
        Label(ans_bs, textvariable=life_bs, font=('calibre', 10, 'normal')).pack(side='top', anchor='e',padx=2, pady=2)
        ans_bs.pack(side='left', padx=5, pady=10)


"""------------------------------additional-functionalities------------------------"""
def weld_options():
    if var2.get():
        Label(frame_weld, text="Class :", font=('calibre', 10, 'bold')).pack(side='left', padx=5, pady=5
                                                                             , anchor='nw')
        option_menu = OptionMenu(frame_weld, clicked, *options)
        option_menu.configure(font=('arial', 12, 'bold'), bg='#fe9d41', activebackground='#41a2fe',
                              highlightthickness=2.5)
        option_menu['menu'].config(bg='#fe9d41', font=('arial', 12, 'bold'))
        clicked.trace('w', bs_class_calc)
        option_menu.pack(side='left', padx=5, pady=5, anchor='nw')

        Label(frame_weld, text="Kf (according to API)", font=('calibre', 10, 'bold')).pack(
            side='left', padx=5, pady=5, anchor='nw')
        Entry(frame_weld, highlightbackground="black",highlightthickness='1',relief='sunken',textvariable=kf, font=('calibre', 10, 'normal')).pack(side='left', padx=5, pady=5,
                                                                                anchor='nw')

    else:
        for widget in frame_weld.winfo_children():
            if not isinstance(widget, Checkbutton):
                widget.destroy()

def addn_sa_temp():
    if var1.get():

        Label(row2, text="Stress Amplitude II (MPa) :", font=('calibre', 12, 'bold')).pack(side='left', padx=10,
                                                                                           pady=10, anchor='nw')
        Entry(row2, highlightbackground="black",highlightthickness='1', textvariable=sa2, font=('calibre', 12, 'normal'), width=15).pack(side='left', padx=10, pady=10,
                                                                                     anchor='nw')

        Label(row2, text=f"Nodal Temperature II (K):", font=('calibre', 12, 'bold')).pack(
            side='left', padx=10, pady=10, anchor='nw')
        Entry(row2, highlightbackground="black",highlightthickness='1', textvariable=temp_val2, font=('calibre', 12, 'normal'), width=15).pack(side='left', padx=10,
                                                                                           pady=10, anchor='nw')
    else:
        for widget in row2.winfo_children():
            widget.destroy()

'''------------------------frames for displaying ui--------------------------'''

row1 = Frame(root)
Label(row1, text="Stress Amplitude (MPa) :", font=('calibre', 12, 'bold')).pack(side='left', padx=10, pady=10,
                                                                                anchor='nw')
Entry(row1, highlightbackground="#4A4A4A",highlightthickness='1', textvariable=sa, font=('calibre', 12, 'normal'), width=15).pack(side='left', padx=10, pady=10, anchor='nw')

Label(row1, text=f"Nodal Temperature (K):", font=('calibre', 12, 'bold')).pack(side='left', padx=10,
                                                                               pady=10, anchor='nw')
Entry(row1, highlightbackground="#4A4A4A",highlightthickness='1', textvariable=temp_val, font=('calibre', 12, 'normal'), width=15).pack(side='left', padx=10, pady=10,
                                                                                  anchor='nw')

Label(row1, text="Process Time (in days):", font=('calibre', 12, 'bold')).pack(side='left', padx=10, pady=10,
                                                                               anchor='nw')
Entry(row1, highlightbackground="#4A4A4A",highlightthickness='1', textvariable=t, font=('calibre', 12, 'normal'), width=13).pack(side='left', padx=10, pady=10, anchor='nw')
row1.pack(anchor='nw')

Checkbutton(row1, fg='black', bg='#41a2fe', relief='raised', font=('calibre', 12, 'bold'), padx=5, pady=5,
            text='Multiple Range', variable=var1, anchor='e', command=addn_sa_temp).pack(side='right',
                                                                                         padx=80, pady=10,
                                                                                         anchor='ne')

# frame if multiple sa values selected
row2 = Frame(root)
row2.pack(anchor='nw')

''' ------------------------------frame for middle of the body-------------------------------------------------'''
midframe = Frame(root)
# --------------------------------frame for radio buttons-----------------------------------------
frame2 = Frame(midframe)
label2 = Label(frame2, text="Choose one of the options below(Et value):", font=('calibre', 12, 'bold'))
label2.pack(anchor=W)

Lb1 = Radiobutton(frame2, text='Carbon, Low Alloy, Series 4xx', value=1, variable=selectionvar,
                  command=lambda: sigma_fnc(frame3))
Lb2 = Radiobutton(frame2,
                  text='Series 3xx High Alloy Steels, Nickel-Chromium-Iron Alloy, Nickel-Iron-Chromium Alloy, and Nickel-Copper Alloy',
                  value=2, variable=selectionvar, command=lambda: sigma_fnc(frame3))
Lb3 = Radiobutton(frame2, text='Wrought 70-30 Copper-Nickel', value=3, variable=selectionvar,
                  command=lambda: sigma_fnc(frame3))
Lb4 = Radiobutton(frame2, text='Nickel-Chromium-Molybdenum-Iron, Alloys X, G, C-4, and C-276', value=4,
                  variable=selectionvar, command=lambda: sigma_fnc(frame3))
Lb5 = Radiobutton(frame2, text='High strength bolting', value=5, variable=selectionvar,
                  command=lambda: sigma_fnc(frame3))

Lb1.pack(anchor=W)
Lb2.pack(anchor=W)
Lb3.pack(anchor=W)
Lb4.pack(anchor=W)
Lb5.pack(anchor=W)

# -------------------declaration: frame3 for additional sigma values & frame2 for selected radio button -------------------
frame3 = Frame(midframe)
frame3.grid(row=0, column=2, padx=20, pady=20, sticky='n')

frame2.grid(row=0, column=0, padx=20, pady=20, sticky='nw')
midframe.pack(padx=20, anchor='nw')

''''---------------------------- middle frame done: over ---------------------'''

# ------------------------------- frame if weld selected-------------------
frame_weld = Frame(root)
Checkbutton(frame_weld, fg='#000000', bg='#41a2fe', relief='raised', font=('calibre', 12, 'bold'), text='Weld',
            variable=var2, command=weld_options, padx=5, pady=5).pack(side='left', padx=5, pady=5,
                                                                      anchor='nw')
frame_weld.pack(padx=5, pady=5, anchor='nw')

#submit button
btn = Button(root, borderwidth=5, bg='#00d870', activebackground='#fcaa4d', text='SUBMIT', font=('calibre', 12, 'bold'),
             command=fatigue_calc, padx=20, pady=5)
btn.pack(padx=20, pady=10)

''' ---------------frames for displaying answers --------------------'''
#frame for displaying the graph wrt metal chosen by user
graph_frame = Frame(root,highlightbackground="#4A4A4A",highlightthickness='1')
graph_frame.pack(padx=20, pady=20, side='left', anchor='ne')
#frame for displaying the ans after required calculations
answer_frame = Frame(root,highlightbackground="#4A4A4A",highlightthickness='1')
answer_frame.pack(padx=20, pady=20, side='left', anchor='nw')
#frame for displaying the BS std graph
graph_frame2 = Frame(root,highlightbackground="#4A4A4A",highlightthickness='1')
graph_frame2.pack(padx=20, pady=20, side='left', anchor='nw')

root.mainloop()
#end of the root frame
