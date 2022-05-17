#!/usr/bin/env python3
''' computes a grid of LVG models '''

import numpy
import io
import os
import time
import copy
import multiprocessing
from subprocess import Popen, PIPE

###################### BEGIN INPUT ######################

UTILIZED_NUMBER_OF_CORES = 28 # None - all available CPU cores will be used

RESULTED_GRID_FILENAME = "grid_HNCO.txt"
SLOWLY_INCREASE_BEAMING = True
PREVIOUS_BEAMING = None # should be a number or None
UTILIZE_PREVIOUS_SOLUTIONS = True

HDENSITIES = numpy.array([5.6e7])
GAS_TEMPERATURES = numpy.arange(50, 110, 10)
#GAS_TEMPERATURES = [80]
SPECIFIC_COLUMN_DENSITIES = numpy.array([2.5e12])
#SPECIFIC_COLUMN_DENSITIES = numpy.array([3.5e12])
COLLISION_PARTNERS_FRACTIONS = [0.2, 1.0]  # partners for HNCO: He , H2
#### ***************
INITIAL_SOLUTION = 2 # 0 - from file; 1 - optically thin; 2 - LTE
INITIAL_SOLUTION_FILENAME = "Populations/Pn.txt"
FILE_WITH_LAMDA_DATA = "hnco.dat"
MAXIMUM_DpopDt_OR_popDiff = 1.e-10
MAXIMUM_NUMBER_OF_ITERATIONS = 5000000
BEAMING = [1.0]
LINE_WIDTH = 0.1
LIST_OF_TRANSITIONS = [1, 24]
#### Dust parameters, dust emission is computed as:
#### J = DUST_DILLUTION_FACTORS * (1 - exp(DUST_OPTICAL_DEPTHS_AT_FREQS0 * (nu/DUST_FREQS0)^DUST_P)) * planck_function(DUST_TEMPERATURES,nu)
DUST_TEMPERATURES = numpy.arange(50, 310, 10) #[K]
DUST_DILLUTION_FACTORS = [1.0]
DUST_OPTICAL_DEPTHS_AT_FREQS0 = [1.0]
DUST_FREQS0 = [337.08e9] # [Hz]
DUST_P = [1.7]
#### HII region emission parameters, HII region emission is computed as:
#### J = HII_DILLUTION_FACTORS * (1-exp(tauHII)) * planck_function(ELECTRON_TEMPERATURES,nu),
#### tauHII = (TURNOVER_FREQS/nu)^2
ELECTRON_TEMPERATURES = [8000.0] # [K]
HII_DILLUTION_FACTORS = [0.0]
def nut_EM(EMeasure, Telec):
    return 1.e9 * (0.082 * Telec**(-1.35) * EMeasure)**0.476
#TURNOVER_FREQS = [1.e8, 1.e9, 1.e10, 1.e11] # [Hz]
TURNOVER_FREQS = [1.e11] # [Hz]
#
CMB_TEMPERATURE = 2.728 # [K]
# whether there will be the dust within the maser region; 1 - yes, 0 - no
INNER_DUST_INCLUDED = 0
#
DUST_TEMPERATURE_EQUAL_TO_GAS_TEMPERATURE = False

###################### END INPUT ######################

LIST_OF_TRANSITIONS.sort()
LIST_OF_TRANSITIONS = numpy.array(LIST_OF_TRANSITIONS)
NUMBER_OF_SPECTRAL_LINES = len(LIST_OF_TRANSITIONS)
TRANSITIONS_FREQS = []

HDENSITIES = numpy.array(HDENSITIES)
GAS_TEMPERATURES = numpy.array(GAS_TEMPERATURES)
SPECIFIC_COLUMN_DENSITIES = numpy.array(SPECIFIC_COLUMN_DENSITIES)
BEAMING = numpy.array(BEAMING)
DUST_TEMPERATURES = numpy.array(DUST_TEMPERATURES)
DUST_DILLUTION_FACTORS = numpy.array(DUST_DILLUTION_FACTORS)
DUST_OPTICAL_DEPTHS_AT_FREQS0 = numpy.array(DUST_OPTICAL_DEPTHS_AT_FREQS0)
DUST_FREQS0 = numpy.array(DUST_FREQS0)
DUST_P = numpy.array(DUST_P)
ELECTRON_TEMPERATURES = numpy.array(ELECTRON_TEMPERATURES)
HII_DILLUTION_FACTORS = numpy.array(HII_DILLUTION_FACTORS)
TURNOVER_FREQS = numpy.array(TURNOVER_FREQS)

def get_transition_frequencies(filename):
    ''' finds transition frequencies given their ids '''
    fin = open(filename, 'r')
    for i in range(0, 6):
        line = fin.readline()
    number_of_levels = int(line.split()[0])
    fin.readline()
    for i in range(0, number_of_levels):
        fin.readline()
    fin.readline()
    line = fin.readline()
    number_of_rad_transitions = int(line.split()[0])
    fin.readline()
    for i in range(0, number_of_rad_transitions):
        line = fin.readline().split()
        trans_id = int(line[0])
        for j in LIST_OF_TRANSITIONS:
            if j[1] == trans_id:
                TRANSITIONS_FREQS.append(float(line[4]))
    fin.close()
    if len(TRANSITIONS_FREQS) < NUMBER_OF_SPECTRAL_LINES:
        exit("I didn't find all transitions from LIST_OF_TRANSITIONS in the given LAMDA file\n")

class input_parameters:
    ''' contains the set of parameters for a given model '''
    def __init__(self, ident, nH, Tg, N_dV, Td, Wd, tauDust, freqDust, pDust, \
                 Te, Whii, turnFreq, beamH, prev_NdV=None, prev_beam=None):
        ''' constructor '''
        self.ident = "_" + str(ident)
        self.nH = nH
        self.Tg = Tg
        self.N_dV = N_dV
        self.Td = Td
        self.Wd = Wd
        self.tauDust = tauDust
        self.freqDust = freqDust
        self.pDust = pDust
        self.Te = Te
        self.Whii = Whii
        self.turnFreq = turnFreq
        self.beamH = beamH
        self.prevBeam = prev_beam
        self.prevNdV = prev_NdV
        self.init_sol_file = ""
        self.final_sol_file = ""
        if UTILIZE_PREVIOUS_SOLUTIONS:
            if self.__allowed_to_write_pops_file():
                self.final_sol_file = INITIAL_SOLUTION_FILENAME + str(N_dV) + "_" + str(nH) + "_" + str(Tg) + "_" + str(Td) + "_" + str(beamH)
            if prev_NdV is None and prev_beam is None and INITIAL_SOLUTION == 0:
                self.init_sol_file = INITIAL_SOLUTION_FILENAME
            if prev_NdV is None and not prev_beam is None:
                self.init_sol_file = INITIAL_SOLUTION_FILENAME + str(N_dV) + "_" + str(nH) + "_" + str(Tg) + "_" + str(Td) + "_" + str(prev_beam)
            if not prev_NdV is None:
                self.init_sol_file = INITIAL_SOLUTION_FILENAME + str(prev_NdV) + "_" + str(nH) + "_" + str(Tg) + "_" + str(Td) + "_" + str(beamH)
        else:
            if INITIAL_SOLUTION == 0:
                self.init_sol_file = INITIAL_SOLUTION_FILENAME
        if self.init_sol_file != "":
            if not os.path.isfile(self.init_sol_file):
                print("can't find the input file with initial solution from previous calculations: " + self.init_sol_file + " my pars: " + \
                    INITIAL_SOLUTION_FILENAME + str(N_dV) + "_" + str(nH) + "_" + str(Tg) + "_" + str(Td) + "_" + str(beamH)
                )
                if INITIAL_SOLUTION == 0:
                    self.init_sol_file = INITIAL_SOLUTION_FILENAME
                else:
                    self.init_sol_file = ""

    def __allowed_to_write_pops_file(self):
        ''' checks if the given model will produce the file with populations that will be used by other models '''
        if (self.Wd == DUST_DILLUTION_FACTORS[int(DUST_DILLUTION_FACTORS.shape[0]/2)] and \
            self.tauDust == DUST_OPTICAL_DEPTHS_AT_FREQS0[int(DUST_OPTICAL_DEPTHS_AT_FREQS0.shape[0]/2)] and \
            self.freqDust == DUST_FREQS0[int(DUST_FREQS0.shape[0]/2)] and \
            self.pDust == DUST_P[int(DUST_P.shape[0]/2)] and \
            self.Te == ELECTRON_TEMPERATURES[int(ELECTRON_TEMPERATURES.shape[0]/2)] and \
            self.Whii == HII_DILLUTION_FACTORS[int(HII_DILLUTION_FACTORS.shape[0]/2)] and \
            self.turnFreq == TURNOVER_FREQS[int(TURNOVER_FREQS.shape[0]/2)]
        ):
            return True
        else:
            return False


def prepare_input(pars, in_pops_file="", out_pops_file=""):
    ''' prepares input for a give model parameteres set '''
    ortho_para = numpy.min( [3.0, 9.0 * numpy.exp(-170.6 / pars.Tg)] )
    cin = ""
    cin += "# Gas kinetic temperature, K\n"
    cin += str(pars.Tg) + "\n"
    cin += "# Molecular hydrogen density, nH_2, cm^-3\n"
    cin += str(pars.nH) + '\n'
    cin += "# Mean molecular weight per H2 molecule in AMU,         dust-to-gas mass ratio (both are not important if the dust inside the maser region is absent)\n"
    cin += "2.8                                                     0.01\n"
    cin += "# Number of molecular species\n"
    cin += str(1) + '\n'
    cin += "# Specific column density, cm^-3 s\n"
    cin += str(pars.N_dV) + '\n'
    cin += "# Molecular abundance (wrt H2)\n"
    cin += str(1.e-7) + '\n'
    cin += "# Abundances of collision partners with respect to the total number of H2 molecules, the order is the same as in input LAMDA datafile:\n"
    '''for icoll in range(0, len(COLLISION_PARTNERS_FRACTIONS)):
        cin += str(COLLISION_PARTNERS_FRACTIONS[icoll]) + "\n"
    '''
    cin += str(COLLISION_PARTNERS_FRACTIONS[0]) + "\n"
    cin += str(COLLISION_PARTNERS_FRACTIONS[1] / (1.+ortho_para)) + "\n"
    cin += str(COLLISION_PARTNERS_FRACTIONS[1] - COLLISION_PARTNERS_FRACTIONS[1] / (1.+ortho_para)) + "\n"
    cin += "\n"


    cin += "# Initial solution:\n"
    cin += "# 0 - from file; 1 - optically thin; 2 - LTE\n"
    if UTILIZE_PREVIOUS_SOLUTIONS and in_pops_file != "":
        cin += "0" + "\n"
    else:
        cin += str(INITIAL_SOLUTION) + "\n"
    cin += "# Name of file with initial solution (can be empty):\n"
    cin += in_pops_file + "\n"
    cin += "# Name of file with final solution (can be empty):\n"
    cin += out_pops_file + "\n"
    cin += "# Name of LAMDA file with with molecular data:\n"
    cin += FILE_WITH_LAMDA_DATA + "\n"
    cin += "# Stopping criteria:\n"
    cin += "# maximum length of Dn/Dt vector where		| maximum number of\n"
    cin += "# Dn/Dt - time derivative of populations    | iterations\n"
    cin += str(MAXIMUM_DpopDt_OR_popDiff) + "\t" + str(MAXIMUM_NUMBER_OF_ITERATIONS) + "\n"
    cin += "#\n"
    cin += "# beamH, beaming for optical depth parallel to the line of sight\n"
    cin += "# this is eps^-1 = D(ln r) / D(ln V) quantity given in e.g. Sobolev et al. 1997, Cragg et al. 2005\n"
    cin += str(pars.beamH) + "\n"
    cin += "#\n"
    cin += "# line width [km/s] and line profile shape; this is used to account for line overlapping; if it is <=0 overlapping will not be taken into account; r - rectangular line profile, g - Gaussian line profile\n"
    cin += str(LINE_WIDTH) + " r\n"

    cin += "# Take mean intensity from file? 0 - no, 1 - yes\n"
    cin += str(0) + "\n"
    cin += "# name of file with mean intensity (it can be empty):\n"
    cin += "cloudy_spectrum.dat\n"
    cin += "# dust emission: J = Wd * (1 - exp(tau0 * (nu/nu0)^p)) * planck_function(Td,nu) (see Sobolev et al. 1997, see van der Walt 2014 for optically thin case)\n"
    cin += "# where Wd - dillution factor, nu0 - frequency [Hz], nu - radiative transition frequency, p - some number, Td - dust temperature [K]\n"
    cin += "# Wd	tau0	nu0 [Hz]	p	Td [K]		inner_dust\n"
    cin += str(pars.Wd) + " " + str(pars.tauDust) + " " + str(pars.freqDust) + " " + str(pars.pDust) + " " + str(pars.Td) + " " + str(INNER_DUST_INCLUDED) + "\n"
    cin += "# inner dust temperature [K],           inner dust mass absorption coefficient at frequency nu0 [cm^2/g]\n"
    cin += "220                                     3.41\n"
    cin += "# HII region emission: J = WHII * (1-exp(tauHII)) * planck_function(Te,nu), tauHII = (turnFreq/nu)^2 (see van der Walt 2014 and Appendix A from Sobolev et al. 1997)\n"
    cin += "# where WHII - dillution factor, turnFreq - turnover frequency of HII region, Te - electron temperature [K]\n"
    cin += "# WHII	turnFreq [Hz]		Te [K]		HII region is centered on the line-of-sight (1) or not (0)\n"
    cin += str(pars.Whii) + " " + str(nut_EM(pars.turnFreq,pars.Te)) + " " + str(pars.Te) + " 0\n"
    cin += "# Cosmic microwave background temperature in [K] (can be zero) :\n"
    cin += str(CMB_TEMPERATURE) + "\nEOF\n"
    return cin

def compute_model(pars):
    ''' computes a single LVG model with a given parameters set '''
    input_pops_file = pars.init_sol_file
    if SLOWLY_INCREASE_BEAMING and pars.beamH > 1 and pars.prevNdV is None:
        max_beaming = int(pars.beamH)
        first_beaming = 1
        if UTILIZE_PREVIOUS_SOLUTIONS and not pars.prevBeam is None:
            first_beaming = pars.prevBeam
        temp_beams = numpy.arange(first_beaming, max_beaming, 1)
        temp_pars = copy.deepcopy(pars)
        for temp_beam in temp_beams:
            temp_pars.beamH = temp_beam
            cin = prepare_input(temp_pars, input_pops_file, pars.init_sol_file + temp_pars.ident)
            p = Popen(['./LVG_LRT.exe', '1'], stdout=PIPE, stdin=PIPE, stderr=PIPE, universal_newlines=True)
            out, err = p.communicate(input=cin)
            input_pops_file = pars.init_sol_file + temp_pars.ident

    cin = prepare_input(pars, input_pops_file, pars.final_sol_file)
    p = Popen(['./LVG_LRT.exe', '1'], stdout=PIPE, stdin=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = p.communicate(input=cin)

    if SLOWLY_INCREASE_BEAMING and pars.beamH > 1 and pars.prevNdV is None:
        os.remove(input_pops_file)

    in_pars = numpy.array([pars.nH, pars.Tg, pars.N_dV, pars.Td, pars.Wd, pars.tauDust, pars.freqDust, \
                           pars.pDust, pars.Te, pars.Whii, pars.turnFreq, pars.beamH])
    in_pars = "".join(["{} ".format(in_pars[x]) for x in range(len(in_pars))])
    out_errs = ""

    fin = io.StringIO(err)
    while True:
        line = fin.readline()
        if len(line) == 0:
            break
        if line[0] == "#":
            out_errs += " " + line.rstrip("\n")
    fin.close()

    if p.returncode == 0:
        res = ""
        iline = 0
        fin = io.StringIO(out)
        line = fin.readline()
        while iline < NUMBER_OF_SPECTRAL_LINES:
            line = fin.readline()
            line = line.split()
            if len(line) == 0:
                exit("I didn't find all transitions given in LIST_OF_TRANSITIONS")
            if int(line[0]) == LIST_OF_TRANSITIONS[iline][0] and int(line[1]) == LIST_OF_TRANSITIONS[iline][1]:
                res += " " + line[4] + " " + line[5] + " " + line[6]
                iline = iline + 1
        fin.close()
        return in_pars + res + out_errs
    return "# process has failed, return=" + str(p.returncode) + ", pars: " + in_pars + "\t, cerr: " + out_errs

def make_grid(output_filename="grid_results.txt"):
    ''' computes a grid of models '''
    tasks = []
    ids = 0
    for N_dVi in SPECIFIC_COLUMN_DENSITIES:
        for nHi in HDENSITIES:
            for Tgi in GAS_TEMPERATURES:
                for Wdi in DUST_DILLUTION_FACTORS:
                    if (Wdi == 0.0 and INNER_DUST_INCLUDED == 0):
                        dust_temp_array = numpy.array([DUST_TEMPERATURES[0]])
                        tau_dust_array = numpy.array([DUST_OPTICAL_DEPTHS_AT_FREQS0[0]])
                        freq_dust_array = numpy.array([DUST_FREQS0[0]])
                        dust_p_array = numpy.array([DUST_P[0]])
                    else:
                        dust_temp_array = DUST_TEMPERATURES
                        tau_dust_array = DUST_OPTICAL_DEPTHS_AT_FREQS0
                        freq_dust_array = DUST_FREQS0
                        dust_p_array = DUST_P
                    if (Wdi == 0.0):
                        tau_dust_array = numpy.array([DUST_OPTICAL_DEPTHS_AT_FREQS0[0]])
                    if (DUST_TEMPERATURE_EQUAL_TO_GAS_TEMPERATURE):
                        dust_temp_array = numpy.array([Tgi])
                    for Tdi in dust_temp_array:
                        for tauDusti in tau_dust_array:
                            for freqDusti in freq_dust_array:
                                for pDusti in dust_p_array:
                                    for Whiii in HII_DILLUTION_FACTORS:
                                        if (Whiii == 0.0):
                                            te_array = numpy.array([ELECTRON_TEMPERATURES[0]])
                                            turnFreqi_array = numpy.array([TURNOVER_FREQS[0]])
                                        else:
                                            te_array = ELECTRON_TEMPERATURES
                                            turnFreqi_array = TURNOVER_FREQS
                                        for Tei in te_array:
                                            for turnFreqi in turnFreqi_array:
                                                for beamHi in BEAMING:
                                                    tasks.append(input_parameters( ids, \
                                                        nHi, Tgi, N_dVi, Tdi, Wdi, tauDusti, freqDusti, \
                                                        pDusti, Tei, Whiii, turnFreqi, beamHi
                                                    ))
                                                    ids += 1
    print("there will be " + str(ids) + " models")
    pool = multiprocessing.Pool(processes=UTILIZED_NUMBER_OF_CORES)
    results = []
    do_work = pool.map_async(compute_model, tasks, callback=results.append)
    do_work.wait() # Wait on the results

    results = numpy.array(results[0])
    fout = open(output_filename, 'w')
    fout.write("# nH Tg N_dV Td Wd tauDust freqDust pDust Te Whii turnFreq beamH, ")
    for i in TRANSITIONS_FREQS:
        fout.write(str(i) + " ")
    fout.write("\n")
    for i in range(len(results)):
        fout.write(results[i] + "\n")
    fout.close()

def make_grid_saving_populations(output_filename="grid_results.txt"):
    ''' computes a grid of models, uses already computed populations as an initial solution '''
    print("there will be no more than " + str(len(HDENSITIES) * len(GAS_TEMPERATURES) * len(SPECIFIC_COLUMN_DENSITIES) * \
      len(BEAMING) * len(DUST_TEMPERATURES) * len(DUST_DILLUTION_FACTORS) * \
      len(DUST_OPTICAL_DEPTHS_AT_FREQS0) * len(DUST_FREQS0) * len(DUST_P) * \
      len(ELECTRON_TEMPERATURES) * len(HII_DILLUTION_FACTORS) * len(TURNOVER_FREQS) \
      ) + " models")

    fout = open(output_filename, 'w')
    fout.write("# nH Tg N_dV Td Wd tauDust freqDust pDust Te Whii turnFreq beamH, ")
    for i in TRANSITIONS_FREQS:
        fout.write(str(i) + " ")
    fout.write("\n")
    fout.close()
    
    #results = numpy.empty(shape=(0))
    ids = 0
    prev_beam = PREVIOUS_BEAMING
    for beamHi in BEAMING:
        Previous_N_dVi = None
        for N_dVi in SPECIFIC_COLUMN_DENSITIES:
            tasks = []
            temp_results = []
            for nHi in HDENSITIES:
                if N_dVi < 1.5e9 * nHi:
                    for Tgi in GAS_TEMPERATURES:
                        for Wdi in DUST_DILLUTION_FACTORS:
                            if (Wdi == 0.0 and INNER_DUST_INCLUDED == 0):
                                dust_temp_array = numpy.array([DUST_TEMPERATURES[0]])
                                tau_dust_array = numpy.array([DUST_OPTICAL_DEPTHS_AT_FREQS0[0]])
                                freq_dust_array = numpy.array([DUST_FREQS0[0]])
                                dust_p_array = numpy.array([DUST_P[0]])
                            else:
                                dust_temp_array = DUST_TEMPERATURES
                                tau_dust_array = DUST_OPTICAL_DEPTHS_AT_FREQS0
                                freq_dust_array = DUST_FREQS0
                                dust_p_array = DUST_P
                            if (Wdi == 0.0):
                                tau_dust_array = numpy.array([DUST_OPTICAL_DEPTHS_AT_FREQS0[0]])
                            if (DUST_TEMPERATURE_EQUAL_TO_GAS_TEMPERATURE):
                                dust_temp_array = numpy.array([Tgi])
                            for Tdi in dust_temp_array:
                                for tauDusti in tau_dust_array:
                                    for freqDusti in freq_dust_array:
                                        for pDusti in dust_p_array:
                                            for Whiii in HII_DILLUTION_FACTORS:
                                                if (Whiii == 0.0):
                                                    te_array = numpy.array([ELECTRON_TEMPERATURES[0]])
                                                    turnFreqi_array = numpy.array([TURNOVER_FREQS[0]])
                                                else:
                                                    te_array = ELECTRON_TEMPERATURES
                                                    turnFreqi_array = TURNOVER_FREQS
                                                for Tei in te_array:
                                                    for turnFreqi in turnFreqi_array:
                                                        tasks.append(input_parameters( ids, \
                                                            nHi, Tgi, N_dVi, Tdi, Wdi, tauDusti, freqDusti, \
                                                            pDusti, Tei, Whiii, turnFreqi, beamHi, \
                                                            Previous_N_dVi, prev_beam \
                                                        ))
                                                        ids += 1
            if len(tasks) > 0:
                pool = multiprocessing.Pool(processes=UTILIZED_NUMBER_OF_CORES)
                do_work = pool.map_async(compute_model, tasks, callback=temp_results.append)
                do_work.wait() # Wait on the results
                cur_result = numpy.array(temp_results[0])
                #results = numpy.append(results, cur_result)
                fout = open(output_filename, 'a')
                for result_i in cur_result:
                    fout.write(result_i + "\n")
                fout.close()
            Previous_N_dVi = N_dVi
        prev_beam = beamHi

if __name__ == "__main__":
    if not os.path.isfile("LVG_LRT.exe"):
        exit("can't find LVG_LRT.exe program")
    if not os.path.isfile(FILE_WITH_LAMDA_DATA):
        exit("can't find LAMDA input file")
    if INITIAL_SOLUTION == 0 and not UTILIZE_PREVIOUS_SOLUTIONS:
        if not os.path.isfile(INITIAL_SOLUTION_FILENAME):
            exit("can't find the input file with initial solution")
    start_time = time.time()
    get_transition_frequencies(FILE_WITH_LAMDA_DATA)
    if UTILIZE_PREVIOUS_SOLUTIONS:
        make_grid_saving_populations(RESULTED_GRID_FILENAME)
    else:
        make_grid(RESULTED_GRID_FILENAME)
    print("elapsed time = " + str(time.time() - start_time))
