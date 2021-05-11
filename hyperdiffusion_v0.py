import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.integrate as integrate
import math
import os
from timeit import default_timer as timer

start = timer()

rho_pos_array = [10, 15]  # Nr of radial pos. between 0 and 0.8
time_step_array = [0.0001]  # Timestep in [s]
hyper_diff_array = [50]  # Value if the uses hyperdiff

t_end = 0.01  # End simulation time [s} (t_start is 0)
t_save = 0.001
rho_tor_cut_off = 0.8

# 'Choose hyperdiffusion'
# '0 = No hyperdiffusion'
# '1 = Hyperdiffusion added as extra D and V'
# '2 = Hyperdiffusion added as extra D and source'
d_hyperdiff_chooser = 1


# 'Choose source feature'
# '0 = proper integration'
# '1 = Discretization of integral '
choose_source = 1

# 'Choose diffusion'
# '0 = Constant D'
# '1 = D prop (1/n dn/dr)^2'
# '2 = D prop (1/n dn/dr)^2 + 0.5  m^2/s'
# '3 = IFSPPL +0.1'
choose_d = 3

choose_convergence_factor = 0.0001  # Average over all spatail steps, abs value
max_iterations = 10  # Number of max iterations

# 'Define values to create initial profiles (Should fix so values are taken from file)'
b_pos = 0.96  # In rho_tor_norm
b_height = 4  # Value at the top of pedestal [m^-3]
b_sol = 0.1  # Value at the sol [m^-3]
b_width = 0.02  # Widrh of pedestal in rho_tor_norm
b_slope = 0.05  #


# 'Define torus geometry'
r_major = 2.5
r_minor = 1


# 'Create the density profiles:'
def f_ped(r, b_pos_def, b_height_def, b_sol_def, b_width_def, b_slope_def):
    """
    Create a density profile using an mtanh profile
    :param r: position in rho_norm where the density is wanted
    :type r: numpy float scalar or array
    :param b_pos_def: position of density pedestal (-)
    :type b_pos_def: numpy float
    :param b_height_def: height of density pedestal (m^-3)
    :type b_height_def: numpy float
    :param b_sol_def: sol value for density pedestal (m^-3)
    :type b_sol_def: numpy float
    :param b_width_def: width of density pedestal (-)
    :type b_width_def: numpy float
    :param b_slope_def: slope of density pedestal (?)
    :type b_slope_def: numpy float
    :return: value of density at r
    :rtype: numpy float scalar or array (sime dimensions as r)
    """
    return (b_height_def - b_sol_def) / 2 * (mtanh((b_pos_def - r) / (2 * b_width_def), b_slope_def) + 1) + b_sol_def


def mtanh(x, b_slope_def):
    """
    modified tanh function
    :param x: position to evaluate function
    :type x: numpy float scalar or array
    :param b_slope_def: interior slope
    :type b_slope_def: numpy float
    :return: value of the mtanh function evaluated at x
    :rtype: numpy float scalar or array (same dimensions as x)

    See https://pdfs.semanticscholar.org/5dc9/029eb9614a0128ae7c3f16ae6c4e54be4ac5.pdf
    for the mtanh definition
    """
    return ((1 + b_slope_def * x) * np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x))


def diff(dndr_n, h_diff_d):  # 'Returns a D_eff'
    def d_function(x):
        return 50 * x ** 2

    a_ln_critical = 0.2
    d_dummy_org = 0
    d_dummy = 0

    def d_ifs_ppl(a_ln):
        if a_ln < a_ln_critical:
            d_ifs = 0
        else:
            d_ifs = 100 * min(a_ln - a_ln_critical, np.sqrt(abs(a_ln - a_ln_critical)))
        return d_ifs

    if choose_d == 0:
        d_dummy_org = 5
    elif choose_d == 1:
        d_dummy_org = d_function(-dndr_n)
    elif choose_d == 2:
        d_dummy_org = d_function(-dndr_n) + 0.5
    elif choose_d == 3:
        d_dummy_org = 0.1 + d_ifs_ppl(-dndr_n)

    if d_hyperdiff_chooser == 1 or d_hyperdiff_chooser == 2:
        d_dummy = d_dummy_org + h_diff_d  # 'Adds hyperdiff, if chosen'
    elif d_hyperdiff_chooser == 0:
        d_dummy = d_dummy_org

    return d_dummy, d_dummy_org  # 'Sends back original + hyperdiff D and original D'


def pinch(dndr_n, h_diff_v):  # 'Need gradient to calculate V_hyperdiff'
    pinch_dummy_org = -0.5
    pinch_dummy = 0

    if d_hyperdiff_chooser == 1:
        pinch_dummy = pinch_dummy_org + h_diff_v * dndr_n  # 'Add Hyper pinch if chosen'
    elif d_hyperdiff_chooser == 0 or d_hyperdiff_chooser == 2:
        pinch_dummy = pinch_dummy_org

    return pinch_dummy, pinch_dummy_org  # 'Sends back original + hyperdiff V and original V'


def source(position, h):
    integral_value = 0

    def source_function(x):
        return 1000 * 10 ** 18 * x * (x - 1) ** 2 / (10 ** 19)

    def source_function_with_vprime(x):
        return 4 * math.pi ** 2 * r_major * x * source_function(x)

    if choose_source == 0:  # '0 = proper integration'
        if position == 0:
            integral_value = integrate.quad(source_function_with_vprime, 0, h / 2)[0]
        elif position > (rho_tor_cut_off - 3 * h / 2):
            integral_value = integrate.quad(source_function_with_vprime, position - h / 2, position + h)[0]
        else:
            integral_value = integrate.quad(source_function_with_vprime, position - h / 2, position + h / 2)[0]
    elif choose_source == 1:  # '1 = discretization of integral (Like for the time derivaviate)'
        if position == 0:
            integral_value = 2 * math.pi ** 2 * r_major * (h / 2) ** 2 * source_function(position)
        elif position > (rho_tor_cut_off - 3 * h / 2):
            integral_value = 2 * math.pi ** 2 * r_major * ((position + h) ** 2 - (position - h / 2) ** 2) \
                             * source_function(position)

        else:
            integral_value = 2 * math.pi ** 2 * r_major * ((position + h / 2) ** 2 - (position - h / 2) ** 2) \
                             * source_function(position)
    return integral_value


def d_hyper_source(radial_nr, dens_vector, h_diff_h_source, h):
    d_hdiff_dummy_s = 0

    if d_hyperdiff_chooser == 2:  # 'Only d_hyperdiff_chooser == 2 adds a hypersource'
        if radial_nr == 0:
            d_hdiff_dummy_s = dens_vector[radial_nr] * (-h_diff_h_source * vprime(h / 2) / h) \
                              + dens_vector[radial_nr + 1] * (h_diff_h_source * vprime(h / 2) / h)
        else:
            d_hdiff_dummy_s = dens_vector[radial_nr - 1] * (h_diff_h_source * vprime(radial_nr * h - h / 2) / h) + \
                              dens_vector[radial_nr] * (-h_diff_h_source * vprime(radial_nr * h + h / 2)
                                                        / h - h_diff_h_source * vprime(radial_nr * h - h / 2) / h) + \
                              dens_vector[radial_nr + 1] * (h_diff_h_source * vprime(radial_nr * h + h / 2) / h)
    return d_hdiff_dummy_s


def vprime(position):  # calculates vprim. Needs position in r [m]
    return 4 * math.pi ** 2 * r_major * position


def delta_v_func(position, h):  # Calculates delta_V. Needs position in r [m]
    if position == 0:
        delta_v_value = (2 * math.pi ** 2 * r_major) * (h / 2) ** 2  # Volume element from zero to h/2
    elif position > (rho_tor_cut_off - 3 * h / 2):
        delta_v_value = (2 * math.pi ** 2 * r_major) * (
                (position + h) ** 2 - (position - h / 2) ** 2)  # Volume element from pos-h/2 to pos+h
    else:
        delta_v_value = (2 * math.pi ** 2 * r_major) * (
                (position + h / 2) ** 2 - (position - h / 2) ** 2)  # Volume element from pos-h/2 to pos+h/2
    return delta_v_value


def norm_dens_der(array_of_density, h):
    array_of_derivatives = np.zeros(len(array_of_density) - 1)
    x_norm = 0
    while x_norm < len(array_of_derivatives):
        array_of_derivatives[x_norm] = ((array_of_density[x_norm + 1] - array_of_density[x_norm]) / h) / (
                (array_of_density[x_norm + 1] + array_of_density[x_norm]) / 2)
        x_norm += 1
    return array_of_derivatives


def check_convergence(new_solution, old_solution):  # Check abs error
    total_error = 0
    a = 0
    while a < (len(new_solution) - 1):  # The last point in the solution are fixed, dont take it into account
        total_error = total_error + np.abs(new_solution[a] - old_solution[a])
        a += 1
    norm_error = total_error / len(new_solution)
    if norm_error > choose_convergence_factor:
        conv = False
    else:
        conv = True
    return conv


def plotter(x_grid_plotter, data_plotter, nr_splits_plotter, name, t_step, spatial_nr, h_diff_plot):
    counter_for_plot = 1
    plt.figure()
    plt.plot(x_grid_plotter, data_plotter[0, :], label='t=0, tstart')
    while counter_for_plot < nr_splits_plotter - 1:
        plt.plot(x_grid_plotter, data_plotter[round(counter_for_plot * t_end / ((nr_splits_plotter - 1) * t_save)), :],
                 label='t=' + str(counter_for_plot * t_end / (nr_splits_plotter - 1)))
        counter_for_plot = counter_for_plot + 1
    plt.plot(x_grid_plotter, data_plotter[round(t_end / t_save), :], label='t=' + str(t_end) + ', tend')
    plt.title(name + '. Spatial nr: ' + str(spatial_nr) + ", t_step: " + str(t_step) + ", h_diff: " + str(h_diff_plot))
    plt.legend()
    plt.savefig(name + ".png")
    plt.show()
    return


def implicit_d_and_v_solver(ne_matrix_solver, x_grid, h, convergnce_nr, d_matrix, v_matrix, t_step, h_diff):
    # 'Create matrix to store for use in solvers'
    dens_matrix_lin = np.zeros((len(x_grid), len(x_grid)))
    value_array_lin = np.zeros(len(x_grid))
    t = 1  # 'Start to calculate the first timestep, i.e. t=1*timestep'
    old_profile = ne_matrix_solver[t - 1, 0:len(x_grid)]
    old_time_profile = ne_matrix_solver[t - 1, 0:len(x_grid)]
    convergence_counter = 0
    while t < (t_end / t_step + 1):
        norm_derivates_one_timestep = norm_dens_der(old_profile, h)
        # 'Add the zeroth entry'
        dndr_n_plus = norm_derivates_one_timestep[0]
        v_prime_plus = vprime(h / 2)
        delta_v = delta_v_func(0, h)
        dens_matrix_lin[0, 0] = -delta_v / t_step - diff(dndr_n_plus, h_diff)[0] * v_prime_plus / h - \
                                pinch(dndr_n_plus, h_diff)[0] * v_prime_plus / 2
        dens_matrix_lin[0, 1] = diff(dndr_n_plus, h_diff)[0] * v_prime_plus / h - \
                                pinch(dndr_n_plus, h_diff)[0] * v_prime_plus / 2
        value_array_lin[0] = -delta_v * old_time_profile[0] / t_step - source(0, h) + \
                             d_hyper_source(0, old_time_profile[:], h_diff, h)
        x = 1  # start from 2nd gridpoint
        while x < len(x_grid) - 1:
            delta_v = delta_v_func(x * h, h)
            dndr_n_minus = norm_derivates_one_timestep[x - 1]
            dndr_n_plus = norm_derivates_one_timestep[x]
            v_prime_minus = vprime(x * h - h / 2)
            v_prime_plus = vprime(x * h + h / 2)

            dens_matrix_lin[x, x - 1] = diff(dndr_n_minus, h_diff)[0] * v_prime_minus / h + \
                                        pinch(dndr_n_minus, h_diff)[0] * v_prime_minus / 2
            dens_matrix_lin[x, x] = -delta_v / t_step - diff(dndr_n_plus, h_diff)[0] * v_prime_plus / h - \
                                    pinch(dndr_n_plus, h_diff)[0] * v_prime_plus / 2 \
                                    - diff(dndr_n_minus, h_diff)[0] * v_prime_minus / h + pinch(dndr_n_minus, h_diff)[
                                        0] * v_prime_minus / 2
            dens_matrix_lin[x, x + 1] = diff(dndr_n_plus, h_diff)[0] * v_prime_plus / h - \
                                        pinch(dndr_n_plus, h_diff)[0] * v_prime_plus / 2
            value_array_lin[x] = -delta_v * old_time_profile[x] / t_step - source(x * h, h) + \
                                 d_hyper_source(x, old_time_profile[:], h_diff, h)
            x += 1

        # 'Add the last entry seperately'
        dens_matrix_lin[x, x] = 1  # 'Fixed value boundary'
        value_array_lin[x] = old_time_profile[x]
        # 'Solve it'
        solution = np.linalg.solve(dens_matrix_lin, value_array_lin)

        # 'Check convergence'
        if convergence_counter == 0:
            convergence = False
        else:
            convergence = check_convergence(solution, old_profile)

        if (convergence is True) or (convergence_counter > max_iterations):
            old_time_profile = solution
            convergence_counter = 0
            if (t % (t_save / t_step)) == 0:
                ne_matrix_solver[round(t * t_step / t_save), 0:len(x_grid)] = solution
                convergnce_nr[round(t * t_step / t_save)] = convergence_counter
                x_d = 0
                while x_d < len(x_grid) - 1:
                    d_matrix[round(t * t_step / t_save), x_d, 0] = diff(norm_derivates_one_timestep[x_d], h_diff)[
                        0]  # 'Sparar D_orginal+D_hyperdiff'
                    v_matrix[round(t * t_step / t_save), x_d, 0] = pinch(norm_derivates_one_timestep[x_d], h_diff)[
                        0]  # 'Sparar V_orginal+V_hyperdiff'
                    d_matrix[round(t * t_step / t_save), x_d, 1] = diff(norm_derivates_one_timestep[x_d], h_diff)[
                        1]  # 'Sparar D_orginal'
                    v_matrix[round(t * t_step / t_save), x_d, 1] = pinch(norm_derivates_one_timestep[x_d], h_diff)[
                        1]  # 'Sparar V_orginal'
                    x_d += 1
            t += 1
        else:
            old_profile = solution
            convergence_counter += 1
    return ne_matrix_solver, convergnce_nr, d_matrix, v_matrix


def initiliazation(spatial_nr, t_save_init):
    x_grid = np.linspace(0, rho_tor_cut_off, spatial_nr)
    h = rho_tor_cut_off / (spatial_nr - 1)
    ne_matrix = np.zeros((round(t_end / t_save_init + 1), len(x_grid) + len(x_grid_outer)))
    # Fills the outer part
    ne_matrix[0, 0:len(x_grid)] = n_initialize_values_spline(x_grid)
    ne_matrix[0, (len(x_grid)):(len(x_grid) + len(x_grid_outer))] = n_initialize_values_spline(x_grid_outer)
    t = 1
    while t < (t_end / t_save_init + 1):
        x = 0
        while x < (len(x_grid_outer)):
            ne_matrix[t, len(x_grid) + x] = ne_matrix[0, len(x_grid) + x]
            x += 1
        t += 1

    convergnce_nr = np.zeros(round(t_end / t_save_init + 1))
    d_matrix = np.zeros((round(t_end / t_save_init + 1), len(x_grid) + len(x_grid_outer), 2))
    v_matrix = np.zeros((round(t_end / t_save_init + 1), len(x_grid) + len(x_grid_outer), 2))
    return x_grid, h, ne_matrix, convergnce_nr, d_matrix, v_matrix


# 'Initializerar datan'
solver_name = 'Implicit D+V seperation'
nr_grid_outer = 21
h_outer = (1 - rho_tor_cut_off) / (nr_grid_outer - 1)  # 'Radiell distance between gridpoints outer'
x_grid_outer = np.linspace(rho_tor_cut_off + h_outer, 1, nr_grid_outer - 1)

# 'Dummy radial grid, kallar initial profile'
x_init = np.linspace(0, 1, 1001)  # 'Determines the initalization grid_size'
value_init = f_ped(x_init, b_pos, b_height, b_sol, b_width, b_slope)
n_initialize_values_spline = interp.InterpolatedUnivariateSpline(x_init, value_init)

# 'Creates all dictionaries'
x_grid_dict = {}
h_dict = {}
ne_matrix_dict = {}
convergnce_nr_dict = {}
d_matrix_dict = {}
v_matrix_dict = {}
x_grid_total = {}

# 'Creates grids'
for i in range(len(rho_pos_array)):  # loop over rho_pos
    for j in range(len(time_step_array)):  # loop over tme_steps
        for k in range(len(hyper_diff_array)):  # loop over hyper_diffs
            x_grid_dict[i, j, k], h_dict[i, j, k], ne_matrix_dict[i, j, k], convergnce_nr_dict[i, j, k], d_matrix_dict[
                i, j, k], v_matrix_dict[i, j, k] \
                = initiliazation(rho_pos_array[i], t_save)
            x_grid_total[i, j, k] = np.append(x_grid_dict[i, j, k], x_grid_outer)

# 'Call the solvers here'
for i in range(len(rho_pos_array)):  # loop over rho_pos
    for j in range(len(time_step_array)):  # loop over tme_steps
        for k in range(len(hyper_diff_array)):  # loop over hyper_diffs
            Data = implicit_d_and_v_solver(ne_matrix_dict[i, j, k], x_grid_dict[i, j, k], h_dict[i, j, k],
                                           convergnce_nr_dict[i, j, k], d_matrix_dict[i, j, k],
                                           v_matrix_dict[i, j, k], time_step_array[j],
                                           hyper_diff_array[k])  # 'Returns nematrix at 0 and solvername at 1'
            ne_matrix_dict[i, j, k] = Data[0]
            convergnce_nr_dict[i, j, k] = Data[1]
            d_matrix_dict[i, j, k] = Data[2]
            print("Simulation for rho_nr: " + str(rho_pos_array[i]) + ", t_step: " + str(
                time_step_array[j]) + ", hyper_diff: " + str(hyper_diff_array[k]) + ", is done")

# 'Plotting'
for i in range(len(rho_pos_array)):  # loop over rho_pos
    for j in range(len(time_step_array)):  # loop over tme_steps
        for k in range(len(hyper_diff_array)):  # loop over hyper_diffs
            plotter(x_grid_total[i, j, k], ne_matrix_dict[i, j, k], 6, 'Density', time_step_array[j], rho_pos_array[i],
                    hyper_diff_array[k])

plt.show()

# 'Save the data'
print("Saving data")
for i in range(len(rho_pos_array)):  # loop over rho_pos
    for j in range(len(time_step_array)):  # loop over tme_steps
        for k in range(len(hyper_diff_array)):  # loop over hyper_diffs
            if os.path.exists('data_rho' + str(rho_pos_array[i]) + '_tstep' + str(time_step_array[j]) + '_hdiff' + str(
                    hyper_diff_array[k]) + '_hdiffmetod' + str(d_hyperdiff_chooser) + '_Dmethod' + str(
                    choose_d) + '_tend' + str(t_end) + '.txt'):
                print('rho: ' + str(rho_pos_array[i]) + ', tstep: ' + str(time_step_array[j]) + ', hdiff: ' + str(
                    hyper_diff_array[k]) + ', hdiffmetod: ' + str(d_hyperdiff_chooser) + ', Dmethod: ' + str(
                    choose_d) + ', tend: ' + str(t_end))
                print("file exist, removing it")
                os.remove('data_rho' + str(rho_pos_array[i]) + '_tstep' + str(time_step_array[j]) + '_hdiff' + str(
                    hyper_diff_array[k]) + '_hdiffmetod' + str(d_hyperdiff_chooser) + '_Dmethod' + str(
                    choose_d) + '_tend' + str(t_end) + '.txt')
            else:
                print('rho: ' + str(rho_pos_array[i]) + ', tstep: ' + str(time_step_array[j]) + ', hdiff: ' + str(
                    hyper_diff_array[k]) + ', hdiffmetod: ' + str(d_hyperdiff_chooser) + ', Dmethod: ' + str(
                    choose_d) + ', tend: ' + str(t_end))
                print("file doesn't exist")
            file = open('data_rho' + str(rho_pos_array[i]) + '_tstep' + str(time_step_array[j]) + '_hdiff' + str(
                hyper_diff_array[k]) + '_hdiffmetod' + str(d_hyperdiff_chooser) + '_Dmethod' + str(
                choose_d) + '_tend' + str(t_end) + '.txt', 'a')
            counter = 0
            while counter < len(x_grid_total[i, j, k]):
                file.write(str('%2.9f' % x_grid_total[i, j, k][counter]) + "\t" + str(
                    '%2.9f' % ne_matrix_dict[i, j, k][-1, counter]) + "\n")
                counter = counter + 1

end = timer()
print("Time elapsed v0: " + str(end - start))
