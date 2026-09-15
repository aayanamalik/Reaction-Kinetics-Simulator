import matplotlib.pyplot as plt
import numpy as np
from time import time
from scipy.integrate import solve_ivp

# Constants
R = 8.31446261815324 # J K-1 mol-1

def get_rate_constant(A_factor, Ea, T):
    # Returns the rate constant using the Arrhenius equation.
    return A_factor * np.exp(-Ea/(R*T))
    
def get_net_rate(concs, k_forward, k_reverse, reactant_orders, product_orders):
    """
    Calculates the net rate of a reaction

    Args:
        concs = an array containing the concentrations of each substance
        k_forward, k_reverse = k-values for forward and backward specifically
        reactant_orders = power of reactants e.g. [A]^a, [B]^b, ... to be used in rate equation

    Returns:
        The net rate
    """
    
    forward_rate = k_forward
    reverse_rate = k_reverse

    """
    Rate equation calculation
    E.g. for Haber Process: forward rate = [N2]^1 * [H2]^3
                            reverse rate = [NH3]^2
                                net rate = forward rate - reverse rate
    """
    
    for i in range(len(concs)):
        if reactant_orders[i] != 0:
            forward_rate *= concs[i]**reactant_orders[i]
        if product_orders[i] != 0:
            reverse_rate *= concs[i]**product_orders[i]
    
    return forward_rate - reverse_rate

def get_temperature_change(
    net_rate, volume, delta_H, UA, current_temp, surroundings_temp, density, Cp
):
    """
    Calculates the temperature change per second (dT/dt) of a reaction.
    
    Args:
        net_rate: rate of reaction (mol dm^-3 s^-1)
        volume: reaction volume (dm^3)
        delta_H: enthalpy change (J mol^-1)
        UA: overall heat transfer coefficient * area (W K^-1)
        current_temp: current reactor temperature (K)
        surroundings_temp: ambient/cooling temperature (K)
        density: mixture density (kg m^-3)
        Cp: specific heat capacity (J kg^-1 K^-1)
    """
    # Heat rates in Watts (J/s)
    # (mol dm^-3 s^-1) * dm^3 * (J mol^-1) = J/s
    Q_generated = net_rate * volume * -delta_H
    
    # (W K^-1) * K = W (J/s)
    Q_removed = UA * (current_temp - surroundings_temp)
    
    # Mass calculation requires volume in m^3 (1 dm^3 = 0.001 m^3)
    volume_m3 = volume / 1000.0
    mass = volume_m3 * density  # kg
    thermal_inertia = mass * Cp  # J/K
    
    delta_T = (Q_generated - Q_removed) / thermal_inertia  # K/s
    return delta_T

def simple_reaction():
    # Constants
    A0 = 10 # 100 moles of initial reactant
    k = 0.01 # Rate constant is 0.1s^-1
    S = 500 # Seconds over which reaction takes place
    
    t = np.arange(S)
    
    reactant_over_time = A0*np.exp(-k*t)
    product_over_time = A0*(1 - np.exp(-k*t))
    
    plt.plot(reactant_over_time, c = 'r', label = 'Reactant')
    plt.plot(product_over_time, c = 'b', label = 'Product')
    plt.legend()
    
    plt.ylabel('Concentration (moldm-3)')
    plt.xlabel('Time (s)')
    
    plt.show()

def k_comparison():
    A0 = 100 # 100 moles of initial reactant
    k_values = np.array([0.001,0.01,0.1,1]) # Rate constant
    S = 1000 # Seconds over which reaction takes place
    
    t = np.arange(S) # An array: [0,1,2,...,S-2,S-1]
    
    reactant_over_time = np.zeros((4,S)) # Array with four rows for each k value
    for i in range(k_values.size): # Cycles i for each k value
        reactant_over_time[i] = A0*np.exp(-(k_values[i])*t) # Variable i selects each row and calculates designated k value
        plt.plot(reactant_over_time[i], label = f'k = {k_values[i]}') # Graph each row, label: k = e.g. 0.1
    
    
    plt.legend()
    plt.ylabel('Concentration of reactant (moldm-3)')
    plt.xlabel('Time (s)')
    
    plt.show()

def temp_comparison():
    # Define Values
    A0 = 10 # 10 moles of initial reactant
    S = 5000 # Seconds over which reaction takes place
    Ea = 50000 # Activation energy in joules
    A_factor = 1e6 # Pre-exponential factor

    # Variables for while loop
    ask_for_temp = True
    count = 1
    temp_values = []

    # Loop to ask for temperature values
    while ask_for_temp:
        temp = input(f'Enter value {count}: ')
        if temp == 'd':
            temp_values = np.array([273.15,283.15,293.15,298.15,303.15])
            ask_for_temp = False
        elif temp == 'f':
            ask_for_temp = False
        else:
            temp_values.append(float(temp))
        
        count+=1

    temp_values = np.array(temp_values)
    
    t = np.arange(S) # 1-dimensional array with values from 0 to S-1
    
    k_values = A_factor*np.exp(-Ea/(R*temp_values)) # Find k values from temperature values using Arrhenius

    reactant_over_time = np.zeros((temp_values.size,S)) # Empty array to take in values of reactant amount for each temperature value
    for i in range(temp_values.size):
        reactant_over_time[i] = A0*np.exp(-(k_values[i])*t) # Iterate through each temperature (row) in reactants_over_time
        plt.plot(reactant_over_time[i], label = str(f'{temp_values[i]}K')) # Plot values

    plt.legend()
    plt.ylabel('Concentration of reactant (moldm-3)')
    plt.xlabel('Time (s)')
    
    plt.show()

    print('Temperature vs Rate Constant:')
    for temp, k in zip(temp_values, k_values):
        print(f'Temperature: {temp:.0f}K | K: {k:.6f}')
    
    print(f'\nConcentration after 1200 seconds (20 minutes) from {A0} moldm-3:')
    for i in range(temp_values.size):
        print(f'Temperature: {temp_values[i]}K | Concentration: {reactant_over_time[i][1200]:.2f} moldm-3')

def consecutive_reaction():
    # This models the reaction A->B->C
    
    # Define constants
    A0 = 10 # mol
    k1 = 0.05 # s^-1
    k2 = 0.01 # s^-1
    S = 300 # Seconds over which reaction takes place

    names_of_reactants = ['A', 'B', 'C']
    t = np.arange(S) # An array: [0,1,2,...,S-2,S-1] to calculate concentration at each second

    reactants_over_time = np.zeros((3,S))
    reactants_over_time[0] = A0*np.exp(-k1*t)                                     # A
    reactants_over_time[1] = k1*A0*(np.exp(-k1*t)-np.exp(-k2*t))/(k2-k1)          # B
    reactants_over_time[2] = A0 - reactants_over_time[0] - reactants_over_time[1] # C

    for i in range(3):
        plt.plot(reactants_over_time[i], label = f'Concentration of {names_of_reactants[i]}')
        
    plt.legend()
    plt.ylabel('Concentration (moldm-3)')
    plt.xlabel('Time (s)')
    plt.title('A modelling of A->B->C')
    
    plt.show()
    
    print('Final Concentrations:')
    for i in range(3):
        print(f'{names_of_reactants[i]}: {reactants_over_time[i,-1]:.2f} moldm-3')
    
    print(f'\nMax concentration of B: {reactants_over_time[1].max():.2f} moldm-3')

def equilibria():
    # This models the reaction A<->B analytically

    # Define constants
    A0 = 10 # mol
    B0 = 0 # mol
    k_forward = 0.05 # s^-1
    k_reverse = 0.01 # s^-1
    S = 240 # Seconds over which reaction takes place
    
    
    t = np.arange(S) # An array: [0,1,2,...,S-2,S-1] to calculate concentration at each second

    reactant_over_time = k_reverse*(A0+B0)/(k_forward+k_reverse) + (A0*k_forward - B0*k_reverse)*np.exp((-k_forward-k_reverse)*t)/(k_forward+k_reverse)
    product_over_time = A0 + B0 - reactant_over_time

    plt.plot(reactant_over_time, label = 'Reactant')
    plt.plot(product_over_time, label = 'Product')

    plt.legend()
    plt.ylabel('Concentration (moldm-3)')
    plt.xlabel('Time (s)')
    plt.title('A modelling of A<->B')

    plt.show()
    
    print('Final Concentrations:')
    
    print(f'A: {reactant_over_time[-1]:.2f} moldm-3')
    print(f'B: {product_over_time[-1]:.2f} moldm-3')
    print(f'\nThis is a 1:{(product_over_time[-1] / reactant_over_time[-1]):.0f} ratio')
    print(f'Forward k value = {k_forward}')
    print(f'Reverse k value = {k_reverse}')

def haber_process_explicit():
    # This models the Haber Process under non-isothermal conditions
    
    # Time constants
    start_time = 0  # s
    delta_time = 0.01  # s
    end_time = 50000 # s
    iterations = int(end_time / delta_time) # Iterations calculated so for loop can be used (more efficient)
    
    # Chemistry constants
    A0 = 1  # N2 mol
    B0 = 3  # H2 mol
    C0 = 0  # NH3 mol
    T0 = 450 # Starting temperature in K
    forward_Ea = 110000 # J/mol assuming iron catalyst
    reverse_Ea = 180000 # J/mol assuming iron catalyst
    A_factor = 5e4 # Pre-exponential factor
    reactant_orders = [1,3,0] # Exponential powers for the rate equation: [N2]^1, [H2]^3
    product_orders = [0,0,2]  # [NH3]^2
    coeffs = [-1,-3,2] # 1 N2 molecule lost, 3 H2 molecules lost, and 2 NH3 molecules gained per unit of the reaction
    
    # Thermodynamic and reaction chamber constants
    # These assume a pressure of 200 atm and partial equilibrium (50% H2, 25% N2, 25% NH3)
    delta_H = -92400 # J per mol of N2
    Cp = 2000 # J/kg/C Specific heat capacity (currently a constant)
    density = 60.0 # kg/m3 (currently a constant)
    volume = 1 # dm3
    surroundings_temp = 300 # K
    UA = 1.0 # W/K
    
    # Data storage
    time_values = []
    A_conc_over_time = []
    B_conc_over_time = []
    C_conc_over_time = []
    temp_over_time = []
    
    # Initialise starting values
    current_time = start_time
    current_A_conc = A0 / volume
    current_B_conc = B0 / volume
    current_C_conc = C0 / volume
    current_T = T0
    
        
    for _ in range(iterations):
        
        # Add initial values to x coords and y coords. This also repeats once the improved Euler's method calculation below has been completed.
        time_values.append(current_time)
        A_conc_over_time.append(current_A_conc)
        B_conc_over_time.append(current_B_conc)
        C_conc_over_time.append(current_C_conc)
        temp_over_time.append(current_T)
        
        # New k values
        k_forward = get_rate_constant(A_factor, forward_Ea, current_T)
        k_reverse = get_rate_constant(A_factor, reverse_Ea, current_T)

        # These two are required to make predictions of the new concs & temp
        current_rate = get_net_rate([current_A_conc,current_B_conc,current_C_conc], k_forward, k_reverse, reactant_orders, product_orders)
        delta_T = get_temperature_change(current_rate, volume, delta_H, UA, current_T, surroundings_temp, density, Cp)
        
        # Predictions made using current concentrations and rates
        pred_A = current_A_conc + current_rate * coeffs[0] * delta_time
        pred_B = current_B_conc + current_rate * coeffs[1] * delta_time
        pred_C = current_C_conc + current_rate * coeffs[2] * delta_time
        pred_T = current_T + delta_T * delta_time
        
        # Further predictions made using predicted concentrations
        k_forward_pred = get_rate_constant(A_factor, forward_Ea, pred_T)
        k_reverse_pred = get_rate_constant(A_factor, reverse_Ea, pred_T)
        pred_rate = get_net_rate([pred_A,pred_B,pred_C], k_forward_pred, k_reverse_pred, reactant_orders, product_orders)
        delta_T_pred = get_temperature_change(pred_rate, volume, delta_H, UA, pred_T, surroundings_temp, density, Cp)
        
        # Concentrations increase/decrease by the average of the previous conc and the predicted conc
        # This is Heun's method (improved Euler's method)
        current_A_conc += 0.5 * (current_rate + pred_rate) * coeffs[0] * delta_time
        current_B_conc += 0.5 * (current_rate + pred_rate) * coeffs[1] * delta_time
        current_C_conc += 0.5 * (current_rate + pred_rate) * coeffs[2] * delta_time
        current_T += 0.5 * (delta_T + delta_T_pred) * delta_time
        
        current_time += delta_time

    plt.subplot(2,1,1) # Two rows, one column, this is the first graph

    plt.plot(time_values, A_conc_over_time, label = "Concentration of N2")
    plt.plot(time_values, B_conc_over_time, label = "Concentration of H2")
    plt.plot(time_values, C_conc_over_time, label = "Concentration of NH3")
    
    plt.legend()
    plt.ylabel("Concentration (moldm-3)")
    plt.title(f"A modelling of the Haber Process starting at {T0}K")
    
    plt.subplot(2,1,2)
    
    plt.plot(time_values, temp_over_time)
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    
    plt.show()

    print(f"Final Temperature: {temp_over_time[-1]:.2f}K")
    print("\nFinal Concentrations:")
    print(f"N2: {A_conc_over_time[-1]:.2f} moldm-3")
    print(f"H2: {B_conc_over_time[-1]:.2f} moldm-3")
    print(f"NH3: {C_conc_over_time[-1]:.2f} moldm-3\n")


def haber_process_implicit():

    start_time = 0 # seconds
    end_time   = 400 # seconds

    # Chemistry constants
    A0 = 1  # N2 mol
    B0 = 3  # H2 mol
    C0 = 0  # NH3 mol
    T0 = 450 # Starting temperature in K
    forward_Ea = 110000 # J/mol assuming iron catalyst
    reverse_Ea = 180000 # J/mol assuming iron catalyst
    A_factor = 5e4 # Pre-exponential factor
    reactant_orders = [1,3,0] # Exponential powers for the rate equation: [N2]^1, [H2]^3
    product_orders = [0,0,2]  # [NH3]^2
    coeffs = [-1,-3,2] # 1 N2 molecule lost, 3 H2 molecules lost, and 2 NH3 molecules gained per unit of the reaction
    
    # Thermodynamic and reaction chamber constants
    # These assume a pressure of 200 atm and partial equilibrium (50% H2, 25% N2, 25% NH3)
    delta_H = -92400 # J per mol of N2
    Cp = 2000 # J/kg/C Specific heat capacity (currently a constant)
    density = 60.0 # kg/m3 (currently a constant)
    volume = 1 # dm3
    surroundings_temp = 650 # K
    UA = 2.0 # W/K

    def rhs(t, y):
        concs, T = y[0:3], y[3] # y = [N2, H2, NH3, T]

        k_forward = get_rate_constant(A_factor, forward_Ea, T)
        k_reverse = get_rate_constant(A_factor, reverse_Ea, T)

        rate = get_net_rate(concs, k_forward, k_reverse, reactant_orders, product_orders)

        dT_dt = get_temperature_change(rate, volume, delta_H, UA, T, surroundings_temp, density, Cp)

        return [rate * coeffs[0], # d[N2]/dt
                rate * coeffs[1], # d[H2]/dt
                rate * coeffs[2], # d[NH3]/dt
                dT_dt]            # dT/dt

    y0 = [A0/volume, B0/volume, C0/volume, T0]

    sol = solve_ivp(
        rhs,
        (start_time, end_time),
        y0,
        method='Radau',
        t_eval=np.linspace(start_time, end_time, 1000),
        rtol=1e-6,                       # relative tolerance
        atol=[1e-9, 1e-9, 1e-9, 1e-6]    # absolute tolerance (backup for rtol)
    )

    time_values      = sol.t
    A_conc_over_time = sol.y[0]
    B_conc_over_time = sol.y[1]
    C_conc_over_time = sol.y[2]
    temp_over_time   = sol.y[3]

    plt.subplot(2,1,1) # Two rows, one column, this is the first graph
    
    plt.plot(time_values, A_conc_over_time, label = "Concentration of N2")
    plt.plot(time_values, B_conc_over_time, label = "Concentration of H2")
    plt.plot(time_values, C_conc_over_time, label = "Concentration of NH3")
    
    plt.legend()
    plt.ylabel("Concentration (moldm-3)")
    plt.title(f"A modelling of the Haber Process starting at {T0}K\nwith a cooling jacket at {surroundings_temp}K")
    
    plt.subplot(2,1,2)
    
    plt.plot(time_values, temp_over_time)
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    
    plt.show()

    print(f"Final Temperature: {temp_over_time[-1]:.2f}K")
    print("\nFinal Concentrations:")
    print(f"N2: {A_conc_over_time[-1]:.2f} moldm-3")
    print(f"H2: {B_conc_over_time[-1]:.2f} moldm-3")
    print(f"NH3: {C_conc_over_time[-1]:.2f} moldm-3\n")

    

def menu():

    menu_options = { # Menu containing the different subroutines
    '1': simple_reaction,
    '2': k_comparison,
    '3': temp_comparison,
    '4': consecutive_reaction,
    '5': equilibria,
    '6': haber_process_explicit,
    '7': haber_process_implicit
    }

    start = True

    if start:
        print("1: Simple Reaction")
        print("2: K-value Comparison")
        print("3: Temperature Comparison")
        print("4: Consecutive Reaction")
        print("5: Equilibria")
        print("6: Non-Isothermal Equilibria (Implicit)")
        print("7: Non-Isothermal Equilibria (Explicit)")
        print("0: Exit")
    
        choice = input("Enter: ")
    
        if choice in menu_options:
            
            menu_options[choice]()
            
        elif choice == '0':
            print("Program exited!")
            # break
        else:
            print("Please enter a valid option.")

menu()