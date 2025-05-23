import itertools
import numpy as np

# Parameters
T = 24  # number of time periods
r = 0.05  # interest rate
g = 0.02  # OD flow growth rate
s = 0.1  # savings rate for concurrent station construction
v = 24  # value of time
mu_m = 6e6  # maintenance cost per station per year
mu_o = 1e6  # operating cost per train per year
L = 2  # link length (miles)
V_rail = 45  # train speed (mph)
V_car = 30  # car speed (mph)
taccess = 12  # access time to rail station (min)
tparking = 12  # parking time for car (min)
B = 1e9  # total budget

# Station build costs
station_costs = {'C': 80e6, 'D': 80e6, 'E': 80e6, 'F': 80e6}
station_order = ['C', 'D', 'E', 'F']

# Car OD distances and flows
car_dists = {
    ('A', 'B'): L, 
    ('A', 'C'): (5/4 + (15/16)**(1/2)) * L,
    ('A', 'D'): (3/2 + (15/4)**(1/2)) * L,
    ('A', 'E'): (5/4 + (15/16)**(1/2)) * L,
    ('A', 'F'): (3/2 + (15/4)**(1/2)) * L,
    ('B', 'C'): (1/4 + (15/16)**(1/2)) * L,
    ('B', 'D'): (1/2 + (15/4)**(1/2)) * L,
    ('B', 'E'): (1/4 + (15/16)**(1/2)) * L,
    ('B', 'F'): (1/2 + (15/4)**(1/2)) * L,
    ('C', 'D'): (1/4 + (15/16)**(1/2)) * L,
    ('E', 'F'): (1/4 + (15/16)**(1/2)) * L
}
rail_dists = {
    ('A', 'B'): L,
    ('A', 'C'): 2*L, 
    ('A', 'D'): 3*L, 
    ('A', 'E'): 2*L, 
    ('A', 'F'): 3*L,
    ('B', 'C'): L, 
    ('B', 'D'): 2*L, 
    ('B', 'E'): L, 
    ('B', 'F'): 2*L, 
    ('C', 'D'): L, 
    ('E', 'F'): L
}
OD_flows_0 = {
    ('A', 'B'): 44000,
    ('A', 'C'): 41000, 
    ('A', 'D'): 13000, 
    ('A', 'E'): 46000, 
    ('A', 'F'): 36000,
    ('B', 'C'): 33000, 
    ('B', 'D'): 10000, 
    ('B', 'E'): 37000,  
    ('B', 'F'): 29000, 
    ('C', 'D'): 9000, 
    ('E', 'F'): 30000
}

# Generate valid construction schedules
def valid_schedules():
    all_scheds = itertools.product(range(T), repeat=4)
    for sched in all_scheds:
        tC, tD, tE, tF = sched
        if tD >= tC and tF >= tE and tC != tE:
            yield dict(zip(station_order, sched))

# Compute cost for a given schedule
def compute_cost(schedule):
    NPC = 0
    y = {t: {station: 0 for station in station_order} for t in range(T)}
    for station, t in schedule.items():
        for future_t in range(t, T):
            y[future_t][station] = 1

    total_spent = 0
    for t in range(T):
        # Construction costs
        c_new = sum((y[t][station] - y[t - 1][station]) * station_costs[station] if t > 0 else y[t][station] * station_costs[station] for station in station_order)
        sCD = y[t]['C'] * y[t]['D'] * (station_costs['C'] + station_costs['D']) * s
        sEF = y[t]['E'] * y[t]['F'] * (station_costs['E'] + station_costs['F']) * s
        sCDF = y[t]['C'] * y[t]['D'] * y[t]['F'] * (station_costs['F']) * s
        sDEF = y[t]['D'] * y[t]['E'] * y[t]['F'] * (station_costs['D']) * s
        Cc = c_new - (sCD + sEF)

        # Maintenance cost
        Cm = mu_m * sum(y[t].values()) * L

        # Operating cost
        k = sum(y[t].values()) - y[t]['C'] * y[t]['E'] - y[t]['D'] * y[t]['F']
        R = (2 + 2 * k) * L / V_rail
        q = sum(f0 * ((1 + g) ** t) for f0 in OD_flows_0.values())
        hAB = np.sqrt(R * mu_o / (q * v)) if q > 0 else 1
        N = round(R / hAB)
        Co = mu_o * N

        # User cost
        Cu = 0
        for (i, j), f0 in OD_flows_0.items():
            if (i == 'A') or ((i,j) == ('B', 'A')):
                h = hAB
            else:
                h = (y[t]['C']*y[t]['E'])*(2*hAB) + ((y[t]['C']-y[t]['E'])**2)*(hAB)
            ft = f0 * ((1 + g) ** t)
            served = 1 if (i not in station_order or y[t][i]) and (j not in station_order or y[t][j]) else 0
            t_car = car_dists[(i, j)] / V_car * 60 + tparking
            t_rail = taccess + h/2 + rail_dists[(i, j)] / V_rail * 60
            Cu += v * ft * (t_rail * served + t_car * (1 - served))

        # Total cost for period
        CS = Cc + Cm + Co
        total_spent += CS
        if total_spent > (t + 1) / T * B:
            return float('inf')
        NPC += (CS + Cu) / ((1 + r) ** (t + 1))
    return NPC

# Find the optimal schedule
best_sched = None
best_cost = float('inf')
for sched in valid_schedules():
    cost = compute_cost(sched)
    if cost < best_cost:
        best_cost = cost
        best_sched = sched

print(best_sched)
print("Minimum NPC: " + str(round(best_cost, 2)))