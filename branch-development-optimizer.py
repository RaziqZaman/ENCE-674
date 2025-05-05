import itertools
import numpy as np

# Parameters
T = 24  # number of time periods
r = 0.05  # interest rate
g = 0.02  # OD flow growth rate
s = 0.1  # savings rate for concurrent station construction
v = 25  # value of time
mu_m = 1000000  # maintenance cost per station per year
mu_o = 2000000  # operating cost per train per year
L = 1  # max one-way distance for a branch line
V_rail = 45  # train speed (km/h)
V_car = 30  # car speed (km/h)
taccess = 5  # access time to rail station (min)
tparking = 10  # parking time for car (min)
B = 5e7  # total budget

# Station build costs
station_costs = {'C': 100000, 'D': 120000, 'E': 110000, 'F': 130000}
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
    ('A', 'B'): 2*L,
    ('A', 'C'): L, 
    ('A', 'D'): L, 
    ('A', 'E'): L, 
    ('A', 'F'): L,
    ('B', 'C'): L, 
    ('B', 'D'): L, 
    ('B', 'E'): L, 
    ('B', 'F'): L, 
    ('C', 'D'): L, 
    ('E', 'F'): L
}
OD_flows_0 = {
    ('A', 'B'): 100000,
    ('A', 'C'): 120000, ('A', 'D'): 140000, ('A', 'E'): 100000, ('A', 'F'): 80000,
    ('B', 'C'): 110000, ('B', 'D'): 130000, ('B', 'E'): 90000,  ('B', 'F'): 70000
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
        Cc = c_new - (sCD + sEF)

        # Maintenance cost
        Cm = mu_m * sum(y[t].values())

        # User cost
        Cu = 0
        for (i, j), f0 in OD_flows_0.items():
            ft = f0 * ((1 + g) ** t)
            served = 1 if (i not in station_order or y[t][i]) and (j not in station_order or y[t][j]) else 0
            t_car = car_dists[(i, j)] / V_car * 60 + tparking
            t_rail = taccess + rail_dists[(i, j)] / V_rail * 60
            Cu += v * ft * (t_rail * served + t_car * (1 - served))

        # Operating cost
        k = sum(y[t].values()) - y[t]['C'] * y[t]['E'] - y[t]['D'] * y[t]['F']
        R = (2 + 2 * k) * L / V_rail
        q = 0.5 * sum(f0 * ((1 + g) ** t) for f0 in OD_flows_0.values())
        hAB = np.sqrt(R * mu_o / (q * v)) if q > 0 else 1
        N = R / hAB
        Co = mu_o * N

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