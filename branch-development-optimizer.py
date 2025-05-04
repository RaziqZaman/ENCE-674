import itertools
from math import sqrt

# Parameters (same as before)
stations = ['C', 'D', 'E', 'F']
T = 3
r = 0.03
g = 0.05
savings_rate = 0.2
mu_m = 5000
mu_o = 10000
v_user = 15
B = 5e6
taccess = 5
tparking = 5
V_rail = 40
V_car = 30
L = 10
od_flows_0 = {('C','D'): 500, ('E','F'): 600}
d_rail = {('C','D'): 10, ('E','F'): 12}
d_car = d_rail.copy()
c = {'C': 1e6, 'D': 1.2e6, 'E': 1.1e6, 'F': 1.3e6}

def present_value(amount, t):
    return amount / ((1 + r)**t)

def build_schedule_cost(schedule):
    # schedule: dict {station: period (1..T) or 0 if not built}
    total_npc = 0
    construction_costs = [0]*T
    maintenance_costs = [0]*T
    operating_costs = [0]*T
    user_costs = [0]*T

    for t in range(1, T+1):
        built_now = [s for s in stations if schedule[s] == t]
        built_so_far = [s for s in stations if 0 < schedule[s] <= t]

        # Construction cost
        const_cost = sum(c[s] for s in built_now)
        if 'C' in built_now and 'D' in built_now:
            const_cost -= (c['C'] + c['D']) * savings_rate
        if 'E' in built_now and 'F' in built_now:
            const_cost -= (c['E'] + c['F']) * savings_rate
        construction_costs[t-1] = const_cost

        # Maintenance cost
        maintenance_costs[t-1] = mu_m * len(built_so_far)

        # Operating cost
        y_CE = 1 if 'C' in built_so_far and 'E' in built_so_far else 0
        y_DF = 1 if 'D' in built_so_far and 'F' in built_so_far else 0
        k = len(built_so_far) - y_CE - y_DF
        R = (2 + 2 * k) * L / V_rail
        q = 0.5 * sum((1 + g)**t * od_flows_0[od] for od in od_flows_0)
        h_AB = sqrt(R * mu_o / (q * v_user)) if q > 0 else 1
        N = R / h_AB
        operating_costs[t-1] = mu_o * N

        # User cost
        user_c = 0
        for (i, j), f0 in od_flows_0.items():
            y_j = 1 if j in built_so_far else 0
            f_ij = f0 * (1 + g)**t
            h_ij = 2 * sqrt(R * mu_o / (f_ij * v_user)) if f_ij > 0 else 1
            t_rail = taccess + h_ij/2 + d_rail[(i, j)] / V_rail
            t_car = d_car[(i, j)] / V_car + tparking
            user_c += f_ij * (t_rail * y_j + t_car * (1 - y_j))
        user_costs[t-1] = v_user * user_c

    # Budget check
    cumulative_cost = 0
    feasible = True
    for t in range(1, T+1):
        yearly = construction_costs[t-1] + maintenance_costs[t-1] + operating_costs[t-1]
        cumulative_cost += yearly
        if cumulative_cost > (t / T) * B:
            feasible = False
            break

    if not feasible:
        return float('inf')  # Infeasible due to budget

    # Total NPC
    for t in range(1, T+1):
        total = (construction_costs[t-1] + maintenance_costs[t-1] +
                 operating_costs[t-1] + user_costs[t-1])
        total_npc += present_value(total, t)

    return total_npc

# Exhaustive search
best_schedule = None
best_cost = float('inf')

# Enumerate all (T+1)^4 = 256 schedules
for values in itertools.product(range(T+1), repeat=4):  # 0 means not built
    schedule = dict(zip(stations, values))
    cost = build_schedule_cost(schedule)
    if cost < best_cost:
        best_cost = cost
        best_schedule = schedule

# Output
print("\nBest Schedule:")
for s in stations:
    t = best_schedule[s]
    if t == 0:
        print(f"Station {s}: not built")
    else:
        print(f"Station {s}: build in period {t}")
print(f"Total Net Present Cost: ${best_cost:,.2f}")
