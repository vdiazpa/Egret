#def build_custom_uc_model(data):
#    uc_form = ucgen.UCFormulation(
#        reserve_vars = 'CA_power_avail_vars',
#        ramping_limits = 'CA_ramping_limits',
#        production_costs = 'basic_production_costs_envelope',
#       uptime_downtime = 'DEKT_UT_DT',
#        startup_costs = 'CA_startup_costs',
#       network_constraints = 'ptdf_power_flow')

from pyomo.environ import *
from egret.data.data_utils import _read_from_file

file_path = "/egret/examples/unit_commitment/tiny_uc_tc_solution.json"

model_data_dict = _read_from_file(file_path, file_type="json")
system = model_data_dict['system']
elements = model_data_dict['elements']

# Replace with your actual JSON filename
model_data_dict = _read_from_file(file_path, file_type="json")

system = model_data_dict.get("system", {})
elements = model_data_dict.get("elements", {})

num_periods = len(system['time_keys'])

m = ConcreteModel()

# Define Model Sets
m.TimePeriods        = RangeSet(1, num_periods)
m.ThermalGenerators  = Set(initialize=elements['generator'].keys())
m.TransmissionLines  = Set(initialize=elements['branch'].keys())
m.Buses              = Set(initialize=elements['bus'].keys())

#Create dicts for Model Parameters
thermal_limit_dict = { br: br_data['rating_long_term'] for br, br_data in elements['branch'].items()}

def time_mapper(values):
    return {t+1: val for t, val in enumerate(values)}

pmax_dict = {}
for g_name, g_data in elements['generator'].items():
    pmax = g_data.get('p_max', None)
    if isinstance(pmax, dict) and pmax.get('data_type') == 'time_series':
        for t, val in time_mapper(pmax['values']).items():
            pmax_dict[(g_name, t)] = val
    elif isinstance(pmax, (int, float)):
        for t in range(1, num_periods + 1):
            pmax_dict[(g_name, t)] = pmax

p0_dict = { gen: gen_data['initial_p_output'] for gen, gen_data in elements['generator'].items()}

# Define Model Parameters
m.PowerGeneratedT0   = Param(m.ThermalGenerators, initialize=p0_dict)
m.ThermalLimit       = Param(m.TransmissionLines, initialize=thermal_limit_dict)
m.MaximumPowerOutput = Param(m.ThermalGenerators, m.TimePeriods, initialize=pmax_dict)

#### Define Variables & Bounds
#1 - Power
def power_bounds_rule(m, g, t):
    return (0, m.MaximumPowerOutput[g,t])
m.PowerGenerated        = Var(m.ThermalGenerators, m.TimePeriods, within=NonNegativeReals, bounds=power_bounds_rule) 
m.MaximumPowerAvailable = Var(m.ThermalGenerators, m.TimePeriods, within=NonNegativeReals)

#2 - Binaries
m.UnitOn    = Var(m.ThermalGenerators, m.TimePeriods, within=Binary)
m.UnitStart = Var(m.ThermalGenerators, m.TimePeriods, within=Binary)
m.UnitStop  = Var(m.ThermalGenerators, m.TimePeriods, within=Binary)

# 3 - LineFlow
def line_bounds_rule(m, l, t):
    return (-m.ThermalLimit[l], m.ThermalLimit[l])
m.LinePower = Var(m.TransmissionLines, m.TimePeriods, bounds=line_bounds_rule)

# 4 - Voltage Angles
m.Theta = Var(m.Buses ,m.TimePeriods, bounds = (-180,180)) #radians or degrees? 

#### Model Constraints

ref_bus = system['reference_bus']
def fix_reference_bus_rule(m, t):
    return m.Angle[ref_bus, t] == system['reference_bus_angle']
m.ReferenceBusAngle = Constraint(m.TimePeriods, rule=fix_reference_bus_rule)

def enforce_max_capacity_rule(m, g, t):
    return m.MaximumPowerAvailable[g,t] <= m.MaximumPowerOutput[g,t]*m.UnitOn[g,t]
m.EnforceMaxCapacity = Constraint(m.ThermalGenerators, m.TimePeriods, rule=enforce_max_capacity_rule)

# def enforce_max_available_ramp_up_rates_rule(m, g, t):
#        # 4 cases, split by (t-1, t) unit status (RHS is defined as the delta from m.PowerGenerated[g, t-1])
#        # (0, 0) - unit staying off:   RHS = maximum generator output (degenerate upper bound due to unit being off)
#        # (0, 1) - unit switching on:  RHS = startup ramp limit
#        # (1, 0) - unit switching off: RHS = standard ramp limit minus startup ramp limit plus maximum power output (degenerate upper bound due to unit off)
#        # (1, 1) - unit staying on:    RHS = standard ramp limit plus power generated in previous time period
#         if _ramp_up_not_needed(m,g,t):
#             return Constraint.Skip
#         if t == m.InitialTime:
#             return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedT0[g] + \
#                                                    m.ScaledNominalRampUpLimit[g,t] * m.UnitOnT0[g] + \
#                                                    m.ScaledStartupRampLimit[g,t] * (m.UnitOn[g, t] - m.UnitOnT0[g]) + \
#                                                    m.MaximumPowerOutput[g,t] * (1 - m.UnitOn[g, t])
#         else:
#             return m.MaximumPowerAvailable[g, t] <= m.PowerGenerated[g, t-1] + \
#                                                   m.ScaledNominalRampUpLimit[g,t] * m.UnitOn[g, t-1] + \
#                                                   m.ScaledStartupRampLimit[g,t] * (m.UnitOn[g, t] - m.UnitOn[g, t-1]) + \
#                                                   m.MaximumPowerOutput[g,t] * (1 - m.UnitOn[g, t])
    
#     model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_up_rates_rule)
    
    
#     # the following constraint encodes Constraint 20 defined in Carrion and Arroyo.
    
#     def enforce_ramp_down_limits_rule(m, g, t):
#         # 4 cases, split by (t-1, t) unit status:
#         # (0, 0) - unit staying off:   RHS = maximum generator output (degenerate upper bound)
#         # (0, 1) - unit switching on:  RHS = standard ramp-down limit minus shutdown ramp limit plus maximum generator output - this is the strangest case.
#         #NOTE: This may never be physically true, but if a generator has ScaledShutdownRampLimit >> MaximumPowerOutput, this constraint causes problems
#         # (1, 0) - unit switching off: RHS = shutdown ramp limit
#         # (1, 1) - unit staying on:    RHS = standard ramp-down limit
#         if _ramp_down_not_needed(m,g,t):
#             return Constraint.Skip
#         if t == m.InitialTime:
#             return m.PowerGeneratedT0[g] - m.PowerGenerated[g, t] <= \
#                  m.ScaledNominalRampDownLimit[g,t] * m.UnitOn[g, t] + \
#                  m.ScaledShutdownRampLimitT0[g]  * (m.UnitOnT0[g] - m.UnitOn[g, t]) + \
#                  m.MaximumPowerOutput[g,t] * (1 - m.UnitOnT0[g])
#         else:
#             return m.PowerGenerated[g, t-1] - m.PowerGenerated[g, t] <= \
#                  m.ScaledNominalRampDownLimit[g,t]  * m.UnitOn[g, t] + \
#                  m.ScaledShutdownRampLimit[g,t-1]  * (m.UnitOn[g, t-1] - m.UnitOn[g, t]) + \
#                  m.MaximumPowerOutput[g,t] * (1 - m.UnitOn[g, t-1])
    
#     model.EnforceScaledNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)