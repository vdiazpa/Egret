from pyomo.environ import *
from egret.data.data_utils import _read_from_file

# enter subhorizon 
subh_n = 2
#file_path = "examples/unit_commitment/tiny_uc_tc_solution.json"
file_path = "examples/unit_commitment/tiny_example.json"

model_data_dict = _read_from_file(file_path, file_type="json")

system      = model_data_dict.get("system", {})
elements    = model_data_dict.get("elements", {})
num_periods = len(system['time_keys'])

p_maxx = {}
p_t0   = {}

for key in elements['generator'].keys():
    p_maxx[key] = elements['generator'][key]['p_max']
    p_t0[key]   = elements['generator'][key]['initial_p_output']


def build_model(s_e):

    m = ConcreteModel()

    # Define Model Sets
    m.TimePeriods        = s_e
    m.ThermalGenerators  = Set(initialize=elements['generator'].keys())
    m.TransmissionLines  = Set(initialize=elements['branch'].keys())
    m.Buses              = Set(initialize=elements['bus'].keys())

    # Define Model Parameters
    m.PowerGeneratedT0   = Param(m.ThermalGenerators, initialize=p_t0)
    m.MaximumPowerOutput = Param(m.ThermalGenerators, initialize=p_maxx)

    #### Define Variables & Bounds
    #1 - Power
    def power_bounds_rule(m, g, t):
        return (0, m.MaximumPowerOutput[g])
    m.PowerGenerated        = Var(m.ThermalGenerators, m.TimePeriods, within=NonNegativeReals, bounds=power_bounds_rule) 
    #m.MaximumPowerAvailable = Var(m.ThermalGenerators, m.TimePeriods, within=NonNegativeReals)

    #2 - Binaries
    m.UnitOn    = Var(m.ThermalGenerators, m.TimePeriods, within=Binary)
    m.UnitStart = Var(m.ThermalGenerators, m.TimePeriods, within=Binary)
    m.UnitStop  = Var(m.ThermalGenerators, m.TimePeriods, within=Binary)

    #### Model Constraints
    m.coupling_constraints = ConstraintList()
    m.logical_constraints  = ConstraintList()

    def enforce_max_capacity_rule(m, g, t):
        return m.PowerGenerated[g,t] <= m.MaximumPowerOutput[g]*m.UnitOn[g,t]
    m.EnforceMaxCapacity = Constraint(m.ThermalGenerators, m.TimePeriods, rule=enforce_max_capacity_rule)

    for g in m.ThermalGenerators:
        for t in m.TimePeriods: 
            if t in m.TimePeriods and t-1 in m.TimePeriods: 
                m.logical_constraints.add(m.UnitStart[g,t] - m.UnitStop[g,t] == m.UnitOn[g,t] - m.UnitOn[g,t-1])
            else:
                m.coupling_constraints.add(m.UnitStart[g,t] - m.UnitStop[g,t] == m.UnitOn[g,t] - m.UnitOn[g,t])    

    m.pprint()

    return m

subproblems = []
for t in range(1,num_periods,subh_n):
    t_j = [i for i in range(t,t+subh_n)]
    model = build_model(t_j)
    subproblems.append(model)



#m.ThermalLimit       = Param(m.TransmissionLines, initialize=thermal_limit_dict)

# # 3 - LineFlow
# def line_bounds_rule(m, l, t):
#    return (-m.ThermalLimit[l], m.ThermalLimit[l])
# m.LinePower = Var(m.TransmissionLines, m.TimePeriods, bounds=line_bounds_rule)
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