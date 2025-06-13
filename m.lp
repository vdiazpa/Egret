\* Source Pyomo model name=unknown *\

min 
Objective:
+100.0 UnitOn(GEN1_0_t_1)
+100.0 UnitOn(GEN1_0_t_2)
+150.0 UnitOn(GEN1_1_t_1)
+150.0 UnitOn(GEN1_1_t_2)

s.t.

c_e_coupling_constraints(1)_:
-1 UnitOn(GEN1_0_t_1)
+1 UnitStart(GEN1_0_t_1)
-1 UnitStop(GEN1_0_t_1)
= -1

c_e_coupling_constraints(2)_:
-1 UnitOn(GEN1_1_t_1)
+1 UnitStart(GEN1_1_t_1)
-1 UnitStop(GEN1_1_t_1)
= -1

c_e_logical_constraints(1)_:
+1 UnitOn(GEN1_0_t_1)
-1 UnitOn(GEN1_0_t_2)
+1 UnitStart(GEN1_0_t_2)
-1 UnitStop(GEN1_0_t_2)
= 0

c_e_logical_constraints(2)_:
+1 UnitOn(GEN1_1_t_1)
-1 UnitOn(GEN1_1_t_2)
+1 UnitStart(GEN1_1_t_2)
-1 UnitStop(GEN1_1_t_2)
= 0

c_u_EnforceMaxCapacity(GEN1_0_t_1)_:
-455.0 UnitOn(GEN1_0_t_1)
+1 PowerGenerated(GEN1_0_t_1)
<= 0

c_u_EnforceMaxCapacity(GEN1_0_t_2)_:
-455.0 UnitOn(GEN1_0_t_2)
+1 PowerGenerated(GEN1_0_t_2)
<= 0

c_u_EnforceMaxCapacity(GEN1_1_t_1)_:
-555.0 UnitOn(GEN1_1_t_1)
+1 PowerGenerated(GEN1_1_t_1)
<= 0

c_u_EnforceMaxCapacity(GEN1_1_t_2)_:
-555.0 UnitOn(GEN1_1_t_2)
+1 PowerGenerated(GEN1_1_t_2)
<= 0

bounds
   0 <= UnitOn(GEN1_0_t_1) <= 1
   0 <= UnitOn(GEN1_0_t_2) <= 1
   0 <= UnitOn(GEN1_1_t_1) <= 1
   0 <= UnitOn(GEN1_1_t_2) <= 1
   0 <= UnitStart(GEN1_0_t_1) <= 1
   0 <= UnitStart(GEN1_0_t_2) <= 1
   0 <= UnitStart(GEN1_1_t_1) <= 1
   0 <= UnitStart(GEN1_1_t_2) <= 1
   0 <= UnitStop(GEN1_0_t_1) <= 1
   0 <= UnitStop(GEN1_0_t_2) <= 1
   0 <= UnitStop(GEN1_1_t_1) <= 1
   0 <= UnitStop(GEN1_1_t_2) <= 1
   0 <= PowerGenerated(GEN1_0_t_1) <= 455.0
   0 <= PowerGenerated(GEN1_0_t_2) <= 455.0
   0 <= PowerGenerated(GEN1_1_t_1) <= 555.0
   0 <= PowerGenerated(GEN1_1_t_2) <= 555.0
binary
  UnitOn(GEN1_0_t_1)
  UnitOn(GEN1_0_t_2)
  UnitOn(GEN1_1_t_1)
  UnitOn(GEN1_1_t_2)
  UnitStart(GEN1_0_t_1)
  UnitStart(GEN1_0_t_2)
  UnitStart(GEN1_1_t_1)
  UnitStart(GEN1_1_t_2)
  UnitStop(GEN1_0_t_1)
  UnitStop(GEN1_0_t_2)
  UnitStop(GEN1_1_t_1)
  UnitStop(GEN1_1_t_2)
end
