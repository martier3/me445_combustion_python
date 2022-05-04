import cantera as ct
import csv
import matplotlib.pyplot as plt

gas = ct.Solution("gr130.yami")
gas.TP = 300, ct.one_atm*10
gas.set_equivalence_ratio(1, "CH4, "O2:1.0, N2:3.7")

f = ct.FreeFlame(gas)
f.set_refine_criteria(ration=3, slope=0.06, curve=0.5)
f.transport_model = "Mix"

f.solve(1, auto=True)