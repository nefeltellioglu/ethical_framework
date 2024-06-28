import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import numpy as np
import pandas as pd

from ethical_sir import SIRParams, SIRInitialConditions, SIRSolution, optimal_initial_conditions, loss_clinical_burden, loss_equity_of_burden, loss_equity_of_vaccination, sir_vacc



# ================================

# TODO This still needs to be fixed to include the values from conmat.
params = SIRParams(0.40, 0.35, 0.35, 0.30, 0.2)
pop_size_1 = 900
pop_size_2 = 100
ts = np.arange(0, 100, 1 / 24)

# tmp_ic = optimal_initial_conditions(params, pop_size_1, pop_size_2, 0, 0)
# tmp_sol = sir_vacc(params, tmp_ic["opt_init_cond"], ts)

foo = []
# Iterate from 0.1 up to 0.45 in steps of 0.025. then do the same for a and b
for a in [0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45]:
    for b in [0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45]:
        tmp_ic = optimal_initial_conditions(params, ts, pop_size_1, pop_size_2, a, b)
        tmp_sol = sir_vacc(params, tmp_ic["opt_init_cond"], ts)
        foo.append(
            {
                "a": a,
                "b": b,
                "total_infections_1": tmp_sol.total_infections()["inf_in_1"],
                "total_infections_2": tmp_sol.total_infections()["inf_in_2"],
                "total_vaccinated_1": tmp_sol.total_vaccinated()["vacc_1"],
                "total_vaccinated_2": tmp_sol.total_vaccinated()["vacc_2"],
                "loss_clinical_burden": loss_clinical_burden(tmp_sol),
                "loss_equity_of_burden": loss_equity_of_burden(tmp_sol),
                "loss_equity_of_vaccination": loss_equity_of_vaccination(tmp_sol),
            }
        )

df = pd.DataFrame(foo)
df.to_csv("scratch-demo-1.csv")

# Do a matplotlib plot.
# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_clinical_burden` is the colour.

fig_demo_1_heatmap = "scratch-demo-1-heatmap-CB.png"
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_clinical_burden"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Clinical Burden Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_1_heatmap)



# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_equity_of_burden` is the colour.

fig_demo_2_heatmap = "scratch-demo-1-heatmap-EB.png"
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_equity_of_burden"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Equity of Burden Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_2_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_equity_of_vaccination` is the colour.

fig_demo_3_heatmap = "scratch-demo-1-heatmap-EV.png"
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_equity_of_vaccination"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Equity of Vaccination Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_3_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_clinical_burden + loss_equity_of_burden + loss_equity_of_vaccination` is the colour.

fig_demo_4_heatmap = "scratch-demo-1-heatmap-ALL.png"
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_clinical_burden"] + df["loss_equity_of_burden"] + df["loss_equity_of_vaccination"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Combined Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_4_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `total_vaccine_<X>` is the colour.

fig_demo_5_1_heatmap = "scratch-demo-1-heatmap-TV1.png"
fig_demo_5_2_heatmap = "scratch-demo-1-heatmap-TV2.png"

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_vaccinated_1"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Vaccinations 1")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_5_1_heatmap)

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_vaccinated_2"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Vaccinations 2")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_5_2_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `total_infections_<X>` is the colour.

fig_demo_6_1_heatmap = "scratch-demo-1-heatmap-TI1.png"
fig_demo_6_2_heatmap = "scratch-demo-1-heatmap-TI2.png"

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_infections_1"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Infections 1")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_6_1_heatmap)

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_infections_2"], cmap="viridis", s=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Infections 2")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_6_2_heatmap)
