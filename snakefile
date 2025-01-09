
rule all:
    input:
        "out/grid_database-2024-10-14_manuscript.pkl",
        "out/grid_database-2024-10-28_limited_vaccine.pkl",
        "out/2024-10-14_manuscript/hm_inf_vacc.png",
        "out/2024-10-28_limited_vaccine/hm_inf_vacc.png",
        "out/grid_database-2024-12-02_limited_low_R0.pkl",
        "out/grid_database-2024-12-02_limited_high_R0.pkl",
        "out/grid_database-2024-12-02_unlimited_low_R0.pkl",
        "out/grid_database-2024-12-02_unlimited_high_R0.pkl",
        "out/2024-12-02_limited_low_R0/hm_inf_vacc.png",
        "out/2024-12-02_limited_high_R0/hm_inf_vacc.png",
        "out/2024-12-02_unlimited_high_R0/hm_inf_vacc.png",
        "out/2024-12-02_unlimited_low_R0/hm_inf_vacc.png",
        "out/2024-10-14_manuscript/trajectories.png",
        "out/2024-10-28_limited_vaccine/trajectories.png",
        "out/2024-12-02_limited_low_R0/trajectories.png",
        "out/2024-12-02_limited_high_R0/trajectories.png",
        "out/2024-12-02_unlimited_high_R0/trajectories.png",
        "out/2024-12-02_unlimited_low_R0/trajectories.png",
	
        # New plots added by AEZ on 2025-01-08
        "out/2024-10-14_manuscript/glamorous-trajectories.png",
        "out/2024-10-14_manuscript/glamorous-loss_surfaces.png",

rule make_grid_database_ode:
    input:
        "ethics/model.py",
        py = "create-grid-database.py",
        config = "config/config-{config_date_name}.json"
    output:
        "out/grid_database-{config_date_name}.pkl"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        python {input.py} {input.config}
        """

rule plot_all_optimization_results_2:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "plot-all-optimization-results-2.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/hm_inf_vacc.png",
        "out/{config_date_name}/hm_all_burdens.png",
	"out/{config_date_name}/hm_total_vacc_across_all.png"
        
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_selected_trajectories:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "plot-selected-trajectories.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/trajectories.png",
        "out/{config_date_name}/example-optimisation-Burden.png",
	"out/{config_date_name}/example-optimisation-Individual_loss.png",
	"out/{config_date_name}/example-optimisation-Normalized_individual_loss.png",
	"out/{config_date_name}/example-optimisation-Aggregated_loss.png"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        python {input.py} {input.config}
        """

rule plot_glamorous_trajectories_and_loss:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "plot-glamorous-trajectories-and-loss-surfaces.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/glamorous-trajectories.png",
        "out/{config_date_name}/glamorous-loss_surfaces.png",
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript"
    shell:
        """
        python {input.py} {input.config}
        """
