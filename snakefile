
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

        # New plots added by AEZ on 2025-01-09
        "out/2024-10-14_manuscript/ab-heatmap-data.csv",
        "out/2024-10-14_manuscript/ab-heatmap-vaccination-and-clinical-burden.png",
        "out/2024-10-14_manuscript/ab-heatmap-group-vaccination.png",
        "out/2024-10-14_manuscript/ab-heatmap-clinical-burden.png"

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


rule make_ab_heatmap_data:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "create-ab-heatmap-data.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/ab-heatmap-data.csv"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_ab_heatmaps_vaccination_and_clinical_burden:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/{config_date_name}/ab-heatmap-data.csv",
        py = "plot-ab-heatmaps-vaccination-and-clinical-burden.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/ab-heatmap-vaccination-and-clinical-burden.png"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_ab_heatmaps_group_vaccination:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/{config_date_name}/ab-heatmap-data.csv",
        py = "plot-ab-heatmaps-group-vaccination.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/ab-heatmap-group-vaccination.png"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_ab_heatmaps_clinical_burden:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/{config_date_name}/ab-heatmap-data.csv",
        py = "plot-ab-heatmaps-clinical-burden-at-optimal.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/ab-heatmap-clinical-burden.png"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript"
    shell:
        """
        python {input.py} {input.config}
        """
