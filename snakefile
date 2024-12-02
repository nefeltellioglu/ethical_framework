
rule all:
    input:
        "out/grid_database-2024-10-14_manuscript.pkl",
        "out/grid_database-2024-10-28_limited_vaccine.pkl",
        "out/2024-10-14_manuscript/example-optimisation-results.png",
        "out/2024-10-14_manuscript/hm_vac_2_perc_across_all.png",
        "out/2024-10-14_manuscript/vacc-vs-inf-group-1.png",
        "out/2024-10-28_limited_vaccine/example-optimisation-results.png",
        "out/2024-10-28_limited_vaccine/hm_vac_2_perc_across_all.png",
        "out/2024-10-28_limited_vaccine/vacc-vs-inf-group-1.png",
        "out/grid_database-2024-12-02_limited_low_R0.pkl",
        "out/grid_database-2024-12-02_limited_high_R0.pkl",
        "out/grid_database-2024-12-02_unlimited_low_R0.pkl",
        "out/grid_database-2024-12-02_unlimited_high_R0.pkl",
        "out/2024-12-02_limited_low_R0/example-optimisation-results.png",
        "out/2024-12-02_limited_low_R0/hm_vac_2_perc_across_all.png",
        "out/2024-12-02_limited_low_R0/vacc-vs-inf-group-1.png",
        "out/2024-12-02_limited_high_R0/example-optimisation-results.png",
        "out/2024-12-02_limited_high_R0/hm_vac_2_perc_across_all.png",
        "out/2024-12-02_limited_high_R0/vacc-vs-inf-group-1.png",
        "out/2024-12-02_unlimited_high_R0/example-optimisation-results.png",
        "out/2024-12-02_unlimited_high_R0/hm_vac_2_perc_across_all.png",
        "out/2024-12-02_unlimited_high_R0/vacc-vs-inf-group-1.png",
        "out/2024-12-02_unlimited_low_R0/example-optimisation-results.png",
        "out/2024-12-02_unlimited_low_R0/hm_vac_2_perc_across_all.png",
        "out/2024-12-02_unlimited_low_R0/vacc-vs-inf-group-1.png",
        "out/2024-10-14_manuscript/tractectories.png",
        "out/2024-10-28_limited_vaccine/tractectories.png",
        "out/2024-12-02_limited_low_R0/tractectories.png",
        "out/2024-12-02_limited_high_R0/tractectories.png",
        "out/2024-12-02_unlimited_high_R0/tractectories.png",
        "out/2024-12-02_unlimited_low_R0/tractectories.png"

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


rule plot_example_optimisation_result:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        config = "config/config-{config_date_name}.json",
        py = "plot-example-optimisation-result.py"
    output:
        "out/{config_date_name}/example-optimisation-results.png"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_all_optimization_results:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "plot-all-optimization-results.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/hm_vac_2_perc_across_all.png",
        "out/{config_date_name}/hm_vac_2_perc_across_all.svg",
        "out/{config_date_name}/hm_vac_1_perc_across_all.png",
        "out/{config_date_name}/hm_vac_1_perc_across_all.svg",
        "out/{config_date_name}/hm_inf_1_perc_across_all.png",
        "out/{config_date_name}/hm_inf_1_perc_across_all.svg",
        "out/{config_date_name}/hm_inf_2_perc_across_all.png",
        "out/{config_date_name}/hm_inf_2_perc_across_all.svg",
	"out/{config_date_name}/hm_cli_burden_across_all.png",
        "out/{config_date_name}/hm_cli_burden_across_all.svg",
        "out/{config_date_name}/hm_total_vacc_across_all.png",
        "out/{config_date_name}/hm_total_vacc_across_all.svg"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_grid_infection_outcomes:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "plot-grid-outcomes-infections.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/vacc-vs-inf-group-1.png",
        "out/{config_date_name}/vacc-vs-inf-group-1.svg",
        "out/{config_date_name}/vacc-vs-inf-group-2.png",
        "out/{config_date_name}/vacc-vs-inf-group-2.svg"
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
        "out/{config_date_name}/tractectories.png",
        "out/{config_date_name}/tractectories.svg",
        
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        python {input.py} {input.config}
        """
        

