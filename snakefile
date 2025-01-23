
rule all:
    input:
        # Unlimited vaccination scenario
        "out/grid_database-2024-10-14_manuscript.pkl",
        "out/2024-10-14_manuscript/glamorous-trajectories.png",
        "out/2024-10-14_manuscript/glamorous-loss_surfaces.png",
        "out/2024-10-14_manuscript/ab-heatmap-data.csv",
        "out/2024-10-14_manuscript/ab-heatmap-vaccination-and-clinical-burden.png",
        # Limited vaccination scenario
        "out/grid_database-2024-10-28_limited_vaccine.pkl",
        "out/2024-10-28_limited_vaccine/glamorous-trajectories.png",
        "out/2024-10-28_limited_vaccine/glamorous-loss_surfaces.png",
        "out/2024-10-28_limited_vaccine/ab-heatmap-data.csv",
        "out/2024-10-28_limited_vaccine/ab-heatmap-vaccination-and-clinical-burden.png",

rule make_grid_database_ode:
    input:
        "ethics/model.py",
        py = "create-grid-database.py",
        config = "config/config-{config_date_name}.json"
    output:
        "out/grid_database-{config_date_name}.pkl"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_glamorous_trajectories_and_loss:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "plot_trajectories_and_loss_surfaces.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/glamorous-trajectories.png",
        "out/{config_date_name}/glamorous-loss_surfaces.png",
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine"
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
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine"
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
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine"
    shell:
        """
        python {input.py} {input.config}
        """
