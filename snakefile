
rule all:
    input:
        # Simulation database files
        "out/grid_database-2024-10-14_manuscript.pkl",
        "out/grid_database-2024-10-28_limited_vaccine.pkl",
        "out/grid_database-2024-12-02_limited_low_R0.pkl",
        "out/grid_database-2024-12-02_limited_high_R0.pkl",
        "out/grid_database-2024-12-02_unlimited_low_R0.pkl",
        "out/grid_database-2024-12-02_unlimited_high_R0.pkl",

        # Plots for the unlimited vaccination model in the main text
        "out/2024-10-14_manuscript/glamorous-trajectories-with-legend.svg",
        "out/2024-10-14_manuscript/glamorous-loss_surfaces.png",
        "out/2024-10-14_manuscript/ab-heatmap-data.csv",
        "out/2024-10-14_manuscript/ab-heatmap-vaccination-and-clinical-burden.png",
        "out/2024-10-14_manuscript/ab-heatmap-group-vaccination.png",
        "out/2024-10-14_manuscript/ab-heatmap-clinical-burden.png",
        "out/manuscript/2024-10-14_manuscript-trajectories-heatmaps.png",
        "out/manuscript/2024-10-14_manuscript-triangles.png",

        # Plots for the limited vaccination model in the main text
        "out/2024-10-28_limited_vaccine/glamorous-trajectories-with-legend.svg",
        "out/2024-10-28_limited_vaccine/glamorous-loss_surfaces.png",
        "out/2024-10-28_limited_vaccine/ab-heatmap-data.csv",
        "out/2024-10-28_limited_vaccine/ab-heatmap-vaccination-and-clinical-burden.png",
        "out/2024-10-28_limited_vaccine/ab-heatmap-group-vaccination.png",
        "out/2024-10-28_limited_vaccine/ab-heatmap-clinical-burden.png",
        "out/manuscript/2024-10-28_limited_vaccine-trajectories-heatmaps.png",
        "out/manuscript/2024-10-28_limited_vaccine-triangles.png",

        # Additional results included in the supplementary information.
        #   Unlimited vaccine high R0
        "out/2024-12-02_unlimited_high_R0/glamorous-trajectories-with-legend.svg",
        "out/2024-12-02_unlimited_high_R0/glamorous-loss_surfaces.png",
        "out/2024-12-02_unlimited_high_R0/ab-heatmap-data.csv",
        "out/2024-12-02_unlimited_high_R0/ab-heatmap-vaccination-and-clinical-burden.png",
        "out/2024-12-02_unlimited_high_R0/ab-heatmap-group-vaccination.png",
        "out/2024-12-02_unlimited_high_R0/ab-heatmap-clinical-burden.png",
        "out/manuscript/2024-12-02_unlimited_high_R0-trajectories-heatmaps.png",
        "out/manuscript/2024-12-02_unlimited_high_R0-triangles.png",
        #   Limited vaccine high R0
        "out/2024-12-02_limited_high_R0/glamorous-trajectories-with-legend.svg",
        "out/2024-12-02_limited_high_R0/glamorous-loss_surfaces.png",
        "out/2024-12-02_limited_high_R0/ab-heatmap-data.csv",
        "out/2024-12-02_limited_high_R0/ab-heatmap-vaccination-and-clinical-burden.png",
        "out/2024-12-02_limited_high_R0/ab-heatmap-group-vaccination.png",
        "out/2024-12-02_limited_high_R0/ab-heatmap-clinical-burden.png",
        "out/manuscript/2024-12-02_limited_high_R0-trajectories-heatmaps.png",
        "out/manuscript/2024-12-02_limited_high_R0-triangles.png",
        #   Unlimited vaccine low R0
        "out/2024-12-02_unlimited_low_R0/glamorous-trajectories-with-legend.svg",
        "out/2024-12-02_unlimited_low_R0/glamorous-loss_surfaces.png",
        "out/2024-12-02_unlimited_low_R0/ab-heatmap-data.csv",
        "out/2024-12-02_unlimited_low_R0/ab-heatmap-vaccination-and-clinical-burden.png",
        "out/2024-12-02_unlimited_low_R0/ab-heatmap-group-vaccination.png",
        "out/2024-12-02_unlimited_low_R0/ab-heatmap-clinical-burden.png",
        "out/manuscript/2024-12-02_unlimited_low_R0-trajectories-heatmaps.png",
        "out/manuscript/2024-12-02_unlimited_low_R0-triangles.png",
        #   Limited vaccine low R0
        "out/2024-12-02_limited_low_R0/glamorous-trajectories-with-legend.svg",
        "out/2024-12-02_limited_low_R0/glamorous-loss_surfaces.png",
        "out/2024-12-02_limited_low_R0/ab-heatmap-data.csv",
        "out/2024-12-02_limited_low_R0/ab-heatmap-vaccination-and-clinical-burden.png",
        "out/2024-12-02_limited_low_R0/ab-heatmap-group-vaccination.png",
        "out/2024-12-02_limited_low_R0/ab-heatmap-clinical-burden.png",
        "out/manuscript/2024-12-02_limited_low_R0-trajectories-heatmaps.png",
        "out/manuscript/2024-12-02_limited_low_R0-triangles.png",


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
        python {input.py} {input.config}
"""

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
        #python {input.py} {input.config}
"""

rule plot_glamorous_trajectories_and_loss:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        "out/grid_database-{config_date_name}.pkl",
        py = "plot-trajectories-and-loss-surfaces.py",
        config = "config/config-{config_date_name}.json",
    output:
        "out/{config_date_name}/glamorous-trajectories-with-legend.svg",
        "out/{config_date_name}/glamorous-loss_surfaces.png",
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
        #config_date_name = "2024-10-14_manuscript"
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
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
        #config_date_name = "2024-10-14_manuscript"
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
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
        #config_date_name = "2024-10-14_manuscript"
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
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
        #config_date_name = "2024-10-14_manuscript"
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
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
        #config_date_name = "2024-10-14_manuscript"
    shell:
        """
        python {input.py} {input.config}
        """

# NOTE This may look a little weird, but because of the dependency,
# you should be able to edit the SVG component file and it will
# recognise that it needs to update the image. This provides an easy
# way to edit the trajectory legend position and then recompile the
# final versions with a command rather than needing to combine things
# manually in inkscape.
rule plot_manuscript_figure_1:
    input:
        traj = "out/{config_date_name}/glamorous-trajectories-with-legend.svg",
    output:
        "out/manuscript/{config_date_name}-trajectories-heatmaps.{ext}"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        convert \\( {input.traj} out/{wildcards.config_date_name}/glamorous-loss_surfaces_global.svg -append \\) -background white -gravity center {output}
        """

rule plot_manuscript_figure_2:
    input:
        tris = "out/{config_date_name}/ab-heatmap-vaccination-and-clinical-burden.png",
    output:
        out = "out/manuscript/{config_date_name}-triangles.png"
    wildcard_constraints:
        config_date_name = "2024-10-14_manuscript|2024-10-28_limited_vaccine|2024-12-02_limited_low_R0|2024-12-02_limited_high_R0|2024-12-02_unlimited_low_R0|2024-12-02_unlimited_high_R0"
    shell:
        """
        cp {input.tris} {output.out}
        """
