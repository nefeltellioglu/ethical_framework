
rule all:
    input:
        "out/grid_database.pkl",
        "out/vacc-vs-inf-group-1.png",
        #"out/example-optimisation-result.png"
        "out/vacc-vs-inf-group-1.png",
        "out/hm_vac_2_across_all_perc.png"


rule make_grid_database_ode:
    input:
        "ethics/model.py",
        py = "create-grid-database.py",
        config = "config/config-2024-10-28_limited_vaccine.json"
    output:
        "out/grid_database.pkl"
    shell:
        """
        python {input.py} {input.config}
        """


rule plot_example_optimisation_result:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        db = "out/grid_database.pkl",
        py = "plot-example-optimisation-result.py"
    """
    output:
        "out/example-optimisation-result.png",
        "out/example-optimisation-result.svg"

    """
    shell:
        """
        python {input.py}
        """
rule plot_all_optimization_results:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        db = "out/grid_database.pkl",
        py = "plot-all-optimization-results.py"
    output:
        "out/hm_vac_2_across_all_perc.png",
        "out/hm_vac_2_across_all_perc.svg",
        "out/hm_vac_1_across_all_perc.png",
        "out/hm_vac_1_across_all_perc.svg",
        "out/hm_inf_1_across_all_perc.png",
        "out/hm_inf_1_across_all_perc.svg",
        "out/hm_inf_2_across_all_perc.png",
        "out/hm_inf_2_across_all_perc.svg"

    shell:
        """
        python {input.py}
        """

rule plot_grid_infection_outcomes:
    input:
        "ethics/model.py",
        "ethics/optimisation.py",
        db = "out/grid_database.pkl",
        py = "plot-grid-outcomes-infections.py"
    output:
        "out/vacc-vs-inf-group-1.png",
        "out/vacc-vs-inf-group-1.svg",
        "out/vacc-vs-inf-group-2.png",
        "out/vacc-vs-inf-group-2.svg"
    shell:
        """
        python {input.py}
        """
