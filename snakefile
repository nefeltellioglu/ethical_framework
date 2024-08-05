
rule all:
    input:
        "out/grid_database.pkl",
        "out/vacc-vs-inf-group-1.png"



rule make_grid_database:
    input:
        py = "create-grid-database.py"
    output:
        "out/grid_database.pkl"
    shell:
        """
        python {input.py}
        """


rule plot_grid_infection_outcomes:
    input:
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
