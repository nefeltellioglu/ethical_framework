# Ethical Framework

## Usage

To run the computations

```
$ python scratch-fancy-model.py
```

The results of this are output to `scratch-fancy-ODE.csv`.

To lint the code and keep documentation up to date:

```
$ bash housekeeping.sh
```

### Environment

```
conda env create -f environment.yaml
```

### Snakemake pipeline

There is a `snakefile` to configure building the grid-search database
and to generate some basic plots. To run the whole pipeline on a
single core use the following command:

```
snakemake -c1 all
```

## Configuration

To avoid magic numbers appearing thoughout the code, there is a single
file `config/config-YYY-MM-DD.json` which should be used as the single
source of truth regarding parameter values. If you want to change one
of the parameters, it should simply be a matter of adjusting this
file.

## `grid-search-opt` branch goals

- TODO Double check that the percentages in the `BurdenParams` class
  can safely be renamed to proportions as this is what they appear to
  be.
- TODO Work out if we should be thinking about the Pareto front?
- TODO Use grid search as the optimisation strategy
  + DONE Create a database of the model solutions to store the heavy
    compute so various optimisation question can be answered quickly
    be querying the database. See SQL example below.
  + DONE Rename `ethical_sir_fancy.py` to `ethics/model.py` or
    whatever allows us to import just the modelling code as
    `ethics.model`. Then there can be a `ethics.optimisation` for the
    optimisation.
  + DONE Configure the database construction with a file. This will
    provide a nicer alternative to having hard-coded parameters in the
    code. This should live in `config/`.
  + DONE Write an optimisation method that searches the database for
    the best initial condition.
  + DONE Generate plots to check that the optimal vaccination strategy
    (i.e. the initial condition) is (practically) identifiable.
  + TODO Implement a normalisation strategy so that all of the
    objectives are on a comparable scale. We could do this by looking
    at the vertices of the a/b simplex and then scaling each loss term
    to take values form 0 to 1.
  + TODO Generate plots to check that we are using a suitable number
    of replicates to get the average value.

The `create-grid-database.py` script iterates over a large combination
of parameters and initial conditions and computes multiple simulations
per pair. This information can be queried to work out what is the
optimal initial condition (i.e. vaccination scheme) to minimise the
loss with an arbitrary loss function. To create the database needed to
do this run the `create-grid-database.py` script.

### Database creation example

This script would create a database containing the model simulation
results. We aren't using a database, but it might be helpful to have
this here as a model for how the data is stored.

```sql
CREATE TABLE model_parameters (
id INTEGER PRIMARY KEY,
parameter1 REAL,
parameter2 REAL
);

CREATE TABLE initial_condition (
id INTEGER PRIMARY KEY,
condition1 REAL,
condition2 REAL
);

CREATE TABLE configuration (
id INTEGER PRIMARY KEY,
model_parameters_id INTEGER,
initial_condition_id INTEGER,
);

CREATE TABLE outcome (
id INTEGER PRIMARY KEY,
configuration_id INTEGER,
seed INTEGER,
outcome1 REAL,
outcome2 REAL
);

CREATE TABLE burden_parameters (
id INTEGER PRIMARY KEY,
parameter1 REAL,
parameter2 REAL
);
```

Note that we have split the `outcome` and the `configuration` table
because when using a stochastic model the outcomes are random
variables and this allows us to reference a shared configuration
rather than replicating that data in each outcome. It adds a bit of
complexity for the deterministic model case but should simplify things
in the stochastic case.

### Optimisation

The optimisation function takes an identifier for model and burden
parameters, and the database. It returns the identifier of the optimal
initial condition.

```python
def optimal_initial_condition(model_param_id: int,
                              burden_param_id: int,
							  db: dict) -> int:
```
