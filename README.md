# Ethical Framework

## Usage

To run the computations

```
$ python transmission_model.py
```

To lint the code and keep documentation up to date:

```
$ bash housekeeping.sh
```

### Environment

```
conda env create -f environment.yaml
```

## `grid-search-opt` branch goals

- TODO Use grid search as the optimisation strategy
  + TODO Create a database of the model solutions to store the heavy
    compute so various optimisation question can be answered quickly
    be querying the database. See SQL example below.
  + TODO Configure the database construction with a file.
  + TODO Write an optimisation method that searches the database for
    the best initial condition.
  + TODO Generate plots to check that we are using a suitable number
    of replicates to get the average value.
  + TODO Generate plots to check that the optimal vaccination strategy
    (i.e. the initial condition) is (practically) identifiable.

The `create-grid-database.py` script iterates over a large combination
of parameters and initial conditions and computes multiple simulations
per pair. This information can be queried to work out what is the
optimal initial condition (i.e. vaccination scheme) to minimise the
loss with an arbitrary loss function. To create the database needed to
do this run the following script.

```
$ python create-grid-database.py
```

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

CREATE TABLE loss_function_configuration (
id INTEGER PRIMARY KEY,
parameter1 REAL,
parameter2 REAL
);
```
