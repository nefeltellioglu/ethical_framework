

def optimal_initial_condition(a: float,
                              b: float,
                              model_param_id: int,
                              burden_param_id: int,
			      db: dict) -> int:
    """
    Find the optimal initial condition for a given model and
    burden parameters given the database and the a and b parameters
    for the objective function.
    """
    # 1. Get all the configurations that correspond to the desired model parameters.
    # 2. Get all the outcomes that correspond to the configurations from the previous step.
    # 3. Compute the loss associated with each of the outcomes given the burden parameters.
    # 4. Return the identifier of the optimal initial condition.
    raise NotImplementedError
