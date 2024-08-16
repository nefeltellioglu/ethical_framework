import ethics.model as em

def optimal_initial_condition(a: float,
                              b: float,
                              model_param_id: int,
                              burden_param_id: int,
			      db: dict) -> int:
    """
    Find the optimal initial condition for a given model and
    burden parameters given the database and the a and b parameters
    for the objective function.

    NOTE The "id" field in the database is use as a primary key so it
    can be used to select a unique record.
    """
    # 1. Get all the configurations that correspond to the desired model parameters.
    # 2. Get all the outcomes that correspond to the configurations from the previous step.
    # 3. Compute the loss associated with each of the outcomes given the burden parameters.
    # 4. Return the identifier of the optimal initial condition.

    configs = [c for c in db["configurations"] if c["model_parameter_id"] == model_param_id]
    config_ids = [c["id"] for c in configs]
    ocs = [o for o in db["outcomes"] if o["configuration_id"] in config_ids]
    bp = [bp for bp in db["burden_parameters"] if bp["id"] == burden_param_id][0]["parameters"]

    best_ic, best_loss = None, float("inf")
    for oc in ocs:
        tmp_config = [c for c in configs if c["id"] == oc["configuration_id"]][0]
        tmp_ic = [ic for ic in db["initial_conditions"] if ic["id"] == tmp_config["initial_condition_id"]][0]
        tmp_loss = em.loss(oc["outcome"], tmp_ic, bp, a, b)
        if tmp_loss < best_loss:
            best_loss = tmp_loss
            best_ic = tmp_ic["id"]

    return best_ic
