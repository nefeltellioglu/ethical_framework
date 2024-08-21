import ethics.model as em
import typing as tp


def optimal_initial_condition(
    a: float, b: float, model_param_id: int, burden_param_id: int, db: dict
) -> tp.Tuple[int, float]:
    """
    Find the optimal initial condition for a given model and
    burden parameters given the database and the a and b parameters
    for the loss function.

    Args:
    a: The a parameter for the loss function.
    b: The b parameter for the loss function.
    model_param_id: The identifier of the model parameters in the database.
    burden_param_id: The identifier of the burden parameters in the database.
    db: The database containing the model, burden, and outcome data.

    NOTE The "id" field in the database is use as a primary key so it
    can be used to select a unique record.
    """
    # 1. Get all the configurations that correspond to the desired model parameters.
    # 2. Get all the outcomes that correspond to the configurations from the previous step.
    # 3. Compute the loss associated with each of the outcomes given the burden parameters.
    # 4. Return the identifier of the optimal initial condition.

    configs = [
        c for c in db["configurations"] if c["model_parameters_id"] == model_param_id
    ]
    config_ids = [c["id"] for c in configs]
    ocs = [o for o in db["outcomes"] if o["configuration_id"] in config_ids]
    bp = [bp for bp in db["burden_parameters"] if bp["id"] == burden_param_id]
    assert len(bp) == 1
    bp = bp[0]["parameters"]

    best_ic, best_loss = None, float("inf")
    for oc in ocs:
        tmp_config = [c for c in configs if c["id"] == oc["configuration_id"]]
        assert len(tmp_config) == 1
        tmp_config = tmp_config[0]
        tmp_ic = [
            ic
            for ic in db["initial_conditions"]
            if ic["id"] == tmp_config["initial_condition_id"]
        ]
        assert len(tmp_ic) == 1
        tmp_ic_id = tmp_ic[0]["id"]
        tmp_ic = tmp_ic[0]["value"]
        tmp_loss = em.loss(oc["outcome"], tmp_ic, bp, a, b)
        if tmp_loss < best_loss:
            best_loss = tmp_loss
            best_ic = tmp_ic_id

    return best_ic, best_loss
