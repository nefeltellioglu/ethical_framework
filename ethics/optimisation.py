import ethics.model as em
import typing as tp


LossTermTransformer = tp.Callable[
    tp.Tuple[float, float, float], tp.Tuple[float, float, float]
]


def optimal_initial_condition(
    a: float,
    b: float,
    model_param_id: int,
    burden_param_id: int,
    db: dict,
    normalise: bool,
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
    normalise: whether to normalise the loss function terms.

    NOTE The "id" field in the database is use as a primary key so it
    can be used to select a unique record.
    """
    if normalise:
        # 1. Get the maximum and minimum values of the loss terms individually.
        # 2. Construct a transformation that normalises them to a range [0, 1].
        # 3. Get the optimal initial condition with normalised loss terms.

        _, max_tcb = extreme_initial_condition(
            0.0, 0.0, model_param_id, burden_param_id, db, minimise=False
        )
        _, min_tcb = extreme_initial_condition(
            0.0, 0.0, model_param_id, burden_param_id, db, minimise=True
        )
        _, max_ecb = extreme_initial_condition(
            1.0, 0.0, model_param_id, burden_param_id, db, minimise=False
        )
        _, min_ecb = extreme_initial_condition(
            1.0, 0.0, model_param_id, burden_param_id, db, minimise=True
        )
        _, max_evb = extreme_initial_condition(
            0.0, 1.0, model_param_id, burden_param_id, db, minimise=False
        )
        _, min_evb = extreme_initial_condition(
            0.0, 1.0, model_param_id, burden_param_id, db, minimise=True
        )

        print("--------------------------------------------------------")
        print("Optimising initial condition with normalised loss terms.")
        print("Normalising loss terms based on the following ranges:")
        print(f"TCB range: [{min_tcb}, {max_tcb}]")
        print(f"ECB range: [{min_ecb}, {max_ecb}]")
        print(f"EVB range: [{min_evb}, {max_evb}]")
        transform = lambda x: (
            (x[0] - min_tcb) / (max_tcb - min_tcb),
            (x[1] - min_ecb) / (max_ecb - min_ecb),
            (x[2] - min_evb) / (max_evb - min_evb),
        )

        best_ic_id, best_loss = extreme_initial_condition(
            a,
            b,
            model_param_id,
            burden_param_id,
            db,
            minimise=True,
            transform=transform,
        )
        print(f"Best loss: {best_loss}")
        print("--------------------------------------------------------")
        return best_ic_id, best_loss
    else:
        best_ic_id, best_loss = extreme_initial_condition(
            a, b, model_param_id, burden_param_id, db, minimise=True
        )
        print(f"Best loss: {best_loss}")
        return best_ic_id, best_loss


def extreme_initial_condition(
    a: float,
    b: float,
    model_param_id: int,
    burden_param_id: int,
    db: dict,
    minimise: bool = True,
    transform: LossTermTransformer = lambda x: x,
) -> tp.Tuple[int, float]:
    """
    Find the most extreme initial condition for a given model and
    burden parameters given the database and the a and b parameters
    for the loss function.

    Args:
    a: The a parameter for the loss function.
    b: The b parameter for the loss function.
    model_param_id: The identifier of the model parameters in the database.
    burden_param_id: The identifier of the burden parameters in the database.
    db: The database containing the model, burden, and outcome data.
    minimise: whether to minimise the loss function (which by default is true).
    transform: a function to transform the loss terms before computing the loss (which by default is the identity function).

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

    best_ic, best_loss = None, float("inf") if minimise else float("-inf")
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
        tmp_loss_tcb, tmp_loss_ecb, tmp_loss_evb = transform(
            em.loss_terms(oc["outcome"], tmp_ic, bp)
        )
        tmp_loss = (1 - a - b) * tmp_loss_tcb + a * tmp_loss_ecb + b * tmp_loss_evb
        if (minimise and (tmp_loss < best_loss)) or (
            (not minimise) and (tmp_loss > best_loss)
        ):
            best_loss = tmp_loss
            best_ic = tmp_ic_id

    return best_ic, best_loss
