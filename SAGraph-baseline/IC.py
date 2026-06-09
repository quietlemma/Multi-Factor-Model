import torch_geometric.datasets as ds
import random
import ndlib.models.epidemics as ep
import ndlib.models.ModelConfig as mc

def _infected_user_count(g, status):
    return sum(
        1 for node, state in status.items()
        if state == 1 and g.nodes[node].get("node_type") == "user"
    )

def IC_multifactor(g, config, seed, rounds=100):
    result = []
    user_nodes = {
        node for node, data in g.nodes(data=True)
        if data.get("node_type") == "user"
    }
    seed = [node for node in seed if node in user_nodes]

    for iter in range(rounds):
        model_temp = ep.IndependentCascadesModel(g)
        config_temp = mc.Configuration()
        config_temp.add_model_initial_configuration('Infected', seed)

        for a, b in g.edges():
            weight = config.config["edges"]['threshold'][(a, b)]
            config_temp.add_edge_configuration('threshold', (a, b), weight)

        model_temp.set_initial_status(config_temp)

        total_no = 0
        for j in range(5):
            model_temp.iteration_bunch(1)
            total_no += _infected_user_count(g, model_temp.status)

        result.append(total_no)

    return result
