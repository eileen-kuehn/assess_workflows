import json
import math
import click
import random
import logging
import hashlib
import time
import pickle

from assess.generators.gnm_importer import CSVTreeBuilder
from assess_workflows.utils.multicoreresult import MulticoreResult
from assess_workflows.utils.utils import do_multicore

from treedistancegenerator.ted_generator import TEDGenerator
from treedistancegenerator.costs.costs import *
from treedistancegenerator.operations.random_operation import RandomOperation
from treedistancegenerator.filter.node_filter import skip_leaf, skip_inner_node, skip_no_node, \
    skip_all_but_attribute_nodes

from assess_workflows.generic.structure import Structure


@click.group()
@click.option("--basepath", "basepath", multiple=False, required=True)
@click.option("--workflow-name", "workflow_name", multiple=False, required=True)
@click.option("--step", "step", default=1, multiple=False)
@click.option("--save", "save", is_flag=True)
@click.option("--use_input", "use_input", is_flag=True,
              help="Use input file specified for current task.")
@click.pass_context
def cli(ctx, basepath, workflow_name, step, save, use_input):
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name, step=step)
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input


class DeleteAttributeTreeEditOperation(DeleteTreeEditOperation):
    def __init__(self, probability):
        self._probability = probability

    def __call__(self, node, mapping_reference=None, tree_reference=None):
        current_node = node.dao()
        del current_node["node_id"]
        samples = random.sample(range(len(node.traffic)), int(math.floor(self._probability * len(node.traffic))))
        for index in sorted(samples, reverse=True):
            current_node["traffic"].pop(index)
        parent = self._valid_parent(node, mapping_reference)
        current_tree_node = tree_reference.add_node(parent=parent, **current_node)
        try:
            current_tree_node.ppid = parent.pid
        except AttributeError:
            pass
        return [current_tree_node]


def _generate_perturbated_tree(kwargs):
    """
    :param kwargs:
    :param filepath: Path to consider
    :param probabilities: List of probabilites
    :param repeat: How often to repeat a single probability
    :param insert_probability: Probability to insert item
    :param delete_probability: Probability to delete item
    :param change_probability: Probability to change item
    :param move_probability: Probability to move item
    :param leaf_nodes_only: Only include leaf nodes?
    :param internal_nodes_only: Only include internal nodes?
    :param attribute_nodes_only: Only include attribute nodes?
    :param cost: True or False
    :return:
    """
    result = MulticoreResult()
    filepath = kwargs.get("filepath", None)
    probabilities = kwargs.get("probabilities", [])
    repeat = kwargs.get("repeat", 1)
    insert_probability = kwargs.get("insert_probability", 0)
    delete_probability = kwargs.get("delete_probability", 0)
    change_probability = kwargs.get("change_probability", 0)
    move_probability = kwargs.get("move_probability", 0)
    leaf_nodes_only = kwargs.get("leaf_nodes_only", False)
    internal_nodes_only = kwargs.get("internal_nodes_only", False)
    attribute_nodes_only = kwargs.get("attribute_nodes_only", False)
    cost = kwargs.get("cost", True)

    tree_builder = CSVTreeBuilder()
    tree = tree_builder.build(filepath)
    if tree is not None:
        result.setdefault(filepath, {})
        result[filepath]["tree"] = tree
        result[filepath].setdefault("perturbated_tree", {})
        for probability in probabilities:
            if attribute_nodes_only:
                ted_generator = TEDGenerator(costs=[],
                                             operation_generator=RandomOperation(
                                                 delete_probability=1,
                                                 delete_operation=DeleteAttributeTreeEditOperation(
                                                     probability=probability)),
                                             probability=1,
                                             skip_node=skip_all_but_attribute_nodes)
            else:
                ted_generator = TEDGenerator(costs=[TreeEditDistanceCost(),
                                                    FanoutWeightedTreeEditDistanceCost(),
                                                    SubtreeWeightedTreeEditDistanceCost(),
                                                    SubtreeHeightWeightedTreeEditDistanceCost(),
                                                    SubtreeWeightedTreeEditDistanceCostWithMove()] if cost else [],
                                             operation_generator=RandomOperation(insert_probability=insert_probability,
                                                                                 delete_probability=delete_probability,
                                                                                 edit_probability=change_probability,
                                                                                 move_probability=move_probability),
                                             probability=probability,
                                             skip_node=skip_leaf if internal_nodes_only else (
                                                 skip_inner_node if leaf_nodes_only else skip_no_node))
            for _ in range(repeat):
                perturbated_tree = ted_generator.generate(tree)
                result[filepath]["perturbated_tree"].setdefault(probability, []).append(perturbated_tree)
                # reload tree
                tree = tree_builder.build(filepath)
    return result


@click.command()
@click.option("--seed", "seed", type=int)
@click.option("--repeat", "repeat", type=int, default=1)
@click.option("--probabilities", "probabilities", type=float, multiple=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.option("--insert_probability", "insert_probability", type=float, default=0)
@click.option("--delete_probability", "delete_probability", type=float, default=0)
@click.option("--change_probability", "change_probability", type=float, default=0)
@click.option("--move_probability", "move_probability", type=float, default=0)
@click.option("--cost", "cost", default=True, type=bool)
@click.option("--leaf_nodes_only", "leaf_nodes_only", default=False, type=bool)
@click.option("--internal_nodes_only", "internal_nodes_only", default=False, type=bool)
@click.option("--attribute_nodes_only", "attribute_nodes_only", default=False, type=bool)
@click.pass_context
def generate_perturbated_tree(ctx, seed, repeat, probabilities, insert_probability, cost,
                              delete_probability, change_probability, move_probability, pcount,
                              leaf_nodes_only, internal_nodes_only, attribute_nodes_only):
    if seed is not None:
        random.seed(seed)
    results = MulticoreResult()
    if ctx.obj.get("use_input"):
        structure = ctx.obj.get("structure", None)
        with open(structure.input_file_path(), "r") as input_file:
            json_data = json.load(input_file)
            samples = json_data["data"]["samples"]
            if pcount > 1:
                data = [{
                    "filepath": item,
                    "repeat": repeat,
                    "probabilities": probabilities,
                    "insert_probability": insert_probability,
                    "delete_probability": delete_probability,
                    "change_probability": change_probability,
                    "move_probability": move_probability,
                    "leaf_nodes_only": leaf_nodes_only,
                    "internal_nodes_only": internal_nodes_only,
                    "attribute_nodes_only": attribute_nodes_only,
                    "cost": cost
                } for sample in samples for item in sample]
                multicore_results = do_multicore(
                    count=pcount,
                    target=_generate_perturbated_tree,
                    data=data
                )
                for result in multicore_results:
                    results += result
            else:
                for sample in samples:
                    for item in sample:
                        results += _generate_perturbated_tree({
                            "filepath": item,
                            "repeat": repeat,
                            "probabilities": probabilities,
                            "insert_probability": insert_probability,
                            "delete_probability": delete_probability,
                            "change_probability": change_probability,
                            "move_probability": move_probability,
                            "leaf_nodes_only": leaf_nodes_only,
                            "internal_nodes_only": internal_nodes_only,
                            "attribute_nodes_only": attribute_nodes_only,
                            "cost": cost
                        })
        if ctx.obj.get("save"):
            # instead of storing all results as one, we split them per base tree
            # a header is used to map all individual stores
            results_header = {}
            for name, result in results.items():
                nick = '%s%02s%s' % (hashlib.sha1(name).hexdigest(), random.getrandbits(8), time.strftime('%H%M%S'))
                with open(structure.intermediate_file_path(file_type="pkl", variant=nick), "w") as output_file:
                    results_header[name] = output_file.name
                    pickle.dump(
                        MulticoreResult({name: result}),
                        output_file
                    )
            with open(structure.intermediate_file_path(file_type="pkl"), "w") as output_file:
                    pickle.dump(results_header, output_file)


cli.add_command(generate_perturbated_tree)

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.WARNING)
    logging.getLogger("EXCEPTION").setLevel(logging.INFO)
    cli(obj={}, auto_envvar_prefix="DISS")
