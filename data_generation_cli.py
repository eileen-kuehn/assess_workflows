import json
import click
import random
import logging
import hashlib
import time
import cPickle as pickle

from assess.generators.gnm_importer import CSVTreeBuilder
from assess_workflows.utils.multicoreresult import MulticoreResult
from assess_workflows.utils.utils import do_multicore

from utility.exceptions import ExceptionFrame
from utility.report import LVL

from treedistancegenerator.ted_generator import TEDGenerator
from treedistancegenerator.costs.costs import *
from treedistancegenerator.operations.random_operation import RandomOperation

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
    cost = kwargs.get("cost", True)

    tree_builder = CSVTreeBuilder()
    tree = tree_builder.build(filepath)
    if tree is not None:
        result.setdefault(filepath, {})
        result[filepath]["tree"] = tree
        result[filepath].setdefault("perturbated_tree", {})
        for probability in probabilities:
            ted_generator = TEDGenerator(costs=[TreeEditDistanceCost(),
                                                FanoutWeightedTreeEditDistanceCost(),
                                                SubtreeWeightedTreeEditDistanceCost(),
                                                SubtreeHeightWeightedTreeEditDistanceCost()] if cost else [],
                                         operation_generator=RandomOperation(insert_probability=insert_probability,
                                                                             delete_probability=delete_probability,
                                                                             edit_probability=change_probability,
                                                                             move_probability=move_probability),
                                         probability=probability)
            for _ in xrange(repeat):
                perturbated_tree = ted_generator.generate(tree)
                result[filepath]["perturbated_tree"].setdefault(probability, []).append(perturbated_tree)
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
@click.pass_context
def generate_perturbated_tree(ctx, seed, repeat, probabilities, insert_probability, cost,
                              delete_probability, change_probability, move_probability, pcount):
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
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix="DISS")
