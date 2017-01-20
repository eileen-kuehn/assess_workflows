from __future__ import print_function
import os
import sys
import json
import math
import click
import logging
import importlib
import cPickle as pickle

from assess.exceptions.exceptions import EventNotSupportedException
from assess.generators.event_generator import NodeGenerator
from assess_workflows.utils.multicoreresult import MulticoreResult
from utility.report import LVL
from utility.exceptions import ExceptionFrame

import assess
from assess.generators.gnm_importer import CSVTreeBuilder, GNMCSVEventStreamer
from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import output_results, determine_version, do_multicore
from assess.decorators.decorator import Decorator


@click.group()
@click.option("--basepath", "basepath", multiple=False, required=True)
@click.option("--workflow-name", "workflow_name", multiple=False, required=True)
@click.option("--step", "step", default=1, multiple=False)
@click.option("--configuration", "configuration", multiple=False,
              help="Location of configuration file")
@click.option("--start", "start", default=0, multiple=False,
              help="Start index of trees to consider for measurements.")
@click.option("--maximum", "maximum", default=float("inf"), multiple=False, metavar="INTEGER",
              help="Maximum number of trees to consider for measurements.")
@click.option("--json", "json", is_flag=True,
              help="Provide JSON output formatting")
@click.option("--save", "save", is_flag=True)
@click.option("--use_input", "use_input", is_flag=True,
              help="Use input file specified for current task.")
@click.pass_context
def cli(ctx, basepath, workflow_name, step, configuration, start, maximum, json, save, use_input):
    ctx.obj["json"] = json or save
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input
    ctx.obj["start"] = start
    ctx.obj["maximum"] = maximum
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name, step=step)
    config_module = ctx.obj["structure"].configuration_module()
    importlib.import_module(config_module)
    configdict = sys.modules[config_module]
    ctx.obj["configurations"] = configdict.configurations


@click.command(short_help="Calculate the distance vector for given data.")
@click.option("--trees", "trees", type=click.Path(), multiple=True,
              help="Path of trees to consider for distance measurement.")
@click.option("--representatives", "representatives", type=click.Path(), multiple=True,
              help="Path of representatives to measure distance to.")
@click.pass_context
def process_as_vector(ctx, trees, representatives):
    results = _init_results()
    results["files"] = trees
    results["prototypes"] = representatives
    tree_paths = _get_input_files(trees, minimum=ctx.obj["start"], maxlen=ctx.obj["maximum"])
    prototype_paths = _get_input_files(representatives)

    # build prototypes
    prototypes = _initialise_prototypes(prototype_paths)

    def path_generator():
        for path in tree_paths:
            yield (path, len(prototypes))

    results["results"] = _process_configurations(
        prototypes=prototypes,
        configurations=ctx.obj["configurations"],
        event_generator=path_generator
    )
    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess.__file__)),
        source="%s (%s)" % (__file__, "process_as_vector")
    )


@click.command()
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def batch_process_from_pkl(ctx, pcount):
    results = _init_results()
    results["distance"] = []
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path(file_type="pkl")
        with open(file_path, "r") as input_file:
            if pcount > 1:
                data = []
                # HEADER
                # ######
                # results are split into a header and data files
                # see data_generation_cli.generate_perturbated_tree
                tree_metadata = pickle.load(input_file)
                results["files"] = tree_metadata.keys()
                for key, pkl_path in tree_metadata.items():
                    # tree is stored in "tree"
                    # distorted trees in "perturbated_tree"
                    data.append({
                        "data_pkl_path": pkl_path,
                        "data_pkl_key": key,
                        "configurations": ctx.obj["configurations"]
                    })
                result_list = (
                    do_multicore(
                        count=pcount,
                        data=data,
                        target=_process_configurations_for_row
                    )
                )
                for result in result_list:
                    results["results"].append(result['results'])
                    results["distance"].append(result['precalculated_costs'])
            else:
                raise NotImplementedError
    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess.__file__)),
        source="%s (%s)" % (__file__, "batch_process_from_pkl")
    )


@click.command()
@click.pass_context
def batch_process_as_vector(ctx):
    results = []

    if ctx.obj.get("use_input", False):
        def path_generator():
            for path in value[1:]:
                yield (path, 1)

        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data")
            for key, values in input_data.items():
                for value in values:
                    results.append(_init_results())
                    if len(value) == 1:
                        # element is file and prototype at the same time
                        value.append(value[0])
                    # first element is tree, second is representative
                    results[-1]["files"] = value[:1]
                    results[-1]["prototypes"] = value[1:]
                    results[-1]["results"] = _process_configurations(
                        prototypes=_initialise_prototypes(value[:1]),
                        configurations=ctx.obj["configurations"],
                        event_generator=path_generator
                    )
                    results[-1]["key"] = key

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess.__file__)),
        source="%s (%s)" % (__file__, "batch_process_as_vector")
    )


@click.command(short_help="Calculate the distance matrix for given trees.")
@click.option("--trees", "trees", type=click.Path(), multiple=True,
              help="Path of trees to consider for pair-wise distance measurement.")
@click.option("--skip_upper", "skip_upper", default=False,
              help="Skip calculations for upper part of matrix.")
@click.option("--skip_diagonal", "skip_diagonal", default=False,
              help="Skip calculations for diagonal of matrix.")
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def process_as_matrix(ctx, trees, skip_upper, skip_diagonal, pcount):
    if len(trees) == 0 and ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            # can be list of lists or flat list
            trees = json.load(input_file).get("data").values()[0]
    results = _init_results()
    results["files"] = results["prototypes"] = trees[:]

    # if we have a flat list, check, otherwise, just take it
    if type(trees[0]) == list:
        tree_paths = trees
        nested = True
    else:
        tree_paths = _get_input_files(trees, minimum=ctx.obj["start"], maxlen=ctx.obj["maximum"])
        nested = False

    if pcount > 1:
        to_process = []
        if nested:
            to_process = tree_paths
        else:
            to_process.append(tree_paths)
        while to_process:
            data = []
            single_tree_paths = to_process.pop(0)
            # prepare blocks of data
            block_size = len(single_tree_paths) / float(pcount) #/ 4.0
            assert block_size > 1, "Blocksize is too small for proper parallelisation: %s" % block_size
            index_value = int(math.ceil(len(single_tree_paths) / block_size))
            for row_idx in range(index_value):
                for col_idx in range(index_value):
                    if skip_upper and col_idx > row_idx:
                        continue
                    row_trees = single_tree_paths[int(row_idx * block_size):min(
                        int((row_idx + 1) * block_size), len(single_tree_paths))]
                    col_trees = single_tree_paths[int(col_idx * block_size):min(
                        int((col_idx + 1) * block_size), len(single_tree_paths))]
                    data.append({
                        "tree_paths": row_trees,
                        "prototype_paths": col_trees,
                        "configurations": ctx.obj["configurations"]
                    })

            result_list = do_multicore(
                count=pcount,
                target=_process_as_matrix,
                data=data
            )
            final_decorators = []
            row_idx = 0
            col_idx = -1
            for result_index, result_entry in enumerate(result_list):
                # calculate the exact position within matrix to help decorators updating their results
                col_idx += 1
                if col_idx >= ((row_idx + 1) if skip_upper else index_value):
                    row_idx += 1
                    col_idx = 0
                current_results = result_entry.get("results", [])
                # each of the results has the same configuration of decorators, so we can get one
                # exemplary list of decorators to process all results
                for decorator_key in current_results[0].get("decorator", {}):
                    for index, current_result in enumerate(current_results):
                        try:
                            # if decorator already exists, we only need to add current data
                            decorator = final_decorators[index][decorator_key]
                            current_decorator = type(decorator)()
                            current_decorator._data = current_result.get("decorator", {})[decorator_key]
                            current_decorator.row_idx = [row_idx]
                            current_decorator.col_idx = [col_idx]
                            decorator += current_decorator
                        except (IndexError, KeyError):
                            # if decorator does not exist, we load it and will later add data
                            decorator = Decorator.from_name(decorator_key)
                            decorator._data = current_result.get("decorator", {})[decorator_key]
                            decorator.row_idx = [row_idx]
                            decorator.col_idx = [col_idx]
                            try:
                                final_decorators[index].setdefault(decorator_key, decorator)
                            except IndexError:
                                final_decorators.append({decorator_key: decorator})
            # format updated data
            finals = result_list[0]
            for index, final in enumerate(finals.get("results", [])):
                for value in final_decorators[index].values():
                    data = value.descriptive_data()
                    final.get("decorator", {})[data.keys()[0]] = data.values()[0]
            results.setdefault("results", []).append(finals["results"])
    else:
        to_process = []
        if nested:
            to_process = tree_paths
        else:
            to_process.append(tree_paths)
        while to_process:
            single_tree_paths = to_process.pop(0)
            # build prototypes
            prototypes = _initialise_prototypes(single_tree_paths)

            def path_generator():
                for tree_index, tree_path in enumerate(single_tree_paths):
                    maxlen = len(single_tree_paths)
                    if skip_upper and skip_diagonal:
                        maxlen = tree_index
                    elif skip_upper:
                        maxlen = tree_index + 1
                    yield (tree_path, maxlen)

            results.setdefault("results", []).append(_process_configurations(
                prototypes=prototypes,
                configurations=ctx.obj["configurations"],
                event_generator=path_generator
            ))

    output_results(
        ctx=ctx,
        results=results,
        version=os.path.dirname(assess.__file__),
        source="%s (%s)" % (__file__, "process_as_matrix")
    )


# TODO: add configurations
def _process_as_matrix(args):
    """
    :param tree_paths:
    :param prototype_paths:
    :param configurations:
    :return:
    """
    results = {}
    configurations = args.get("configurations", [])
    tree_paths = args.get("tree_paths", [])
    prototypes = _initialise_prototypes(args.get("prototype_paths", []))

    def path_generator():
        for tree_path in tree_paths:
            yield (tree_path, len(prototypes))

    if len(tree_paths) == 0 or len(prototypes) == 0:
        # skip stuff where nothing needs to be calculated...
        return results
    results["results"] = _process_configurations(
        prototypes=prototypes,
        configurations=configurations,
        event_generator=path_generator
    )
    return results


def _init_results():
    return {
        "files": None,
        "prototypes": None,
        "results": []
    }


def _process_configurations_for_row(kwargs):
    results = []
    precalculated_costs = {}
    with ExceptionFrame():
        data_pkl_path = kwargs.pop('data_pkl_path')
        data_pkl_key = kwargs.pop('data_pkl_key')
        with open(data_pkl_path) as input_pkl:
            raw_data = pickle.load(input_pkl)[data_pkl_key]
        tree = raw_data['tree']
        prototypes = [ptree for ptrees in raw_data["perturbated_tree"].values() for ptree in ptrees]
        # extract precalculated distances
        for perturbated_tree in prototypes:
            # append precalculated distances
            for cost_model, distance in perturbated_tree.distance.items():
                precalculated_costs.setdefault(cost_model.__class__.__name__, []).append(distance)
        # calculate actual results
        configurations = kwargs.get("configurations", {})
        for configuration in configurations:
            for algorithm_def in configuration.get("algorithms", []):
                for signature_def in configuration.get("signatures", []):
                    signature = signature_def()
                    algorithm = algorithm_def(signature=signature)
                    algorithm.prototypes = prototypes
                    decorator = None
                    streamer = tree
                    result = _perform_calculation(
                        tree=streamer,
                        algorithm=algorithm,
                        decorator_def=configuration.get("decorator", None),
                        maxlen=None
                    )
                    if decorator:
                        decorator.update(result)
                    else:
                        decorator = result
                    if decorator:
                        results.append({
                            "algorithm": "%s" % algorithm,
                            "signature": "%s" % signature,
                            "event_streamer": "%s" % GNMCSVEventStreamer(csv_path=None),
                            "decorator": decorator.descriptive_data()
                        })
    return {'results': results, 'precalculated_costs': precalculated_costs}


def _process_configurations(prototypes, configurations, event_generator):
    """
    Method parses all different possibilities of calculations that need to be performed based
    on given configurations. The calculations itself are performed on the prototypes and the
    trees (as well as info how many of the prototypes should be processed) that are passed via
    a generator. Method returns the collected results.

    :param prototypes: List of prototypes that are used on the algorithm
    :param configurations: List of configurations to parse
    :param event_generator: Generator returning tuples (tree_path, maxlen)
    :return: List of results
    """
    results = []
    for configuration in configurations:
        for algorithm_def in configuration.get("algorithms", []):
            for signature_def in configuration.get("signatures", []):
                for event_streamer in configuration.get("event_streamer", [GNMCSVEventStreamer]):
                    signature = signature_def()
                    algorithm = algorithm_def(signature=signature)
                    algorithm.prototypes = prototypes
                    decorator = None
                    streamer = None
                    for (event, maxlen) in event_generator():
                        streamer = event_streamer(csv_path=event)
                        result = _perform_calculation(
                            tree=streamer,
                            algorithm=algorithm,
                            decorator_def=configuration.get("decorator", None),
                            maxlen=maxlen
                        )
                        if decorator:
                            decorator.update(result)
                        else:
                            decorator = result
                    if decorator:
                        results.append({
                            "algorithm": "%s" % algorithm,
                            "signature": "%s" % signature,
                            "event_streamer": "%s" % streamer if streamer is not None
                            else event_streamer(csv_path=None),
                            "decorator": decorator.descriptive_data()
                        })
    return results


def _perform_calculation(tree, algorithm, decorator_def, maxlen=float("Inf")):
    """
    Method performs the actual calculation given an event stream (tree), the configured algorithm
    to be used as well as the definition for decorators to be created.

    :param tree: Event generator based on a tree
    :param algorithm: Configured algorithm containing prototypes as well as signature
    :param decorator_def: Description of decorators to be created
    :param maxlen: Maximum length that is considered for prototypes
    :return: Decorator containing requested results
    """
    decorator = decorator_def()
    decorator.wrap_algorithm(algorithm=algorithm)
    # starting a new tree not to mix former results with current
    algorithm.start_tree(maxlen=maxlen)
    for event in tree.event_iter():
        try:
            algorithm.add_event(event)
        except EventNotSupportedException:
            pass
    algorithm.finish_tree()
    return decorator


def _get_input_files(file_paths=None, minimum=0, maxlen=float("Inf")):
    """
    Method takes a path to a file. The file might either contain the data directly, or
    it might be a file containing paths to files with actual data.

    :param file_paths: Paths to files to check
    :param minimum: Index to start reading paths from
    :param maxlen: Maximum amount of trees to return
    :return: List of data file paths
    """
    result = []
    for file_path in file_paths:
        # check if file contains list of paths or data
        if len(result) < maxlen:
            with open(file_path, "r") as input_file:
                # check first line of file
                index = 0
                line = next(input_file)
                if "#" in line[0]:
                    # I guess it is a valid CSV file, starting with comment
                    result.append(file_path)
                elif len(line.split(",")) > 4:
                    # found CSV formatted strings
                    result.append(file_path)
                else:
                    if index >= minimum and len(result) < maxlen:
                        result.append(line.strip())
                    for index, line in enumerate(input_file):
                        if index >= minimum and len(result) < maxlen:
                            result.append(line.strip())
                        else:
                            break
    return result


def _initialise_prototypes(prototype_paths):
    """
    Method initialises the prototype trees from given file paths.

    :param prototype_paths: List of paths to prototypes
    :return: List of trees
    """
    prototypes = []
    tree_builder = CSVTreeBuilder()
    for prototype_path in prototype_paths:
        prototypes.append(tree_builder.build(prototype_path))
    return prototypes


cli.add_command(process_as_vector)
cli.add_command(process_as_matrix)
cli.add_command(batch_process_as_vector)
cli.add_command(batch_process_from_pkl)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
