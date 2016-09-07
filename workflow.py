import os
import json
import click
import logging
import datetime
import subprocess

from utility.report import LVL
from utility.exceptions import ExceptionFrame

import assess
from assess.generators.gnm_importer import CSVTreeBuilder, GNMCSVEventStreamer


@click.group()
@click.option("--configuration", "configuration", multiple=False,
              help="Location of configuration file", required=True)
@click.option("--start", "start", default=0, multiple=False,
              help="Start index of trees to consider for measurements.")
@click.option("--maximum", "maximum", default=float("inf"), multiple=False, metavar="INTEGER",
              help="Maximum number of trees to consider for measurements.")
@click.option("--pcount", "pcount", default=1,
              help="Number of processes to start")
@click.option("--hosts", "hosts", default="localhost",
              help="List of hosts to start calculation on")
@click.option("--json", "json", is_flag=True,
              help="Provide JSON output formatting")
@click.pass_context
def cli(ctx, configuration, start, maximum, pcount, hosts, json):
    ctx.obj["json"] = json
    ctx.obj["start"] = start
    ctx.obj["maximum"] = maximum
    configdict = {}
    execfile(configuration, configdict)
    ctx.obj["configurations"] = configdict["configurations"][:]


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
    _output_results(results=results, format_json=ctx.obj["json"])


@click.command(short_help="Calculate the distance matrix for given trees.")
@click.option("--trees", "trees", type=click.Path(), multiple=True,
              help="Path of trees to consider for pair-wise distance measurement.")
@click.option("--skip-upper", "skip_upper", is_flag=True,
              help="Skip calculations for upper part of matrix.")
@click.option("--skip-diagonal", "skip_diagonal", is_flag=True,
              help="Skip calculations for diagonal of matrix.")
@click.pass_context
def process_as_matrix(ctx, trees, skip_upper, skip_diagonal):
    results = _init_results()
    results["files"] = results["prototypes"] = trees
    tree_paths = _get_input_files(trees, minimum=ctx.obj["start"], maxlen=ctx.obj["maximum"])

    # build prototypes
    prototypes = _initialise_prototypes(tree_paths)

    def path_generator():
        for tree_index, tree_path in enumerate(tree_paths):
            maxlen = len(tree_paths)
            if skip_upper and skip_diagonal:
                maxlen = tree_index
            elif skip_upper:
                maxlen = tree_index + 1
            yield (tree_path, maxlen)

    results["results"] = _process_configurations(
        prototypes=prototypes,
        configurations=ctx.obj["configurations"],
        event_generator=path_generator
    )
    _output_results(results=results, format_json=ctx.obj["json"])


def _init_results():
    return {
        "files": None,
        "prototypes": None,
        "results": []
    }


def _output_results(results=None, format_json=False):
    if format_json:
        dump = {
            "meta": {
                "date": "%s" % datetime.datetime.now(),
                "version": subprocess.check_output(
                    ["git", "describe"],
                    cwd=os.path.dirname(assess.__file__)).strip()
            },
            "data": results
        }
        print(json.dumps(dump, indent=2))
    else:
        print(results)


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
    for event in tree:
        algorithm.add_event(event)
    algorithm.finish_tree()
    return decorator


def _get_input_files(file_paths=None, minimum=0, maxlen=float("Inf")):
    """
    Method takes a path to a file. The file might either contain the data directly, or
    it might be a file containing paths to files with actual data.

    :param file_path: Path to file to check
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

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={})
