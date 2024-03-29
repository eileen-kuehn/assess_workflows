import os
import glob
import json
import pickle

import click
import random
import logging
import subprocess
import pandas as pd

from assess.exceptions.exceptions import TreeInvalidatedException
from assess.prototypes.simpleprototypes import Prototype

from assess_workflows.utils.multicoreresult import MulticoreResult
try:
    from dbutils.sqlcommand import SQLCommand
except:
    pass
from gnmutils.exceptions import DataNotInCacheException

from assess.generators.gnm_importer import CSVTreeBuilder
from assess.events.events import ProcessStartEvent, ProcessExitEvent, TrafficEvent

import assess_workflows
from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import output_results, determine_version, do_multicore


@click.group()
@click.option("--basepath", "basepath", multiple=False, required=True)
@click.option("--workflow-name", "workflow_name", multiple=False, required=True)
@click.option("--prototype-cache-preloaded-only", "prototype_cache_preloaded_only", default=False)
@click.option("--step", "step", default=1, multiple=False)
@click.option("--save", "save", is_flag=True)
@click.option("--use_input", "use_input", is_flag=True,
              help="Use input file specified for current task.")
@click.pass_context
def cli(ctx, basepath, workflow_name, step, save, use_input, prototype_cache_preloaded_only):
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name, step=step)
    ctx.obj["json"] = True
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input
    ctx.obj["prototype_cache_preloaded_only"] = prototype_cache_preloaded_only


def _relevant_files_for_context(ctx, path):
    filenames = []
    if ctx.obj.get("prototype_cache_preloaded_only", False):
        filenames.extend([found_file.replace("_prototypes.pkl", ".csv") for found_file in
                          glob.glob("%s/*/*-process_prototypes.pkl" % path)])
    else:
        filenames.extend(glob.glob("%s/*/*-process.csv" % path))
    return filenames


@click.command()
@click.option("--trees", "trees", type=click.Path(), multiple=True,
              help="Path of trees to consider for distance measurement.")
@click.option("--representatives", "representatives", type=click.Path(), multiple=True,
              help="Path of representatives to measure distance to.")
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def index_valid_hdf_trees(ctx, trees, representatives, pcount):
    structure = ctx.obj.get("structure", None)
    results = {}
    paths = [(key, value) for key, values in {
        "trees": trees,
        "representatives": representatives
    }.items() for value in values if values]
    if pcount > 1:
        trees = do_multicore(
            count=pcount,
            target=_valid_hdf_tree,
            data=paths
        )
        for category, tree, name in trees:
            results.setdefault(category, []).append(_write_tree_to_pkl(structure, tree, name))
    else:
        for filename in paths:
            trees = _valid_hdf_tree(filename)
            for category, tree, name in trees:
                results.setdefault(category, []).append(_write_tree_to_pkl(structure, tree, name))
    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_valid_hdf_trees")
    )


def _write_tree_to_pkl(structure, tree, name):
    pkl_filename = structure.intermediate_file_path(variant=str(name), file_type="pkl")
    with open(pkl_filename, "wb") as pkl_file:
        pickle.dump(tree, pkl_file)
    return pkl_filename


def _valid_hdf_tree(args):
    """
    :param args: Tuple from category and file path
    :return: Tuple from category, tree, and tree name
    """
    category, filename = args
    results = []
    df = pd.read_hdf(filename, key="train_events")
    label = df.index.get_level_values(0)[0]
    events = df.index.get_level_values(1).unique()
    for event in events:
        tree_data = df.xs((label, event), level=('label', 'evtNum'))
        tree = Prototype()
        last_node = None
        for node_id, values in tree_data.iterrows():
            # add tme values
            tree_values = dict(values)
            tree_values["tme"] = tree_values["exit_tme"] = 0
            parent_node_id = int(values.pop("motherIndex"))
            name = values.pop("PDG")
            if last_node is not None:
                node = tree.add_node(
                    name,
                    parent_node_id=parent_node_id,
                    node_id=node_id,
                    **tree_values
                )
            else:
                node = tree.add_node(
                    name,
                    node_id=node_id,
                    **tree_values
                )
            last_node = node
        if tree:
            results.append((category, tree, "%s-%s" % (
                os.path.basename(filename).split(".")[0], event),))
    return results


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def index_valid_trees(ctx, paths, pcount):
    """
    Method walks the given paths and reads all trees that are found within. For each tree that
    can successfully be read, it is appended to the results list. This list can be used for
    further processing.

    :param ctx: Click context
    :param paths: The paths to scan for valid tree data
    """
    results = []
    filenames = []
    for path in paths:
        filenames.extend(_relevant_files_for_context(ctx, path))
    if pcount > 1:
        results.extend(do_multicore(
            count=pcount,
            target=_valid_tree,
            data=filenames
        ))
    else:
        for filename in filenames:
            result = _valid_tree(filename)
            if result is not None:
                results.append(result)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_valid_trees"))


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def index_data_by_tme(ctx, paths, pcount):
    results = MulticoreResult()
    filenames = []
    for path in paths:
        filenames.extend(_relevant_files_for_context(ctx, path))
    if pcount > 1:
        result_list = do_multicore(
            count=pcount,
            target=_data_by_tme,
            data=filenames
        )
        for result in result_list:
            results += result
    else:
        for filename in filenames:
            results += _data_by_tme(filename)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_data_by_tme"))


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def index_data_by_uid(ctx, paths, pcount):
    results = MulticoreResult()
    filenames = []
    for path in paths:
        filenames.extend(_relevant_files_for_context(ctx, path))
    if pcount > 1:
        result_list = do_multicore(
            count=pcount,
            target=_data_by_uid,
            data=filenames
        )
        for result in result_list:
            results += result
    else:
        for filename in filenames:
            results += _data_by_uid(filename)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_data_by_uid"))


@click.command()
@click.option("--paths", "paths", multiple=True)
@click.pass_context
def index_data_by_number_of_nodes(ctx, paths):
    results = {}
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            paths = list(paths)
            paths.extend(json.load(input_file)["data"])
    filenames = []
    for path in paths:
        filenames.extend(_relevant_files_for_context(ctx, path))
    for filename in filenames:
        if not os.path.isfile(path):
            count = _line_count(filename)
            results.setdefault(count, []).append(filename)
        else:
            count = _line_count(path)
            results.setdefault(count, []).append(path)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_data_by_number_of_nodes"))


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.pass_context
def index_data_by_number_of_traffic_events(ctx, paths):
    results = {}
    for path in paths:
        for filename in glob.glob("%s/*/*-traffic.csv" % path):
            count = _line_count(filename)
            results.setdefault(count, []).append(filename.replace("-traffic.csv", "-process.csv"))

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_data_by_number_of_traffic_events")
    )


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.pass_context
def index_data_by_number_of_events(ctx, paths):
    results = {}
    filenames = []
    for path in paths:
        filenames.extend(_relevant_files_for_context(ctx, path))
    for filename in filenames:
        basename = os.path.basename(filename)
        db_id = basename.split("-")[0]
        # access process file
        count = _line_count(filename=filename) * 2  # start and finishing of a process
        # access traffic file (if existent)
        traffic_count = _line_count(
            filename=os.path.join(os.path.dirname(filename), "%s-traffic.csv" % db_id))
        results.setdefault(count + traffic_count, []).append(filename)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_data_by_number_of_events")
    )


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.pass_context
def index_data_by_number_of_payloads(ctx, paths):
    results = {}
    for path in paths:
        for filename in glob.glob("%s/*/*-process.csv" % path):
            count = _payload_line_count(filename)
            results.setdefault(count, []).append(filename)
    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_data_by_number_of_payloads")
    )


@click.command()
@click.pass_context
def index_data_by_activity(ctx):
    results = {}
    with SQLCommand(providerName="PostgresDBProvider",
                    connectionString="dbname=gnm user=gnm") as sql_command:
        fields = ["payload_id", "activity", "task_monitor_id", "status_name", "workernode.name", "job.run"]
        sql_results = sql_command.execute("select %s from payload_result "
                                          "inner join workernode on payload_result.workernode_id=workernode.id "
                                          "inner join payload on payload_result.payload_id=payload.id "
                                          "inner join job on job.id=payload.job_id "
                                          "where payload_id!=%%s and"
                                          "(activity=%%s or activity=%%s or activity=%%s or activity=%%s) and "
                                          "(status_name=%%s or status_name=%%s or status_name=%%s or status_name=%%s)" % ",".join(fields),
                                      ["", "reprocessing", "production", "analysis", "analysis-crab3", "SUCCEEDED", "FAILED", "DONE", "ABORTED"])
        for sql_result in sql_results:
            result = dict(zip(fields, sql_result))
            current_result = results.setdefault(result["activity"], {})
            current_result.setdefault(result["status_name"], {}).setdefault(result["task_monitor_id"], []).append(
                os.path.join("/home/fq8360/data/gnm/payloads", os.path.join(
                    os.path.join(result["workernode.name"], result["job.run"]),
                    "%s-process.csv" % result["payload_id"])))
    if ctx.obj.get("save", False):
        output_results(
            ctx=ctx,
            results=results,
            version=determine_version(os.path.dirname(assess_workflows.__file__)),
            source="%s (%s)" % (__file__, "index_data_by_activity")
        )


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def index_process_names(ctx, paths, pcount):
    filenames = []
    result_set = set()
    for path in paths:
        filenames.extend(_relevant_files_for_context(ctx, path))
    if pcount > 1:
        result_list = do_multicore(
            count=pcount,
            target=_process_names,
            data=filenames
        )
        for result in result_list:
            result_set.union(result)
    else:
        for filename in filenames:
            result_set.union(_process_names(filename))

    output_results(
        ctx=ctx,
        results={"process_names": [name for name in result_set]},
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_process_names")
    )


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def index_tree_statistics(ctx, paths, pcount):
    filenames = []
    results = MulticoreResult()
    for path in paths:
        filenames.extend(_relevant_files_for_context(ctx, path))
    if pcount > 1:
        result_list = do_multicore(
            count=pcount,
            target=_tree_statistics,
            data=filenames
        )
        for result in result_list:
            results += result
    else:
        for filename in filenames:
            results += _tree_statistics(filename)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_tree_statistics")
    )


@click.command()
@click.option("--range_width", "range_width", default=.1)
@click.option("--maximum_chain", "maximum_chain", default=100)
@click.pass_context
def squish_index_into_ranges(ctx, range_width, maximum_chain):
    results = {}
    if ctx.obj.get("use_input", False):
        def probability_function(one, two):
            return abs(one-two)/(max(one, two)*range_width)

        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            input_data = json.load(input_file)["data"]
            # number of nodes, events, tmes
            keys = [int(key) for key in input_data.keys()]
            keys.sort()
            # indizes of numbers overlapping within probability_function
            probabilities = [index for index, probability in enumerate(
                probability_function(one, two) for one, two in zip(
                    keys, keys[1:])) if probability <= 1]
            # pick ranges from sequences of overlapping numbers
            while probabilities:
                index = probabilities[0]
                # check all indexes that belong to a chain
                key_chain = [keys[index]]
                merged = {keys[index]: input_data.pop(str(keys[index]))}
                while index in probabilities:
                    probabilities.pop(0)
                    index += 1
                    try:
                        merged[keys[index]] = input_data.pop(str(keys[index]))
                        key_chain.append(keys[index])
                    except KeyError:  # TODO: is this still needed?
                        continue
                # split chain if it is too long
                key_collection = [key_chain]
                while key_collection:
                    current_key = key_collection.pop()
                    if len(current_key) > maximum_chain:
                        if len(current_key) % maximum_chain == 0:
                            split_point = maximum_chain
                        else:
                            split_point = int(len(current_key) / (len(current_key) // maximum_chain + 1))
                        while current_key:
                            key_collection.append(current_key[:split_point])
                            current_key = current_key[split_point:]
                    else:
                        mean = sum(current_key) / len(current_key)
                        merged_files = [merged.pop(key) for key in current_key]
                        input_data[str(mean)] = [value for values in merged_files for value in values]
            results = input_data

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "squish_index_into_ranges")
    )


@click.command()
@click.option("--include_key", "include_key", default="lambda key, value: False")
@click.pass_context
def subset_data(ctx, include_key):
    """
    Example Usage of this method:
        include_key="'lambda key, value: int(key) > 25'"

    :param ctx:
    :param include_key:
    :return:
    """
    results = {}
    if ctx.obj.get("use_input", False):
        include_key = eval(include_key)
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)
            for key, value in input_data.items():
                if isinstance(value, dict):
                    for inner_key, inner_value in value.items():
                        if isinstance(inner_value, dict):
                            for inner_inner_key, inner_inner_value in inner_value.items():
                                if include_key(inner_inner_key, inner_inner_value):
                                    results.setdefault(key, {}).setdefault(inner_key, {})[inner_inner_key] = inner_inner_value
                        else:
                            if include_key(inner_key, inner_value):
                                results.setdefault(key, {})[inner_key] = inner_value
                else:
                    if include_key(key, value):
                        results[key] = value

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "subset_data")
    )


@click.command()
@click.option("--seed", "seed", type=int)
@click.option("--repeat", "repeat", type=int, default=1)
@click.option("--count", "count", type=int, required=True)
@click.option("--minimum", "minimum", type=int, default=0)
@click.option("--skip_key", "skip_key", default=False)
@click.pass_context
def pick_samples(ctx, seed, repeat, count, minimum, skip_key):
    results = {}

    if seed is not None:
        random.seed(seed)
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)
            working_data = input_data
            try:
                if skip_key:
                    working_data = {value for values in input_data.values() for value in values}
                for key, values in working_data.items():
                    try:
                        for _ in range(repeat):
                            results.setdefault(key, []).append(random.sample(values, count))
                    except ValueError:
                        if minimum == 0 or minimum <= len(values):
                            results.setdefault(key, []).append(values)
                        else:
                            continue
            except AttributeError:
                key = "samples"
                for _ in range(repeat):
                    results.setdefault(key, []).append(random.sample(working_data, count))

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "pick_samples")
    )


@click.command()
@click.option("--skip_key", "skip_key", default=False, type=bool)
@click.pass_context
def aggregate_samples(ctx, skip_key):
    """
    Method aggregates nested dictionaries into a flat one. If it already is a flat dictionary,
    than data is kept and written.

    :param ctx:
    :return:
    """
    results = {}

    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)

            for key, values in input_data.items():
                if isinstance(values, list):
                    results.setdefault(key if not skip_key else "samples", []).extend(values)
                else:
                    try:
                        if len(values[0]) > 1:
                            # flattening data
                            results.setdefault(key, []).append(
                                [element for value in values[0] for element in value])
                        else:
                            # data can be kept
                            results.setdefault(key, []).append(values[0])
                    except KeyError:
                        to_check = [values]
                        while to_check:
                            current_item = to_check.pop(0)
                            try:
                                while current_item:
                                    _, value = current_item.popitem()
                                    to_check.append(value)
                            except KeyError:
                                results.setdefault(key, []).append(current_item)
                            except AttributeError:
                                results.setdefault(key, []).extend(current_item)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "aggregate_samples")
    )


def _line_count(filename):
    return int(subprocess.check_output(["wc", "-l", filename]).strip().split()[0]) - 2


def _payload_line_count(filename):
    result = 0
    try:
        grep = subprocess.Popen(["grep", "(condor_starter)", filename], stdout=subprocess.PIPE)
        wc_result = subprocess.check_output(["wc", "-l"], stdin=grep.stdout)
        grep.wait()
        result = wc_result.strip().split()[0]
    except subprocess.CalledProcessError as e:
        print(e)
    return result


def _data_by_tme(filename):
    results = MulticoreResult()
    tree_builder = CSVTreeBuilder()
    try:
        tree = tree_builder.build(filename)
    except (DataNotInCacheException, TreeInvalidatedException):
        pass
    else:
        if tree is not None:
            node = next(tree.node_iter())
            results.setdefault(node.tme, []).append(filename)
    return results


def _data_by_uid(filename):
    results = MulticoreResult()
    tree_builder = CSVTreeBuilder()
    try:
        tree = tree_builder.build(filename)
    except (DataNotInCacheException, TreeInvalidatedException):
        pass
    else:
        if tree is not None:
            uids = set()
            for node in tree.node_iter():
                if node.uid not in uids:
                    uids.add(node.uid)
                    results.setdefault(node.uid, []).append(filename)
    return results


def _valid_tree(filename):
    tree_builder = CSVTreeBuilder()
    try:
        tree = tree_builder.build(filename)
        if tree:
            return filename
    except (DataNotInCacheException, TreeInvalidatedException):
        pass


def _process_names(filename):
    result = set()
    tree_builder = CSVTreeBuilder()
    try:
        tree = tree_builder.build(filename)
    except DataNotInCacheException:
        tree = None
    except TreeInvalidatedException:
        tree = None
    if tree is not None:
        for node in tree.node_iter():
            try:
                if "(" in node.node[0]:
                    result.add(node.name)
            except IndexError:
                pass
    return result


def _tree_statistics(filename):
    result = MulticoreResult()
    tree_builder = CSVTreeBuilder()
    try:
        tree = tree_builder.build(filename)
    except DataNotInCacheException:
        tree = None
    except TreeInvalidatedException:
        tree = None
    if tree is not None:
        for event in tree.event_iter():
            file_dict = result.setdefault(
                filename, {"process": {}, "traffic": {}, "traffic_count": {}})
            if isinstance(event, ProcessStartEvent) or isinstance(event, ProcessExitEvent):
                file_dict["process"][event.tme] = file_dict["process"].get(event.tme, 0) + 1
            elif isinstance(event, TrafficEvent):
                file_dict["traffic"][event.tme] = file_dict["traffic"].get(event.tme, 0) + 1
                file_dict["traffic_count"][event.tme] = file_dict["traffic_count"].get(
                    event.tme, 0) + (event.in_cnt + event.out_cnt)
    return result


cli.add_command(index_valid_trees)
cli.add_command(index_valid_hdf_trees)
cli.add_command(index_data_by_tme)
cli.add_command(index_data_by_uid)
cli.add_command(index_data_by_number_of_nodes)
cli.add_command(index_data_by_number_of_traffic_events)
cli.add_command(index_data_by_number_of_events)
cli.add_command(index_data_by_number_of_payloads)
cli.add_command(index_data_by_activity)
cli.add_command(index_process_names)
cli.add_command(index_tree_statistics)
cli.add_command(squish_index_into_ranges)
cli.add_command(subset_data)
cli.add_command(pick_samples)
cli.add_command(aggregate_samples)

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.WARNING)
    logging.getLogger("EXCEPTION").setLevel(logging.INFO)
    cli(obj={}, auto_envvar_prefix='DISS')
