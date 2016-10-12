import os
import glob
import json
import click
import random
import logging
import subprocess

from assess.exceptions.exceptions import TreeInvalidatedException
from gnmutils.exceptions import DataNotInCacheException
from utility.exceptions import ExceptionFrame
from utility.report import LVL

from assess.generators.gnm_importer import CSVTreeBuilder
from assess.events.events import ProcessStartEvent, ProcessExitEvent, TrafficEvent

import assess_workflows
from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import output_results, determine_version


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
    ctx.obj["json"] = True
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.pass_context
def index_data_by_number_of_nodes(ctx, paths):
    results = {}
    for path in paths:
        for filename in glob.glob("%s/*/*-process.csv" % path):
            count = _line_count(filename)
            results.setdefault(count, []).append(filename)

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
            results.setdefault(count, []).append(filename)

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
    for path in paths:
        for filename in glob.glob("%s/*/*-process.csv" % path):
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
def index_tree_statistics(ctx, paths):
    tree_builder = CSVTreeBuilder()
    results = {}
    for path in paths:
        for filename in glob.glob("%s/*/*-process.csv" % path):
            try:
                tree = tree_builder.build(filename)
            except DataNotInCacheException:
                tree = None
            except TreeInvalidatedException:
                tree = None
            if tree:
                for event in tree.event_iter():
                    file_dict = results.setdefault(
                        filename,
                        {"process": {}, "traffic": {}, "traffic_count": {}})
                    if isinstance(event, ProcessStartEvent) or isinstance(event, ProcessExitEvent):
                        file_dict["process"][event.tme] = file_dict["process"].get(event.tme, 0) + 1
                    elif isinstance(event, TrafficEvent):
                        file_dict["traffic"][event.tme] = file_dict["traffic"].get(event.tme, 0) + 1
                        file_dict["traffic_count"][event.tme] = file_dict["traffic_count"].get(
                            event.tme, 0) + (event.in_cnt + event.out_cnt)
            else:
                print("skipping %s as it is invalid" % filename)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "index_tree_statistics")
    )


@click.command()
@click.option("--range_width", "range_width", default=.1)
@click.pass_context
def squish_index_into_ranges(ctx, range_width):
    results = {}
    if ctx.obj.get("use_input", False):
        def probability_function(one, two):
            return abs(one-two)/(max(one, two)*range_width)

        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)
            keys = [int(key) for key in input_data.keys()]
            keys.sort()
            probabilities = {index: probability for index, probability in
                             enumerate(probability_function(one, two) for one, two in
                                       zip(keys, keys[1:])) if probability <= 1}
            sorted_indexes = [key for key in probabilities]
            sorted_indexes.sort()

            for index in sorted_indexes:
                try:
                    # check all indexes that belong to a chain
                    key_collection = [keys[index]]
                    merged = input_data.pop(str(keys[index]), [])
                    while index in sorted_indexes:
                        index += 1
                        key_collection.append(keys[index])
                        merged.extend(input_data.pop(str(keys[index])))
                    # FIXME: a qualitative splitting should be done for not getting too wide ranges
                    mean = sum(key_collection)/len(key_collection)
                    input_data[str(mean)] = merged
                except KeyError:
                    continue  # current index has already been removed
            results = input_data

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "squish_index_into_ranges")
    )


@click.command()
@click.option("--seed", "seed", type=int)
@click.option("--repeat", "repeat", type=int, default=1)
@click.option("--count", "count", type=int, required=True)
@click.pass_context
def pick_samples(ctx, seed, repeat, count):
    results = {}

    if seed is not None:
        random.seed(seed)
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)
            for key, values in input_data.items():
                try:
                    for _ in xrange(repeat):
                        results.setdefault(key, []).append(random.sample(values, count))
                except ValueError:
                    continue

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "pick_samples")
    )


def _line_count(filename):
    return int(subprocess.check_output(["wc", "-l", filename]).strip().split()[0]) - 2


cli.add_command(index_data_by_number_of_nodes)
cli.add_command(index_data_by_number_of_traffic_events)
cli.add_command(index_data_by_number_of_events)
cli.add_command(index_tree_statistics)
cli.add_command(squish_index_into_ranges)
cli.add_command(pick_samples)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
