import os
import sys
import json
import click
import importlib

import logging
import assess_workflows
from assess.exceptions.exceptions import TreeInvalidatedException
from gnmutils.exceptions import DataNotInCacheException

from utility.exceptions import ExceptionFrame
from utility.report import LVL

from assess.generators.gnm_importer import CSVTreeBuilder

from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import output_results, determine_version


@click.group()
@click.option("--basepath", "basepath", multiple=False, required=True)
@click.option("--workflow-name", "workflow_name", multiple=False, required=True)
@click.option("--step", "step", default=1, multiple=False)
@click.option("--save", "save", is_flag=True)
@click.option("--use_input", "use_input", is_flag=True,
              help="Use input file specified for current task.")
@click.option("--configuration", "configuration", multiple=False,
              help="Location of configuration file")
@click.pass_context
def cli(ctx, basepath, workflow_name, step, save, use_input, configuration):
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name, step=step)
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input
    config_module = ctx.obj["structure"].configuration_module()
    importlib.import_module(config_module)
    configdict = sys.modules[config_module]
    ctx.obj["configurations"] = configdict.configurations


@click.command()
@click.pass_context
def analyse_diamonds(ctx):
    results = {}
    ctx.obj["json"] = True
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        tree_builder = CSVTreeBuilder()
        signature_builders = ctx.obj.get("configurations", [{}])[0].get("signatures", [])

        with open(file_path, "r") as input_file:
            analysis_files = json.load(input_file).get("data", None)
            for node_count, tree_path in analysis_files.items():
                try:
                    tree = tree_builder.build(tree_path[0][0])
                except DataNotInCacheException:
                    tree = None
                except TreeInvalidatedException:
                    tree = None
                if tree is not None:
                    for signature_builder in signature_builders:
                        signature = signature_builder()
                        node_dict = {}
                        for node in tree.node_iter():
                            current_signatures = signature.get_signature(node, node.parent())
                            node_dict.setdefault(current_signatures[0], set()).add(current_signatures[1])
                        diamond_values = [len(signatures) - 1 for signatures in node_dict.values() if len(signatures) > 1]
                        current_result = results.setdefault(node_count, {}).setdefault(
                            signature._signatures[0]._height, {})
                        current_result.setdefault("raw", []).append(diamond_values)
                        current_result.setdefault("identities", []).append(len(node_dict))
                        current_result.setdefault("diamonds", []).append(len(diamond_values))
                        current_result.setdefault("files", []).append(tree_path)
    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "analyse_diamonds")
    )


@click.command()
@click.pass_context
def analyse_metric(ctx):
    results = ""
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            analysis_data = json.load(input_file).get("data", None)

            if analysis_data:
                results += "\n# Metric Analysis\n\n"
                results += "## Files for further reference\n\n"
                for index, file in enumerate(analysis_data["files"]):
                    results += "* [%s]: %s\n" % (index, file)
                results += "\n## Diagonal in data\n\n"
                results += "--> Diagonal should be equal to 0 everywhere\n\n"
                decorator_data = analysis_data["results"][0]["decorator"]
                for key in decorator_data:
                    if key in ["normalized_matrix", "matrix"]:
                        results += "Identified %s, checking now\n\n" % key
                        matrix_data = decorator_data[key]
                        valid = True
                        for row_idx, row_value in enumerate(matrix_data):
                            if row_value[0][row_idx] != 0:
                                valid = False
                                results += "* Identified Problem in Row/Col %s: %s\n" % (row_idx, row_value[0][row_idx])
                        if valid:
                            results += "* No issues found with diagonal\n"
                results += "\n## Symmetry in data\n\n"
                results += "--> Data should be equal when distance is a metric\n\n"
                for key in decorator_data:
                    if key in ["normalized_matrix", "matrix"]:
                        results += "Identified %s, checking now\n\n" % key
                        matrix_data = decorator_data[key]
                        valid = True
                        for row_idx, row_value in enumerate(matrix_data):
                            for col_idx in range(row_idx + 1, len(matrix_data)):
                                if row_value[0][col_idx] != matrix_data[col_idx][0][row_idx]:
                                    valid = False
                                    results += "* Identified Problem for (%s:%s - %s:%s): %s != %s\n" % (row_idx, col_idx, col_idx, row_idx, row_value[0][col_idx], matrix_data[col_idx][0][row_idx])
                        if valid:
                            results += "* No issues found with symmetry\n"
                results += "\n## Checking for Metric vs. Pseudo-Metric\n\n"
                results += "--> For a metric different objects can never have distance 0\n\n"
                for key in decorator_data:
                    if key in ["normalized_matrix", "matrix"]:
                        results += "Identified %s, checking now\n\n" % key
                        matrix_data = decorator_data[key]
                        valid = True
                        for row_idx, row_value in enumerate(matrix_data):
                            for col_idx, value in enumerate(row_value):
                                if row_idx == col_idx:
                                    continue
                                if value == 0:
                                    valid = False
                                    results += "* Identified distance of 0 for (%s:%s): %s == 0\n" % (row_idx, col_idx, value)
                        if valid:
                            results += "* Distance might be a metric\n"
                        else:
                            results += "\n* Distance is no metric, but might be a pseudo-metric\n"

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "analyse_metric"),
        file_type="md"
    )


cli.add_command(analyse_metric)
cli.add_command(analyse_diamonds)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
