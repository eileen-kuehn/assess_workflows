import os
import sys
import json
import click
import importlib

import logging
import assess_workflows
from assess.exceptions.exceptions import TreeInvalidatedException
from assess_workflows.utils.multicoreresult import MulticoreResult
from gnmutils.exceptions import DataNotInCacheException

from utility.exceptions import ExceptionFrame
from utility.report import LVL

from assess.generators.gnm_importer import CSVTreeBuilder

from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import output_results, determine_version, do_multicore


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


def _analyse_compression(kwargs):
    filepath = kwargs.get("filepath", None)
    node_count = kwargs.get("node_count", None)
    signature_builders = kwargs.get("signature_builders", None)
    result = MulticoreResult()
    tree_builder = CSVTreeBuilder()
    try:
        tree = tree_builder.build(filepath)
    except (DataNotInCacheException, TreeInvalidatedException):
        pass
    else:
        if tree is not None:
            for signature_builder in signature_builders:
                signature = signature_builder()
                compression = [set() for _ in range(signature.count)]
                for node in tree.node_iter():
                    identities = signature.get_signature(node, node.parent())
                    for index, identity in enumerate(identities):
                        compression[index].add(identity)
                # write results
                # {node_count: {signature_1: value, signature_2: value}}
                current = result.setdefault(node_count, {})
                for index, single_signature in enumerate(signature._signatures):
                    current.setdefault(repr(single_signature), []).append(len(compression[index]))
    return result


@click.command()
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def analyse_compression(ctx, pcount):
    results = MulticoreResult()
    ctx.obj["json"] = True
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        signature_builders = ctx.obj.get("configurations", [{}])[0].get("signatures", [])

        with open(file_path, "r") as input_file:
            analysis_files = json.load(input_file).get("data", None)
            if pcount > 1:
                data = [{
                    "node_count": node_count,
                    "filepath": tree_path[0],
                    "signature_builders": signature_builders} for node_count, tree_paths in
                        analysis_files.items() for tree_path in tree_paths]
                multicore_results = do_multicore(
                    count=pcount,
                    target=_analyse_compression,
                    data=data
                )
                for result in multicore_results:
                    results += result
            else:
                for node_count, tree_paths in analysis_files.items():
                    for tree_path in tree_paths:
                        results += _analyse_compression({
                            "node_count": node_count,
                            "filepath": tree_path[0],
                            "signature_builders": signature_builders
                        })

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "analyse_compression")
    )


@click.command()
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def analyse_diamond_perturbations(ctx, pcount):
    results = MulticoreResult()
    ctx.obj["json"] = True
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        signature_builders = ctx.obj.get("configurations", [{}])[0].get("signatures", [])

        with open(file_path, "r") as input_file:
            analysis_files = json.load(input_file).get("data", None)
            if pcount > 1:
                # combine data
                data = [{
                            "filepath": path[0],
                            "signature_builders": signature_builders
                        } for paths in analysis_files.values() for path in paths]
                multicore_results = do_multicore(
                    count=pcount,
                    target=_analyse_diamond_perturbation,
                    data=data)
                for result in multicore_results:
                    results += result
            else:
                for tree_paths in analysis_files.values():
                    for tree_path in tree_paths:
                        results += _analyse_diamond_perturbation({
                            "filepath": tree_path[0], "signature_builders": signature_builders})

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "analyse_diamond_perturbation")
    )


@click.command()
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def analyse_diamonds(ctx, pcount):
    results = MulticoreResult()
    ctx.obj["json"] = True
    if ctx.obj.get("use_input", False):
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()
        signature_builders = ctx.obj.get("configurations", [{}])[0].get("signatures", [])

        with open(file_path, "r") as input_file:
            analysis_files = json.load(input_file).get("data", None)
            if pcount > 1:
                data = [{
                    "node_count": node_count,
                    "filepath": tree_path[0],
                    "signature_builders": signature_builders} for node_count, tree_paths in
                        analysis_files.items() for tree_path in tree_paths]
                multicore_results = do_multicore(
                    count=pcount,
                    target=_analyse_diamonds,
                    data=data
                )
                for result in multicore_results:
                    results += result
            else:
                for node_count, tree_paths in analysis_files.items():
                    for tree_path in tree_paths:
                        results += _analyse_diamonds({
                            "node_count": node_count,
                            "filepath": tree_path[0],
                            "signature_builders": signature_builders})

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


def _analyse_diamond_perturbation(kwargs):
    """
    :param kwargs: dict with keys filepath and signature_builders
    :return:
    """
    filepath = kwargs.get("filepath", None)
    signature_builders = kwargs.get("signature_builders", None)
    tree_builder = CSVTreeBuilder()
    perturbation_results = MulticoreResult()
    try:
        tree = tree_builder.build(filepath)
    except (DataNotInCacheException, TreeInvalidatedException):
        pass
    else:
        if tree is not None:
            for signature_builder in signature_builders:
                diamonds = {}
                node_signatures = set()
                signature = signature_builder()
                for node in tree.node_iter():
                    current_signature = signature.get_signature(node, node.parent())
                    node_signatures.add(current_signature[0])
                    diamond = diamonds.setdefault(current_signature[0], {})
                    diamond.setdefault("signatures", set()).add(current_signature[1])
                    diamond.setdefault("nodes", set()).add(node)
                diamonds = {key: diamond for key, diamond in diamonds.items() if
                            len(diamond.get("signatures", set())) > 1}
                diamond_perturbation = {}
                for diamond_key, diamond in diamonds.items():
                    # found a diamond
                    result = diamond_perturbation.setdefault(diamond_key, {"factor": 1, "nested": 0, "nodes": set()})
                    result["factor"] *= len(diamond.get("signatures", set()))
                    for node in diamond.get("nodes"):
                        to_check = set(node.children_list())
                        result["nodes"].add(node)
                        while to_check:
                            child = to_check.pop()
                            child_signatures = signature.get_signature(child, child.parent())
                            if child_signatures[0] not in diamonds:
                                # child is only node, not diamond, so also take care on its children
                                result["nodes"].add(child)
                                to_check.update(child.children_list())
                            else:
                                # take care that the diamond is initialised as nested diamond
                                diamond_perturbation[child_signatures[0]] = {
                                    "factor": result["factor"],
                                    "nested": result["nested"] + 1,
                                    "nodes": set()
                                }
                diamond_count = len(diamond_perturbation)

                perturbations = [(diamond.get("factor", 2) - 1) * len(diamond.get(
                    "nodes", [])) for diamond in diamond_perturbation.values()]
                perturbation_result = perturbation_results.setdefault(
                    signature._signatures[0]._height, {}).setdefault(diamond_count, {})
                perturbation_result.setdefault("perturbations", []).append(sum(perturbations))
                perturbation_result.setdefault("node_counts", []).append(len(node_signatures))
                perturbation_result.setdefault("raw", []).append({key: {
                    "factor": value["factor"],
                    "nested": value["nested"],
                    "nodes": len(value["nodes"])} for key, value in diamond_perturbation.items()})
    return perturbation_results


def _analyse_diamonds(kwargs):
    """
    :param kwargs: dict containing keys node_count, filepath and signature_builders
    :return:
    """
    node_count = kwargs.get("node_count", None)
    filepath = kwargs.get("filepath", None)
    signature_builders = kwargs.get("signature_builders", None)
    result = MulticoreResult()
    tree_builder = CSVTreeBuilder()
    try:
        tree = tree_builder.build(filepath)
    except (DataNotInCacheException, TreeInvalidatedException):
        pass
    else:
        if tree is not None:
            for signature_builder in signature_builders:
                signature = signature_builder()
                node_dict = {}
                for node in tree.node_iter():
                    current_signatures = signature.get_signature(node, node.parent())
                    node_dict.setdefault(current_signatures[0], set()).add(current_signatures[1])
                diamond_values = [len(signatures) - 1 for signatures in node_dict.values() if
                                  len(signatures) > 1]
                current_result = result.setdefault(node_count, {}).setdefault(
                    signature._signatures[0]._height, {})
                current_result.setdefault("raw", []).append(diamond_values)
                current_result.setdefault("identities", []).append(len(node_dict))
                current_result.setdefault("diamonds", []).append(len(diamond_values))
                current_result.setdefault("files", []).append(filepath)
    return result

cli.add_command(analyse_compression)
cli.add_command(analyse_metric)
cli.add_command(analyse_diamonds)
cli.add_command(analyse_diamond_perturbations)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
