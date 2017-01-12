import os
import sys
import json
import click
import importlib

import logging
import assess_workflows
from assess.exceptions.exceptions import TreeInvalidatedException
from assess_workflows.utils.multicoreresult import MulticoreResult
from assess_workflows.utils.statistics import errors, uncorrelated_relative_error
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
    """
    Method returns output file that follows the following format:

    {
        node_count: {
            p_value: {
                "raw": [[diamond levels], ...],
                "identities": [identity_count, ...],
                "diamonds": [diamond_count, ...],
                "files": [file_path, ...]
            }
        }
    }

    :param ctx:
    :param pcount:
    :return:
    """
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
                results += "\n# Metric Analysis\n"
                decorator_data = analysis_data["results"][0]["decorator"]
                for key in decorator_data:
                    if key in ["normalized_matrix", "matrix"]:
                        results += "\n## %s Analysis\n\n" % key
                        matrix_data = decorator_data[key]
                        diagonals = []
                        all_values = []
                        diagonal_issue_counter = 0
                        symmetry_issue_counter = 0
                        metric_issue_counter = 0
                        for row_idx, row_value in enumerate(matrix_data):
                            diagonals.append(row_value[0][row_idx])
                            if row_value[0][row_idx] != 0:
                                diagonal_issue_counter += 1
                            for col_idx in range(row_idx + 1, len(matrix_data)):
                                if row_value[0][col_idx] != matrix_data[col_idx][0][row_idx]:
                                    symmetry_issue_counter += 1
                                if col_idx != row_idx:
                                    all_values.append(
                                        (row_value[0][col_idx], matrix_data[col_idx][0][row_idx],))
                                    if (row_value[0][col_idx] == 0 or
                                            matrix_data[col_idx][0][row_idx] == 0) and col_idx != row_idx:
                                        metric_issue_counter += 1
                        results += "\n### Diagonal in data\n\n"
                        results += "--> Diagonal should be equal to 0 everywhere\n\n"
                        if diagonal_issue_counter == 0:
                            results += "* No issues found with diagonal\n"
                        else:
                            results += "* Identified %s problems in diagonal" % diagonal_issue_counter
                            results += "\n#### Error for Diagonal\n\n"
                            mean, std_error, relative_std_error = errors(diagonals)
                            results += "* absolute error: %s +- %s\n" % (mean, std_error)
                            results += "* relative error: %s +- %s\n" % (mean, relative_std_error)
                        results += "\n### Symmetry in data\n\n"
                        results += "--> Data should be equal when distance is a metric\n\n"
                        if symmetry_issue_counter == 0:
                            results += "* No issues found with symmetry\n"
                        else:
                            results += "* Identified %s problems in symmetry" % symmetry_issue_counter
                            results += "\n#### Error for Symmetry\n\n"
                            mean, relative_std_error = uncorrelated_relative_error(all_values)
                            results += "* uncorrelated relative error: %s +- %s\n" % (mean, relative_std_error)
                        results += "\n### Checking for Metric vs. Pseudo-Metric\n\n"
                        results += "--> For a metric different objects can never have distance 0\n\n"
                        if metric_issue_counter == 0:
                            results += "* No issues found with equality, distance might be a metric\n"
                        else:
                            results += "* Identified %s problems with equality\n" % metric_issue_counter
                            results += "* Distance is no metric, but might be a pseudo-metric\n"
                results += "## Files for further reference\n\n"
                for index, current_file in enumerate(analysis_data["files"]):
                    results += "* [%s]: %s\n" % (index, current_file)

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "analyse_metric"),
        file_type="md"
    )


def _analyse_diamond_perturbation(kwargs):
    """
    {
        p_count: {
            diamond_count: {
                "profile_distortions": [],              # profile distortion based on frequency
                "profile_distortions_signatures": [],   # profile distortion based on set count
                "distance_errors": []                   # distance error based on frequency
                "distance_errors_signatures": []        # distance error based on set count
                "signature_counts": [],                 # nr of signatures in tree
                "node_counts": [],                      # nr of nodes in tree
                "raw": [{
                    "level": diamond_level,
                    "nested": nesting_level,
                    "nodes": node_count,
                    "signatures": signature_count
                }, ...]
            }
        }
    }

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
                node_count = 0
                for node in tree.node_iter():
                    node_count += 1
                    current_signature = signature.get_signature(node, node.parent())
                    node_signatures.add(current_signature[0])
                    diamond = diamonds.setdefault(current_signature[0], {})
                    diamond.setdefault("signatures", set()).add(current_signature[1])
                    diamond.setdefault("nodes", set()).add(node)
                diamonds = {key: diamond for key, diamond in diamonds.items() if
                            len(diamond.get("signatures", set())) > 1}
                diamond_perturbation = {}
                for diamond_key, diamond in diamonds.items():
                    # found a diamond, that represents several diamond nodes
                    result = diamond_perturbation.setdefault(diamond_key, {
                        "nested": 0,
                        "nodes": set(),
                        "signatures": set()
                    })
                    result["level"] = max(0, len(diamond.get("signatures", set())) - 1)
                    for node in diamond.get("nodes"):
                        to_check = set(node.children_list())
                        result["nodes"].add(node)
                        result["signatures"].add(signature.get_signature(node, node.parent)[0])
                        while to_check:
                            child = to_check.pop()
                            result["nodes"].add(child)
                            child_signatures = signature.get_signature(child, child.parent())
                            result["signatures"].add(child_signatures[0])
                            to_check.update(child.children_list())
                            if child_signatures[0] in diamonds:
                                # diamond is a nested diamond, so initialise it here
                                diamond_perturbation[child_signatures[0]] = {
                                    "level": 1,
                                    "nested": result["nested"] + 1,
                                    "nodes": set(),
                                    "signatures": set()
                                }
                diamond_count = len(diamond_perturbation)
                perturbation_result = perturbation_results.setdefault(
                    signature._signatures[0]._height, {}).setdefault(diamond_count, {})
                perturbation_result.setdefault("profile_distortions", []).append(
                    sum([len(diamond.get("nodes", [])) * diamond["level"]
                         for diamond in diamond_perturbation.values()])
                )
                perturbation_result.setdefault("profile_distortions_signatures", []).append(
                    sum([len(diamond.get("signatures", [])) * diamond["level"]
                         for diamond in diamond_perturbation.values()])
                )
                perturbation_result.setdefault("distance_errors", []).append(
                    sum([len(diamond.get("nodes", [])) for diamond in diamond_perturbation.values()])
                )
                perturbation_result.setdefault("distance_errors_signatures", []).append(
                    sum([len(diamond.get("signatures", [])) for diamond in diamond_perturbation.values()])
                )
                perturbation_result.setdefault("signature_counts", []).append(len(node_signatures))
                perturbation_result.setdefault("node_counts", []).append(node_count)
                perturbation_result.setdefault("raw", []).append({key: {
                    "level": value["level"],
                    "nested": value["nested"],
                    "nodes": len(value["nodes"]),
                    "signatures": len(value["signatures"])} for key, value in diamond_perturbation.items()})
    return perturbation_results


def _analyse_diamonds(kwargs):
    """
    Method expects an ensemble signature in configuration were signature at position 0 has length
    n - 1 whereas signature at position 1 has length n (criterium for diamonds). It then builds
    a dictionary for given signatures from position 0 and builds a collection from signatures at
    position 1. The number of signatures that are associated to the different keys is then relevant
    to determine the diamonds. When more than one signature is assigned, then we got a diamond.

    Method creates different fields in output file:

    * raw: contains the levels of the diamonds within a given tree
    * identities: number of identities for the whole tree
    * diamonds: number of diamonds within the tree (independent from level)
    * diamond_nodes: number of nodes that make up the diamonds
    * files: files that were used

    In addition, all of these fields are associated to a given signature_builder. It defines the
    actual height that is analysed. Meaning, the p value that is used to index the output file.

    {
        node_count: {
            p_value: {
                "raw": {
                    "levels": [[diamond level, ...], ...],
                    "nodes": [[diamond nodes, ...], ...]
                }
                "identities": [identity_count, ...],
                "diamonds": [diamond_count, ...],
                "diamond_nodes": [diamond_node_count, ...],
                "node_counts": [node_count, ...],
                "files": [file_path, ...]
            }
        }
    }

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
                current_node_count = 0
                for node in tree.node_iter():
                    current_node_count += 1
                    current_signatures = signature.get_signature(node, node.parent())
                    current_node = node_dict.setdefault(current_signatures[0], {})
                    current_node.setdefault("nodes", set()).add(node)
                    current_node.setdefault("signatures", set()).add(current_signatures[1])
                diamonds = {signature: {
                    "nodes": len(signature_values.get("nodes", set())),
                    "levels": len(signature_values.get("signatures", set())) - 1
                } for signature, signature_values in node_dict.items() if
                            len(signature_values.get("signatures", set())) > 1}
                current_result = result.setdefault(node_count, {}).setdefault(
                    signature._signatures[0]._height, {})
                raw_result = current_result.setdefault("raw", {"levels": [], "nodes": []})
                raw_result["levels"].append([diamond.get("levels", 0) for diamond in diamonds.values()])
                raw_result["nodes"].append([diamond.get("nodes", 0) for diamond in diamonds.values()])
                current_result.setdefault("node_counts", []).append(current_node_count)
                current_result.setdefault("identities", []).append(len(node_dict))
                current_result.setdefault("diamonds", []).append(len(diamonds))
                current_result.setdefault("diamond_nodes", []).append(
                    sum([diamond.get("nodes", 0) for diamond in diamonds.values()]))
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
