import os
import json
import click

import logging
import assess_workflows

from utility.exceptions import ExceptionFrame
from utility.report import LVL

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
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input


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

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
