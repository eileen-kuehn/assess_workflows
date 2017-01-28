import os
import re
import csv
import json

import assess_workflows
import click
import logging

from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import output_results, determine_version
from utility.exceptions import ExceptionFrame
from utility.report import LVL


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
    ctx.obj["json"] = False
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input


@click.command()
@click.option("--symmetric", "symmetric", default=True)
@click.pass_context
def transform_matrix_to_adjacency_list(ctx, symmetric):
    if ctx.obj.get("use_input", False):
        ctx.obj["json"] = True
        results = {}
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)
            files = input_data["files"]
            for result_idx, result in enumerate(input_data["results"][0]):
                decorator = result["decorator"]["normalized_matrix"]
                for row_idx, row in enumerate(decorator):
                    for col_idx, col in enumerate(row[0]):
                        if col_idx == row_idx:
                            continue
                        results.setdefault(files[result_idx][row_idx], {})[files[result_idx][col_idx]] = col
                        if symmetric:
                            results.setdefault(files[result_idx][col_idx], {})[files[result_idx][row_idx]] = col
        output_results(
            ctx=ctx,
            results=results,
            version=determine_version(os.path.dirname(assess_workflows.__file__)),
            source="%s (%s)" % (__file__, "transform_matrix_to_adjacency_list")
        )


@click.command()
@click.pass_context
def transform_matrix_to_csv(ctx):
    if ctx.obj.get("use_input", False):
        results = ""
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)
            files = input_data["files"]
            for result_idx, result in enumerate(input_data["results"][0]):
                decorator = result["decorator"]["normalized_matrix"]
                maximum_index = len(decorator)
                results += ",".join(files[result_idx])
                results += "\n"
                for row_index in xrange(0, maximum_index):
                    row = [0 for _ in xrange(row_index+1)]
                    for col_index in xrange(row_index+1, maximum_index):
                        row.append(decorator[col_index][0][row_index])
                    results += "%s\n" % ",".join([str(item) for item in row])
        output_results(
            ctx=ctx,
            results=results,
            file_type="csv",
            version=determine_version(os.path.dirname(assess_workflows.__file__)),
            source="%s (%s)" % (__file__, "transform_matrix_to_csv")
        )


@click.command()
@click.pass_context
def transform_matrix_to_sql(ctx):
    if ctx.obj.get("use_input", False):
        result = "INSERT INTO object_distances (a, b, d) VALUES\n"
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)

            # contains list of lists
            # each list is a row within the matrix, inside the matrix, everything above the diagonal
            # including diagonal is 0, so results can be skipped here
            data = input_data["results"][0]["decorator"]["normalized_matrix"]
            for row_index in xrange(0, len(data)):
                for column_index in xrange(0, row_index):
                    if len(result) > 46:  # length of INSERT INTO...
                        result += ",\n"
                    result += "(%d,%d,%s)" % (row_index, column_index, data[row_index][column_index])
        output_results(
            ctx=ctx,
            results=result + ";",
            version=determine_version(os.path.dirname(assess_workflows.__file__)),
            source="%s (%s)" % (__file__, "transform_matrix_to_sql")
        )


@click.command()
@click.option("--paths", "paths", type=click.Path(), multiple=True,
              help="Path of mapping results to consider for SQL generation.")
@click.pass_context
def transform_mapping_to_sql(ctx, paths):
    if ctx.obj.get("use_input", False) or len(paths) > 0:
        result = ""
        update_cmd = "update payload_result set payload_id='%s' where id=%s;\n"

        for path in paths:
            print("starting with %s" % path)
            with open(path, "r") as input_file:
                reader = csv.reader(input_file, quotechar="'")
                for row in reader:
                    result += update_cmd % (re.match("([\d-]+)", row[1]).group(), re.match("(\d+)", row[0]).group())

        output_results(
            ctx=ctx,
            results=result + ";",
            version=determine_version(os.path.dirname(assess_workflows.__file__)),
            source="%s (%s)" % (__file__, "transform_mapping_to_sql")
        )

cli.add_command(transform_matrix_to_adjacency_list)
cli.add_command(transform_matrix_to_sql)
cli.add_command(transform_matrix_to_csv)
cli.add_command(transform_mapping_to_sql)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
