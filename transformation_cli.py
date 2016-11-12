import os
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
@click.pass_context
def transform_matrix_to_adjacency_list(ctx):
    if ctx.obj.get("use_input", False):
        ctx.obj["json"] = True
        result = {}
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)
            files = input_data["files"]
            data = input_data["results"][0]["decorator"]["normalized_matrix"]
            for row_idx, row in enumerate(data):
                result[files[row_idx]] = {}
                for col_idx, col in enumerate(row[0]):
                    if col_idx == row_idx:
                        continue
                    result[files[row_idx]][files[col_idx]] = col
        output_results(
            ctx=ctx,
            results=result,
            version=determine_version(os.path.dirname(assess_workflows.__file__)),
            source="%s (%s)" % (__file__, "transform_matrix_to_adjacency_list")
        )


@click.command()
@click.pass_context
def transform_matrix_to_csv(ctx):
    if ctx.obj.get("use_input", False):
        result = ""
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)

            data = input_data["results"][0]["decorator"]["normalized_matrix"]
            maximum_index = len(data)
            result += "%s\n" % ",".join([str(number) for number in range(1, maximum_index+1)])
            for row_index in xrange(0, maximum_index):
                # write first part of matrix: 0 - index
                result += "%s" % ",".join([str(element) for element in data[row_index][0:row_index]])
                # write second part of matrix: 1 - len(data)
                for col_index in xrange(row_index, maximum_index):
                    result += ",%s" % str(data[col_index][row_index])
                result += "\n"
        output_results(
            ctx=ctx,
            results=result,
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

cli.add_command(transform_matrix_to_adjacency_list)
cli.add_command(transform_matrix_to_sql)
cli.add_command(transform_matrix_to_csv)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
