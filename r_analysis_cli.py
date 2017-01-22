import os
import json
import logging

import click
import rpy2
import math

from assess_workflows.generic.structure import Structure
from assess_workflows.utils.statistics import uncorrelated_relative_error, \
    uncorrelated_relative_distance_deviation, uncorrelated_relative_max_distance_deviation, \
    uncorrelated_relative_deviation_and_standard_error
from assess_workflows.utils.utils import output_r_data
from utility.exceptions import ExceptionFrame
from utility.report import LVL

from assess.algorithms.statistics.setstatistics import SetStatistics
from assess.algorithms.statistics.splittedstatistics import SplittedStatistics


@click.group()
@click.option("--basepath", "basepath", multiple=False, required=True)
@click.option("--workflow-name", "workflow_name", multiple=False, required=True)
@click.option("--step", "step", default=1, multiple=False)
@click.option("--configuration", "configuration", multiple=False,
              help="Location of configuration file")
@click.option("--save", "save", is_flag=True)
@click.option("--use_input", "use_input", is_flag=True,
              help="Use input file specified for current task.")
@click.pass_context
def cli(ctx, basepath, workflow_name, step, configuration, save, use_input):
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name, step=step)
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input


@click.command()
@click.pass_context
def analyse_compression(ctx):
    structure = ctx.obj.get("structure")

    if ctx.obj.get("use_input", False):
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            from rpy2.robjects.packages import importr
            from rpy2.robjects.lib.dplyr import DataFrame
            import rpy2.robjects.lib.ggplot2 as ggplot2

            grdevices = importr("grDevices")
            datatable = importr("data.table")
            base = importr("base")

            input_data = json.load(input_file).get("data", None)
            result_dt = None
            for node_count, signature_values in input_data.items():
                for signature, compression_values in signature_values.items():
                    current_dt = datatable.data_table(
                        signature=signature,
                        node_count=base.as_integer(node_count),
                        compression=base.unlist(compression_values)
                    )
                    if result_dt is None:
                        result_dt = current_dt
                    else:
                        result_dt = datatable.rbindlist([result_dt, current_dt])
            # summarize data table
            summarized_values = (DataFrame(result_dt)
                                 .group_by("signature", "node_count")
                                 .summarize(compression_mean="mean(compression)",
                                            relative_compression_mean="mean(compression/node_count)",
                                            compression_stderror="sd(compression)/sqrt(length(compression))",
                                            relative_compression_stderror="sd(compression/node_count)/sqrt(length(compression))"))
            absolute_plot = ggplot2.ggplot(summarized_values) + ggplot2.aes_string(
                x="node_count", y="compression_mean", color="signature") + \
                ggplot2.geom_point() + ggplot2.geom_errorbar(width=.01) + ggplot2.aes_string(
                ymin="compression_mean-compression_stderror",
                ymax="compression_mean+compression_stderror")
            relative_plot = ggplot2.ggplot(summarized_values) + ggplot2.aes_string(
                x="node_count", y="relative_compression_mean", color="signature") + \
                ggplot2.geom_point() + ggplot2.geom_errorbar(width=.01) + ggplot2.aes_string(
                ymin="relative_compression_mean-relative_compression_stderror",
                ymax="relative_compression_mean+relative_compression_stderror")
            absolute_filename = os.path.join(structure.exploratory_path(), "absolute_compression.png")
            relative_filename = os.path.join(structure.exploratory_path(), "relative_compression.png")
            for file_name, plot in {absolute_filename: absolute_plot,
                                    relative_filename: relative_plot}.items():
                grdevices.png(file_name)
                plot.plot()
                grdevices.dev_off()

            if ctx.obj.get("save", False):
                # save model data for further adaptations
                rdata_filename = structure.intermediate_file_path(file_type="RData")
                output_r_data(
                    ctx=ctx, filename=rdata_filename, absolute_plot=absolute_plot,
                    relative_plot=relative_plot, absolute_filename=absolute_filename,
                    relative_filename=relative_filename, summarized_values=summarized_values,
                    result_dt=result_dt
                )


@click.command()
@click.pass_context
def analyse_diamond_perturbations(ctx):
    def _diamond_perturbation_plot(data, mean, stderror):
        import rpy2.robjects.lib.ggplot2 as ggplot2
        return ggplot2.ggplot(data) + ggplot2.aes_string(
            x="diamond_count", y="%s" % mean, color="p_value") + ggplot2.geom_point() + \
            ggplot2.geom_errorbar(width=.01) + \
            ggplot2.aes_string(ymin="%s-%s" % (mean, stderror), ymax="%s+%s" % (mean, stderror))
    structure = ctx.obj.get("structure")

    if ctx.obj.get("use_input", False):
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            from rpy2.robjects.packages import importr
            import rpy2.robjects.lib.ggplot2 as ggplot2
            from rpy2.robjects.lib.dplyr import DataFrame

            grdevices = importr("grDevices")
            datatable = importr("data.table")
            base = importr("base")

            input_data = json.load(input_file).get("data", None)
            result_dt = None
            for p_count, diamond_values in input_data.items():
                for diamond_count, diamond_samples in diamond_values.items():
                    current_dt = datatable.data_table(
                        p_value=p_count,
                        diamond_count=base.as_integer(diamond_count),
                        profile_distortion=base.unlist(diamond_samples.get("profile_distortions")),
                        profile_distortion_signatures=base.unlist(diamond_samples.get("profile_distortions_signatures")),
                        distance_error=base.unlist(diamond_samples.get("distance_errors")),
                        distance_error_signatures=base.unlist(diamond_samples.get("distance_errors_signatures")),
                        node_count=base.unlist(diamond_samples.get("node_counts")),
                        signature_count=base.unlist(diamond_samples.get("signature_counts"))
                    )
                    if result_dt is None:
                        result_dt = current_dt
                    else:
                        result_dt = datatable.rbindlist([result_dt, current_dt])
            # summarize data table
            summarized_values = (DataFrame(result_dt)
                                 .group_by("p_value", "diamond_count")
                                 .summarize(profile_distortion_mean="mean(profile_distortion)",
                                            profile_distortion_signatures_mean="mean(profile_distortion_signatures)",
                                            distance_error_mean="mean(distance_error)",
                                            distance_error_signatures_mean="mean(distance_error_signatures)",
                                            relative_profile_distortion_mean="mean(profile_distortion/node_count)",
                                            relative_profile_distortion_signatures_mean="mean(profile_distortion_signatures/signature_count)",
                                            relative_distance_error_mean="mean(distance_error/node_count)",
                                            relative_distance_error_signatures_mean="mean(distance_error_signatures/signature_count)",
                                            profile_distortion_stderror="sd(profile_distortion)/sqrt(length(profile_distortion))",
                                            profile_distortion_signatures_stderror="sd(profile_distortion_signatures)/sqrt(length(profile_distortion_signatures))",
                                            distance_error_stderror="sd(distance_error)/sqrt(length(distance_error))",
                                            distance_error_signatures_stderror="sd(distance_error_signatures)/sqrt(length(distance_error_signatures))",
                                            relative_profile_distortion_stderror="sd(profile_distortion/node_count)/sqrt(length(profile_distortion))",
                                            relative_profile_distortion_signatures_stderror="sd(profile_distortion_signatures/signature_count)/sqrt(length(profile_distortion_signatures))",
                                            relative_distance_error_stderror="sd(distance_error/node_count)/sqrt(length(distance_error))",
                                            relative_distance_error_signatures_stderror="sd(distance_error_signatures/node_count)/sqrt(length(distance_error_signatures))"))
            absolute_profile_distortion_plot = _diamond_perturbation_plot(
                summarized_values, "profile_distortion_mean", "profile_distortion_stderror")
            absolute_profile_distortion_filename = os.path.join(structure.exploratory_path(), "profile_distortion.png")
            relative_profile_distortion_plot = _diamond_perturbation_plot(
                summarized_values, "relative_profile_distortion_mean", "relative_profile_distortion_stderror")
            relative_profile_distortion_filename = os.path.join(structure.exploratory_path(), "relative_profile_distortion.png")
            absolute_profile_distortion_signatures_plot = _diamond_perturbation_plot(
                summarized_values, "profile_distortion_signatures_mean", "profile_distortion_signatures_stderror")
            absolute_profile_distortion_signatures_filename = os.path.join(structure.exploratory_path(), "profile_distortion_signatures.png")
            relative_profile_distortion_signatures_plot = _diamond_perturbation_plot(
                summarized_values, "relative_profile_distortion_signatures_mean", "relative_profile_distortion_signatures_stderror")
            relative_profile_distortion_signatures_filename = os.path.join(structure.exploratory_path(), "relative_profile_distortion_signatures.png")
            absolute_distance_error_plot = _diamond_perturbation_plot(
                summarized_values, "distance_error_mean", "distance_error_stderror")
            absolute_distance_error_filename = os.path.join(structure.exploratory_path(), "distance_error.png")
            relative_distance_error_plot = _diamond_perturbation_plot(
                summarized_values, "relative_distance_error_mean", "relative_distance_error_stderror")
            relative_distance_error_filename = os.path.join(structure.exploratory_path(), "relative_distance_error.png")
            absolute_distance_error_signatures_plot = _diamond_perturbation_plot(
                summarized_values, "distance_error_signatures_mean", "distance_error_signatures_stderror")
            absolute_distance_error_signatures_filename = os.path.join(structure.exploratory_path(), "distance_error_signatures.png")
            relative_distance_error_signatures_plot = _diamond_perturbation_plot(
                summarized_values, "relative_distance_error_signatures_mean", "relative_distance_error_signatures_stderror")
            relative_distance_error_signatures_filename = os.path.join(structure.exploratory_path(), "relative_distance_error_signatures.png")

            for file_name, plot in {absolute_profile_distortion_filename: absolute_profile_distortion_plot,
                                    relative_profile_distortion_filename: relative_profile_distortion_plot,
                                    absolute_profile_distortion_signatures_filename: absolute_profile_distortion_signatures_plot,
                                    relative_profile_distortion_signatures_filename: relative_profile_distortion_signatures_plot,
                                    absolute_distance_error_filename: absolute_distance_error_plot,
                                    relative_distance_error_filename: relative_distance_error_plot,
                                    absolute_distance_error_signatures_filename: absolute_distance_error_signatures_plot,
                                    relative_distance_error_signatures_filename: relative_distance_error_signatures_plot}.items():
                grdevices.png(file_name)
                plot.plot()
                grdevices.dev_off()

            if ctx.obj.get("save", False):
                # save model data for further adaptations
                rdata_filename = structure.intermediate_file_path(file_type="RData")
                output_r_data(
                    ctx=ctx, filename=rdata_filename,
                    absolute_profile_distortion_plot=absolute_profile_distortion_plot,
                    absolute_profile_distortion_filename=absolute_profile_distortion_filename,
                    relative_profile_distortion_plot=relative_profile_distortion_plot,
                    relative_profile_distortion_filename=relative_profile_distortion_filename,
                    absolute_profile_distortion_signatures_plot=absolute_profile_distortion_signatures_plot,
                    absolute_profile_distortion_signatures_filename=absolute_profile_distortion_signatures_filename,
                    relative_profile_distortion_signatures_plot=relative_profile_distortion_signatures_plot,
                    relative_profile_distortion_signatures_filename=relative_profile_distortion_signatures_filename,
                    absolute_distance_error_plot=absolute_distance_error_plot,
                    absolute_distance_error_filename=absolute_distance_error_filename,
                    relative_distance_error_plot=relative_distance_error_plot,
                    relative_distance_error_filename=relative_distance_error_filename,
                    absolute_distance_error_signatures_plot=absolute_distance_error_signatures_plot,
                    absolute_distance_error_signatures_filename=absolute_distance_error_signatures_filename,
                    relative_distance_error_signatures_plot=relative_distance_error_signatures_plot,
                    relative_distance_error_signatures_filename=relative_distance_error_signatures_filename,
                    summarized_values=summarized_values, result_dt=result_dt
                )


@click.command()
@click.pass_context
def analyse_diamond_level(ctx):
    structure = ctx.obj.get("structure")

    if ctx.obj.get("use_input", False):
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            from rpy2.robjects.packages import importr
            import rpy2.robjects.lib.ggplot2 as ggplot2
            from rpy2 import robjects

            grdevices = importr("grDevices")
            base = importr("base")
            datatable = importr("data.table")
            brewer = importr("RColorBrewer")

            input_data = json.load(input_file).get("data", None)
            result_dt = None
            for node_count, p_value_list in input_data.items():
                for p_value in p_value_list:
                    current_result = datatable.data_table(
                        p_value=base.as_integer(p_value),
                        level=robjects.IntVector([value for elem in p_value_list[p_value].get(
                            "raw", [0]) for value in elem] or [0])
                    )
                    if result_dt is None:
                        result_dt = current_result
                    else:
                        result_dt = datatable.rbindlist([result_dt, current_result])
            robjects.r("""
            complete_count <- function(data) {
                require(data.table)
                result <- setkey(data, p_value, level)[CJ(unique(p_value), unique(level)), .N, by=.EACHI]
            }
            """)
            complete_count = robjects.r["complete_count"]
            tmp_dt = complete_count(result_dt)
            plot = ggplot2.ggplot(tmp_dt) + ggplot2.aes_string(
                x="as.factor(p_value)", y="as.factor(level)", fill="N") + ggplot2.geom_tile(
                color="white", size=.1) + ggplot2.coord_equal() + ggplot2.scale_fill_gradientn(
                trans="log", colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Count")
            plot_filename = os.path.join(structure.exploratory_path(), "diamond_level.png")
            grdevices.png(plot_filename)
            plot.plot()
            grdevices.dev_off()

            if ctx.obj.get("save", False):
                # save model data for further adaptations
                rdata_filename = structure.intermediate_file_path(file_type="RData")
                output_r_data(
                    ctx=ctx, filename=rdata_filename, plot=plot, plot_filename=plot_filename,
                    result_dt=result_dt, tmp_dt=tmp_dt)


@click.command()
@click.pass_context
def analyse_diamonds(ctx):
    def _analyse_diamonds_plot(data, x_value, y_base):
        import rpy2.robjects.lib.ggplot2 as ggplot2
        return ggplot2.ggplot(data) + ggplot2.aes_string(x="%s" % x_value, y="%s_mean" % y_base,
                                                         color="as.factor(p_value)") + \
            ggplot2.geom_point() + ggplot2.geom_errorbar(width=.01) + ggplot2.aes_string(
            ymin="%s_mean-%s_stderror" % (y_base, y_base), ymax="%s_mean+%s_stderror" %
                                                                (y_base, y_base))

    def _analyse_diamonds_pplot(data, x_value, y_base):
        import rpy2.robjects.lib.ggplot2 as ggplot2
        return ggplot2.ggplot(data) + ggplot2.aes_string(x="%s" % x_value, y="%s_mean" % y_base) + \
            ggplot2.geom_ribbon(ggplot2.aes_string(ymin="%s_min" % y_base, ymax="%s_max" % y_base),
                                fill="green", alpha=.4) + ggplot2.geom_point() + \
            ggplot2.geom_errorbar(width=.01) + ggplot2.aes_string(
            ymin="%s_mean-%s_stderror" % (y_base, y_base),
            ymax="%s_mean+%s_stderror" % (y_base, y_base))
    structure = ctx.obj.get("structure")

    if ctx.obj.get("use_input", False):
        file_path = structure.input_file_path()
        with open(file_path, "r") as input_file:
            from rpy2.robjects.packages import importr
            import rpy2.robjects.lib.ggplot2 as ggplot2
            from rpy2.robjects.lib.dplyr import DataFrame
            from rpy2 import robjects

            r = robjects.r
            grdevices = importr("grDevices")
            datatable = importr("data.table")
            base = importr("base")
            utils = importr("utils")

            input_data = json.load(input_file).get("data", None)
            result_dt = None
            for node_count, p_value_list in input_data.items():
                for p_value in p_value_list:
                    current_result = datatable.data_table(
                        node_count=base.as_integer(node_count),
                        p_value=base.as_integer(p_value),
                        current_node_count=base.as_integer(p_value_list[p_value].get("node_counts", [0])),
                        diamonds=base.unlist(p_value_list[p_value].get("diamonds", [0])),
                        diamond_nodes=base.unlist(p_value_list[p_value].get("diamond_nodes", [0])),
                        identity_count=base.unlist(p_value_list[p_value].get("identities", [0]))
                    )
                    if result_dt is None:
                        result_dt = current_result
                    else:
                        result_dt = datatable.rbindlist([result_dt, current_result])
            # summarise data table
            summarized_values = (DataFrame(result_dt)
                                 .group_by("p_value", "node_count")
                                 .summarize(diamond_mean="mean(diamonds)",
                                            diamond_nodes_mean="mean(diamond_nodes)",
                                            relative_diamond_mean="mean(diamonds/identity_count)",
                                            relative_diamond_nodes_mean="mean(diamond_nodes/current_node_count)",
                                            diamond_stderror="sd(diamonds)/sqrt(length(diamonds))",
                                            diamond_nodes_stderror="sd(diamond_nodes)/sqrt(length(diamond_nodes))",
                                            relative_diamond_stderror="sd(diamonds/identity_count)/sqrt(length(diamonds))",
                                            relative_diamond_nodes_stderror="sd(diamond_nodes/current_node_count)/sqrt(length(diamond_nodes))"))
            summarized_pvalues = (DataFrame(result_dt)
                                  .group_by("p_value")
                                  .summarize(diamond_mean="mean(diamonds)",
                                             diamond_nodes_mean="mean(diamond_nodes)",
                                             diamond_min="min(diamonds)",
                                             diamond_nodes_min="min(diamond_nodes)",
                                             diamond_max="max(diamonds)",
                                             diamond_nodes_max="max(diamond_nodes)",
                                             relative_diamond_mean="mean(diamonds/identity_count)",
                                             relative_diamond_nodes_mean="mean(diamond_nodes/current_node_count)",
                                             relative_diamond_min="min(diamonds/identity_count)",
                                             relative_diamond_nodes_min="min(diamond_nodes/current_node_count)",
                                             relative_diamond_max="max(diamonds/identity_count)",
                                             relative_diamond_nodes_max="max(diamond_nodes/current_node_count)",
                                             diamond_stderror="sd(diamonds)/sqrt(length(diamonds))",
                                             diamond_nodes_stderror="sd(diamond_nodes)/sqrt(length(diamond_nodes))",
                                             relative_diamond_stderror="sd(diamonds/identity_count)/sqrt(length(diamonds))",
                                             relative_diamond_nodes_stderror="sd(diamond_nodes/current_node_count)/sqrt(length(diamond_nodes))"))
            absolute_diamonds_plot = _analyse_diamonds_plot(summarized_values, "node_count", "diamond")
            relative_diamonds_plot = _analyse_diamonds_plot(summarized_values, "node_count", "relative_diamond")
            absolute_diamond_nodes_plot = _analyse_diamonds_plot(summarized_values, "node_count", "diamond_nodes")
            relative_diamond_nodes_plot = _analyse_diamonds_plot(summarized_values, "node_count", "relative_diamond_nodes")
            absolute_diamonds_pplot = _analyse_diamonds_pplot(summarized_pvalues, "p_value", "diamond")
            relative_diamonds_pplot = _analyse_diamonds_pplot(summarized_pvalues, "p_value", "relative_diamond")
            absolute_diamond_nodes_pplot = _analyse_diamonds_pplot(summarized_pvalues, "p_value", "diamond_nodes")
            relative_diamond_nodes_pplot = _analyse_diamonds_pplot(summarized_pvalues, "p_value", "relative_diamond_nodes")

            absolute_diamonds_filename = os.path.join(structure.exploratory_path(), "ncount_diamonds.png")
            relative_diamonds_filename = os.path.join(structure.exploratory_path(), "ncount_relative_diamonds.png")
            absolute_diamond_nodes_filename = os.path.join(structure.exploratory_path(), "ncount_diamond_nodes.png")
            relative_diamond_nodes_filename = os.path.join(structure.exploratory_path(), "ncount_relative_diamond_nodes.png")
            absolute_diamonds_pfilename = os.path.join(structure.exploratory_path(), "pcount_diamonds.png")
            relative_diamonds_pfilename = os.path.join(structure.exploratory_path(), "pcount_relative_diamonds.png")
            absolute_diamond_nodes_pfilename = os.path.join(structure.exploratory_path(), "pcount_diamond_nodes.png")
            relative_diamond_nodes_pfilename = os.path.join(structure.exploratory_path(), "pcount_relative_diamond_nodes.png")
            for filename, plot in {absolute_diamonds_filename: absolute_diamonds_plot,
                                   relative_diamonds_filename: relative_diamonds_plot,
                                   absolute_diamond_nodes_filename: absolute_diamond_nodes_plot,
                                   relative_diamond_nodes_filename: relative_diamond_nodes_plot,
                                   absolute_diamonds_pfilename: absolute_diamonds_pplot,
                                   relative_diamonds_pfilename: relative_diamonds_pplot,
                                   absolute_diamond_nodes_pfilename: absolute_diamond_nodes_pplot,
                                   relative_diamond_nodes_pfilename: relative_diamond_nodes_pplot}.items():
                grdevices.png(filename)
                plot.plot()
                grdevices.dev_off()

            if ctx.obj.get("save", False):
                # save model data for further adaptations
                rdata_filename = structure.intermediate_file_path(file_type="RData")
                output_r_data(
                    ctx=ctx, filename=rdata_filename, absolute_diamonds_plot=absolute_diamonds_plot,
                    relative_diamonds_plot=relative_diamonds_plot,
                    absolute_diamond_nodes_plot=absolute_diamond_nodes_plot,
                    relative_diamond_nodes_plot=relative_diamond_nodes_plot,
                    absolute_diamonds_pplot=absolute_diamonds_pplot,
                    relative_diamonds_pplot=relative_diamonds_pplot,
                    absolute_diamond_nodes_pplot=absolute_diamond_nodes_pplot,
                    relative_diamond_nodes_pplot=relative_diamond_nodes_pplot,
                    absolute_diamonds_filename=absolute_diamonds_filename,
                    relative_diamonds_filename=relative_diamonds_filename,
                    absolute_diamond_nodes_filename=absolute_diamond_nodes_filename,
                    relative_diamond_nodes_filename=relative_diamond_nodes_filename,
                    absolute_diamonds_pfilename=absolute_diamonds_pfilename,
                    relative_diamonds_pfilename=relative_diamonds_pfilename,
                    absolute_diamond_nodes_pfilename=absolute_diamond_nodes_pfilename,
                    relative_diamond_nodes_pfilename=relative_diamond_nodes_pfilename,
                    summarized_values=summarized_values,
                    summarized_pvalues=summarized_pvalues, result_dt=result_dt
                )


@click.command()
@click.pass_context
def analyse_attribute_weight(ctx):
    structure = ctx.obj.get("structure")

    if ctx.obj.get("save", False):
        if ctx.obj.get("use_input", False):
            file_path = structure.input_file_path()
            with open(file_path, "r") as input_file:
                from rpy2.robjects.packages import importr
                import rpy2.robjects.lib.ggplot2 as ggplot2
                from rpy2 import robjects

                base = importr("base")
                grdevices = importr("grDevices")
                datatable = importr("data.table")
                brewer = importr("RColorBrewer")

                input_data = json.load(input_file).get("data", None)
                data_trees = input_data.get("files", [])

                # pre-cache values for data.table generation
                result_indexes = []
                weights = []
                statistics = []
                trees = []
                tree_sizes = []
                prototypes = []
                prototype_sizes = []
                decorators = []
                values = []
                # first collect the data for the first data.table before creating the next one
                for result_index, result in enumerate(input_data.get("results", [])):
                    for result_entry in result:
                        try:
                            data_tree_sizes = result_entry.get("decorator", {})["data"]["prototypes"]["original"][0]
                        except KeyError:
                            data_tree_sizes = []
                        algorithm = result_entry.get("algorithm", None)
                        for key in ["SetStatistics", "SplittedStatistics"]:
                            if key in algorithm:
                                statistic = key
                        weight = float(algorithm.split("=")[-1].split(")")[0])
                        for decorator_key in result_entry.get("decorator", {}):
                            if "matrix" not in decorator_key:
                                # skip other decorators
                                continue
                            decorator = result_entry.get("decorator").get(decorator_key, None)
                            for index, ensemble in enumerate(decorator):
                                for tree in ensemble:
                                    tree = [value if value is not None else 0 for value in tree]

                                    result_indexes.extend([result_index for _ in range(len(tree))])
                                    weights.extend([weight for _ in range(len(tree))])
                                    statistics.append([statistic for _ in range(len(tree))])
                                    trees.append([data_trees[result_index][index] for _ in range(len(tree))])
                                    tree_sizes.extend([data_tree_sizes[index] for _ in range(len(tree))])
                                    prototypes.extend(data_trees[result_index])
                                    prototype_sizes.extend(data_tree_sizes)
                                    decorators.extend([decorator_key for _ in range(len(tree))])
                                    values.extend(tree)
                result_dt = datatable.data_table(
                    repetition=base.unlist(result_indexes),
                    weight=base.unlist(weights),
                    statistic=base.unlist(statistics),
                    tree=base.unlist(trees),
                    tree_size=base.unlist(tree_sizes),
                    prototype=base.unlist(prototypes),
                    prototype_size=base.unlist(prototype_sizes),
                    decorator=base.unlist(decorators),
                    value=base.unlist(values)
                )
                result_indexes = []
                weights = []
                statistics = []
                decorators = []
                errors = []
                trees = []
                tree_sizes = []
                prototypes = []
                prototype_sizes = []
                # create the data for the second data.table
                for result_index, result in enumerate(input_data.get("results", [])):
                    for result_entry in result:
                        try:
                            data_tree_sizes = result_entry.get("decorator", {})["data"]["prototypes"]["original"][0]
                        except KeyError:
                            data_tree_sizes = []
                        algorithm = result_entry.get("algorithm", None)
                        for key in ["SetStatistics", "SplittedStatistics"]:
                            if key in algorithm:
                                statistic = key
                        weight = float(algorithm.split("=")[-1].split(")")[0])
                        for decorator_key in result_entry.get("decorator", {}):
                            if "matrix" not in decorator_key:
                                # skip other decorators
                                continue
                            decorator = result_entry.get("decorator").get(decorator_key, None)
                            for index, ensemble in enumerate(decorator):
                                for tree in ensemble:
                                    for column_index in range(index + 1):
                                        left_value = tree[column_index]
                                        right_value = decorator[column_index][0][index]
                                        if left_value is None or right_value is None:
                                            continue
                                        error = uncorrelated_relative_max_distance_deviation([
                                            (left_value, data_tree_sizes[index] * 2,
                                             right_value, data_tree_sizes[column_index] * 2,)
                                        ], None if index != column_index else 0)
                                        result_indexes.append(result_index)
                                        weights.append(weight)
                                        statistics.append(statistic)
                                        decorators.append(decorator_key)
                                        errors.append(error)
                                        trees.append(data_trees[result_index][index])
                                        tree_sizes.append(data_tree_sizes[index])
                                        prototypes.append(data_trees[result_index][column_index])
                                        prototype_sizes.append(data_tree_sizes[column_index])
                calculated_dt = datatable.data_table(
                    weight=base.unlist(weights),
                    statistic=base.unlist(statistics),
                    decorator=base.unlist(decorators),
                    error=base.unlist(errors),
                    tree=base.unlist(trees),
                    tree_size=base.unlist(tree_sizes),
                    prototype=base.unlist(prototypes),
                    prototype_size=base.unlist(prototype_sizes)
                )

            robjects.r("""
                create_cut <- function(dt, error_field, decorator, statistics, diagonal=FALSE) {
                    require(data.table)
                    tmp <- NULL
                    if (diagonal) {
                        tmp <- dt[statistic==statistics & decorator==decorator & tree==prototype, ]
                    } else {
                        tmp <- dt[statistic==statistics & decorator==decorator & tree!=prototype, ]
                    }
                    breaks <- min(30, length(unique(unlist(tmp[, error_field, with=F]))))
                    if (breaks > 1) {
                        tmp$cut <- cut(unlist(tmp[, error_field, with=F]), breaks=breaks, right=F)
                    } else {
                        tmp$cut <- factor(unlist(tmp[, error_field, with=F]), ordered=T)
                    }
                    tmp <- tmp[,.(count=.N), by=list(cut, weight)]
                    setkey(tmp, cut, weight)
                    tmp <- tmp[CJ(factor(levels(tmp[,cut]), levels(tmp$cut), ordered=T), tmp[,weight], unique=T)]
                    tmp
                }
            """)
            create_cut = robjects.r["create_cut"]
            tmp_dt = create_cut(calculated_dt, "error", "matrix", "SplittedStatistics")
            # create a heatmap for our errors
            error_heatmap = ggplot2.ggplot(tmp_dt) + ggplot2.aes_string(x="weight", y="cut", fill="count") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Count")
            error_heatmap_filename = os.path.join(structure.exploratory_path(), "error_heatmap.png")
            # create heatmap for diagonal
            tmp_dt = create_cut(calculated_dt, "error", "matrix", "SplittedStatistics", "TRUE")
            diagonal_error_heatmap = ggplot2.ggplot(tmp_dt) + ggplot2.aes_string(x="weight", y="cut", fill="count") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Count")
            diagonal_error_heatmap_filename = os.path.join(structure.exploratory_path(), "error_heatmap_diagonal.png")
            # heatmap for SetStatistics
            tmp_dt = create_cut(calculated_dt, "error", "matrix", "SetStatistics")
            set_error_heatmap = ggplot2.ggplot(tmp_dt) + ggplot2.aes_string(x="weight", y="cut", fill="count") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Count")
            set_error_heatmap_filename = os.path.join(structure.exploratory_path(), "set_error_heatmap.png")
            # heatmap for diagonal for SetStatistics
            tmp_dt = create_cut(calculated_dt, "error", "matrix", "SetStatistics", "TRUE")
            diagonal_set_error_heatmap = ggplot2.ggplot(tmp_dt) + ggplot2.aes_string(x="weight", y="cut", fill="count") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Count")
            diagonal_set_error_heatmap_filename = os.path.join(structure.exploratory_path(), "set_error_heatmap_diagonal.png")

            robjects.r("""
                create_cut_tree_sizes <- function(dt, error_field, decorator, statistics, selected_weight, diagonal=FALSE) {
                    require(data.table)
                    tmp <- NULL
                    if (diagonal) {
                        tmp <- dt[statistic==statistics & decorator==decorator & weight==selected_weight & tree==prototype, ]
                    } else {
                        tmp <- dt[statistic==statistics & decorator==decorator & weight==selected_weight & tree!=prototype, ]
                    }
                    breaks <- min(30, length(unique(tmp$tree_size)))
                    if (breaks > 1) {
                        tmp$cut <- cut(tmp$tree_size, breaks=breaks, right=F, ordered_result=T)
                        tmp$pcut <- cut(tmp$prototype_size, breaks=breaks, right=F, ordered_result=T)
                    } else {
                        tmp$cut <- factor(tmp$tree_size, ordered=T)
                        tmp$pcut <- factor(tmp$prototype_size, ordered=T)
                    }
                    tmp <- tmp[,.(mean=mean(get(error_field))), by=list(cut, pcut)]
                    setkey(tmp, cut, pcut)
                    tmp <- tmp[CJ(factor(levels(tmp[,cut]), levels(tmp$cut), ordered=T), tmp[,pcut], unique=T)]
                    setkey(tmp, pcut, cut)
                    tmp <- tmp[CJ(factor(levels(tmp[,pcut]), levels(tmp$pcut), ordered=T), tmp[,cut], unique=T)]
                    tmp
                }
            """)
            create_cut_tree_sizes = robjects.r["create_cut_tree_sizes"]
            size_tmp_dt = create_cut_tree_sizes(calculated_dt, "error", "matrix", "SplittedStatistics", 0)
            # create heatmap plot for tree sizes
            tree_size_heatmap = ggplot2.ggplot(size_tmp_dt) + ggplot2.aes_string(x="cut", y="pcut", fill="mean") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Error")
            tree_size_heatmap_filename = os.path.join(structure.exploratory_path(), "tree_size_heatmap.png")
            # tree sizes for for diagonal
            size_tmp_dt = create_cut_tree_sizes(calculated_dt, "error", "matrix", "SplittedStatistics", 0, "TRUE")
            diagonal_tree_size_heatmap = ggplot2.ggplot(size_tmp_dt) + ggplot2.aes_string(x="cut", y="pcut", fill="mean") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Error")
            diagonal_tree_size_heatmap_filename = os.path.join(structure.exploratory_path(), "tree_size_heatmap_diagonal.png")
            # create heatmap for tree sizes for SetStatistics
            # a different value for weight is required here, because only for values != 0, 0.5 or 1
            # we get values for plotting
            size_tmp_dt = create_cut_tree_sizes(calculated_dt, "error", "matrix", "SetStatistics", 0.1)
            set_tree_size_heatmap = ggplot2.ggplot(size_tmp_dt) + ggplot2.aes_string(x="cut", y="pcut", fill="mean") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Error")
            set_tree_size_heatmap_filename = os.path.join(structure.exploratory_path(), "set_tree_size_heatmap.png")
            # heatmap for tree sizes for diagonal for SetStatistics
            size_tmp_dt = create_cut_tree_sizes(calculated_dt, "error", "matrix", "SetStatistics", 0.1, "TRUE")
            diagonal_set_tree_size_heatmap = ggplot2.ggplot(size_tmp_dt) + ggplot2.aes_string(x="cut", y="pcut", fill="mean") + \
                            ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name='Reds'), na_value="white", name="Error")
            diagonal_set_tree_size_heatmap_filename = os.path.join(structure.exploratory_path(), "set_tree_size_heatmap_diagonal.png")

            robjects.r("""
            create_cut_diagonal <- function(dt, decorator, statistics) {
                tmp <- dt[statistic==statistics & decorator==decorator & tree==prototype,]
                breaks <- min(30, length(unique(tmp$tree_size)))
                if (breaks > 1) {
                    tmp$cut <- cut(tmp$tree_size, breaks=30, right=F, ordered_result=T)
                } else {
                    tmp$cut <- factor(tmp$tree_size, ordered=T)
                }
                tmp <- tmp[,.(mean=mean(error)),by=list(weight, cut)]
                tmp
            }
            """)
            create_cut_diagonal = robjects.r["create_cut_diagonal"]
            diagonal_tmp_dt = create_cut_diagonal(calculated_dt, "matrix", "SplittedStatistics")
            # diagonal special plot
            diagonal_heatmap = ggplot2.ggplot(diagonal_tmp_dt) + ggplot2.aes_string(x="weight", y="cut", fill="mean") + \
                ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name="Reds"), na_value="white", name="Error")
            diagonal_heatmap_filename = os.path.join(structure.exploratory_path(), "diagonal.png")
            # set specific diagonal
            diagonal_tmp_dt = create_cut_diagonal(calculated_dt, "matrix", "SetStatistics")
            set_diagonal_heatmap = ggplot2.ggplot(diagonal_tmp_dt) + ggplot2.aes_string(x="weight", y="cut", fill="mean") + \
                ggplot2.geom_tile(color="white", size=.1) + ggplot2.scale_fill_gradientn(
                colours=brewer.brewer_pal(n=9, name="Reds"), na_value="white", name="Error")
            set_diagonal_heatmap_filename = os.path.join(structure.exploratory_path(), "set_diagonal.png")

            # perform the plotting
            for plot, filename in [(error_heatmap, error_heatmap_filename,),
                                   (set_error_heatmap, set_error_heatmap_filename,),
                                   (tree_size_heatmap, tree_size_heatmap_filename,),
                                   (set_tree_size_heatmap, set_tree_size_heatmap_filename,),
                                   (diagonal_error_heatmap, diagonal_error_heatmap_filename,),
                                   (diagonal_set_error_heatmap, diagonal_set_error_heatmap_filename,),
                                   (diagonal_tree_size_heatmap, diagonal_tree_size_heatmap_filename,),
                                   (diagonal_set_tree_size_heatmap, diagonal_set_tree_size_heatmap_filename,),
                                   (diagonal_heatmap, diagonal_heatmap_filename,),
                                   (set_diagonal_heatmap, set_diagonal_heatmap_filename,)]:
                grdevices.png(filename)
                plot.plot()
                grdevices.dev_off()

            # save model data for further adaptations
            rdata_filename = structure.intermediate_file_path(file_type="RData")
            output_r_data(
                ctx=ctx, filename=rdata_filename, result_dt=result_dt, calculated_dt=calculated_dt,
                error_heatmap=error_heatmap, error_heatmap_filename=error_heatmap_filename,
                set_error_heatmap=set_error_heatmap, set_error_heatmap_filename=set_error_heatmap_filename,
                tree_size_heatmap=tree_size_heatmap, tree_size_heatmap_filename=tree_size_heatmap_filename,
                set_tree_size_heatmap=set_tree_size_heatmap, set_tree_size_heatmap_filename=set_tree_size_heatmap_filename,
                diagonal_error_heatmap=diagonal_error_heatmap, diagonal_error_heatmap_filename=diagonal_error_heatmap_filename,
                diagonal_tree_size_heatmap=diagonal_tree_size_heatmap, diagonal_tree_size_heatmap_filename=diagonal_tree_size_heatmap_filename,
                diagonal_set_tree_size_heatmap=diagonal_set_tree_size_heatmap, diagonal_set_tree_size_heatmap_filename=diagonal_set_tree_size_heatmap_filename,
                diagonal_heatmap=diagonal_heatmap, diagonal_heatmap_filename=diagonal_heatmap_filename,
                set_diagonal_heatmap=set_diagonal_heatmap, set_diagonal_heatmap_filename=set_diagonal_heatmap_filename
            )


@click.command()
@click.option("--seed", "seed", type=int, default=None)
@click.pass_context
def analyse_attribute_metric(ctx, seed):
    """
    Method performs specific analysis on SetStatistics and SplittedStatistics based on two different
    measurements. First considers the change of the mean of one specific distribution to change the
    overlap with the input distribution. The second changes the number of samples picked for a given
    distribution.

    Method only produces plots and no further results.

    :param ctx:
    """
    structure = ctx.obj.get("structure")

    if ctx.obj.get("save", False):
        from rpy2.robjects.packages import STAP, importr
        import rpy2.robjects.lib.ggplot2 as ggplot2
        from rpy2.robjects.lib.dplyr import DataFrame
        from rpy2 import robjects

        with open(os.path.join(os.path.join(structure.base_script_path(), "R"),
                               "statistics.R")) as statistics_r_file:
            statistics_r_string = statistics_r_file.read()
        statistics = STAP(statistics_r_string, "statistics")
        with open(os.path.join(os.path.join(structure.base_script_path(), "R"),
                  "utils.R")) as utils_r_file:
            utils_r_string = utils_r_file.read()
        rutils = STAP(utils_r_string, "rutils")

        base = importr("base")
        stats = importr("stats")
        grdevices = importr("grDevices")
        datatable = importr("data.table")

        if seed is not None:
            base.set_seed(seed)

        base_mean = 5
        # generate initial distribution to learn attribute statistics from
        distribution = stats.rnorm(1000, base_mean, 1)
        # generate distribution to test learned statistics from
        validation_distribution = stats.rnorm(1000, base_mean, 1)

        set_statistics = SetStatistics(
            convert=lambda value: int(round(value)), unconvert=lambda value: value)
        splitted_statistics = SplittedStatistics()
        for value in distribution:
            # add value to current statistic
            set_statistics.add(value)
            splitted_statistics.add(value)

        validation_distribution_dt = datatable.data_table(validation_distribution)
        validation_plot = ggplot2.ggplot(validation_distribution_dt) + ggplot2.aes_string(x="V1") + \
                          ggplot2.geom_histogram(binwidth=.5)

        set_values = [entry for value in ([statistic.value] * statistic.count for
                                          statistic in set_statistics) for entry in value]
        set_values_dt = datatable.data_table(base.unlist(set_values))
        set_plot = ggplot2.ggplot(set_values_dt) + ggplot2.aes_string(x="V1") + ggplot2.geom_histogram(binwidth=1)
        # output current results
        splitted_plot = ggplot2.ggplot(datatable.data_table(x=rutils.plot_x_get_range(validation_plot))) + \
                        ggplot2.aes_string(x="x")
        for index, statistic in enumerate(splitted_statistics):
            splitted_plot = splitted_plot + ggplot2.stat_function(
                fun=statistics.stats_dnorm_height,
                args=robjects.FloatVector([statistic.value, math.sqrt(
                    statistic.variance if statistic.variance is not None else 0), statistic.count]))

        # generate plots
        distribution_filename = os.path.join(structure.exploratory_path(), "validation_distribution.png")
        set_filename = os.path.join(structure.exploratory_path(), "set_distribution.png")
        splitted_filename = os.path.join(structure.exploratory_path(), "splitted_distribution.png")

        for filename, plot in {distribution_filename: validation_plot,
                               set_filename: set_plot,
                               splitted_filename: splitted_plot}.items():
            grdevices.png(filename)
            plot.plot()
            grdevices.dev_off()

        current_mean = base_mean - .11
        overlaps, set_distances, splitted_distances = [], [], []
        while len(overlaps) == 0 or overlaps[-1] > 0.1:
            current_mean += .1
            for _ in range(10):
                overlaps.append(statistics.stats_gaussians_overlap(base_mean, current_mean, 1)[0])
                current_distribution = stats.rnorm(1000, mean=current_mean, sd=1)
                set_distance, _ = _distance_and_statistic(set_statistics, current_distribution)
                splitted_distance, _ = _distance_and_statistic(splitted_statistics, current_distribution)
                set_distances.append(set_distance)
                splitted_distances.append(splitted_distance)
        overlap_dt = datatable.data_table(overlaps=base.rep(base.unlist(overlaps), 2),
                                          distances=base.unlist(set_distances + splitted_distances),
                                          type=base.as_factor(base.rep(base.unlist(
                                              ["set", "splitted"]), each=len(overlaps))))
        robjects.globalenv["optimal_overlaps"] = rutils.analysis_attribute_distances_optimal_overlaps
        upper_overlap_plot, lower_overlap_plot, summarized_overlap_dt, ratio_summarized_overlap_dt = \
            _upper_and_lower_plot("overlaps", overlap_dt, "optimal_overlaps", 2000)
        _layout_plots(os.path.join(structure.exploratory_path(), "overlaps.png"),
                      upper_overlap_plot, lower_overlap_plot)

        current_count = 2000
        counts, set_distances, splitted_distances = [], [], []
        while current_count > 0:
            for _ in range(10):
                counts.append(current_count)
                current_distribution = stats.rnorm(current_count, mean=base_mean, sd=1)
                set_distance, _ = _distance_and_statistic(set_statistics, current_distribution)
                splitted_distance, _ = _distance_and_statistic(splitted_statistics, current_distribution)
                set_distances.append(set_distance)
                splitted_distances.append(splitted_distance)
            current_count -= 11
        count_dt = datatable.data_table(counts=base.rep(base.unlist(counts), 2),
                                        distances=base.unlist(set_distances + splitted_distances),
                                        type=base.as_factor(base.rep(base.unlist(
                                            ["set", "splitted"]), each=len(counts))))
        robjects.globalenv["optimal_counts"] = rutils.analysis_attribute_distances_optimal_counts
        upper_count_plot, lower_count_plot, summarized_count_dt, ratio_summarized_count_dt = \
            _upper_and_lower_plot("counts", count_dt, "optimal_counts", 1000)
        _layout_plots(os.path.join(structure.exploratory_path(), "counts.png"),
                      upper_count_plot, lower_count_plot)

        rdata_filename = structure.intermediate_file_path(file_type="RData")
        output_r_data(
            ctx=ctx, filename=rdata_filename, validation_distribution_dt=validation_distribution_dt,
            validation_plot=validation_plot, distribution_filename=distribution_filename,
            set_values_dt=set_values_dt, set_plot=set_plot, set_filename=set_filename,
            splitted_filename=splitted_filename, splitted_plot=splitted_plot,
            overlap_dt=overlap_dt, count_dt=count_dt, summarized_count_dt=summarized_count_dt,
            ratio_summarized_count_dt=ratio_summarized_count_dt,
            summarized_overlap_dt=summarized_overlap_dt, ratio_summarized_overlap_dt=ratio_summarized_overlap_dt,
            upper_overlap_plot=upper_overlap_plot, lower_overlap_plot=lower_overlap_plot,
            upper_count_plot=upper_count_plot, lower_count_plot=lower_count_plot,
            optimal_overlaps=robjects.r["optimal_overlaps"], optimal_counts=robjects.r["optimal_counts"]
        )


def _upper_and_lower_plot(variant, overlap_dt, optimal_rfunction_name, base_distance):
    from rpy2 import robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects.lib.ggplot2 import ggplot2
    from rpy2.robjects.lib.dplyr import DataFrame

    long_format_dt = (DataFrame(overlap_dt)
                      .group_by(variant, "type")
                      .summarize(distance_mean="mean(distances)",
                                 distance_stderror="sd(distances)/sqrt(length(distances))"))
    ratio_format_dt = (DataFrame(long_format_dt)
                       .group_by(variant, "type")
                       .mutate(mean_ratio="distance_mean/%s(%s, %s)" % (optimal_rfunction_name, variant, base_distance),
                               ratio_error="distance_stderror/distance_mean"))

    upper_plot = (ggplot2.ggplot(long_format_dt) + ggplot2.aes_string(
        x=variant, y="distance_mean", colour="type") + ggplot2.geom_point() +
                          ggplot2.geom_errorbar() + ggplot2.aes_string(
        ymin="distance_mean-distance_stderror", ymax="distance_mean+distance_stderror") +
                          ggplot2.stat_function(
                              fun=optimal_rfunction_name, args=robjects.IntVector([base_distance])) +
                          ggplot2.theme(
                              **{"legend.justification": robjects.IntVector([1,1]),
                                 "legend.position": robjects.IntVector([1,1]),
                                 "axis.title.x": ggplot2.element_blank(),
                                 "axis.text.x": ggplot2.element_blank()}))
    lower_plot = (ggplot2.ggplot(ratio_format_dt) + ggplot2.aes_string(x=variant, y="mean_ratio", colour="type") +
                  ggplot2.geom_point() + ggplot2.geom_errorbar(ggplot2.aes_string(
        ymin="mean_ratio-ratio_error", ymax="mean_ratio+ratio_error")) +
                  ggplot2.theme(**{"legend.position": "none"}))
    return upper_plot, lower_plot, long_format_dt, ratio_format_dt


def _layout_plots(filename, first_plot, second_plot):
    from rpy2 import robjects
    from rpy2.robjects.packages import importr
    gridextra = importr("gridExtra")
    grdevices = importr("grDevices")

    robjects.r("""
        adapt_width <- function(plot_1, plot_2) {
            require(ggplot2)
            require(grid)
            require(gridExtra)
            plot_1 <- ggplot_gtable(ggplot_build(plot_1))
            plot_2 <- ggplot_gtable(ggplot_build(plot_2))
            maxWidth = unit.pmax(plot_1$widths[2:3], plot_2$widths[2:3])
            plot_1$widths[2:3] <- maxWidth
            plot_2$widths[2:3] <- maxWidth
            list(plot_1, plot_2)
        }
    """)
    adapt_width = robjects.r["adapt_width"]
    plots = adapt_width(first_plot, second_plot)

    grdevices.png(filename)
    gridextra.grid_arrange(plots[0], plots[1], heights=robjects.IntVector([3, 1]))
    grdevices.dev_off()


def _plot_distances_for_values(structure, distances, values, variant):
    from rpy2.robjects.packages import importr
    from rpy2.robjects.lib.dplyr import DataFrame
    import rpy2.robjects.lib.ggplot2 as ggplot2

    datatable = importr("data.table")
    base = importr("base")
    grdevices = importr("grDevices")

    distance_overlap_dt = datatable.data_table(values=base.unlist(values), distances=base.unlist(distances))
    summarized_values = (DataFrame(distance_overlap_dt)
                          .group_by("values")
                          .summarize(distance_mean="mean(distances)",
                                     distance_stderror="sd(distances)/sqrt(length(distances))"))

    filename = os.path.join(structure.exploratory_path(), "%s_distances.png" % variant)
    grdevices.png(filename)
    (ggplot2.ggplot(summarized_values) + ggplot2.aes_string(x="values", y="distance_mean") +
     ggplot2.geom_point() + ggplot2.geom_errorbar(width=.01) +
     ggplot2.aes_string(ymin="distance_mean-distance_stderror",
                        ymax="distance_mean+distance_stderror")).plot()
    grdevices.dev_off()


def _distance_and_statistic(base_statistic, distribution_values):
    current_statistic = type(base_statistic)()
    try:
        current_statistic._convert = base_statistic._convert
        current_statistic._unconvert = base_statistic._unconvert
    except AttributeError:
        pass
    distance = base_statistic.count()
    for value in distribution_values:
        current_distance = base_statistic.distance(value, count=current_statistic.count(value))
        if current_distance <= 0.5:
            distance -= 1 - current_distance
        else:
            distance += current_distance
        current_statistic.add(value)
    return distance, current_statistic

cli.add_command(analyse_compression)
cli.add_command(analyse_diamonds)
cli.add_command(analyse_diamond_perturbations)
cli.add_command(analyse_diamond_level)
cli.add_command(analyse_attribute_metric)
cli.add_command(analyse_attribute_weight)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
