import os
import logging

import click
import rpy2
import math

from assess_workflows.generic.structure import Structure
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
def analyse_attribute_metric(ctx):
    structure = ctx.obj.get("structure")

    if ctx.obj.get("save", False):
        from rpy2.robjects.packages import STAP, importr
        import rpy2.robjects.lib.ggplot2 as ggplot2
        from rpy2.robjects.lib.dplyr import DataFrame

        with open(os.path.join(os.path.join(structure.base_script_path(), "R"),
                               "statistics.R")) as statistics_r_file:
            statistics_r_string = statistics_r_file.read()
        statistics = STAP(statistics_r_string, "statistics")

        utils = importr("utils")
        base = importr("base")
        stats = importr("stats")
        grdevices = importr("grDevices")
        datatable = importr("data.table")
        fpc = importr("fpc")

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

        # generate plots
        filename = os.path.join(structure.exploratory_path(), "validation_distribution.png")
        grdevices.png(filename)
        (ggplot2.ggplot(validation_distribution_dt) + ggplot2.aes_string(x="V1") + ggplot2.geom_histogram(binwidth=.5)).plot()
        grdevices.dev_off()

        set_values = [entry for value in (
            [statistic.value] * statistic.count for statistic in set_statistics)
                      for entry in value]
        set_values_dt = datatable.data_table(base.unlist(set_values))
        set_filename = os.path.join(structure.exploratory_path(), "set_distribution.png")
        grdevices.png(set_filename)
        (ggplot2.ggplot(set_values_dt) + ggplot2.aes_string(x="V1") + ggplot2.geom_histogram(binwidth=1)).plot()
        grdevices.dev_off()

        # output current results
        splitted_dt = None
        for index, statistic in enumerate(splitted_statistics):
            current_dt = datatable.data_table(
                    value=base.as_factor(index),
                    distribution=stats.rnorm(
                        statistic.count,
                        mean=statistic.value,
                        sd=math.sqrt(statistic.variance if statistic.variance is not None else 0)))
            if splitted_dt is None:
                splitted_dt = current_dt
            else:
                splitted_dt = datatable.rbindlist([splitted_dt, current_dt])

        splitted_filename = os.path.join(structure.exploratory_path(), "splitted_distribution.png")
        grdevices.png(splitted_filename)
        (ggplot2.ggplot(splitted_dt) + ggplot2.aes_string(
            x="distribution", y="..count..", fill="value") +
         ggplot2.geom_density(position="stack", adjust=1)).plot()
        grdevices.dev_off()

        current_mean = base_mean - .1
        overlaps = []
        set_distances = []
        splitted_distances = []
        while len(overlaps) == 0 or overlaps[-1] > 0.1:
            current_mean += .1
            for _ in range(10):
                overlaps.append(statistics.stats_gaussians_overlap(base_mean, current_mean, 1)[0])
                current_distribution = stats.rnorm(1000, mean=current_mean, sd=1)
                set_distance, _ = _distance_and_statistic(set_statistics, current_distribution)
                splitted_distance, _ = _distance_and_statistic(splitted_statistics, current_distribution)
                set_distances.append(set_distance)
                splitted_distances.append(splitted_distance)
        _plot_distances_for_values(structure, set_distances, overlaps, "set_overlap")
        _plot_distances_for_values(structure, splitted_distances, overlaps, "splitted_overlap")

        current_count = 2000
        counts = []
        set_distances = []
        splitted_distances = []
        while current_count > 0:
            for _ in range(10):
                counts.append(current_count)
                current_distribution = stats.rnorm(current_count, mean=base_mean, sd=1)
                set_distance, _ = _distance_and_statistic(set_statistics, current_distribution)
                splitted_distance, _ = _distance_and_statistic(splitted_statistics, current_distribution)
                set_distances.append(set_distance)
                splitted_distances.append(splitted_distance)
            current_count -= 10
        _plot_distances_for_values(structure, set_distances, counts, "set_counts")
        _plot_distances_for_values(structure, splitted_distances, counts, "splitted_counts")


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

cli.add_command(analyse_attribute_metric)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
