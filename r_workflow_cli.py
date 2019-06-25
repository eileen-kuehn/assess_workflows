import os
import click
import logging

from cpy2py import TwinMaster, meta

from assess_workflows.generic.structure import Structure
from utility.exceptions import ExceptionFrame
from utility.report import LVL


@click.group()
@click.option("--basepath", "basepath", multiple=False, required=True)
@click.option("--workflow-name", "workflow_name", multiple=False, required=True)
@click.option("--step", "step", default=1, multiple=False)
@click.option("--use_input", "use_input", is_flag=True,
              help="Use input file specified for current task.")
@click.pass_context
def cli(ctx, basepath, workflow_name, step, use_input):
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name, step=step)
    ctx.obj["use_input"] = use_input


@click.command()
@click.pass_context
def performance_plots(ctx):
    twinterpreter.execute(do_performance_plots, ctx.obj)


def do_performance_plots(ctx_obj):
    from rpy2.robjects.packages import STAP, importr
    from rpy2.robjects.lib.dplyr import DataFrame
    import rpy2.robjects.lib.ggplot2 as ggplot2

    # TODO: check out if and how this works....
    dtplyr = importr("dtplyr")
    base = importr("base")
    grdevices = importr("grDevices")

    structure = ctx_obj.get("structure")

    with open(os.path.join(os.path.join(os.path.join(
            structure.basepath, "scripts"), "R"), "assess.R")) as assess_r_file:
        assess_r_string = assess_r_file.read()
    assess = STAP(assess_r_string, "assess")

    if ctx_obj.get("use_input", False):
        raw_data = assess.assess_read(structure.input_file_path())
        dt = assess.assess_get_performances(raw_data)
        dt_summary = assess.assess_summarise_performances(dt)

        performance_objects = ["signature_performance", "distance_performance", "performance"]
        signature_objects = list(base.unique(dt.rx2("signature")))
        print(signature_objects)

        def create_plot(filename, data, performance_object):
            grdevices.png(file=filename)
            (ggplot2.ggplot(data) +
             ggplot2.aes_string(x="dimension", y="mean.%s" % performance_object) +
             ggplot2.geom_point(ggplot2.aes_string(colour="signature")) +
             ggplot2.geom_errorbar(ggplot2.aes_string(
                 ymax="mean.%s+stderror.%s" % (performance_object, performance_object),
                 ymin="mean.%s-stderror.%s" % (performance_object, performance_object),
             ))).plot()
            grdevices.dev_off()

        for performance_object in performance_objects:
            # FIXME: plotting as EPS currently not working
            # create plot with all signatures inside
            create_plot(
                os.path.join(structure.exploratory_path(), "combined_%s.png" % performance_object),
                dt_summary,
                performance_object
            )

            # create plot for each signature
            for signature_object in signature_objects:
                create_plot(
                    os.path.join(structure.exploratory_path(), "%s_%s.png" %
                                 (signature_object, performance_object)),
                    DataFrame(dt_summary).filter("signature=='%s'" % signature_object),
                    performance_object
                )
    print assess

cli.add_command(performance_plots)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        twinterpreter = TwinMaster(twinterpreter_id='/usr/local/bin/python', executable='/usr/local/bin/python')
        twinterpreter.start()
        cli(obj={}, auto_envvar_prefix='DISS')
