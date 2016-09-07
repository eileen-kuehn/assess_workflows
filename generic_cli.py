"""
This file shall try to standardise my current workflows. It should provide a fixed structure for
contents that are created and used within the workflow.

To initialise the structure, the command init_environment can be used. It will create a structure
looking like this:

* base_folder
    * ./data
        * ./raw
        * ./processed
        * ./conferences <--- those should not be touched because they might already be published
    * ./figures
        * ./conferences <--- those should not be touched because they might already be published
        * ./exploratory
    * ./models
    * ./scripts
        * ./helper
            * ./conferences
        * ./md
        * ./python
        * ./R
    * ./text
    * ./workflows       <--- workflows are attached here...

There are specific workflow for tasks that are done. Each of those workflow gets its own folder
and also substructure. It should be self-contained by also has to be able to rely on the outer
folder structure, so that no raw or processed data needs to be copied. Only a specification of
which files to be used is needed here.

* ./workflows
    * ./<unique_workflow_id>
        * configuration.py <--- contains the configuration to run
        * run_workflow.sh  <--- contains the actual steps to be executed to get the results
        * README.md        <--- steps, thoughts, etc. are described in here
        * ./input          <--- either links or accumulated paths to data are stored in here
        * ./exploratory    <--- exploratory results are stored here
        * ./intermediate   <--- intermediate results are stored here
        * ./final          <--- final results are stored here

But how to integrate the workflows? How to integrate the actual research log? How to integrate
the stuff you are doing exploratory? Still no idea... but a first step taken here
"""
import logging
import datetime
import click
import os
import assess
import assess_workflows

from assess_workflows.generic.structure import Structure

from utility.exceptions import ExceptionFrame
from utility.report import LVL


@click.command()
@click.option("--name", "name")
@click.option("--basepath", "basepath")
def create_workflow(name, basepath):
    click.echo("using basepath %s" % basepath)
    structure = Structure(basepath=basepath)
    if name is None:
        name = "%s_%s" % (datetime.date.today(), "workflow")
    try:
        os.makedirs(structure.input_path(name))
        os.makedirs(structure.exploratory_path(name))
        os.makedirs(structure.intermediate_path(name))
        os.makedirs(structure.final_path(name))
        with open(structure.execution_file_path(name), "w") as run_file:
            run_file.write("#!/bin/bash\n")
            run_file.write("export DISS_WORKFLOW=%s\n" % structure.workflow_path(name))
            run_file.write("export DISS_CONFIGURATION=%s\n" %
                           structure.configuration_file_path(name))
            run_file.write("export DISS_WORKFLOW_LIB=%s\n" %
                           os.path.dirname(assess_workflows.__file__))
            run_file.write("export DISS_ASSESS_LIB=%s\n" %
                           os.path.dirname(assess.__file__))
        with open(structure.configuration_file_path(name), "w") as configuration_file:
            pass
    except OSError:
        click.echo("The workflow has already been created before, skipping here...")
    else:
        click.echo("created workflow structure at %s" % structure.workflow_path(name))


if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        create_workflow(auto_envvar_prefix='DISS')
