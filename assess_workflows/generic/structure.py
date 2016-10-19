import os
import datetime

import assess_workflows


class Structure(object):
    def __init__(self, basepath=None, name=None, step=1):
        self.basepath = basepath or os.environ["DISS_BASEPATH"]
        self.workflows_basepath = os.path.join(self.basepath, "workflows")
        self.name = name
        self.step = step

    @staticmethod
    def workflow_cli():
        return os.path.join(os.path.join(os.path.dirname(
            os.path.realpath(assess_workflows.__file__)), os.pardir), "workflow.py")

    @staticmethod
    def generic_cli():
        return os.path.join(os.path.join(os.path.dirname(
            os.path.realpath(assess_workflows.__file__)), os.pardir), "generic_cli.py")

    @staticmethod
    def data_selection_cli():
        return os.path.join(os.path.join(os.path.dirname(
            os.path.realpath(assess_workflows.__file__)), os.pardir), "data_selection_cli.py")

    @staticmethod
    def r_workflow_cli():
        return os.path.join(os.path.join(os.path.dirname(
            os.path.realpath(assess_workflows.__file__)), os.pardir), "r_workflow_cli.py")

    @staticmethod
    def transformation_cli():
        return os.path.join(os.path.join(os.path.dirname(
            os.path.realpath(assess_workflows.__file__)), os.pardir), "transformation_cli.py")

    def input_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "input")

    def final_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "final")

    def intermediate_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "intermediate")

    def exploratory_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "exploratory")

    def workflow_path(self, name=None):
        return os.path.join(self.workflows_basepath, name or self.name)

    def workflow_config_file_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "workflow_configuration.json")

    def execution_file_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "run_workflow.py")

    def configuration_file_path(self, name=None, step=None):
        return os.path.join(self.workflow_path(name or self.name),
                            "configuration_%d.py" % (step or self.step))

    def intermediate_file_path(self, name=None, step=None, variant=None):
        if variant:
            return os.path.join(self.intermediate_path(name), "output_%d-%s.json" %
                                (step or self.step, variant))
        return os.path.join(self.intermediate_path(name), "output_%d.json" % (step or self.step))

    def input_file_path(self, name=None, step=None):
        return os.path.join(self.input_path(name), "input_%d.json" % (step or self.step))

    def final_file_path(self, name=None, step=None):
        return os.path.join(self.final_path(name), "%s_final_%d.json" %
                            (datetime.date.today(), step or self.step))
