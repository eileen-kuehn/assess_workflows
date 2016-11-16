import os
import datetime

import assess_workflows


class Structure(object):
    def __init__(self, basepath=None, name=None, step=1):
        self.basepath = basepath or os.environ["DISS_BASEPATH"]
        self.workflows_basepath = os.path.join(self.basepath, os.environ.get(
            "DISS_WORKFLOWFOLDER_NAME", "workflows"))
        self.name = name
        self.step = step

    @staticmethod
    def analysis_cli():
        return os.path.join(os.path.join(os.path.dirname(
            os.path.realpath(assess_workflows.__file__)), os.pardir), "analysis_cli.py")

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

    @staticmethod
    def clustering_cli():
        return os.path.join(os.path.join(os.path.dirname(
            os.path.realpath(assess_workflows.__file__)), os.pardir), "clustering_cli.py")

    def input_path(self, name=None):
        return self._get_folder_path(self.workflow_path(name or self.name), "input")

    def final_path(self, name=None):
        return self._get_folder_path(self.workflow_path(name or self.name), "final")

    def intermediate_path(self, name=None):
        return self._get_folder_path(self.workflow_path(name or self.name), "intermediate")

    def exploratory_path(self, name=None):
        return self._get_folder_path(self.workflow_path(name or self.name), "exploratory")

    def workflow_path(self, name=None):
        return os.path.join(self.workflows_basepath, name or self.name)

    def workflow_config_file_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "workflow_configuration.json")

    def execution_file_path(self, name=None):
        return os.path.join(self.workflow_path(name or self.name), "run_workflow.py")

    def configuration_module(self, step=None):
        return "configuration_%d" % (step or self.step)

    def configuration_file_path(self, name=None, step=None):
        return os.path.join(self.workflow_path(name or self.name),
                            "configuration_%d.py" % (step or self.step))

    def intermediate_file_path(self, name=None, step=None, variant=None, file_type=None):
        if file_type is None:
            file_type = "json"
        if variant:
            return os.path.join(self.intermediate_path(name), "output_%d-%s.%s" %
                                (step or self.step, variant, file_type))
        return os.path.join(self.intermediate_path(name), "output_%d.%s" % (step or self.step, file_type))

    def input_file_path(self, name=None, step=None):
        return os.path.join(self.input_path(name), "input_%d.json" % (step or self.step))

    def final_file_path(self, name=None, step=None, file_type="json"):
        return os.path.join(self.final_path(name), "%s_final_%d.%s" %
                            (datetime.date.today(), step or self.step, file_type))

    def _get_folder_path(self, workflow_name, folder_name):
        folder_path = os.path.join(workflow_name, folder_name)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        return folder_path
