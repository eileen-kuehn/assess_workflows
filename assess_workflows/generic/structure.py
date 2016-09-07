import os


class Structure(object):
    def __init__(self, basepath=None, name=None):
        self.basepath = basepath or os.environ["DISS_BASEPATH"]
        self.workflows_basepath = os.path.join(self.basepath, "workflows")

    def input_path(self, name):
        return os.path.join(self.workflow_path(name), "input")

    def final_path(self, name):
        return os.path.join(self.workflow_path(name), "final")

    def intermediate_path(self, name):
        return os.path.join(self.workflow_path(name), "intermediate")

    def exploratory_path(self, name):
        return os.path.join(self.workflow_path(name), "exploratory")

    def workflow_path(self, name):
        return os.path.join(self.workflows_basepath, name)

    def execution_file_path(self, name):
        return os.path.join(self.workflow_path(name), "run_workflow.sh")

    def configuration_file_path(self, name):
        return os.path.join(self.workflow_path(name), "configuration.py")
