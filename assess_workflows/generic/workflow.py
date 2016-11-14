import os
import sys
import shlex
import subprocess

from assess_workflows.generic.structure import Structure


class Task(object):
    def __init__(self, cli_path=None, cmd=None, save=True, use_input=False, env={}, **kwargs):
        self.cli_path = cli_path
        self.cmd = cmd
        self.save = save
        self.use_input = use_input
        self.others = kwargs if kwargs else {}
        self.env = env

    def build_subprocess(self, **kwargs):
        command = "python %s --step %s %s %s %s %s" % (
            self.cli_path or Structure.workflow_cli(),
            kwargs.get("index", 0),
            ("--save" if self.save else ""),
            ("--use_input" if self.use_input else ""),
            self.cmd,
            " ".join(
                ["--%s %s" % (key, value) for (key, values) in self.others.items() for
                    value in (values if isinstance(values, list) else [values])]
            )
        )
        return shlex.split(command)


class Workflow(object):
    def __init__(self):
        self._tasks = []

    def add_task(self, cli_path=None, cmd=None, save=None, use_input=None, **kwargs):
        self._tasks.append(Task(cli_path=cli_path, cmd=cmd, save=save, use_input=use_input, **kwargs))

    def finalise(self):
        """
        Convenience function that encapsulates the copy process to finalise the last inermediate
        file that was created.
        """
        current_step = len(self._tasks) + 1
        self.add_task(
            cli_path=Structure.generic_cli(),
            cmd="finalise",
            from_step=current_step - 1
        )

    def prepare_intermediate_as_input(self):
        """
        Convenience function that encapsulates the copy process from intermediate results from
        last step to prepare for next step.
        """
        current_step = len(self._tasks) + 1
        self.add_task(
            cli_path=Structure.generic_cli(),
            cmd="intermediate_as_input",
            from_step=current_step - 1,
            to_step=current_step + 1
        )

    def execute(self, environment_variables=None, start=0, end=None):
        """
        interval [start, end[

        :param environment_variables:
        :param start: Start from index start
        :param end: Consider steps up to end (exclusively)
        :return:
        """
        print("Using environment %s" % environment_variables)
        environment = os.environ
        environment.update(environment_variables)
        environment["PYTHONPATH"] = os.path.abspath(sys.path[0]) + \
            (':' + os.environ['PYTHONPATH'] if 'PYTHONPATH' in os.environ else '')
        for index, task in enumerate(self._tasks[start:end], start=start):
            print("starting workflow %d of %d" % (index+1, len(self._tasks)))
            current_environment = environment.copy()
            current_environment.update(task.env)
            subprocess.check_call(task.build_subprocess(index=index+1), env=current_environment)
