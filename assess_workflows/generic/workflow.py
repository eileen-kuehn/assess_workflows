import os
import sys
import shlex
import subprocess

from assess_workflows.generic.structure import Structure


class Task(object):
    def __init__(self, cli_path=None, cmd=None, save=True, use_input=False, env={}, name=None, **kwargs):
        self.cli_path = cli_path
        self.cmd = cmd
        self.save = save
        self.use_input = use_input
        self.others = kwargs if kwargs else {}
        self.env = env
        self.name = name

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
        self._names = {}
        self._tasks = []

    def add_task(self, cli_path=None, cmd=None, save=None, use_input=None, name=None, **kwargs):
        if name is None:
            name = len(self._tasks) + 1
        if name not in self._names:
            self._names[name] = len(self._tasks) + 1
            self._tasks.append(
                Task(cli_path=cli_path, cmd=cmd, save=save, use_input=use_input, name=name, **kwargs))
        else:
            raise ValueError

    def finalise(self, file_type="json", reference=None, name=None):
        """
        Convenience function that encapsulates the copy process to finalise the last inermediate
        file that was created.
        """
        current_step = len(self._tasks) + 1
        if reference is not None:
            target = self._names.get(reference, None)
        else:
            target = None
        self.add_task(
            cli_path=Structure.generic_cli(),
            cmd="finalise",
            from_step=target or (current_step - 1),
            file_type=file_type,
            name="Finalising %s" % (name or reference or current_step)
        )

    def prepare_intermediate_as_input(self, from_step=None, reference=None, name=None,
                                      file_type="json"):
        """
        Convenience function that encapsulates the copy process from intermediate results from
        last step to prepare for next step.

        :param from_step: The step to copy data from
        """
        current_step = len(self._tasks) + 1
        if reference is not None:
            target = self._names.get(reference, None)
        else:
            target = None
        from_step = target or from_step or (current_step - 1)
        self.add_task(
            cli_path=Structure.generic_cli(),
            cmd="intermediate_as_input",
            from_step=target or from_step or (current_step - 1),
            to_step=current_step + 1,
            name="Preparing %s for %s" % ((name or reference or from_step), current_step + 1),
            file_type=file_type
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
            print("starting task %d of %d (%s)" % (index+1, len(self._tasks), task.name))
            current_environment = environment.copy()
            current_environment.update(task.env)
            subprocess.check_call(task.build_subprocess(index=index+1), env=current_environment)
