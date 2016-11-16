from __future__ import print_function
import json
import datetime
import subprocess
import multiprocessing

from evenmoreutils.files import smart_open


def determine_version(path):
    try:
        version = subprocess.check_output(["git", "describe"], cwd=path).strip()
    except subprocess.CalledProcessError:
        try:
            version = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=path).strip()
        except subprocess.CalledProcessError:
            version = "unversioned"
    return version


def output_results(ctx, results=None, version=None, source=None, variant=None, file_type=None):
    """

    :param ctx:
    :param results:
    :param version:
    :param source:
    :param variant: variant to use, when a workflow supports more than one output file
    :return:
    """
    format_json = ctx.obj.get("json", False)
    file_name = None if not ctx.obj.get("save", False) else \
        ctx.obj.get("structure").intermediate_file_path(variant=variant, file_type=file_type)
    with smart_open(filename=file_name) as output_channel:
        if format_json:
            dump = {
                "meta": {
                    "date": "%s" % datetime.datetime.now(),
                    "source": "%s" % source if source else "unknown",
                    "version": "%s" % version if version else "unknown",
                },
                "data": results
            }
            print(json.dumps(dump, indent=2), file=output_channel)
        else:
            print("# date: %s" % datetime.datetime.now(), file=output_channel)
            print("# source: %s" % source if source else "unknown", file=output_channel)
            print("# version: %s" % version if version else "unkown", file=output_channel)
            print(results, file=output_channel)


def do_multicore(count=1, target=None, data=None):
    pool = multiprocessing.Pool(processes=count)
    result_list = pool.map(target, data)
    pool.close()
    pool.join()
    return [result for result in result_list if result is not None]
