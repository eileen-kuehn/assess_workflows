import contextlib
import json
import datetime
import os
import subprocess
import multiprocessing
import sys
import pandas as pd
import numpy as np


@contextlib.contextmanager
def smart_open(filename=None):
    if filename.endswith(".h5"):
        handle = pd.HDFStore(filename)
    elif filename and filename != '-':
        handle = open(filename, 'w')
    else:
        handle = sys.stdout

    try:
        yield handle
    finally:
        if handle is not sys.stdout:
            handle.close()


def determine_version(path):
    try:
        version = subprocess.check_output(["git", "describe"], cwd=path).strip()
    except subprocess.CalledProcessError:
        try:
            version = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=path).strip()
        except subprocess.CalledProcessError:
            version = "unversioned"
    return version


def output_results(ctx, results=None, version=None, source=None, variant=None, file_type=None,
                   comment_sign="#"):
    """

    :param ctx:
    :param results:
    :param version:
    :param source:
    :param variant: variant to use, when a workflow supports more than one output file
    :return:
    """
    format_json = ctx.obj.get("json", False)
    format_hdf = ctx.obj.get("hdf", False)
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
        elif format_hdf:
            tree_references = [os.path.basename(name).split("-")[-1].split(
                ".")[0] for name in results.get("files", [])]
            # check for representatives
            prototype_references = [os.path.basename(name).split("-")[-1].split(
                ".")[0] for name in results.get("prototypes", [])]
            for index, result in enumerate(results.get("results")):
                data = []
                if not isinstance(result, dict):
                    result = result[0]
                # assuming one matrix decorator per result
                matrix = result.get("decorator").get("matrix")
                for row in matrix:
                    data.append(row[0])
                output_channel.put("df_%d" % index, pd.DataFrame(
                    np.array(data),
                    columns=prototype_references if prototype_references else tree_references,
                    index=tree_references
                ))
                output_channel.get_storer("df_%d" % index).attrs.meta = {
                    "algorithm": result.get("algorithm", None),
                    "signature": result.get("signature", None),
                    "event_streamer": result.get("event_streamer", None),
                    "date": "%s" % datetime.datetime.now(),
                    "source": "%s" % source if source else "unknown",
                    "version": "%s" % version if version else "unknown"
                }
        else:
            print("%s date: %s" % (comment_sign, datetime.datetime.now()), file=output_channel)
            print("%s source: %s" % (comment_sign, source if source else "unknown"), file=output_channel)
            print("%s version: %s" % (comment_sign, version if version else "unkown"), file=output_channel)
            print(results, file=output_channel)


def output_r_data(ctx, filename, **kwargs):
    from rpy2 import robjects
    for key, value in kwargs.items():
        robjects.globalenv[key] = value
    robjects.r.save(*list(kwargs.keys()), file=filename)


def do_multicore(count=1, target=None, data=None):
    pool = multiprocessing.Pool(processes=count)
    result_list = pool.map(target, data)
    results = [result for result in result_list if result is not None]
    pool.close()
    pool.join()
    return results
