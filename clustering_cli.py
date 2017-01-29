import json
import os

import time
import click
import logging

import assess_workflows
from assess_workflows.utils.utils import output_results, determine_version
from dengraph.graphs import graph_io
from utility.exceptions import ExceptionFrame
from utility.report import LVL

from dengraph.quality.silhouette import silhouette_score
from dengraph.quality.calinski_harabasz import calinski_harabasz_score
from dengraph.quality.davies_bouldin import davies_bouldin_score
from dengraph.graphs.adjacency_graph import AdjacencyGraph
from dengraph.dengraph import DenGraphIO

from assess.clustering.clustering import Clustering
from assess.clustering.clusterdistance import ClusterDistance
from assess.generators.gnm_importer import CSVTreeBuilder
from assess.events.events import ProcessStartEvent, ProcessExitEvent, TrafficEvent

from assess_workflows.generic.structure import Structure


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
    ctx.obj["json"] = True
    ctx.obj["save"] = save
    ctx.obj["use_input"] = use_input
    configdict = {}
    execfile(configuration or ctx.obj["structure"].configuration_file_path(), configdict)
    ctx.obj["configurations"] = configdict["configurations"][:]


@click.command()
@click.option("--eta", "eta", type=int, default=5, multiple=True)
@click.option("--epsilon", "epsilon", type=float, default=.1, multiple=True)
@click.pass_context
def perform_precalculated_clustering(ctx, eta, epsilon):
    results = {}

    if ctx.obj.get("use_input", False):
        configuration = ctx.obj.get("configurations", None)[0]
        signature = configuration.get("signatures", [None])[0]
        distance = configuration.get("distances", [None])[0]
        # FIXME: I might also need cache statistics
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path(file_type="csv")  # expecting csv file
        tree_builder = CSVTreeBuilder()

        def header_to_cache(file_path):
            tree = tree_builder.build(file_path)
            # FIXME: does this work as expected?
            tree_index = tree.to_index(
                signature=signature,
                start_support=distance.supported.get(ProcessStartEvent, False),
                exit_support=distance.supported.get(ProcessExitEvent, False),
                traffic_support=distance.supported.get(TrafficEvent, False)
            )
            tree_index.key = file_path
            return tree_index

        with open(file_path) as csv_file:
            # load the graph from precalculated csv distance values
            graph = graph_io.csv_graph_reader(
                (ln for ln in csv_file if ln[0] != '#' and ln != '\n'),
                nodes_header=header_to_cache,
                symmetric=True)

        for single_eta in eta:
            for single_epsilon in epsilon:
                start = time.time()
                clustering = DenGraphIO(
                    base_graph=graph,
                    cluster_distance=single_epsilon,
                    core_neighbours=single_eta
                )
                end = time.time()
                cluster_distance = ClusterDistance(distance=distance)
                clustering.graph.distance = cluster_distance
                print("---> performed clustering with eta %s and epsilon %s in %s" % (eta, epsilon, end - start))
                current_result = results.setdefault("data", []).append({})
                current_result.setdefault("meta", {})["algorithm"] = clustering.__class__.__name__
                current_result.setdefault("meta", {})["eta"] = eta
                current_result.setdefault("meta", {})["epsilon"] = epsilon
                current_result["duration"] = end - start
                for cluster_idx, cluster in enumerate(clustering):
                    current_result.setdefault("clusters", []).append([node.key for node in cluster])  # TODO: determine CR
                    print("[cluster %s] %s" % (cluster_idx, len(cluster)))
                print("[noise] %s" % len(clustering.noise))
                for noise in clustering.noise:
                    current_result.setdefault("noise", []).append(noise.key)
                # for score in [silhouette_score, calinski_harabasz_score, davies_bouldin_score]:
                for score in [silhouette_score]:
                    try:
                        the_score = score(clustering.clusters, clustering.graph)
                    except ValueError:
                        the_score = None
                    current_result.setdefault("scores", {})[score.__name__] = the_score
                    print("Got a %s of %s" % (score.__name__, the_score))

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "perform_precalculated_clustering")
    )


@click.command()
@click.option("--eta", "eta", type=int, default=5)
@click.option("--epsilon", "epsilon", type=float, default=.1)
@click.pass_context
def perform_clustering(ctx, eta, epsilon):
    results = {}

    if ctx.obj.get("use_input", False):
        configuration = ctx.obj.get("configurations", None)[0]
        signature = configuration.get("signatures", [None])[0]
        distance = configuration.get("distances", [None])[0]
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path()

        tree_builder = CSVTreeBuilder()
        clustering = Clustering(distance=distance, cluster_distance=epsilon, core_neighbours=eta)
        with open(file_path, "r") as input_file:
            input_data = json.load(input_file).get("data", None)

            for sample in input_data.get("samples", []):
                tree = tree_builder.build(sample[0])
                # convert tree to index
                tree_index = tree.to_index(
                    signature=signature,
                    start_support=distance.supported.get(ProcessStartEvent, False),
                    exit_support=distance.supported.get(ProcessExitEvent, False),
                    traffic_support=distance.supported.get(TrafficEvent, False)
                )
                clustering[sample[0]] = tree_index
        print("---> performed clustering with eta %s and epsilon %s" % (eta, epsilon))
        results.setdefault("meta", {})["algorithm"] = clustering.clusterer.__class__.__name__
        results.setdefault("meta", {})["eta"] = eta
        results.setdefault("meta", {})["epsilon"] = epsilon
        for cluster in clustering:
            results.setdefault("clusters", []).append([node.key for node in cluster])  # TODO: determine CR
        for noise in clustering.clusterer.noise:
            results.setdefault("noise", []).append(noise.key)
        for score in [silhouette_score, calinski_harabasz_score, davies_bouldin_score]:
            try:
                the_score = score(clustering.clusterer.clusters, clustering.clusterer.graph)
            except ValueError:
                the_score = None
            results.setdefault("scores", {})[score.__name__] = the_score

    output_results(
        ctx=ctx,
        results=results,
        version=determine_version(os.path.dirname(assess_workflows.__file__)),
        source="%s (%s)" % (__file__, "perform_clustering")
    )

cli.add_command(perform_clustering)
cli.add_command(perform_precalculated_clustering)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
