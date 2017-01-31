import json
import os

import time
import click
import logging

import assess_workflows
from assess.exceptions.exceptions import EventNotSupportedException
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
from assess.generators.gnm_importer import CSVTreeBuilder, GNMCSVEventStreamer
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
        distance = configuration.get("distances", [None])[0]
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path(file_type="csv")  # expecting csv file

        graph = _create_graph(ctx, file_path)
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
                print("---> performed clustering with eta %s and epsilon %s in %s" % (single_eta, single_epsilon, end - start))
                results.setdefault("results", []).append({})
                current_result = results["results"][-1]
                current_result.setdefault("meta", {})["algorithm"] = clustering.__class__.__name__
                current_result.setdefault("meta", {})["eta"] = single_eta
                current_result.setdefault("meta", {})["epsilon"] = single_epsilon
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


@click.command()
@click.option("--eta", "eta", type=int, default=5)
@click.option("--epsilon", "epsilon", type=float, default=.1)
@click.pass_context
def perform_classification(ctx, eta, epsilon):
    """
    Method performs a classification. Before the actual classification can be tested, a clustering
    is applied. Those clusters are following used for classification.

    :param ctx:
    :return:
    """
    if ctx.obj.get("use_input", False):
        configuration = ctx.obj.get("configurations", None)[0]
        distance_cls = configuration.get("distances", [None])[0]
        structure = ctx.obj.get("structure", None)
        file_path = structure.input_file_path(file_type="csv")  # expecting csv file

        graph = _create_graph(ctx, file_path)
        clustering = DenGraphIO(
            base_graph=graph,
            cluster_distance=epsilon,
            core_neighbours=eta
        )
        cluster_distance = ClusterDistance(distance=distance_cls())
        clustering.graph.distance = cluster_distance
        # calculate CRs from clusters
        prototype_caches = []
        cluster_names = []
        for cluster in clustering:
            cluster_names.append(cluster[0].key)
            prototype_caches.append(cluster_distance.mean(list(cluster)), prototype=cluster_names[-1])

        results = []
        decorator_def = configuration.get("decorator", None)
        for algorithm_def in configuration.get("algorithms", []):
            for signature_def in configuration.get("signatures", []):
                for event_streamer in configuration.get("event_streamer", [GNMCSVEventStreamer]):
                    signature = signature_def()
                    algorithm = algorithm_def(signature=signature)
                    algorithm.cluster_representatives(
                        signature_prototypes=prototype_caches, prototypes=cluster_names)
                    decorator = decorator_def()
                    decorator.wrap_algorithm(algorithm=algorithm)
                    # starting a new tree not to mix former results with current
                    for node in clustering.graph:
                        tree = event_streamer(csv_path=node.key)
                        algorithm.start_tree()
                        for event in tree.event_iter():
                            try:
                                algorithm.add_event(event)
                            except EventNotSupportedException:
                                pass
                        algorithm.finish_tree()
                    if decorator:
                        results.append({
                            "algorithm": "%s" % algorithm,
                            "signature": "%s" % signature,
                            "event_streamer": "%s" % tree if tree is not None
                            else event_streamer(csv_path=None),
                            "decorator": decorator.descriptive_data()
                        })
        output_results(
            ctx=ctx,
            results=results,
            version=determine_version(os.path.dirname(assess_workflows.__file__)),
            source="%s (%s)" % (__file__, "perform_classification")
        )


def _create_graph(ctx, file_path):
    configuration = ctx.obj.get("configurations", None)[0]
    signature = configuration.get("signatures", [None])[0]
    distance_builder = configuration.get("distances", [None])[0]
    statistics_cls = configuration.get("statistics", [None])[0]
    tree_builder = CSVTreeBuilder()
    distance = distance_builder()

    def header_to_cache(tree_path):
        tree = tree_builder.build(tree_path)
        tree_index = tree.to_index(
            signature=signature(),
            start_support=distance.supported.get(ProcessStartEvent, False),
            exit_support=distance.supported.get(ProcessExitEvent, False),
            traffic_support=distance.supported.get(TrafficEvent, False),
            statistics_cls=statistics_cls
        )
        tree_index.key = tree_path
        return tree_index

    with open(file_path) as csv_file:
        # load the graph from precalculated csv distance values
        graph = graph_io.csv_graph_reader(
            (ln for ln in csv_file if ln[0] != '#' and ln != '\n'),
            nodes_header=header_to_cache,
            symmetric=True
        )
        return graph

cli.add_command(perform_clustering)
cli.add_command(perform_precalculated_clustering)
cli.add_command(perform_classification)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix='DISS')
