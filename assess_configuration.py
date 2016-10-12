from assess.algorithms.incrementaldistancealgorithm import IncrementalDistanceAlgorithm
from assess.algorithms.signatures.signatures import *
from assess.decorators.anomalydecorator import AnomalyDecorator

from assess.algorithms.distances.simpledistance import SimpleDistance
from assess.generators.gnm_importer import GNMCSVEventStreamer, EventStreamer, EventStreamBranchRelabeler


def build_decorator():
    anomaly = AnomalyDecorator(percentage=.15)
    # build decorator chain
    return anomaly

signature_class = ParentChildByNameTopologySignature

def chance_generator(range_values, repetition_values):
    splittings = (item for sublist in map(lambda x,y: [x] * y, range_values, repetition_values) for item in sublist)
    current = 0
    yield current
    for splitting in splittings:
        current += splitting
        yield current

def branch_relabeller(signature, prune_chance, **kwargs):
    csv_event_streamer = GNMCSVEventStreamer(**kwargs)
    event_stream_relabeler = EventStreamBranchRelabeler(
        signature=signature_class(),
        chance=prune_chance,
        streamer=csv_event_streamer,
        seed=1234)
    return EventStreamer(streamer=event_stream_relabeler)


configurations = [{
    "algorithms": [
        lambda **kwargs: IncrementalDistanceAlgorithm(distance=SimpleDistance, **kwargs)
    ], "signatures": [
        lambda: signature_class()
    ], "decorator": build_decorator,
    "event_streamer": (
        lambda **kwargs: branch_relabeller(
            signature=signature_class(), 
            prune_chance=chance, 
            **kwargs) for chance in chance_generator([.01, .05, .1], [20, 4, 6]) for _ in range(2)
    )
}]

