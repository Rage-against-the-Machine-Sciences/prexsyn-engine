import numpy
import numpy.typing
import prexsyn_engine.chemspace
import prexsyn_engine.descriptor
import typing

class MultiThreadedDetokenizer:
    def __init__(self, chemical_space: prexsyn_engine.chemspace.ChemicalSpace, token_def: prexsyn_engine.descriptor.TokenDef = ..., max_outcomes_per_reaction: typing.SupportsInt | typing.SupportsIndex | None = ...) -> None: ...
    def __call__(self, batch_size: typing.SupportsInt | typing.SupportsIndex, tokens: typing.Annotated[numpy.typing.ArrayLike, numpy.int64]) -> list[prexsyn_engine.chemspace.Synthesis]: ...

def detokenize(tokens: typing.Annotated[numpy.typing.ArrayLike, numpy.int64], chemical_space: prexsyn_engine.chemspace.ChemicalSpace, token_def: prexsyn_engine.descriptor.TokenDef = ..., max_outcomes_per_reaction: typing.SupportsInt | typing.SupportsIndex | None = ...) -> prexsyn_engine.chemspace.Synthesis: ...
