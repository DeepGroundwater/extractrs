"""channel_attrs — pure-function API for the Phase-A channel-corridor pipeline.

Import these directly for use in notebooks, scripts, or downstream pipelines:

    from pipelines.channel_attrs import (
        build_corridors,
        build_scaled_corridors,
        order_bankfull_width_m,
        scaled_half_width_m,
        sample_along_lines,
        weighted_transfer,
        assemble,
    )

The ``main()`` entry-points in each submodule handle I/O; the functions
exported here are pure (data in → data out, no filesystem side-effects) and
can be exercised with in-memory data in unit tests or notebooks.
"""

from .corridors import (
    build_corridors,
    build_scaled_corridors,
    order_bankfull_width_m,
    scaled_half_width_m,
)
from .wtd_sample import sample_along_lines
from .transfer import weighted_transfer
from .assemble import assemble

__all__ = [
    "build_corridors",
    "build_scaled_corridors",
    "order_bankfull_width_m",
    "scaled_half_width_m",
    "sample_along_lines",
    "weighted_transfer",
    "assemble",
]
