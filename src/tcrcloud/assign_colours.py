"""Maintainer-only utility to (re)generate V-gene colour palettes.

This script is not part of the `TCRcloud` CLI (it isn't wired into any
subparser in `tcrcloud.TCRcloud`) - it is a standalone tool maintainers run
by hand to (re)generate the hard-coded V-gene colour dictionaries that live
in `tcrcloud.colours` and the per-species colour modules (e.g.
`tcrcloud.colours_mus_musculus`), using the Glasbey algorithm to produce a
set of visually distinct colours.

Because it's a dev-only tool, the `glasbey` dependency it needs is optional
(install with ``pip install glasbey`` or ``pip install -e ".[dev]"``) and is
imported lazily below so the rest of `tcrcloud` doesn't require it.
"""

import importlib

try:
    from glasbey import Glasbey
except ImportError:
    # Deferred: only required if `assign_glassbey_colours_to_dict` is actually
    # called. This keeps the module importable (e.g. for unit tests) even
    # when the optional `glasbey` dependency isn't installed.
    Glasbey = None

# Safety cap on how many times we'll ask Glasbey for a bigger palette before
# giving up. Without this, a pathological case (e.g. Glasbey consistently
# generating mostly white/grey colours) would loop forever.
MAX_PALETTE_ATTEMPTS = 10


def rgb_to_hex(rgb_tuple):
    """Convert an RGB tuple (r, g, b) with int values 0-255 to hex string."""
    r, g, b = rgb_tuple
    return f"#{r:02X}{g:02X}{b:02X}"


def is_white_or_grey(rgb_tuple, white_thresh=240, grey_threshold=30):
    """Return True for white or grey-like colours on white paper."""
    r, g, b = rgb_tuple
    if r >= white_thresh and g >= white_thresh and b >= white_thresh:
        return True
    max_c = max(r, g, b)
    min_c = min(r, g, b)
    if max_c - min_c <= grey_threshold:
        return True
    return False


def assign_glassbey_colours_to_dict(
    target_dict,
    no_black=True,
    overwrite=True,
    white_thresh=240,
    grey_threshold=30,
):
    """Assign a Glasbey-generated palette to a dict's keys as hex colours.

    Args:
        target_dict: Mapping (e.g. `tcrcloud.colours.TRAV`) whose keys will
            each be assigned a generated hex colour. Mutated in place.
        no_black: Passed through to `Glasbey`; excludes black/near-black
            colours from the generated palette.
        overwrite: If True (default), replace every key's existing colour.
            If False, only fill in colours for keys that don't already have
            a truthy value (use this to add new genes without disturbing
            an already-curated palette).
        white_thresh: Passed through to `is_white_or_grey` when filtering
            out near-white candidate colours.
        grey_threshold: Passed through to `is_white_or_grey` when filtering
            out near-grey candidate colours.

    Raises:
        ImportError: if the optional `glasbey` dependency isn't installed.
        RuntimeError: if a sufficiently large, distinct, non-white/grey
            palette could not be generated after `MAX_PALETTE_ATTEMPTS`
            attempts.
    """
    if Glasbey is None:
        raise ImportError(
            "The 'glasbey' package is required to generate colours. "
            "Install it with: pip install glasbey"
        )

    if not isinstance(target_dict, dict):
        raise TypeError("target_dict must be a dict")

    keys = list(target_dict.keys())
    n = len(keys)
    if n == 0:
        return target_dict

    gb = Glasbey(no_black=no_black)
    safe_colours = []
    batch_size = max(2 * n, 64)

    attempts = 0
    while len(safe_colours) < n:
        attempts += 1
        if attempts > MAX_PALETTE_ATTEMPTS:
            raise RuntimeError(
                f"Could not generate {n} sufficiently distinct (non-white/grey) "
                f"colours after {MAX_PALETTE_ATTEMPTS} attempts."
            )

        palette = gb.generate_palette(size=len(safe_colours) + batch_size)
        rgb_palette = gb.convert_palette_to_rgb(palette)

        safe_colours = [
            c
            for c in rgb_palette
            if not is_white_or_grey(c, white_thresh, grey_threshold)
        ]

        if len(safe_colours) < n:
            batch_size *= 2

    for key, colour_rgb in zip(keys, safe_colours[:n]):
        hex_colour = rgb_to_hex(colour_rgb)
        if overwrite or not target_dict.get(key):
            target_dict[key] = hex_colour

    return target_dict


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="Assign Glasbey colours to V-gene dictionaries"
    )
    parser.add_argument(
        "--source",
        default="tcrcloud.colours",
        help="Python module path containing dicts",
    )
    parser.add_argument(
        "--output", default="trav_colours_glassbey.json", help="Output JSON file"
    )
    parser.add_argument(
        "--grey-threshold", type=int, default=30, help="Threshold for grey filtering"
    )
    parser.add_argument(
        "--exclude-black",
        dest="no_black",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude black and near-black colours from the generated palette "
        "(pass --no-exclude-black to allow them). Default: enabled.",
    )
    parser.add_argument(
        "--overwrite",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Overwrite every key's existing colour. Pass --no-overwrite to "
        "only fill in colours for keys that don't already have one. "
        "Default: enabled.",
    )
    args = parser.parse_args()

    module = importlib.import_module(args.source)
    dicts = {
        "TRAV": getattr(module, "TRAV"),
        "TRBV": getattr(module, "TRBV"),
        "TRGV": getattr(module, "TRGV"),
        "TRDV": getattr(module, "TRDV"),
    }

    if hasattr(module, "IGHV"):
        dicts["IGHV"] = getattr(module, "IGHV")
    if hasattr(module, "IGKV"):
        dicts["IGKV"] = getattr(module, "IGKV")
    if hasattr(module, "IGLV"):
        dicts["IGLV"] = getattr(module, "IGLV")

    print(f"Assigning Glasbey colours to dicts from {args.source}...")
    if args.overwrite:
        print(
            "Warning: --overwrite is enabled (default); existing colours will be "
            "replaced. Pass --no-overwrite to only fill in colours for keys that "
            "don't already have one."
        )

    for name, d in dicts.items():
        assign_glassbey_colours_to_dict(
            d,
            no_black=args.no_black,
            overwrite=args.overwrite,
            grey_threshold=args.grey_threshold,
        )
        print(f"  {name}: {len(d)} keys (generated)")

    if "TRAV" in dicts:
        print("Sample TRAV values (first 20):")
        for i, (k, v) in enumerate(dicts["TRAV"].items()):
            if i >= 20:
                break
            print(f"{k}: {v}")

    out_path = args.output
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(dicts, f, indent=2)

    print(f"Wrote updated colours to {out_path}")
