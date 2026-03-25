from glasbey import Glasbey
import importlib


def rgb_to_hex(rgb_tuple):
    """Convert an RGB tuple (r, g, b) with int values 0-255 to hex string."""
    r, g, b = rgb_tuple
    return f"#{r:02X}{g:02X}{b:02X}"


def is_white_or_grey(rgb_tuple, white_thresh=240, grey_threshold=30):
    """Return True for white or grey-like colors on white paper."""
    r, g, b = rgb_tuple
    if r >= white_thresh and g >= white_thresh and b >= white_thresh:
        return True
    max_c = max(r, g, b)
    min_c = min(r, g, b)
    if max_c - min_c <= grey_threshold:
        return True
    return False


def assign_glassbey_colors_to_dict(target_dict, no_black=True, overwrite=True):
    """Assign a Glasbey-generated palette to a dict's keys as hex colors."""
    if not isinstance(target_dict, dict):
        raise TypeError("target_dict must be a dict")

    keys = list(target_dict.keys())
    n = len(keys)
    if n == 0:
        return target_dict

    gb = Glasbey(no_black=no_black)
    safe_colors = []
    batch_size = max(2 * n, 64)

    while len(safe_colors) < n:
        palette = gb.generate_palette(size=len(safe_colors) + batch_size)
        rgb_palette = gb.convert_palette_to_rgb(palette)

        safe_colors = [c for c in rgb_palette if not is_white_or_grey(c)]

        if len(safe_colors) < n:
            batch_size *= 2

    for key, color_rgb in zip(keys, safe_colors[:n]):
        hex_color = rgb_to_hex(color_rgb)
        if overwrite or not target_dict.get(key):
            target_dict[key] = hex_color

    return target_dict


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description="Assign Glasbey colors to V-gene dictionaries")
    parser.add_argument("--source", default="tcrcloud.colours", help="Python module path containing dicts")
    parser.add_argument("--output", default="trav_colors_glassbey.json", help="Output JSON file")
    parser.add_argument("--grey-threshold", type=int, default=30, help="Threshold for grey filtering")
    parser.add_argument("--no-black", action="store_true", default=True, help="Exclude black and near-black")
    parser.add_argument("--overwrite", action="store_true", default=True, help="Overwrite existing values")
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

    print(f"Assigning Glasbey colors to dicts from {args.source}...")

    # set grey threshold runtime override
    def override_is_white_or_grey(rgb_tuple):
        r, g, b = rgb_tuple
        if r >= 240 and g >= 240 and b >= 240:
            return True
        return max(r, g, b) - min(r, g, b) <= args.grey_threshold

    # monkey patch helper for this run
    is_white_or_grey = override_is_white_or_grey

    for name, d in dicts.items():
        assign_glassbey_colors_to_dict(d, no_black=args.no_black, overwrite=args.overwrite)
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

    print(f"Wrote updated colors to {out_path}")
