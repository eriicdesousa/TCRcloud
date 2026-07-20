"""Unit tests for tcrcloud.assign_colours."""

import pytest

import tcrcloud.assign_colours as assign_colours

# ---------------------------------------------------------------------------
# rgb_to_hex
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "rgb, expected",
    [
        ((0, 0, 0), "#000000"),
        ((255, 255, 255), "#FFFFFF"),
        ((255, 0, 128), "#FF0080"),
        ((1, 2, 3), "#010203"),
    ],
)
def test_rgb_to_hex(rgb, expected):
    assert assign_colours.rgb_to_hex(rgb) == expected


# ---------------------------------------------------------------------------
# is_white_or_grey
# ---------------------------------------------------------------------------


def test_is_white_or_grey_detects_white():
    assert assign_colours.is_white_or_grey((250, 250, 250)) is True


def test_is_white_or_grey_detects_grey():
    assert assign_colours.is_white_or_grey((100, 105, 110)) is True


def test_is_white_or_grey_allows_saturated_colour():
    assert assign_colours.is_white_or_grey((200, 20, 20)) is False


def test_is_white_or_grey_respects_custom_thresholds():
    # (200, 210, 232) has a max-min spread of 32, so it isn't grey under the
    # default threshold (30)...
    assert assign_colours.is_white_or_grey((200, 210, 232)) is False
    # ...but is considered grey with a looser custom threshold.
    assert assign_colours.is_white_or_grey((200, 210, 232), grey_threshold=35) is True


# ---------------------------------------------------------------------------
# assign_glassbey_colours_to_dict
# ---------------------------------------------------------------------------


class _FakeGlasbey:
    """Stand-in for `glasbey.Glasbey` that yields deterministic RGB colours."""

    def __init__(self, no_black=True):
        self.no_black = no_black

    def generate_palette(self, size):
        # Return a list of placeholder "palette entries"; the actual values
        # don't matter since convert_palette_to_rgb below produces the real
        # RGB tuples used by the code under test.
        return list(range(size))

    def convert_palette_to_rgb(self, palette):
        # Produce clearly distinct, saturated RGB tuples - none of which are
        # white/grey - so every generated colour passes the safety filter.
        return [
            ((10 * i) % 256, (50 + 5 * i) % 256, (200 - 3 * i) % 256) for i in palette
        ]


class _AllGreyGlasbey(_FakeGlasbey):
    """A fake Glasbey that only ever produces grey colours, to exercise the
    "couldn't generate enough safe colours" failure path."""

    def convert_palette_to_rgb(self, palette):
        return [(128, 128, 128) for _ in palette]


def test_assign_glassbey_colours_raises_without_glasbey_installed(monkeypatch):
    monkeypatch.setattr(assign_colours, "Glasbey", None)

    with pytest.raises(ImportError, match="glasbey"):
        assign_colours.assign_glassbey_colours_to_dict({"TRAV1-1": None})


def test_assign_glassbey_colours_raises_for_non_dict(monkeypatch):
    monkeypatch.setattr(assign_colours, "Glasbey", _FakeGlasbey)

    with pytest.raises(TypeError, match="target_dict must be a dict"):
        assign_colours.assign_glassbey_colours_to_dict(["not", "a", "dict"])


def test_assign_glassbey_colours_noop_for_empty_dict(monkeypatch):
    monkeypatch.setattr(assign_colours, "Glasbey", _FakeGlasbey)

    result = assign_colours.assign_glassbey_colours_to_dict({})

    assert result == {}


def test_assign_glassbey_colours_overwrite_true_replaces_existing(monkeypatch):
    monkeypatch.setattr(assign_colours, "Glasbey", _FakeGlasbey)

    target = {"TRAV1-1": "#000000", "TRAV1-2": "#111111"}
    result = assign_colours.assign_glassbey_colours_to_dict(target, overwrite=True)

    assert result is target
    assert all(v != "#000000" and v != "#111111" for v in target.values())
    assert all(v.startswith("#") and len(v) == 7 for v in target.values())


def test_assign_glassbey_colours_overwrite_false_keeps_existing(monkeypatch):
    monkeypatch.setattr(assign_colours, "Glasbey", _FakeGlasbey)

    target = {"TRAV1-1": "#ABCDEF", "TRAV1-2": None, "TRAV1-3": ""}
    assign_colours.assign_glassbey_colours_to_dict(target, overwrite=False)

    # Existing truthy value is preserved...
    assert target["TRAV1-1"] == "#ABCDEF"
    # ...but falsy values are filled in with generated colours.
    assert target["TRAV1-2"] not in (None, "")
    assert target["TRAV1-3"] not in (None, "")


def test_assign_glassbey_colours_passes_grey_threshold_through(monkeypatch):
    """Regression test: grey/white filtering must be configurable via a
    parameter rather than requiring callers to monkeypatch a module global."""

    monkeypatch.setattr(assign_colours, "Glasbey", _FakeGlasbey)

    calls = []
    real_is_white_or_grey = assign_colours.is_white_or_grey

    def spy(rgb_tuple, white_thresh=240, grey_threshold=30):
        calls.append(grey_threshold)
        return real_is_white_or_grey(rgb_tuple, white_thresh, grey_threshold)

    monkeypatch.setattr(assign_colours, "is_white_or_grey", spy)

    assign_colours.assign_glassbey_colours_to_dict({"TRAV1-1": None}, grey_threshold=99)

    assert calls and all(threshold == 99 for threshold in calls)


def test_assign_glassbey_colours_raises_runtime_error_when_palette_exhausted(
    monkeypatch,
):
    monkeypatch.setattr(assign_colours, "Glasbey", _AllGreyGlasbey)
    monkeypatch.setattr(assign_colours, "MAX_PALETTE_ATTEMPTS", 2)

    with pytest.raises(RuntimeError, match="after 2 attempts"):
        assign_colours.assign_glassbey_colours_to_dict(
            {"TRAV1-1": None, "TRAV1-2": None}
        )
