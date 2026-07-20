"""Unit tests for tcrcloud.colours species-palette lookups."""

import pytest

import tcrcloud.colours as colours

# ---------------------------------------------------------------------------
# _normalize_species
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        (None, "homo_sapiens"),
        ("", "homo_sapiens"),
        ("Homo_Sapiens", "homo_sapiens"),
        ("macaca fascicularis", "macaca_fascicularis"),
        ("Macaca-Fascicularis", "macaca_fascicularis"),
        ("  mus_musculus  ", "mus_musculus"),
    ],
)
def test_normalize_species(raw, expected):
    assert colours._normalize_species(raw) == expected


# ---------------------------------------------------------------------------
# get_vgene_palette / get_color_from_vcall
# ---------------------------------------------------------------------------


def test_get_vgene_palette_returns_species_specific_palette_for_macaca_fascicularis():
    # Regression test: SPECIES_VGENES used to key macaca fascicularis with a
    # space ("macaca fascicularis"), which never matched the normalized,
    # underscore-separated species name produced by `_normalize_species`
    # (and advertised by the CLI as "macaca_fascicularis"). That mismatch
    # silently fell back to the homo_sapiens palette.
    palette = colours.get_vgene_palette("TRBV", species="macaca_fascicularis")
    assert palette == colours.MACACA_FASCICULARIS["TRBV"]
    assert palette != colours.HOMO_SAPIENS["TRBV"]


@pytest.mark.parametrize(
    "species_input",
    ["macaca_fascicularis", "macaca fascicularis", "Macaca-Fascicularis"],
)
def test_get_vgene_palette_accepts_species_name_variants(species_input):
    assert (
        colours.get_vgene_palette("TRBV", species=species_input)
        == colours.MACACA_FASCICULARIS["TRBV"]
    )


def test_get_vgene_palette_falls_back_to_homo_sapiens_for_unknown_species():
    palette = colours.get_vgene_palette("TRBV", species="unknown_species")
    assert palette == colours.HOMO_SAPIENS["TRBV"]


def test_get_color_from_vcall_uses_macaca_fascicularis_palette():
    vcall, expected_color = next(iter(colours.MACACA_FASCICULARIS["TRBV"].items()))
    assert (
        colours.get_color_from_vcall(vcall, species="macaca_fascicularis")
        == expected_color
    )


def test_get_color_from_vcall_returns_default_for_unknown_vcall():
    assert (
        colours.get_color_from_vcall(
            "TRBVUNKNOWN", species="homo_sapiens", default="grey"
        )
        == "grey"
    )


def test_get_color_from_vcall_returns_default_for_empty_or_short_vcall():
    assert colours.get_color_from_vcall("", default="grey") == "grey"
    assert colours.get_color_from_vcall("TR", default="grey") == "grey"
