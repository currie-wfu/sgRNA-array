"""Smoke tests for the Flask UI: routes return expected status codes and content."""

from __future__ import annotations

import pytest

from webapp.app import create_app


@pytest.fixture()
def app():
    a = create_app()
    a.config.update(TESTING=True)
    return a


@pytest.fixture()
def client(app):
    return app.test_client()


def test_index_renders(client):
    response = client.get("/")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert "sgRNA Array Designer" in body
    assert "Array size" in body
    assert 'name="array_size"' in body


def test_index_with_size_param_preselects(client):
    response = client.get("/?size=8")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    # The 8-sgRNA option should be selected; check via the rendered HTML.
    assert '<option value="8" selected' in body
    # There should be 8 rows of crRNA inputs.
    assert body.count('name="crrna_') == 8


def test_index_with_invalid_size_falls_back_to_default(client):
    response = client.get("/?size=99")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert body.count('name="crrna_') == 4


def test_design_post_valid_array_returns_results(client):
    form_data = {
        "array_size": "4",
        "gene_label_1": "Shh",
        "crrna_1": "GGTGACGCGTGTGTACCTGG",
        "gene_label_2": "Pax6",
        "crrna_2": "ACGTACGTACGTACGTAAAA",
        "gene_label_3": "Tbx5",
        "crrna_3": "AACCAATTAACCAATTGCAT",
        "gene_label_4": "Sox2",
        "crrna_4": "GCATGCATGCATGCATAACA",
    }
    response = client.post("/design", data=form_data)
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert "Array built" in body
    assert "Download FASTA" in body
    assert "Shh" in body
    assert "ACGG" in body  # B5
    assert "GAGC" in body  # B3
    assert "ACCA" in body  # stitcher overhang for stop-at-4


def test_design_post_invalid_crrna_re_renders_form_with_errors(client):
    form_data = {
        "array_size": "4",
        "gene_label_1": "Bad",
        "crrna_1": "TOOSHORT",  # invalid length
        "gene_label_2": "Ok",
        "crrna_2": "ACGTACGTACGTACGTAAAA",
        "gene_label_3": "Ok",
        "crrna_3": "AACCAATTAACCAATTGCAT",
        "gene_label_4": "Ok",
        "crrna_4": "GCATGCATGCATGCATAACA",
    }
    response = client.post("/design", data=form_data)
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert "Array built" not in body
    assert "must be exactly 20 nt" in body or "exactly 20" in body
    assert "row-error" in body  # the invalid row gets the error styling
    # Form values should be preserved so the user doesn't have to retype.
    assert "Bad" in body
    assert "ACGTACGTACGTACGTAAAA" in body


def test_design_rejects_invalid_array_size(client):
    response = client.post(
        "/design",
        data={"array_size": "7", "crrna_1": "AAAAAAAAAAAAAAAAAAAA"},
    )
    assert response.status_code == 400


def test_export_fasta_for_built_array(client):
    form_data = {
        "array_size": "4",
        "gene_label_1": "Shh",
        "crrna_1": "GGTGACGCGTGTGTACCTGG",
        "gene_label_2": "Pax6",
        "crrna_2": "ACGTACGTACGTACGTAAAA",
        "gene_label_3": "Tbx5",
        "crrna_3": "AACCAATTAACCAATTGCAT",
        "gene_label_4": "Sox2",
        "crrna_4": "GCATGCATGCATGCATAACA",
    }
    design_response = client.post("/design", data=form_data)
    body = design_response.get_data(as_text=True)
    # Extract a design_id from the download links in the results page.
    import re

    match = re.search(r"/export/([\w\-_]+)/fasta", body)
    assert match is not None, "design_id not found in results page"
    design_id = match.group(1)

    fasta_response = client.get(f"/export/{design_id}/fasta")
    assert fasta_response.status_code == 200
    text = fasta_response.get_data(as_text=True)
    assert text.startswith(">")
    assert "pos1_Shh" in text


def test_export_unknown_design_id_returns_404(client):
    response = client.get("/export/nonexistent/fasta")
    assert response.status_code == 404


def test_export_unknown_format_returns_404(client):
    form_data = {
        "array_size": "4",
        "gene_label_1": "Shh",
        "crrna_1": "GGTGACGCGTGTGTACCTGG",
        "gene_label_2": "Pax6",
        "crrna_2": "ACGTACGTACGTACGTAAAA",
        "gene_label_3": "Tbx5",
        "crrna_3": "AACCAATTAACCAATTGCAT",
        "gene_label_4": "Sox2",
        "crrna_4": "GCATGCATGCATGCATAACA",
    }
    design_response = client.post("/design", data=form_data)
    body = design_response.get_data(as_text=True)
    import re

    design_id = re.search(r"/export/([\w\-_]+)/fasta", body).group(1)
    response = client.get(f"/export/{design_id}/bogus")
    assert response.status_code == 404


def test_assembly_map_svg_renders(client):
    """The assembly map SVG should be embedded in the results page."""
    form_data = {
        "array_size": "4",
        "gene_label_1": "Shh",
        "crrna_1": "GGTGACGCGTGTGTACCTGG",
        "gene_label_2": "Pax6",
        "crrna_2": "ACGTACGTACGTACGTAAAA",
        "gene_label_3": "Tbx5",
        "crrna_3": "AACCAATTAACCAATTGCAT",
        "gene_label_4": "Sox2",
        "crrna_4": "GCATGCATGCATGCATAACA",
    }
    response = client.post("/design", data=form_data)
    body = response.get_data(as_text=True)
    assert "<svg" in body
    assert "Assembly map" in body
    # All four gene labels should appear in the SVG.
    for label in ("Shh", "Pax6", "Tbx5", "Sox2"):
        assert label in body
