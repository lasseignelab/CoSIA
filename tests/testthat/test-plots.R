# plotSpeciesGEx and plotTissueGEx tests
# These run offline — VST data is injected as a list column to match the real
# ExperimentHub format and catch the separate_rows vs unnest regression.

make_plot_cosian <- function() {
    obj <- make_cosian(map_species = "m_musculus", map_tissues = "heart")

    mock_gex <- data.frame(
        Anatomical_entity_name = c("heart", "heart"),
        Ensembl_ID             = c("ENSG00000008710", "ENSMUSG00000032855"),
        Sample_size            = c(3L, 3L),
        Experiment_ID          = c("EXP1", "EXP2"),
        Anatomical_entity_ID   = c("UBERON:0000948", "UBERON:0000948"),
        Species                = c("Homo_sapiens", "Mus_musculus"),
        Scaled_Median_VST      = c(1.0, 1.0),
        mad                    = c(0.1, 0.1),
        stringsAsFactors       = FALSE
    )
    mock_gex$VST <- list(c(10.5, 11.0, 12.0), c(11.0, 11.5, 11.8))

    obj@gex <- mock_gex
    obj@converted_id <- data.frame(
        h_sapiens_ensembl_id  = "ENSG00000008710",
        h_sapiens_symbol      = "PKD1",
        m_musculus_ensembl_id = "ENSMUSG00000032855",
        m_musculus_symbol     = "Pkd1",
        stringsAsFactors      = FALSE
    )
    obj
}

# --- plotSpeciesGEx ---

test_that("plotSpeciesGEx returns a plotly object", {
    obj <- make_plot_cosian()
    fig <- plotSpeciesGEx(obj, "heart", "ENSG00000008710")
    expect_true(inherits(fig, "plotly"))
})

test_that("plotSpeciesGEx VST values are non-NA (list-column unnest)", {
    obj <- make_plot_cosian()
    fig <- plotly::plotly_build(plotSpeciesGEx(obj, "heart", "ENSG00000008710"))

    all_y <- unlist(lapply(fig$x$data, function(tr) tr$y))
    expect_false(
        all(is.na(all_y)),
        label = "all VST values are NA — separate_rows regression likely"
    )
    expect_true(all(all_y > 0, na.rm = TRUE))
})


# --- plotTissueGEx ---

make_tissue_cosian <- function() {
    obj <- make_cosian(map_species = "h_sapiens", map_tissues = c("heart", "brain"))

    mock_gex <- data.frame(
        Anatomical_entity_name = c("heart", "brain"),
        Ensembl_ID             = c("ENSG00000008710", "ENSG00000008710"),
        Sample_size            = c(3L, 4L),
        Experiment_ID          = c("EXP1", "EXP2"),
        Anatomical_entity_ID   = c("UBERON:0000948", "UBERON:0000955"),
        Species                = c("Homo_sapiens", "Homo_sapiens"),
        Scaled_Median_VST      = c(1.0, 1.0),
        mad                    = c(0.1, 0.1),
        stringsAsFactors       = FALSE
    )
    mock_gex$VST <- list(c(10.5, 11.0, 12.0), c(9.5, 10.0, 10.5, 11.0))

    obj@gex <- mock_gex
    obj@converted_id <- data.frame(
        h_sapiens_ensembl_id = "ENSG00000008710",
        h_sapiens_symbol     = "PKD1",
        stringsAsFactors     = FALSE
    )
    obj
}

test_that("plotTissueGEx returns a plotly object", {
    obj <- make_tissue_cosian()
    fig <- plotTissueGEx(obj, "h_sapiens", "ENSG00000008710")
    expect_true(inherits(fig, "plotly"))
})

test_that("plotTissueGEx VST values are non-NA (list-column unnest)", {
    obj <- make_tissue_cosian()
    fig <- plotly::plotly_build(plotTissueGEx(obj, "h_sapiens", "ENSG00000008710"))

    all_y <- unlist(lapply(fig$x$data, function(tr) tr$y))
    expect_false(
        all(is.na(all_y)),
        label = "all VST values are NA — separate_rows regression likely"
    )
    expect_true(all(all_y > 0, na.rm = TRUE))
})

