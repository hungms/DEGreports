# Tests for detect_organism() (HGNC/MGI gene-symbol-based organism
# auto-detection) and validate_msigdb_collections() (the multi-collection
# GSEA argument validator).

test_that("detect_organism identifies human gene symbols (all-caps HGNC convention)", {
    expect_equal(detect_organism(c("TP53", "BRCA1", "HLA-A", "IL6")), "human")
})

test_that("detect_organism identifies mouse gene symbols (title-case MGI convention)", {
    expect_equal(detect_organism(c("Trp53", "Brca1", "Il6", "Gm1234")), "mouse")
})

test_that("detect_organism drops NA/empty entries before classification", {
    expect_equal(detect_organism(c("TP53", NA, "", "BRCA1", "IL6")), "human")
})

test_that("detect_organism errors when no usable gene symbols are given", {
    expect_error(detect_organism(c(NA, "")), "no usable")
    expect_error(detect_organism(character(0)), "no usable")
})

test_that("detect_organism errors (rather than guesses) on ambiguous/mixed input", {
    # An even 50/50 split of human- and mouse-style symbols: neither
    # convention's share strictly exceeds the other's, so detection must
    # refuse to guess.
    expect_error(
        detect_organism(c("TP53", "BRCA1", "Trp53", "Brca1")),
        "Could not confidently detect organism"
    )
})

test_that("detect_organism respects a custom min_fraction threshold", {
    genes <- c("TP53", "BRCA1", "Trp53") # 2/3 human = 0.667
    expect_equal(detect_organism(genes, min_fraction = 0.5), "human")
    expect_error(detect_organism(genes, min_fraction = 0.9), "Could not confidently detect organism")
})

test_that("validate_msigdb_collections accepts each supported collection and rejects unknown ones", {
    for (coll in c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "PID")) {
        expect_equal(validate_msigdb_collections(coll), coll)
    }
    expect_equal(
        validate_msigdb_collections(c("HALLMARK", "GOBP", "REACTOME")),
        c("HALLMARK", "GOBP", "REACTOME")
    )
    expect_error(validate_msigdb_collections("NOT_A_REAL_COLLECTION"), "Unsupported MSigDB collection")
    expect_error(validate_msigdb_collections("KEGG_MEDICUS"), "Unsupported MSigDB collection")
    expect_error(validate_msigdb_collections(character(0)), "non-empty character vector")
    expect_equal(validate_msigdb_collections("KEGG_LEGACY"), "KEGG")
})

test_that(".match_msigdbr_subcollection resolves KEGG against the installed msigdbr's names", {
    # msigdbr >= 10 (MSigDB 2023.2+) renamed `CP:KEGG` to `CP:KEGG_LEGACY`,
    # so the alias must resolve differently depending on what is installed.
    local_mocked_bindings(.msigdbr_subcollections = function(cat) c("CP:KEGG", "CP:REACTOME"))
    expect_equal(.match_msigdbr_subcollection("C2", c("CP:KEGG_LEGACY", "CP:KEGG")), "CP:KEGG")
})

test_that(".match_msigdbr_subcollection prefers the earliest available candidate", {
    local_mocked_bindings(
        .msigdbr_subcollections = function(cat) c("CP:KEGG", "CP:KEGG_LEGACY")
    )
    expect_equal(.match_msigdbr_subcollection("C2", c("CP:KEGG_LEGACY", "CP:KEGG")), "CP:KEGG_LEGACY")
})

test_that(".match_msigdbr_subcollection passes through collections without a subcollection", {
    expect_null(.match_msigdbr_subcollection("H", NULL))
    expect_null(.match_msigdbr_subcollection("H", character(0)))
})

test_that(".match_msigdbr_subcollection falls back to the first candidate when the collection table is unreadable", {
    local_mocked_bindings(.msigdbr_subcollections = function(cat) NULL)
    expect_equal(.match_msigdbr_subcollection("C2", c("CP:KEGG_LEGACY", "CP:KEGG")), "CP:KEGG_LEGACY")
})

test_that(".match_msigdbr_subcollection errors informatively when no candidate exists", {
    local_mocked_bindings(.msigdbr_subcollections = function(cat) c("CP:BIOCARTA"))
    expect_error(
        .match_msigdbr_subcollection("C2", c("CP:KEGG_LEGACY", "CP:KEGG")),
        "Try updating the msigdbr package"
    )
})

test_that("import_msigdbr and import_msigdbr_term2gene fetch a non-empty REACTOME gene set", {
    skip_if_not_installed("msigdbr")
    gene_sets <- import_msigdbr("human", "REACTOME")
    expect_gt(length(gene_sets), 0)
    expect_true(all(startsWith(names(gene_sets), "REACTOME_")))

    term2gene <- import_msigdbr_term2gene("human", "REACTOME")
    expect_true(all(c("term", "gene") %in% colnames(term2gene)))
    expect_gt(nrow(term2gene), 0)
    expect_true(all(startsWith(unique(term2gene$term), "REACTOME_")))
})

test_that("import_msigdbr and import_msigdbr_term2gene fetch a non-empty KEGG gene set", {
    skip_if_not_installed("msigdbr")
    gene_sets <- import_msigdbr("human", "KEGG")
    expect_gt(length(gene_sets), 0)
    # NOTE: msigdbr's `CP:KEGG_LEGACY` subcollection still names its gene
    # sets with the older bare "KEGG_" prefix.
    expect_true(all(startsWith(names(gene_sets), "KEGG_")))

    term2gene <- import_msigdbr_term2gene("human", "KEGG")
    expect_true(all(c("term", "gene") %in% colnames(term2gene)))
    expect_gt(nrow(term2gene), 0)
    expect_true(all(startsWith(unique(term2gene$term), "KEGG_")))
})

