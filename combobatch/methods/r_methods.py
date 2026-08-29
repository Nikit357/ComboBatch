"""The 20 harmonization methods that run through R.

Two calling conventions, both ported from the donor benchmark:

* **In-memory** via ``importr`` — used where the R function takes a matrix directly.
* **File round-trip** — the matrix is written into a ``TemporaryDirectory``, an R script
  reads and writes files there, and Python reads the result back. Brittle-looking, but
  it isolates R's memory per call, which matters over thousands of jobs.

Systematic changes from the donor:

* Every argument that was frozen inside an interpolated R string (``polyFit(dis, 9)``,
  ``RUVIII(k=5)``, ``log_target=FALSE``, ...) is now a declared hyperparameter, threaded
  through :func:`combobatch.methods.rinterop.r_arglist`, which is the injection boundary.
* ``target_group`` is resolved by the caller, never defaulted to a dataset-specific name.
* ``37_fabatch`` is repaired: it called an API that does not exist and silently returned
  *uncorrected* data (see :func:`normalize_fabatch`).
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pandas as pd

from combobatch.methods.base import assert_rnaseq_only, drop_na_genes
from combobatch.methods.python_methods import DEFAULT_CONTROL_GENES
from combobatch.methods.rinterop import (
    as_numpy,
    py2rpy,
    r_arglist,
    r_gc,
    r_literal,
    require_r_packages,
)


def _ro():
    """Return the rpy2 robjects module, imported lazily."""
    import rpy2.robjects as ro

    return ro


def _importr(name: str):
    """Import an R package through rpy2."""
    from rpy2.robjects.packages import importr

    return importr(name)


# ── In-memory methods ─────────────────────────────────────────────────────────────


def normalize_limma(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, **kw
) -> pd.DataFrame:
    """
    ``limma::removeBatchEffect`` (Ritchie et al. 2015).

    Returns
    -------
    Batch-corrected expression matrix.
    """
    require_r_packages("03_limma", ("limma",))
    ro = _ro()
    limma, base = _importr("limma"), _importr("base")

    batches = ann_df.loc[exp_df.index, batch_col].astype(str).values
    r_mat = base.as_matrix(py2rpy(exp_df.T))
    result = as_numpy(limma.removeBatchEffect(r_mat, batch=ro.StrVector(batches)))
    out = pd.DataFrame(result.T, index=exp_df.index, columns=exp_df.columns)
    r_gc()
    return out


def normalize_sva(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, bio_col: str, **kw
) -> pd.DataFrame:
    """
    Surrogate Variable Analysis (Leek & Storey 2012), with biology protected in the model.

    Returns
    -------
    SVA-corrected expression matrix, or the input unchanged when SVA finds no
    surrogate variables to remove.
    """
    require_r_packages("04_sva", ("sva", "limma"))
    from rpy2.rinterface_lib.embedded import RRuntimeError

    ro = _ro()
    sva, limma, base = _importr("sva"), _importr("limma"), _importr("base")

    bio = ann_df.loc[exp_df.index, bio_col].fillna("Unknown").astype(str).values
    frame = ro.DataFrame({"r_bio": ro.StrVector(bio)})
    formula = "~r_bio" if len(set(bio)) >= 2 else "~1"
    mod = ro.r["model.matrix"](ro.Formula(formula), data=frame)
    mod0 = ro.r["model.matrix"](ro.Formula("~1"), data=frame)

    exp_in, _ = drop_na_genes(exp_df)
    r_mat = base.as_matrix(py2rpy(exp_in.T))

    n_sv = sva.num_sv(r_mat, mod)
    if int(n_sv[0]) == 0:
        return exp_df.copy()

    try:
        sv_obj = sva.sva(r_mat, mod, mod0, n_sv=n_sv)
    except RRuntimeError:
        return exp_df.copy()

    surrogates = as_numpy(sv_obj.rx2("sv"))
    if surrogates.ndim == 1:
        surrogates = surrogates.reshape(-1, 1)
    elif surrogates.shape[0] != len(exp_in) and surrogates.shape[1] == len(exp_in):
        surrogates = surrogates.T

    try:
        result = limma.removeBatchEffect(
            r_mat,
            covariates=base.as_matrix(
                py2rpy(pd.DataFrame(surrogates, index=exp_in.index))
            ),
        )
    except RRuntimeError:
        return exp_df.copy()

    out = pd.DataFrame(as_numpy(result).T, index=exp_in.index, columns=exp_in.columns)
    r_gc()
    return out


def normalize_combat(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, bio_col: str, **kw
) -> pd.DataFrame:
    """
    ComBat empirical-Bayes correction (Johnson et al. 2007).

    Returns
    -------
    ComBat-corrected expression matrix.
    """
    require_r_packages("05_combat", ("sva",))
    ro = _ro()
    sva, base = _importr("sva"), _importr("base")

    batches = ann_df.loc[exp_df.index, batch_col].astype(str).values
    bio = ann_df.loc[exp_df.index, bio_col].fillna("Unknown").astype(str).values
    r_mat = base.as_matrix(py2rpy(exp_df.T))

    if len(set(bio)) >= 2:
        mod = ro.r["model.matrix"](
            ro.Formula("~r_bio"), data=ro.DataFrame({"r_bio": ro.StrVector(bio)})
        )
    else:
        mod = ro.rinterface.NULL

    result = sva.ComBat(dat=r_mat, batch=ro.StrVector(batches), mod=mod)
    out = pd.DataFrame(as_numpy(result).T, index=exp_df.index, columns=exp_df.columns)
    r_gc()
    return out


def normalize_combat_seq(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, bio_col: str, **kw
) -> pd.DataFrame:
    """
    ComBat-seq (Zhang et al. 2020), on integer-rounded counts.

    Returns
    -------
    ComBat-seq-corrected expression matrix.
    """
    require_r_packages("06_combat_seq", ("sva",))
    ro = _ro()
    sva, base = _importr("sva"), _importr("base")

    batches = ann_df.loc[exp_df.index, batch_col].astype(str).values
    bio = ann_df.loc[exp_df.index, bio_col].fillna("Unknown").astype(str).values
    counts = np.round(exp_df.values).astype(int).clip(0)
    r_mat = base.as_matrix(
        py2rpy(pd.DataFrame(counts.T, index=exp_df.columns, columns=exp_df.index))
    )
    group = ro.StrVector(bio) if len(set(bio)) >= 2 else ro.rinterface.NULL

    result = sva.ComBat_seq(counts=r_mat, batch=ro.StrVector(batches), group=group)
    out = pd.DataFrame(as_numpy(result).T, index=exp_df.index, columns=exp_df.columns)
    r_gc()
    return out


def normalize_ruv(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    k: int = 2,
    control_genes: tuple[str, ...] = DEFAULT_CONTROL_GENES,
    **kw,
) -> pd.DataFrame:
    """
    RUVg (Risso et al. 2014), using housekeeping genes as negative controls.

    Returns
    -------
    RUV-corrected expression matrix, or the input unchanged when too few control
    genes are present to estimate unwanted variation.
    """
    require_r_packages("09_ruv", ("RUVSeq",))
    ro = _ro()
    ruvseq, base = _importr("RUVSeq"), _importr("base")

    present = [gene for gene in control_genes if gene in exp_df.columns]
    if len(present) < 3:
        return exp_df.copy()

    counts = np.round(exp_df.values).astype(int).clip(0)
    r_counts = base.as_matrix(
        py2rpy(pd.DataFrame(counts.T, index=exp_df.columns, columns=exp_df.index))
    )
    columns = list(exp_df.columns)
    ctrl_idx = ro.IntVector([columns.index(gene) + 1 for gene in present])

    ruv_obj = ruvseq.RUVg(r_counts, cIdx=ctrl_idx, k=k)
    normalized = as_numpy(ruv_obj.rx2("normalizedCounts"))
    out = pd.DataFrame(normalized.T, index=exp_df.index, columns=exp_df.columns)
    r_gc()
    return out


def normalize_mnn(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, k: int = 20, **kw
) -> pd.DataFrame:
    """
    fastMNN (Haghverdi et al. 2018), inverse-projected from the corrected embedding.

    Returns
    -------
    MNN-corrected expression matrix.
    """
    require_r_packages("10_mnn", ("batchelor",))
    from sklearn.decomposition import PCA

    ro = _ro()
    batchelor, base = _importr("batchelor"), _importr("base")

    batches = ann_df.loc[exp_df.index, batch_col].astype(str)
    exp_in, _ = drop_na_genes(exp_df)
    exp_in = exp_in.clip(lower=0.0)
    r_mat = base.as_matrix(py2rpy(exp_in.T))

    result = batchelor.fastMNN(r_mat, batch=ro.StrVector(batches.values), k=k)
    embedding = as_numpy(ro.r["as.matrix"](ro.r["reducedDim"](result, "corrected")))

    pca = PCA(n_components=embedding.shape[1]).fit(exp_in.values)
    recovered = embedding @ pca.components_ + pca.mean_
    out = pd.DataFrame(recovered, index=exp_in.index, columns=exp_in.columns)
    r_gc()
    return out


def normalize_qsmooth(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, **kw
) -> pd.DataFrame:
    """
    Smooth quantile normalization (Hicks et al. 2018).

    Returns
    -------
    qsmooth-normalized expression matrix.
    """
    require_r_packages("14_qsmooth", ("qsmooth",))
    ro = _ro()
    qsmooth, base = _importr("qsmooth"), _importr("base")

    groups = ann_df.loc[exp_df.index, batch_col].astype(str).values
    exp_in, _ = drop_na_genes(exp_df)
    r_mat = base.as_matrix(py2rpy(exp_in.T))

    qs_obj = qsmooth.qsmooth(r_mat, group_factor=ro.StrVector(groups))
    result = as_numpy(qsmooth.qsmoothData(qs_obj))
    out = pd.DataFrame(result.T, index=exp_in.index, columns=exp_in.columns)
    r_gc()
    return out


def normalize_tmm(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    method: str = "TMM",
    prior_count: float = 1.0,
    **kw,
) -> pd.DataFrame:
    """
    TMM normalization via edgeR (Robinson & Oshlack 2010). RNA-seq only.

    Returns
    -------
    log2-CPM after TMM scaling.
    """
    assert_rnaseq_only(ann_df.loc[exp_df.index], batch_col, "22_tmm")
    require_r_packages("22_tmm", ("edgeR",))
    edger, base = _importr("edgeR"), _importr("base")

    counts = np.round(exp_df.values).astype(int).clip(0)
    r_mat = base.as_matrix(
        py2rpy(pd.DataFrame(counts.T, index=exp_df.columns, columns=exp_df.index))
    )
    dge = edger.calcNormFactors(edger.DGEList(counts=r_mat), method=method)
    result = as_numpy(edger.cpm(dge, log=True, prior_count=prior_count))
    out = pd.DataFrame(result.T, index=exp_df.index, columns=exp_df.columns)
    r_gc()
    return out


def normalize_vst(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    blind: bool = True,
    min_genes_for_vst: int = 1000,
    **kw,
) -> pd.DataFrame:
    """
    DESeq2 variance-stabilizing transformation (Love et al. 2014). RNA-seq only.

    ``vst()`` needs at least ``min_genes_for_vst`` genes for its subsetting step;
    below that DESeq2's exact transformation is used instead.

    Returns
    -------
    Variance-stabilized expression matrix.
    """
    assert_rnaseq_only(ann_df.loc[exp_df.index], batch_col, "23_vst")
    require_r_packages("23_vst", ("DESeq2",))
    ro = _ro()
    deseq2, base = _importr("DESeq2"), _importr("base")

    counts = np.round(exp_df.values).astype(int).clip(0)
    r_mat = base.as_matrix(
        py2rpy(pd.DataFrame(counts.T, index=exp_df.columns, columns=exp_df.index))
    )
    col_data = ro.DataFrame({"sample": ro.StrVector(list(exp_df.index))})
    dds = deseq2.DESeqDataSetFromMatrix(
        countData=r_mat, colData=col_data, design=ro.Formula("~1")
    )
    vsd = (
        deseq2.vst(dds, blind=blind)
        if exp_df.shape[1] >= min_genes_for_vst
        else deseq2.varianceStabilizingTransformation(dds, blind=blind)
    )
    result = as_numpy(ro.r["assay"](vsd))
    out = pd.DataFrame(result.T, index=exp_df.index, columns=exp_df.columns)
    r_gc()
    return out


# ── File round-trip methods ───────────────────────────────────────────────────────


def _run_r(script: str) -> None:
    """Execute an R script string."""
    _ro().r(script)


def _read_result(path: str, method_name: str, hint: str = "") -> pd.DataFrame:
    """Read an R-produced CSV, raising a clear error when the script wrote nothing."""
    if not os.path.exists(path):
        raise RuntimeError(f"{method_name} produced no output. {hint}".strip())
    return pd.read_csv(path, index_col=0)


def normalize_fsqn_r(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    target_group: str | None = None,
    **kw,
) -> pd.DataFrame:
    """
    FSQN via the original R package (Franks et al. 2018).

    Each non-reference batch is normalized independently against the reference batch,
    in samples x genes orientation, which is what FSQN expects. Reference samples pass
    through unchanged. With no reference present this degenerates to per-gene quantile
    normalization over all samples.

    Returns
    -------
    FSQN-normalized expression matrix.
    """
    require_r_packages("16_fsqn_r", ("FSQN",))
    exp_in, _ = drop_na_genes(exp_df)
    groups = ann_df.loc[exp_in.index, batch_col].astype(str)
    use_reference = target_group is not None and target_group in groups.values

    with tempfile.TemporaryDirectory() as td:
        ref_path = os.path.join(td, "ref.csv")
        out_path = os.path.join(td, "out_0000.csv")

        if not use_reference:
            query_path = os.path.join(td, "query.csv")
            exp_in.to_csv(ref_path)
            exp_in.to_csv(query_path)
            _run_r(f"""
                library(FSQN)
                ref   <- as.matrix(read.csv({r_literal(ref_path)}, row.names=1,
                                            check.names=FALSE))
                query <- as.matrix(read.csv({r_literal(query_path)}, row.names=1,
                                            check.names=FALSE))
                out   <- quantileNormalizeByFeature(query, ref)
                if (is.null(out)) stop("quantileNormalizeByFeature returned NULL")
                write.csv(out, {r_literal(out_path)}, row.names=TRUE)
                """)
            result = _read_result(out_path, "16_fsqn_r")
            r_gc()
            return result.reindex(exp_df.index)

        reference = exp_in.loc[groups == target_group]
        reference.to_csv(ref_path)
        non_reference = [b for b in groups.unique() if b != target_group]

        blocks = []
        for position, batch in enumerate(non_reference):
            query_path = os.path.join(td, f"query_{position:04d}.csv")
            batch_out = os.path.join(td, f"out_{position:04d}.csv")
            exp_in.loc[groups == batch].to_csv(query_path)
            blocks.append(f"""
                query <- as.matrix(read.csv({r_literal(query_path)}, row.names=1,
                                            check.names=FALSE))
                if (ncol(query) != n_ref_genes) stop("FSQN gene count mismatch")
                out <- quantileNormalizeByFeature(query, ref)
                if (is.null(out)) stop("quantileNormalizeByFeature returned NULL")
                write.csv(out, {r_literal(batch_out)}, row.names=TRUE)
                """)

        _run_r(f"""
            library(FSQN)
            ref         <- as.matrix(read.csv({r_literal(ref_path)}, row.names=1,
                                              check.names=FALSE))
            n_ref_genes <- ncol(ref)
            {''.join(blocks)}
            """)

        parts = [reference]
        for position, batch in enumerate(non_reference):
            parts.append(
                _read_result(
                    os.path.join(td, f"out_{position:04d}.csv"),
                    "16_fsqn_r",
                    f"(batch {batch!r})",
                )
            )

    r_gc()
    return pd.concat(parts).reindex(exp_df.index)


def normalize_tdm(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    target_group: str | None = None,
    log_target: bool = False,
    **kw,
) -> pd.DataFrame:
    """
    Training Distribution Matching (Thompson et al. 2016).

    Returns
    -------
    TDM-normalized expression matrix.
    """
    require_r_packages("19_tdm", ("TDM",))
    groups = ann_df.loc[exp_df.index, batch_col].astype(str)
    exp_in, _ = drop_na_genes(exp_df)
    target_mask = (
        groups == target_group
        if target_group is not None and target_group in groups.values
        else pd.Series(True, index=exp_in.index)
    )

    with tempfile.TemporaryDirectory() as td:
        ref_path = os.path.join(td, "ref.tsv")
        query_path = os.path.join(td, "query.tsv")
        out_path = os.path.join(td, "out.tsv")

        for frame, path in ((exp_in.loc[target_mask], ref_path), (exp_in, query_path)):
            table = frame.T.copy()
            table.index.name = "gene"
            table.reset_index().to_csv(path, sep="\t", index=False)

        _run_r(f"""
            library(TDM)
            result <- tdm_transform(file={r_literal(query_path)},
                                    ref_file={r_literal(ref_path)},
                                    log_target={r_literal(log_target)})
            write.table(result, {r_literal(out_path)}, sep="\\t",
                        row.names=FALSE, quote=FALSE)
            """)
        if not os.path.exists(out_path):
            raise RuntimeError("19_tdm produced no output")
        result = pd.read_csv(out_path, sep="\t")

    result = result.drop(columns=[result.columns[0]])
    r_gc()
    return pd.DataFrame(
        result.values.T.astype(float), index=exp_in.index, columns=exp_in.columns
    ).reindex(exp_df.index)


def normalize_harmonizr(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    algorithm: str = "ComBat",
    **kw,
) -> pd.DataFrame:
    """
    HarmonizR (Voss et al. 2022): NA-aware ComBat/limma over matrix blocks.

    Returns
    -------
    HarmonizR-corrected expression matrix.
    """
    require_r_packages("21_harmonizr", ("HarmonizR",))
    batches = ann_df.loc[exp_df.index, batch_col].astype(str)
    batch_ints = batches.map({b: i + 1 for i, b in enumerate(batches.unique())})

    with tempfile.TemporaryDirectory() as td:
        data_path = os.path.join(td, "data.tsv")
        desc_path = os.path.join(td, "description.csv")
        out_base = os.path.join(td, "out")

        exp_df.T.to_csv(data_path, sep="\t")
        pd.DataFrame(
            {"ID": exp_df.index, "batch": batches.values, "sample": batch_ints.values}
        ).to_csv(desc_path, index=False)

        _run_r(f"""
            library(HarmonizR)
            harmonizR(data_as_input={r_literal(data_path)},
                      description_as_input={r_literal(desc_path)},
                      algorithm={r_literal(algorithm)},
                      output_file={r_literal(out_base)},
                      plot=FALSE, verbosity=0)
            """)
        out_path = out_base + ".tsv"
        if not os.path.exists(out_path):
            raise RuntimeError("21_harmonizr produced no output")
        result = pd.read_csv(out_path, sep="\t", index_col=0)

    r_gc()
    return result.T.reindex(exp_df.index)


def normalize_dwd(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    target_group: str | None = None,
    min_batch_size: int = 5,
    **kw,
) -> pd.DataFrame:
    """
    DWD per-batch projection via DWDLargeR (Qing & Marron 2018).

    Note ``genDWD`` is called without a ``penalty`` argument: DWDLargeR removed it, and
    passing it made every batch fall through to the error handler uncorrected.

    Returns
    -------
    DWD-corrected expression matrix.
    """
    require_r_packages("27_dwd", ("DWDLargeR",))
    exp_in, _ = drop_na_genes(exp_df)
    groups = ann_df.loc[exp_in.index, batch_col].astype(str)
    reference = (
        target_group
        if target_group is not None and target_group in groups.values
        else groups.iloc[0]
    )

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.T.to_csv(exp_path)
        groups.rename("batch").to_csv(ann_path)

        _run_r(f"""
            library(DWDLargeR)
            exp_mat   <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1,
                                            check.names=FALSE))
            ann_df    <- read.csv({r_literal(ann_path)}, row.names=1)
            batches   <- ann_df[colnames(exp_mat), "batch"]
            ref_label <- {r_literal(str(reference))}

            ref_mask <- batches == ref_label
            if (sum(ref_mask, na.rm=TRUE) == 0) ref_mask <- rep(TRUE, ncol(exp_mat))
            X_ref  <- exp_mat[, ref_mask, drop=FALSE]
            result <- exp_mat

            for (b in setdiff(unique(batches), ref_label)) {{
                batch_mask <- batches == b
                if (sum(batch_mask) < {int(min_batch_size)}) next
                X_batch    <- exp_mat[, batch_mask, drop=FALSE]
                X_combined <- cbind(X_ref, X_batch)
                y <- c(rep(1, ncol(X_ref)), rep(-1, ncol(X_batch)))
                tryCatch({{
                    sol <- genDWD(X = X_combined, y = y)
                    w   <- as.numeric(if (!is.null(sol$beta)) sol$beta else sol$w)
                    shift <- mean(as.numeric(t(X_ref) %*% w)) -
                             mean(as.numeric(t(X_batch) %*% w))
                    result[, batch_mask] <- X_batch + shift * w
                }}, error = function(e) {{
                    message("DWD failed for batch '", b, "': ", conditionMessage(e))
                }})
            }}
            write.csv(result, {r_literal(out_path)})
            """)
        result = _read_result(out_path, "27_dwd", "Install DWDLargeR.")

    r_gc()
    return result.T.reindex(exp_df.index)


def normalize_npn(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    npn_func: str = "truncation",
    **kw,
) -> pd.DataFrame:
    """
    Nonparanormal transform per batch via ``huge::huge.npn`` (Liu et al. 2009).

    Every gene's marginal becomes standard normal within each batch, which removes
    between-batch mean and variance shifts without needing a reference.

    Returns
    -------
    NPN-normalized expression matrix.
    """
    require_r_packages("28_npn", ("huge",))
    exp_in, _ = drop_na_genes(exp_df)

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.to_csv(exp_path)
        ann_df.loc[exp_in.index, [batch_col]].to_csv(ann_path)

        _run_r(f"""
            library(huge)
            exp_mat <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1))
            ann_df  <- read.csv({r_literal(ann_path)}, row.names=1)
            batches <- as.character(ann_df[rownames(exp_mat), {r_literal(batch_col)}])
            result  <- exp_mat
            for (b in unique(batches)) {{
                mask <- batches == b
                result[mask, ] <- huge.npn(exp_mat[mask, , drop=FALSE],
                                           npn.func={r_literal(npn_func)},
                                           verbose=FALSE)
            }}
            write.csv(as.data.frame(result), {r_literal(out_path)})
            """)
        result = _read_result(out_path, "28_npn", "Install huge.")

    r_gc()
    return result.reindex(exp_df.index)


def normalize_combat_ref(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    bio_col: str,
    target_group: str | None = None,
    **kw,
) -> pd.DataFrame:
    """
    M-ComBat: ComBat anchored to a reference batch rather than the global mean.

    Returns
    -------
    Reference-anchored ComBat expression matrix.
    """
    require_r_packages("29_combat_ref", ("sva",))
    exp_in, _ = drop_na_genes(exp_df)
    groups = ann_df.loc[exp_in.index, batch_col].astype(str)
    use_reference = target_group is not None and target_group in groups.values

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.T.to_csv(exp_path)
        ann_df.loc[exp_in.index, [batch_col, bio_col]].fillna("Unknown").to_csv(
            ann_path
        )

        ref_arg = f", ref.batch={r_literal(target_group)}" if use_reference else ""
        _run_r(f"""
            library(sva)
            exp_mat <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1,
                                          check.names=FALSE))
            ann_df  <- read.csv({r_literal(ann_path)}, row.names=1)
            batches <- as.character(ann_df[colnames(exp_mat), {r_literal(batch_col)}])
            bio     <- as.character(ann_df[colnames(exp_mat), {r_literal(bio_col)}])
            mod     <- model.matrix(~ bio)
            result  <- ComBat(dat=exp_mat, batch=batches, mod=mod{ref_arg})
            write.csv(as.data.frame(result), {r_literal(out_path)})
            """)
        result = _read_result(out_path, "29_combat_ref", "Check the sva package.")

    r_gc()
    return result.T.reindex(exp_df.index)


def normalize_ruv3prps(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    bio_col: str,
    k_factors: int = 5,
    min_cell_size: int = 2,
    **kw,
) -> pd.DataFrame:
    """
    RUV-III with pseudo-replicate pseudo-samples (Molania et al. 2022).

    Each (biology x batch) cell with enough samples contributes its mean as a
    pseudo-replicate, giving RUV-III the replicate structure it needs.

    Returns
    -------
    RUV-III-PRPS corrected expression matrix.
    """
    require_r_packages("31_ruv3prps", ("ruv",))
    exp_in, _ = drop_na_genes(exp_df)

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.to_csv(exp_path)
        ann_df.loc[exp_in.index, [batch_col, bio_col]].to_csv(ann_path)

        _run_r(f"""
            library(ruv)
            Y       <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1))
            ann_df  <- read.csv({r_literal(ann_path)}, row.names=1)
            batches <- as.character(ann_df[rownames(Y), {r_literal(batch_col)}])
            bio     <- as.character(ann_df[rownames(Y), {r_literal(bio_col)}])
            cells   <- paste(bio, batches, sep="__")
            valid   <- names(which(table(cells) >= {int(min_cell_size)}))

            if (length(valid) == 0) {{
                write.csv(as.data.frame(Y), {r_literal(out_path)})
            }} else {{
                pseudo <- do.call(rbind, lapply(valid, function(cell) {{
                    colMeans(Y[cells == cell, , drop=FALSE])
                }}))
                rownames(pseudo) <- paste0("pseudo__", valid)
                Y_aug  <- rbind(Y, pseudo)
                n_orig <- nrow(Y)
                M <- matrix(0L, nrow=nrow(Y_aug), ncol=length(valid))
                for (j in seq_along(valid)) {{
                    M[c(which(cells == valid[j]), n_orig + j), j] <- 1L
                }}
                k <- max(1L, min({int(k_factors)}, length(valid) - 1L))
                fit  <- RUVIII(Y=Y_aug, M=M, ctl=rep(TRUE, ncol(Y_aug)), k=k)
                newY <- if (is.list(fit)) {{
                    if (!is.null(fit$newY)) fit$newY else fit$corrY
                }} else fit
                result <- newY[seq_len(n_orig), , drop=FALSE]
                rownames(result) <- rownames(Y)
                write.csv(as.data.frame(result), {r_literal(out_path)})
            }}
            """)
        result = _read_result(out_path, "31_ruv3prps", "Check the ruv package.")

    r_gc()
    return result.reindex(exp_df.index)


def normalize_amdbnorm(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    target_group: str | None = None,
    poly_degree: int = 9,
    n_dist_points: int = 500,
    **kw,
) -> pd.DataFrame:
    """
    AMDBNorm: fit a polynomial to the reference batch distribution and correct toward it.

    Returns
    -------
    AMDBNorm-corrected expression matrix.
    """
    require_r_packages("33_amdbnorm", ("AMDBNorm", "DBNorm"))
    exp_in, _ = drop_na_genes(exp_df)
    groups = ann_df.loc[exp_in.index, batch_col].astype(str)
    reference = (
        target_group
        if target_group is not None and target_group in groups.values
        else groups.iloc[0]
    )

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.T.to_csv(exp_path)
        groups.rename("batch").to_csv(ann_path)

        _run_r(f"""
            library(AMDBNorm); library(DBNorm); library(reshape2)
            exp_mat   <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1,
                                            check.names=FALSE))
            batches   <- read.csv({r_literal(ann_path)}, row.names=1)[
                             colnames(exp_mat), "batch"]
            ref_label <- {r_literal(str(reference))}
            ref_mat   <- exp_mat[, batches == ref_label, drop=FALSE]

            dis <- DBNorm::genDistData(na.omit(melt(ref_mat)[, 3]),
                                       {int(n_dist_points)})
            fit <- DBNorm::polyFit(dis, {int(poly_degree)})

            result <- ref_mat
            for (b in setdiff(unique(batches), ref_label)) {{
                corrected <- AMDBNorm(data=exp_mat[, batches == b, drop=FALSE],
                                      ref_batch=ref_mat, fit=fit)
                result    <- cbind(result, corrected)
            }}
            write.csv(as.data.frame(result[, colnames(exp_mat)]), {r_literal(out_path)})
            """)
        result = _read_result(
            out_path, "33_amdbnorm", "Install AMDBNorm and DBNorm from GitHub."
        )

    r_gc()
    return result.T.reindex(exp_df.index)


def normalize_arsyn(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    bio_col: str,
    norm: str = "n",
    logtransf: bool = False,
    **kw,
) -> pd.DataFrame:
    """
    ARSyN: ANOVA-ASCA removal of the batch component (NOISeq).

    Returns
    -------
    ARSyN-corrected expression matrix.
    """
    require_r_packages("34_arsyn", ("NOISeq",))
    exp_in, _ = drop_na_genes(exp_df)

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.T.to_csv(exp_path)
        ann_df.loc[exp_in.index, [batch_col, bio_col]].fillna("Unknown").to_csv(
            ann_path
        )

        _run_r(f"""
            library(NOISeq)
            exp_mat <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1))
            ann_df  <- read.csv({r_literal(ann_path)}, row.names=1)
            factors <- data.frame(
                Batch     = as.character(ann_df[colnames(exp_mat),
                                                {r_literal(batch_col)}]),
                Condition = as.character(ann_df[colnames(exp_mat),
                                                {r_literal(bio_col)}]),
                row.names = colnames(exp_mat))
            mydata   <- readData(data=exp_mat, factors=factors)
            myresult <- ARSyNseq(mydata, factor="Batch", norm={r_literal(norm)},
                                 logtransf={r_literal(logtransf)})
            write.csv(as.data.frame(assayData(myresult)$exprs), {r_literal(out_path)})
            """)
        result = _read_result(out_path, "34_arsyn", "Check the NOISeq package.")

    r_gc()
    return result.T.reindex(exp_df.index)


def normalize_fabatch(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    bio_col: str,
    target_col: str | None = None,
    **kw,
) -> pd.DataFrame:
    """
    FAbatch factor-adjustment batch correction, via ``bapred::fabatch``.

    Repaired from the donor, which called ``fabatch(xtr=, ytr=, type="among")`` and read
    ``$adj.data``. None of those exist in bapred: R raised "unused arguments", the
    surrounding ``tryCatch`` caught it, and the *uncorrected* matrix was written out and
    scored as though it had been corrected. The real signature is
    ``fabatch(x, y, batch)`` returning ``$xadj``, and there is no fallback here — a
    failure is a failure.

    ``bapred`` structurally requires a **two-level** target, so this raises
    ``NotImplementedError`` when the target column has any other number of levels. That
    is an honest skip; auto-binarizing would silently change the question being asked.

    Parameters
    ----------
    target_col
        Annotation column supplying the binary target. Defaults to the biology column.

    Returns
    -------
    FAbatch-corrected expression matrix.
    """
    require_r_packages("37_fabatch", ("bapred",))
    exp_in, _ = drop_na_genes(exp_df)
    resolved_target = target_col or bio_col

    if resolved_target not in ann_df.columns:
        raise NotImplementedError(
            f"37_fabatch: target column {resolved_target!r} is not in the annotation."
        )

    target = ann_df.loc[exp_in.index, resolved_target].fillna("Unknown").astype(str)
    n_levels = target.nunique()
    if n_levels != 2:
        raise NotImplementedError(
            f"37_fabatch needs a two-level target, but {resolved_target!r} has "
            f"{n_levels} level(s). Point --fabatch-target-col at a binary column, or "
            "exclude this method from the run."
        )

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.T.to_csv(exp_path)
        frame = ann_df.loc[exp_in.index, [batch_col]].copy()
        frame["__target__"] = target
        frame.to_csv(ann_path)

        _run_r(f"""
            library(bapred)
            exp_mat <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1,
                                          check.names=FALSE))
            ann_df  <- read.csv({r_literal(ann_path)}, row.names=1)
            batches <- as.factor(ann_df[colnames(exp_mat), {r_literal(batch_col)}])
            y       <- as.factor(ann_df[colnames(exp_mat), "__target__"])

            res    <- fabatch(x = t(exp_mat), y = y, batch = batches)
            result <- t(res$xadj)
            write.csv(as.data.frame(result), {r_literal(out_path)})
            """)
        result = _read_result(
            out_path, "37_fabatch", "Install bapred from CRAN (BiocManager::install)."
        )

    r_gc()
    return result.T.reindex(exp_df.index)


def normalize_harman(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    bio_col: str,
    limit: float = 0.1,
    **kw,
) -> pd.DataFrame:
    """
    Harman: PCA-based correction with an overcorrection constraint.

    Returns
    -------
    Harman-corrected expression matrix.
    """
    require_r_packages("38_harman", ("Harman",))
    exp_in, _ = drop_na_genes(exp_df)

    with tempfile.TemporaryDirectory() as td:
        exp_path = os.path.join(td, "exp.csv")
        ann_path = os.path.join(td, "ann.csv")
        out_path = os.path.join(td, "out.csv")

        exp_in.T.to_csv(exp_path)
        ann_df.loc[exp_in.index, [batch_col, bio_col]].fillna("Unknown").to_csv(
            ann_path
        )

        _run_r(f"""
            library(Harman)
            exp_mat <- as.matrix(read.csv({r_literal(exp_path)}, row.names=1))
            ann_df  <- read.csv({r_literal(ann_path)}, row.names=1)
            batches <- as.character(ann_df[colnames(exp_mat), {r_literal(batch_col)}])
            bio     <- as.character(ann_df[colnames(exp_mat), {r_literal(bio_col)}])
            pc      <- harman(exp_mat, expt=bio, batch=batches,
                              {r_arglist({"limit": limit})})
            write.csv(as.data.frame(reconstructData(pc)), {r_literal(out_path)})
            """)
        result = _read_result(out_path, "38_harman", "Check the Harman package.")

    r_gc()
    return result.T.reindex(exp_df.index)
