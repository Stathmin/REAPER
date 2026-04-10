# satMiner integration (optional)
#
# This module connects RepOrtR RepeatExplorer2/seqclust outputs to satMiner-like steps
# via Python3 adapters. Heavy downstream steps (RepeatMasker/DeconSeq) are intentionally
# not executed unless explicitly enabled/targeted.
#
# Enable by setting:
#   global:
#     satminer:
#       enabled: true

import os


def _satminer_enabled():
    return bool(config.get("global", {}).get("satminer", {}).get("enabled", False))


rule satminer_select_contigs_re2:
    """satMiner step 1c (RE2): select contigs until half total coverage per cluster."""
    input:
        clusters_dir="projects/{project}/samples/{sample}/tarean/seqclust/clustering/clusters",
    output:
        concat_fasta="projects/{project}/samples/{sample}/satminer/contigs_info_concat.fasta",
        contig_list="projects/{project}/samples/{sample}/satminer/selected_contigs.list",
        token="projects/{project}/samples/{sample}/satminer/CONTIG_LIST_READY",
    threads: 1
    conda:
        "../envs/reportr.yaml"
    run:
        if not _satminer_enabled():
            raise ValueError("satMiner rules are disabled. Set global.satminer.enabled: true to run.")
        shell(
            r"""
            set -euo pipefail
            python3 workflows/smk_scripts/satminer_select_contigs_re2.py \
              --clusters-dir "{input.clusters_dir}" \
              --out-fasta "{output.concat_fasta}" \
              --out-list "{output.contig_list}"
            touch "{output.token}"
            """
        )


rule satminer_extract_selected_contigs:
    """satMiner step 1c: extract selected contigs to FASTA."""
    input:
        concat_fasta="projects/{project}/samples/{sample}/satminer/contigs_info_concat.fasta",
        contig_list="projects/{project}/samples/{sample}/satminer/selected_contigs.list",
        token="projects/{project}/samples/{sample}/satminer/CONTIG_LIST_READY",
    output:
        selected_fasta="projects/{project}/samples/{sample}/satminer/selected_contigs.fasta",
        token="projects/{project}/samples/{sample}/satminer/SELECTED_CONTIGS_READY",
    threads: 1
    conda:
        "../envs/reportr.yaml"
    run:
        if not _satminer_enabled():
            raise ValueError("satMiner rules are disabled. Set global.satminer.enabled: true to run.")
        shell(
            r"""
            set -euo pipefail
            python3 workflows/smk_scripts/satminer_extract_selected_contigs.py \
              --fasta "{input.concat_fasta}" \
              --list "{input.contig_list}" \
              --out "{output.selected_fasta}"
            touch "{output.token}"
            """
        )


rule satminer_check_tools:
    """Fail fast if heavy satMiner tools are requested but missing."""
    output:
        touch("logs/satminer_toolcheck.ok")
    conda:
        "../envs/reportr.yaml"
    run:
        if not _satminer_enabled():
            raise ValueError("satMiner rules are disabled. Set global.satminer.enabled: true to run.")
        shell(
            r"""
            set -euo pipefail
            missing=0
            # satDNA analysis after RepeatExplorer2: RepeatMasker-first pipeline.
            # DeconSeq is intentionally not part of this integration.
            for tool in RepeatMasker calcDivergenceFromAlign.pl; do
              if ! command -v "$tool" >/dev/null 2>&1; then
                echo "MISSING: $tool" >&2
                missing=1
              fi
            done
            if [ "$missing" -ne 0 ]; then
              echo "Install into reportr env with: bash scripts/admin/install_satminer_deps.sh reportr" >&2
              exit 1
            fi
            """
        )


rule satminer_repeatmasker_satdna:
    """satDNA analysis step after RepeatExplorer2: run RepeatMasker on selected contigs (config-gated)."""
    input:
        selected_fasta="projects/{project}/samples/{sample}/satminer/selected_contigs.fasta",
        selected_token="projects/{project}/samples/{sample}/satminer/SELECTED_CONTIGS_READY",
        toolcheck="logs/satminer_toolcheck.ok",
    output:
        align="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.align",
        out="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.out",
        divsum="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.align.divsum",
        token="projects/{project}/samples/{sample}/satminer/SATDNA_ANALYSIS_DONE",
    params:
        lib=lambda w: get_param(w.project, "satminer", "repeatmasker_library", sample_id=w.sample),  # type: ignore[name-defined]
        threads=lambda w: int(get_param(w.project, "satminer", "repeatmasker_threads", sample_id=w.sample)),  # type: ignore[name-defined]
    threads:
        lambda w: int(get_param(w.project, "satminer", "repeatmasker_threads", sample_id=w.sample))  # type: ignore[name-defined]
    resources:
        repeatmasker_slots=1
    conda:
        "../envs/reportr.yaml"
    run:
        if not _satminer_enabled():
            raise ValueError("satMiner rules are disabled. Set global.satminer.enabled: true to run.")

        lib = str(params.lib or "").strip()
        if not lib:
            # Explicitly skip if library is not configured (empty in global; per-project opt-in).
            shell(
                r"""
                set -euo pipefail
                mkdir -p "$(dirname "{output.align}")"
                : > "{output.align}"
                echo "SKIPPED: repeatmasker_library not configured" > "{output.align}"
                echo "SKIPPED: repeatmasker_library not configured" >&2
                : > "{output.out}"
                : > "{output.divsum}"
                touch "{output.token}"
                """
            )
            return

        shell(
            r"""
            set -euo pipefail
            outdir="$(dirname "{output.align}")"
            mkdir -p "$outdir"

            # RepeatMasker writes outputs next to input by default; run in outdir.
            cp -f "{input.selected_fasta}" "$outdir/selected_contigs.fasta"
            cd "$outdir"

            RepeatMasker -pa "{params.threads}" -a -nolow -no_is -lib "{params.lib}" "selected_contigs.fasta"

            # Tokenized output path expected by the rule:
            test -f "selected_contigs.fasta.align"
            test -f "selected_contigs.fasta.out"
            calcDivergenceFromAlign.pl -s "selected_contigs.fasta.align.divsum" "selected_contigs.fasta.align"
            touch "{output.token}"
            """
        )


rule satminer_family_divergence:
    """satMiner 2c helper: remap subfamilies to families and compute divergence summary (optional)."""
    input:
        align="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.align",
        token="projects/{project}/samples/{sample}/satminer/SATDNA_ANALYSIS_DONE",
    output:
        fam_align="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.align.fam",
        fam_divsum="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.align.fam.divsum",
        token="projects/{project}/samples/{sample}/satminer/FAMILY_DIVERGENCE_DONE",
    params:
        pattern=lambda w: get_param(w.project, "satminer", "family_pattern_file", sample_id=w.sample),  # type: ignore[name-defined]
    conda:
        "../envs/reportr.yaml"
    run:
        if not _satminer_enabled():
            raise ValueError("satMiner rules are disabled. Set global.satminer.enabled: true to run.")
        pat = str(params.pattern or "").strip()
        if not pat:
            shell(
                r"""
                set -euo pipefail
                mkdir -p "$(dirname "{output.fam_align}")"
                : > "{output.fam_align}"
                echo "SKIPPED: family_pattern_file not configured" > "{output.fam_align}"
                echo "SKIPPED: family_pattern_file not configured" >&2
                : > "{output.fam_divsum}"
                touch "{output.token}"
                """
            )
            return
        shell(
            r"""
            set -euo pipefail
            python3 workflows/smk_scripts/satminer_subfam2fam.py \
              --align "{input.align}" \
              --pattern "{params.pattern}" \
              --out "{output.fam_align}"
            calcDivergenceFromAlign.pl -s "{output.fam_divsum}" "{output.fam_align}"
            touch "{output.token}"
            """
        )


rule satminer_rm_getseq:
    """satMiner 2b helper: extract RM hit sequences into FASTA (optional)."""
    input:
        fasta="projects/{project}/samples/{sample}/satminer/selected_contigs.fasta",
        rm_out="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.out",
        token="projects/{project}/samples/{sample}/satminer/SATDNA_ANALYSIS_DONE",
    output:
        hits_fasta="projects/{project}/samples/{sample}/satminer/repeatmasker/selected_contigs.fasta.out.fas",
        token="projects/{project}/samples/{sample}/satminer/RM_GETSEQ_DONE",
    params:
        min_len=lambda w: int(get_param(w.project, "satminer", "rm_getseq_min_len", sample_id=w.sample)),  # type: ignore[name-defined]
        enabled=lambda w: bool(get_param(w.project, "satminer", "rm_getseq_enabled", sample_id=w.sample)),  # type: ignore[name-defined]
    conda:
        "../envs/reportr.yaml"
    run:
        if not _satminer_enabled():
            raise ValueError("satMiner rules are disabled. Set global.satminer.enabled: true to run.")
        if not bool(params.enabled):
            shell(
                r"""
                set -euo pipefail
                mkdir -p "$(dirname "{output.hits_fasta}")"
                : > "{output.hits_fasta}"
                echo "SKIPPED: rm_getseq_enabled is false" > "{output.hits_fasta}"
                echo "SKIPPED: rm_getseq_enabled is false" >&2
                touch "{output.token}"
                """
            )
            return
        shell(
            r"""
            set -euo pipefail
            python3 workflows/smk_scripts/satminer_rm_getseq.py \
              --fasta "{input.fasta}" \
              --rm-out "{input.rm_out}" \
              --min-len "{params.min_len}" \
              --out "{output.hits_fasta}"
            touch "{output.token}"
            """
        )

