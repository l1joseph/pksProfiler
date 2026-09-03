#!/bin/bash
# Build clb_all_protein.hmm from NCBI colibactin protein sequences
# Run with: conda activate hmm_build && bash scripts/build_clb_protein_hmm.sh

set -euo pipefail

WORKDIR="$(dirname "$0")/../ref/hmm/clb_protein_build"
OUTHMM="$(dirname "$0")/../ref/hmm/clb_all_protein.hmm"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

GENES=(A B C D E F G H I J K L M N O P Q R S)

echo "=== Fetching colibactin proteins from NCBI ==="
for g in "${GENES[@]}"; do
    gene="clb${g}"
    fa="${gene}.fasta"
    if [[ -s "${fa}" ]]; then
        echo "  ${gene}: already fetched ($(grep -c '^>' "${fa}") seqs)"
        continue
    fi
    echo "  Fetching ${gene}..."
    # Use broad query: gene name + colibactin keyword
    esearch -db protein -query "${gene}[Gene Name] colibactin" \
        | efetch -format fasta > "${fa}" 2>/dev/null
    n=$(grep -c '^>' "${fa}" 2>/dev/null || echo 0)
    # Fallback: product name search if gene name returns nothing
    if [[ "${n}" -lt 2 ]]; then
        echo "    Gene name search got ${n} hits, trying product search..."
        esearch -db protein -query "colibactin ${gene}" \
            | efetch -format fasta > "${fa}" 2>/dev/null
        n=$(grep -c '^>' "${fa}" 2>/dev/null || echo 0)
    fi
    echo "    ${n} sequences fetched"
    sleep 1
done

echo ""
echo "=== Sequence counts ==="
for g in "${GENES[@]}"; do
    gene="clb${g}"
    n=$(grep -c '^>' "${gene}.fasta" 2>/dev/null || echo 0)
    echo "  ${gene}: ${n}"
done

echo ""
echo "=== Aligning with MAFFT ==="
for g in "${GENES[@]}"; do
    gene="clb${g}"
    fa="${gene}.fasta"
    msa="${gene}.msa"
    n=$(grep -c '^>' "${fa}" 2>/dev/null || echo 0)
    if [[ "${n}" -lt 1 ]]; then
        echo "  SKIP ${gene}: no sequences"
        continue
    fi
    if [[ -s "${msa}" ]]; then
        echo "  ${gene}: already aligned"
        continue
    fi
    echo "  Aligning ${gene} (${n} seqs)..."
    if [[ "${n}" -eq 1 ]]; then
        # Single sequence — copy as-is (no alignment needed)
        cp "${fa}" "${msa}"
    else
        mafft --auto --quiet "${fa}" > "${msa}"
    fi
done

echo ""
echo "=== Building per-gene HMMs ==="
for g in "${GENES[@]}"; do
    gene="clb${g}"
    msa="${gene}.msa"
    hmm="${gene}.hmm"
    if [[ ! -s "${msa}" ]]; then
        echo "  SKIP ${gene}: no alignment"
        continue
    fi
    if [[ -s "${hmm}" ]]; then
        echo "  ${gene}: already built"
        continue
    fi
    echo "  Building ${gene}.hmm..."
    hmmbuild -n "${gene}" "${hmm}" "${msa}" > "${gene}.hmmbuild.log" 2>&1
done

echo ""
echo "=== Concatenating HMMs ==="
cat clb{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S}.hmm > clb_all_protein.hmm 2>/dev/null || {
    # Concatenate whatever was built
    cat clb*.hmm > clb_all_protein.hmm
}
n_models=$(grep -c '^HMMER3' clb_all_protein.hmm)
echo "  ${n_models} HMM models in concatenated file"

echo ""
echo "=== Pressing (hmmpress) ==="
hmmpress -f clb_all_protein.hmm

echo ""
echo "=== Validation (hmmstat) ==="
hmmstat clb_all_protein.hmm

echo ""
echo "=== Copying to ref/hmm/ ==="
cp clb_all_protein.hmm "${OUTHMM}"
cp clb_all_protein.hmm.h3f "${OUTHMM}.h3f"
cp clb_all_protein.hmm.h3i "${OUTHMM}.h3i"
cp clb_all_protein.hmm.h3m "${OUTHMM}.h3m"
cp clb_all_protein.hmm.h3p "${OUTHMM}.h3p"

echo ""
echo "Done: ${OUTHMM}"
echo "Models: ${n_models}"
