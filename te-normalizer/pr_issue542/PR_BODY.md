# Add missing TE name aliases reported in #542

Follow-up to #542. Of the 73 names reported in that thread, 12 were already
resolvable in `TE_Sequence_Ontology.txt`. This PR adds aliases for 21 of the
remaining 61 — those that are structural variants of superfamilies already present in
the table, plus `TIR/Novosib`, whose target row already exists.

## What changed

- **Only the third column (alias list) of 11 existing rows is extended.**
  No new SO entry line is added, and no SO name or SO ID is modified.
- Therefore `bin/gff2bed.pl` needs **no** corresponding update: as noted in the
  header of `TE_Sequence_Ontology.txt`, that file only has to change when a *new*
  SO entry line is introduced.
- Total: 21 aliases.

## Deliberately kept out: clade-level names (30)

In #542 you noted that clade-level classifications are *"not the convention of
RepeatMasker library naming"* and suggested stripping them from the library before
feeding it to EDTA2. This PR therefore leaves them untouched: names that carry a clade
segment (`LTR/Copia/Ale`, `LTR/Gypsy/Athila`, ...) continue to behave exactly as they do
today.

If you would rather have the table discard the clade segment automatically, that change
affects 30 names across two rows (`Copia_LTR_retrotransposon`,
`Gypsy_LTR_retrotransposon`). I can open it as a separate PR so the two decisions stay
independent.

## Aliases added, grouped by target term

| Target SO name | SO ID | Aliases added |
|---|---|---|
| `Novosib_TIR_transposon` | SO:0000182 | 1: `TIR/Novosib` |
| `repeat_fragment` | SO:0001050 | 1: `Unspecified/NA` |
| `MITE` | SO:0000338 | 1: `MITE/NA` |
| `CACTA_TIR_transposon` | SO:0002285 | 3: `DNA/CACTA/MITE`, `DNA/CACTA/nMITE`, `DNA/CACTA/unknown` |
| `hAT_TIR_transposon` | SO:0002279 | 3: `DNA/hAT/MITE`, `DNA/hAT/nMITE`, `DNA/hAT/unknown` |
| `Mutator_TIR_transposon` | SO:0002280 | 3: `DNA/Mutator/MITE`, `DNA/Mutator/nMITE`, `DNA/Mutator/unknown` |
| `PIF_Harbinger_TIR_transposon` | SO:0002284 | 3: `DNA/Harbinger/MITE`, `DNA/Harbinger/nMITE`, `DNA/Harbinger/unknown` |
| `Tc1_Mariner_TIR_transposon` | SO:0002278 | 3: `DNA/TcMar/MITE`, `DNA/TcMar/nMITE`, `DNA/TcMar/unknown` |
| `P_TIR_transposon` | SO:0001535 | 1: `DNA/P/nMITE` |
| `LTR_retrotransposon` | SO:0000186 | 1: `LTR/mixture` |
| `Penelope_retrotransposon` | SO:0002290 | 1: `PLE/NA` |

## Deliberately **not** included (6 names)

- `ARTEFACT` — per your reply in #542 (2026-03-02): *"ARTEFACT is ambiguous to me, and I am unsure which TE classification fits it better than repeat_fragment, so the current default is probably good for now."* The existing `repeat_fragment` fallback is therefore left as-is.
- `LTR/PIF_Harbinger` — `PIF_Harbinger` is a TIR superfamily; the `LTR/` prefix looks like a misclassification on the reporting side. Encoding it as an alias would legitimise a wrong call and skew LTR vs TIR counts.
- `mixture` — a TEsorter marker for chimeric elements rather than a taxonomic name. Mapping it to `repeat_fragment` would silently flatten every chimeric call into the fallback class and hide its real affiliation.
- `mixture/NA` — same as `mixture`.
- `nLTR` — ambiguous between *non-LTR* and *nested LTR*, which sit in completely different branches of the ontology. Worth confirming with the reporter before assigning a target.
- `nLTR/NA` — same ambiguity as `nLTR`.

`TIR/Novosib` is included because a `Novosib_TIR_transposon` row already exists; only
the `TIR/` spelling was missing.

Note that `LTR/mixture` is included while the bare `mixture` is not: `LTR/mixture`
still names the class (an LTR retrotransposon whose clade is indeterminate), so what
remains after stripping is well defined, whereas `mixture` on its own carries no class
information to fall back on.

## Not included here: 4 judgement calls

`nMITE`, `nMITE/NA`, `Unknown//ClassI` and `Unknown//ClassII` have defensible targets,
but assigning them rests on a judgement (e.g. whether the `n` prefix marks a
non-autonomous copy, or which class an unspecified element belongs to) rather than on a
mechanical strip. They are left out of this PR; I can add them, or explain the reasoning
for each, if you find it useful.

## Verification

- Every added alias was checked to be absent from its target row beforehand
  (case-insensitive), so no duplicate is introduced.
- After patching, all 21 names resolve to an SO term via the alias index, and the
  record count is unchanged (98) — i.e. no SO entry was added.
- Tab separators, comma-without-space alias formatting, LF line endings and the absence
  of a BOM match the original file outside the edited cells.

A per-name breakdown of all 61 reported names (what resolves, what does not, and why) is
available on request.
