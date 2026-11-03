# GOAL
For each candidate projection to query, give it a query ID and classify it into different orthology relation categories between each reference and query CNEEs

# Filtering criteria
A lone ref/query edge (a ref CNEE with exactly one candidate projection, to a
query locus no other ref CNEE also lands on) is never filtered -- it is
always reported as one2one, with its own sequence_identity/gap_fraction
exposed as output columns (see Output) rather than gated on them, so
filtering can happen post-hoc if desired.

Anything with more than one candidate projection -- one2many, many2one, or
complex/many2many by initial graph shape -- goes through relative-score
filtering ("Modified relative scoring algo" in `spurious_hit_filtering.md`;
`stola.spurious_filter`, wired into
`stola.orthology.classify_orthology` as `filter_multi_projections`).
It works per query locus (the merged query span of the eligible
projections landing there), scoring every chain aligned at that locus once:
- Soft depth filter: of the chains overlapping the locus, only the top 20
  (`--soft-depth-limit`) by chain score that actually have aligned bases
  there are scored; a projection on any other chain is pruned as
  repeat-driven (`depth`).
- If only one chain aligns there, its projection is kept only if its
  surrounding flanks have real aligned-block coverage from ANY chain
  (`no_flank_support` otherwise).
- Otherwise, each chain's local score = sequence identity over the locus x
  aligned fraction over the locus plus 10kb flanks each side; a projection
  is kept if its chain's score is > 0.5 (`--relative-score-threshold`) of
  the best chain's at that locus (`low_score` otherwise). Several
  chains can pass at one locus, so many2one calls can survive filtering.

# Orthology categories:
## Ref CNEE categories
- one2one: CNEE single copy in both query and ref; passes filtering criteria
- one2many: single copy in ref, multicopy in query; query projections pass filtering criteria
- many2one: multiple ref CNEEs align to same query locus; ref alignments pass filtering criteria
- complex: complex cases where multiple ref CNEEs overlap multiple the same query loci
- unaligned: ref CNEE has NO chain overlaps
- probable loss: ref CNEE overlaps a chain but only in chain gap AND chain gap not a contig boundary
- assembly gap loss: ref CNEE overlaps a chain but only in chain gap and gap is a contig boundary
- uncertain loss: ref CNEE has only pruned projections

## Query CNEE categories
- one2one: CNEE single copy in both query and ref; passes filtering criteria
- one2many: single copy in ref, multicopy in query; query projections pass filtering criteria
- many2one: multiple ref CNEEs align to same query locus; ref alignments pass filtering criteria
- complex: complex cases where multiple ref CNEEs overlap multiple the same query loci
- pruned: query projection fails relative-score filtering

## Loss sub-categorization (implemented)
Whether a chain gap is a contig boundary is inferred from the query FASTA
alone (no AGP): scan for N-runs, and see whether one (or a small handful of)
length(s) clearly dominate the runs of at least `--min-gap-length` bases
(default 10; shorter runs are ambiguous bases, never gaps) -- see
`assembly_gaps.py`. Four outcomes:
- **high confidence** (one dominant length): only that length counts as a
  boundary marker.
- **medium confidence** (a small set of 2-3 lengths together dominate):
  any of those lengths counts as a boundary marker -- some assemblies
  legitimately use more than one fixed-length gap type (e.g. contig-gap vs
  scaffold-gap).
- **estimated** (at least `MIN_RUNS_FOR_INFERENCE` qualifying N-runs, but no
  dominant length): the assembly records estimated gap sizes, so every
  qualifying N-run counts as a boundary marker. (A hard-masked query genome would also land here,
  and over-call assembly_gap_loss.)
- **none** (no qualifying N-runs, or fewer than `MIN_RUNS_FOR_INFERENCE`):
  nothing counts as a boundary marker.

A lost ref whose chain gap's query-side span overlaps a boundary marker is
`assembly_gap_loss`; every other lost ref is `probable_loss` -- including
one at the very end of its chain (no aligned block on one side, so there is
no query-side gap span to check).

This sub-categorization (plus per-ref pruning counts and pre-pruning
category) is independent of multi-projection filtering and always runs (the
reference and query genomes are required inputs), written to the ref-CNEE TSV described in
output_files.md. It does not change out_orthology's own "lost"/"unaligned"
labels.

# Implementation
- Create network graph where:
  - each node = query or ref CNEE 
- Draw an edge between ref node to query node (i.e. each alignment)
  - Store edge attributes:
    - sequence_identity
    - fraction_gap
- Determine each ref CNEE's initial (pre-filtering) graph shape from its
  connected component: a lone ref/query edge needs no filtering at all; a
  ref with multiple query edges, multiple refs sharing one query edge, or a
  larger many-to-many component are all "multi-projection" and go through
  the relative-score filtering described above (see Filtering criteria) before
  the graph is re-clustered from the surviving projections.
- for each query node (in the post-filtering graph):
  - get connected ref nodes
  - use info re: connected nodes to assign orthology category and ID for that query node (i.e CNEE):
    - if connected to single ref node that has no other connections, then: "one2one"
      - set query ID as "q_[ref CNEE name]"
    - if connected single ref node and:
      - ref node connects multiple query nodes AND those query nodes only connect this reference
      - then: "one2many"
        - set query ID as "q_[ref DNEE name]_[n]"
    - if connects to multiple ref nodes and those ref nodes only connect this query node, then: "many2one"
    - for now, if none of these criteria met, then: "complex"
  - a ref left with no surviving query edge at all (but that did have some
    projection before filtering) becomes "uncertain_loss"
- for each reference node with 0 chain-overlap candidates:
  - set query CNEE ID to "NA"
  - set orthology category as: "missing"
- for each reference node with chain-overlap candidates but no aligned query bases:
  - set query CNEE ID to "NA"
  - set orthology category as: "lost"

# Output
Tab separated list:
- query CNEE w/ ID
- partner ref CNEE(s) as comma separated list
- orthology category

# Orthology flowchart

"Solid" match for ref vs query projection = sequence identity > 0.7 AND gap fraction < 0.1.
If genome FASTAs are not provided, sequence identity is unavailable and the
solid-match decision uses gap fraction only.

```mermaid
flowchart TD

A["For each ref CNEE"] --> B{"Does it overlap a chain?"}
B -- "No" --> C["Missing"]
B -- "Yes" --> D{"Does it have aligned query bases?"}
D -- "No" --> E["Lost"]
D -- "Yes" --> F{"Does it have >=1 solid query projection?"}
F -- "No" --> O["Uncertain Loss"]
F -- "Yes" --> G{"Exactly 1 solid query?"}
G -- "Yes" --> H{"Does the 1 query match other refs?"}
H -- "No" --> I["one2one"]
H -- "Yes" --> J{"Do those other refs match other queries?"}
J -- "No" --> K["many2one"]
J -- "Yes" --> L["complex"]
G -- "No" --> M{"Do all multiple query matches match only this ref?"}
M -- "Yes" --> N["one2many"]
M -- "No" --> L
```

# Orthology flowchart v3
Per-locus relative-score filtering (see Filtering criteria and
`spurious_hit_filtering.md`) runs on every ref CNEE in an ambiguous
(one2many/many2one/complex) component before the graph is re-clustered and
final categories are assigned from the surviving edges; a lone ref/query
edge is never filtered.

```mermaid
flowchart TD

A["For each ref CNEE"] --> B{"Does it overlap a chain?"}
B -- "No" --> C["Missing"]
B -- "Yes" --> D{"Does it have aligned query bases?"}
D -- "No" --> E["Lost"]
D -- "Yes" --> F{"Exactly 1 query projection pre-filtering?"}
F -- "Yes" --> G{"Does the 1 query match other refs?"}
G -- "No" --> H1["one2one"]
G -- "Yes" --> I{"Filter, then: do surviving other refs match other queries?"}
I -- "No surviving edges at all" --> O["Uncertain Loss"]
I -- "No" --> I1["many2one"]
I -- "Yes" --> I2["complex"]
F -- "No" --> J["Filter, then: do all surviving queries match only this ref?"]
J -- "No surviving edges at all" --> O
J -- "No" --> I2
J -- "Yes, more than one survives" --> K["one2many"]
J -- "Yes, exactly one survives" --> H1
```
