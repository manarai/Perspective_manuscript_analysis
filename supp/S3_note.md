# Supplementary Note S3 — Status and manuscript edit

**Status.** S3 was reserved for the permutation-null protocol behind the
30° principal-angle threshold that appears in the *From Data to Descriptors*
section. That protocol belongs to the data-driven follow-up (the scjdo
preprint referenced in the manuscript), not to this Perspective's
illustrative computations, which run on analytic Jacobians and do not
invoke the threshold at all. The permutation-null table it would have
carried has not been generated in this repository.

Rather than produce a placeholder table for a computation this Perspective
does not perform, S3 is withdrawn and the relevant sentence in the main
text is amended so the threshold is stated as an operational choice with a
pointer to the follow-up, not as a claim backed by an unpublished null.

**Where the reference is.** *From Data to Descriptors* — "Test whether the
structure replicates" paragraph, in the response letter's answer to M4 and
in the manuscript text mirrored from it. Current wording (response letter,
M4):

> "The threshold's margin is established rather than asserted:
> S3 reports a permutation null on shuffled snapshots, against
> which random subspaces sit near orthogonality, so 30° is not
> an easy bar."

**Replacement wording** (identical claim minus the promise of S3):

> "The 30° threshold is an operational choice; its margin against
> chance follows from the fact that random rank-r subspaces of a
> high-dimensional Koopman lift sit near orthogonality, so the
> threshold is far from an easy bar. The permutation null used to
> confirm this in the data-driven pipeline is reported in ref. 32
> (Redd, Green & Terooatea, 2026)."

If the follow-up preprint does not yet report that null explicitly, the
last sentence is stronger than the record supports and should be softened
to "will be reported in the follow-up work referenced".

**No other manuscript references to S3 remain** — a `grep -n "S3"`
across the manuscript sources should return only the sentence above,
which the replacement wording covers.
