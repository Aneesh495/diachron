# Note 3: Cognacy and loans

Cognacy is a judgment, not an observation, and the platform treats it as such. A
cognate judgment carries a confidence and an annotator, and conflicting judgments
for one lexeme are retained rather than resolved by fiat. Orthographic similarity
is never equated with cognacy; the transliteration provenance is kept precisely so
that a form's spelling is not mistaken for evidence of shared descent.

Loans are the sharp case. A borrowed form resembles an inherited one but does not
testify to common descent, so folding it into the inherited character would create
false shared presence and distort the tree. The transform annotates loans and, in
the borrowing sensitivity configuration, removes them from the inherited
characters entirely; the borrowing-aware model instead routes them through a
separate contact channel on a bounded candidate network. Partial cognates, which
share a morpheme but not the whole form, are kept as their own classes rather than
merged into a full cognate, so partial sharing is not overclaimed.

The honest limit here is that the platform consumes cognate judgments; it does not
produce them. A dataset whose cognate coding is systematically wrong will mislead
inference no matter how correct the likelihood is. This is why the evidence trail
back to the annotator and source matters: it makes the judgment, and therefore its
potential errors, visible to the researcher rather than buried in a binary matrix.
