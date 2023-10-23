# Note 11: Historical case studies

Five case studies exercise the parts of the platform that pure statistics do not
reach. Each is framed as a question the platform can answer with uncertainty, not
as a settled fact. On the labeled-synthetic corpora the ground truth is known,
which lets each case be a controlled test; on real CLDF collections the same
machinery runs without the ground truth.

A documented loan set. When a loan is injected between geographic neighbors, the
borrowing-aware model's intensity on the corresponding candidate edge rises above
its prior, while the tree-only model instead shows reduced support for the clade
that the loan's shared presence would otherwise suggest. The case shows the
borrowing channel absorbing signal that would otherwise corrupt the tree.

A contested calibration. Running the dated model under the full, half, and absent
calibration and comparing the root-age posteriors makes the calibration dependence
concrete for a single family rather than in the aggregate: the posterior widens and
shifts as the calibration weakens.

An extinct-language tip. A dated tip placed at a positive age tests the time
machinery: the tip's age enters as a constraint and the lineage above it is dated
relative to it, which the sampled-ancestor-free time model handles because the tip
is a genuine node, not a contemporary leaf.

A sparse lexical variety. A language with most concepts missing contributes mostly
non-informative missing cells; the case confirms that such a variety is placed
with wide uncertainty rather than being confidently but wrongly resolved.

Incompatible catalog classifications. Where the posterior tree disagrees with a
catalog classification, the comparison view surfaces the disagreement rather than
hiding it, so a catalog grouping that the data do not support is visible as such.

## Audited annotation: Brittonic fruit terms

The pinned IE-CoR release records the Brittonic group in cognate class 4364
as loan-related to class 227 for the same meaning. The class entry attributes
this analysis to a borrowing from Latin and points to Albert Deshayes,
*Dictionnaire étymologique du breton* (2003), p. 258, alongside the Welsh
dictionary entry. The donor class serves as a reference lineage in the
wordlist. This audit follows the published judgment and does not date the
transfer. [IE-CoR class 4364](https://iecor.clld.org/cognatesets/4364),
[IE-CoR class 227](https://iecor.clld.org/cognatesets/227).

The frozen sample has seven eligible recipient-class varieties and twenty
eligible varieties with the source class. The two classes enter as separate
binary characters. Removing loan-tagged forms removes the recipient class
and changes some absences to missing data where only loan-tagged
classification had justified an absence. The source class retains its
twenty presences. The original lexeme-to-class records and the directed
loan-source annotation remain independently inspectable.

The proximity screen proposes no directed edge between these source and
recipient sets. This is a concrete limit of that contact network. Moreover,
a copy of one binary state can only connect observations under a shared
class. Directly testing this historical annotation requires an explicit
cross-class borrowing observation model. The current campaign leaves that
comparison open.

The numeric audit is in [iecor-loan-fruit.json](../results/iecor-loan-fruit.json).
