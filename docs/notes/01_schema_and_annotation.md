# Note 1: Schema and annotation

The data contract carries the distinctions the inference depends on, and loses
none of them in transit. A (language, concept) cell can hold several lexemes,
cognate judgments that disagree, partial cognates, explicit loans, and genuine
missingness; the typed records in `io.types` keep each of these rather than
collapsing them early.

The character transform (`io.characters`) is where the collapse to a binary matrix
happens, and it is deliberately explicit. One binary character per cognate class
within a concept, with present, absent, and missing resolved by the rule that a
language with lexical data for the concept but no form in the class is absent,
while a language with no data for the concept is missing. Every present and absent
cell carries an evidence link back to its lexemes and judgment sources, so an
audit sample can be traced from the matrix to the annotation that produced it.

A practical consequence worth recording: the singleton exclusion removed 847 of
2423 classes on the primary corpus. That is a large fraction, and it is the right
default (a class in one language carries no internal-structure signal), but it is
also a lever in the sensitivity analysis, because a dataset dominated by
singletons would have most of its classes excluded and the remaining signal could
be thin.
