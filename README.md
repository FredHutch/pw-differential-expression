# Differential Expression (PubWeb)
Workflow for performing differential expression analysis from gene abundance data

## Input Data

The input data for this workflow must consist of two files, a counts table and a
manifest table. The counts table contains the number of reads per gene across
all specimens, while the manifest contains annotations on those specimens which
are used to identify which groups of specimens will be compared.

In addition to the input files, the user must also specify a column of the
manifest which will be used to identify which specimens should be compared.

### Counts Table Format

The counts table must be formatted as either a TSV or CSV, with the file
extension either '.tsv[.gz]' or '.csv[.gz]', as appropriate.

The first column in the counts table must contain the gene ID. Every specimen
defined in the manifest must correspond to a column in the counts table.
There may be additional columns in the counts table, but they will be ignored
in this analysis.

### Manifest Table Format

The manifest table must be formatted as a CSV, with the first column containing
the specimen names. Every additional column defines the annotations used for
those specimens. It is perfectly acceptable to leave a value empty for a particular
piece of metadata, in which case the specimen will be ignored in any analysis
which uses that feature for comparisons.

### Defining Comparisons

Comparisons can be made between specimens using either continuous or categorical
variables. A column in the manifest is treated as a continuous variable when
all of the values are numeric, and it is treated as a categorical variable when
those values are strings.

The user can define the exact comparisons which are tested using a combination
of the parameters:

 - `comp_col`: Column used for comparison
 - `comp_ref`: For categorical comparisons, the label used for the reference group
 - `group_cols`: Comma separated list indicating any columns which contain batch information which will be normalized before testing the comparison group

 Note: Interaction terms are not currently supported.

Consider the following example:

---
| specimen | batch | genotype | age |
| - | - | - | - |
| A1 | A | wt | 16 |
| A2 | A | wt | 11 |
| A3 | A | nod2 | 14 |
| A4 | A | nod2 | 18 |
| B1 | B | wt | 10 |
| B2 | B | wt | 8 |
| B3 | B | nod2 | 6 |
| B4 | B | nod2 | 12 |
---

#### Example 1:

To simply compare the `wt` vs. `nod2` samples (while indicating that `wt` is
the reference group), the user would specify:

 - `comp_col`: `genotype`
 - `comp_ref`: `wt`

In the results, the fold change will indicate the average value of `nod2 / wt`
for all genes, while using a `comp_ref`: `nod2` would yield the inverse
ratio `wt / nod2`.

Note: Any time a categorical comparison is used, the reference value must
be specified to ensure complete clarity.

#### Example 2:

To test for changes in gene expression as a function of age, the user
would specify:

 - `comp_col`: `age`

Note that `comp_ref` is not needed for continuous (numeric) variables. In this
case, the results will indicate the average change in gene expression for 
every unit change in `age` across the dataset.

#### Example 3:

To test for changes in gene expression as a function of age while controlling
for any differences correlated with genotype, the user would specify:

 - `comp_col`: `age`
 - `group_cols`: `genotype`

#### Example 4:

To test for changes in gene expression as a function of genotype while
controlling for any differences correlated with age, the user would specify:

 - `comp_col`: `genotype`
 - `comp_ref`: `wt`
 - `group_cols`: `age`

#### Example 5:

In some cases, the user may want to control for multiple factors, such
as the `batch` indicated in the example above. This can be specified by
the user by listing those columns as a comma-separated list.

For example, to test for changes in gene expression as a function of
genotype while controlling for any differences correlated with age, while
also controlling for batch effects, the user would specify:

 - `comp_col`: `genotype`
 - `comp_ref`: `wt`
 - `group_cols`: `age,group`

