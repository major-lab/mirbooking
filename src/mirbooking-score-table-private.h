/*
 * Reusable constants across score table implementations
 */

#define R 1.987203611e-3
#define T 310.15

#define SEED_OFFSET 1
#define SEED_LENGTH 7

/*
 * Forward rate constant for the [E] + [S] -> [ES] reaction, which is
 * indistinguishable from seed types or supplementary bindings (Wee et al. 2012
 * and Salomon et al. 2015). We picked the value from the latest publication
 * for the seed-only setup.
 *
 * Reference:
 * Liang Meng Wee et al., “Argonaute Divides Its RNA Guide into Domains
 * with Distinct Functions and RNA-Binding Properties,” Cell 151, no. 5
 * (November 21, 2012): 1055–67, https://doi.org/10.1016/j.cell.2012.10.036.
 *
 * William E. Salomon et al., “Single-Molecule Imaging Reveals That Argonaute
 * Reshapes the Binding Properties of Its Nucleic Acid Guides,” Cell 162, no. 1
 * (July 2, 2015): 84–95, https://doi.org/10.1016/j.cell.2015.06.029.
 */
#define KF   2.4e-4 // pM^-1s^-1
#define KCAT 3.6e-2 // s^-1

/*
 * Half-life of a messenger RNA in presence of microRNA is ~2 hours.
 *
 * Reference:
 * Nadya Morozova et al., “Kinetic Signatures of MicroRNA Modes of Action,” RNA
 * 18, no. 9 (September 2012): 1635–55, https://doi.org/10.1261/rna.032284.112.
 */
#define KDEG 9.627e-5 // s^-1

/*
 * For the duplex: 'CUACCUC&GAGGUAG', ViennaRNA reports a binding energy of
 * -9.37 kcal/mol.
 *
 * Wee et al. measured a dissociation constant for a mouse Ago2 protein
 * carrying a guide miRNA with only the seed pairing of 26±2 pM, which
 * correspond to a free energy of -15.02 kcal/mol.
 *
 * On the other hand, Salomon et al. instead measured 15±2 pm, which correspond
 * to -15.36 kcal/mol.
 *
 * In addition, the seed setup experiment in both experiments had 'A', which
 * shoudld account for a -0.56 kcal/mol additional contribution (Schirle et al. 2015).
 *
 * Using RNAup, we folded the reporter with the seed-only sequence which
 * yielded a 0.466 kcal/mol penalty to "open" a 17 nucleotides window around
 * the seed.
 *
 * We thus the latest value and impute the -5.90 kcal/mol gap to AGO2 entropic
 * contribution.
 *
 * Reference:
 * Liang Meng Wee et al., “Argonaute Divides Its RNA Guide into
 * Domains with Distinct Functions and RNA-Binding Properties,” Cell 151, no. 5
 * (November 21, 2012): 1055–67, https://doi.org/10.1016/j.cell.2012.10.036.
 *
 * Nicole T Schirle et al., “Water-Mediated Recognition of T1-Adenosine Anchors
 * Argonaute2 to MicroRNA Targets,” ed. Phillip D Zamore, ELife 4 (September
 * 11, 2015): e07646, https://doi.org/10.7554/eLife.07646.
 */
#define AGO2_SCORE (-5.43f)

/*
 * AGO2 has a slight preference for sites starting with 'A' at position t1.
 *
 * Reference:
 * Nicole T Schirle et al., “Water-Mediated Recognition of T1-Adenosine
 * Anchors Argonaute2 to MicroRNA Targets,” ed. Phillip D Zamore, ELife 4
 * (September 11, 2015): e07646, https://doi.org/10.7554/eLife.07646.
 */
#define T1_ADENOSINE_SCORE (-0.56f)

/*
 * We allow a 'G' nucleation bulge at position t5.
 *
 * There's a penalty of 1.2 kcal/mol for this non-canonical motif that was
 * taken from the authors free energy estimates.
 *
 * Reference:
 * Sung Wook Chi, Gregory J. Hannon, and Robert B. Darnell, “An Alternative
 * Mode of MicroRNA Target Recognition,” Nature Structural & Molecular
 * Biology 19, no. 3 (March 2012): 321–27, https://doi.org/10.1038/nsmb.2230.
 */
#define G_BULGED_SEED_SCORE (1.2f)
