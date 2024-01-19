(* ::Package:: *)

(* # Single Higgs production using gluon fusion
 *
 *)

(* To start off, for interactive development it is convenient
 * to reduce the width of Mathematica formatted output to 65
 * characters, and force it to clear the history (the `Out[]`
 * variables) so that old expressions will not linger on in
 * the memory.
 *)

(*SetOptions[$Output, PageWidth -> 65];
$HistoryLength = 2;*)

(* Load *alibrary* and the [[QCD Feynman rules]] from `amodel-qcd-HHH.m`.
 *)

Get["PathToAlibrary"];
Get[$Apath <> "/amodel-qcd-HHH.m"];

(* Define the number of loops and the kinematics
 *)

NLOOPS = 1;
ggHsprules = {sp[q1, q1] -> 0, sp[q2, q2] -> 0, sp[q1, q2] -> s12/2};

(* ## Diagrams
 *
 * Generate the diagrams with two incoming gluons
 * and one outgoing Higgs boson using QGraf.
 * We immediately require the conservation of 4-momentum.
 *)

diagrams = Diagrams[{"g", "g"}, {"H"}, NLOOPS] /. q -> q1 + q2;
Print["Loaded ", diagrams//Length, " diagrams"];

(* ## Tensor projectors
 *
 * Because the gluon propagator has open Lorentz
 * indices corresponding to the incoming and outgoing particles,
 * we need to project it to scalar values somehow.
 *
 * If $H^{\mu\nu}$ is the amplitude, and $q$ is the momentum
 * flowing through it, there are two tensor structures that it
 * can be made of:
 *
 * $$ H^{\mu\nu} = A g^{\mu\nu} + B q^\mu q^\nu. $$
 *
 * Because of the Ward identities,
 *
 * $$ q_\mu H^{\mu\nu} = H^{\mu\nu} q_\nu = 0, $$
 *
 * this form is further constrained to just
 *
 * $$ H^{\mu\nu} = H ( g^{\mu\nu} - q^\mu q^\nu / q^2 ). $$
 *
 * Then to get the scalar $H$ from $H^{\mu\nu}$ we can construct
 * a projector as
 *
 * $$ H = H^{\mu\nu} g_{\mu\nu} / (d-2). $$
 *
 * Here it is in Mathematica notation.
 *)
 
projector = delta[adj[-1],adj[-3]]*(delta[lor[-1], lor[-3]] - momentum[q1, lor[-3]]*momentum[q2, lor[-1]]/sp[q1, q2]) / (d-2);
 
 (* Apply Feynman rules to get amplitudes out of diagrams.
 *)

amplitudes = diagrams // Map[Amplitude[#] * projector&];

(* ## Zero integral cleanup
 * 
 * Cleanup scaleless integrals. Some of these show up as propagators
 * with zero momentum, which means that a part of the graph is
 * disconnected from the rest, and thus scaleless. We can set
 * these to zero immediately.
 *)

amplitudes2 = amplitudes /. den[0] -> 0 /. momentum[0,_] -> 0;
Print["Non-zero amplitudes: ", amplitudes2//Count[Except[0]], " of ", amplitudes2//Length];

(* In this particular example there is a set of diagrams that are
 * zero by the color factors. For example, those with subdiagrams
 * where the gluons fuse into a single gluon. In principle we
 * could try to skip these during diagram generation, but we
 * don\[CloseCurlyQuote]t need to. Let's compute color factors instead, and see
 * what turns to zero.
 *)

amplitudes3 = amplitudes2 // RunThroughForm[{ "#call colorsum\n" }];
Print["Non-zero amplitudes: ", amplitudes3//Count[Except[0]], " of ", amplitudes3//Length];

(* ## Integral family construction
 *
 * Next we want to define the integral families onto which we
 * shall map the integrals, and which will be used in the IBP
 * reduction later.
 *
 * To this end, start with the set of denominators per diagram.
 * We also extract the set of numerators or constant denominators that might occur.
 *)

loopmomenta = diagrams // CaseUnion[l1];
externalmomenta = diagrams // CaseUnion[q1|q2];
Print["External momenta: ", externalmomenta];
Print["Loop momenta: ", loopmomenta];
FailUnless[Length[loopmomenta] === NLOOPS];

constantdenominators =(NormalizeDens[amplitudes3] // Map[CaseUnion[_den] /* Select[
  FreeQ[Alternatives @@ loopmomenta]]] // Flatten // Union) \
  /. den[mom_, masssq_] :> den[mom, masssq] -> 1 / (sp[mom, mom] - masssq) //. ggHsprules;

denominatorsets = amplitudes3 // NormalizeDens // Map[
  CaseUnion[_den] /* Select[NotFreeQ[Alternatives@@loopmomenta]]
];
Print["Unique denominator sets: ", denominatorsets // DeleteCases[{}] // Union // Length];

(* In principle we could define the integral families by the
 * denominator sets above, one family per denominator set. This
 * is not very efficient though, as there are symmetries between
 * those families. It\[CloseCurlyQuote]s best to first eliminate denominator sets
 * that are symmetric to others.
 *
 * The symmetries manifest most directly in the Feynman parameter
 * space, as permutations of the parameters. In the momenta space
 * this corresponds to loop momenta shifts, and we would like
 * to have a set of momenta shifts that would make symmetric
 * families explicitly identical, or identical to subsets of
 * bigger families, so we could test if a family is symmetric by
 * just asking if the set of denominators a subset of another
 * family.
 *
 * The tool to compute this momenta mapping is [Feynson], and
 * the interface to it is [[SymmetryMaps]].
 *
 * [feynson]: https://github.com/magv/feynson
 *)
momentamaps = SymmetryMaps[denominatorsets, loopmomenta];
Print["Found ", momentamaps // DeleteCases[{}] // Length, " momenta mappings"];

symmetrizeddenominatorsets =
  MapThread[ReplaceAll, {denominatorsets, momentamaps}] //
  NormalizeDens;
  
(* Then, the set of unique supersets of the denominator sets is
 * the set of integral families we need.
 *)

{denominatorsupersets, supersetindices} =
  UniqueSupertopologyMapping[symmetrizeddenominatorsets];
Print["Total integral families: ", denominatorsupersets//Length];

(* Let us then construct the IBP basis objects for each unique
 * denominator superset. These objects are just associations
 * storing denominators, and maps from scalar products into the
 * denominators.
 *
 * Also in the case when the denominator set is not sufficient
 * to form the full linear basis of scalar products, we want to
 * complete it; [[CompleteIBPBasis]] will do this for us.
 *)

bases = denominatorsupersets //
  MapIndexed[CompleteIBPBasis[
    First[#2], #1 // NormalizeDens // Sort, loopmomenta, externalmomenta, ggHsprules]&];
  
(* ## Amplitude conversion
 * 
 * OK, now that we have the IBP bases, we can convert the
 * amplitudes to them.
 *
 * One practical thing to start with is to identify the set of
 * sectors (integral family subsets) that correspond to scaleless
 * integrals. This is also done with [Feynson].
 *)

zerosectors = ZeroSectors[bases];

(* Next, just call FORM to do all the tensor summation and
 * conversion to IBP families. Here, we use our list of numerators
 * to simplify the ampiltude beforehand.
 *)

amplitudesB =
  MapThread[ReplaceAll, {amplitudes3, momentamaps}] /. constantdenominators //
  # * BID^supersetindices & //
  RunThroughForm[{
    "#call contractmomenta\n",
    "#call sort(after-contractmomenta)\n",
    "#call chaincolorT\n",
    "#call chaingammachain\n",
    "#call flavorsumwithcharge\n",
    "#call colorsum\n",
    "#call sort(after-colorsum)\n",
    "#call polarizationsum\n",
    "#call spinsum\n",
    "#call diractrace\n",
    "#call contractmomenta\n",
    FormCallToB[bases],
    "id mt1^2 = mt2;\n",
    FormCallZeroSectors[zerosectors]
  }] //
  MapWithProgress[FasterFactor];

FailUnless[FreeQ[amplitudesB, l1]];

(* Mathematica and alibrary tend to be lazy when it comes to applying
 * the kinematic rules. So we force it to repeatedly replace 
 * all occuring scalar products until none are left.
 * We also require our expression to be devoid of any explicit Higgs masses,
 * since they can be expressed using our invariant s12.
 *)
 
amplitudesB = amplitudesB //. Union[ggHsprules, {mH2 -> s12}];

FailUnless[FreeQ[amplitudesB, sp]];
FailUnless[FreeQ[amplitudesB, mH2]];

(* ## IBP reduction
 * 
 * Next, lets do the IBP reduction.
 *
 * We'll use [[KiraIBP]], a simple interface to IBP with [Kira].
 * It is probably too simplistic to work automatically for larger
 * examples, but for this problem it\[CloseCurlyQuote]s ideal.
 *
 * [kira]: https://kira.hepforge.org/
 *)

amplitudesBibp = amplitudesB // KiraIBP[bases, ReplaceByOne->mt2];

(* Here we use ReplaceByOne to speed up the IBP reduction. Note
 * that Kira can restore the powers of mt2 after the reduction
 * is over, but this is fairly slow at the moment. For this
 * reason we'll skip that, and update the rest of the amplitude
 * to this convention too.
 *)

amplitudesBibp = amplitudesBibp // ReplaceAll[mt2->1];

(* The full amplitude is just the sum of the diagram amplitudes.
 *)

fullAmplitude = amplitudesBibp // Apply[Plus] // Bracket[#, _B, Factor]&;

(* A good correctness check is to see if there is any Xi dependence
 * left. None should remain.
 *)

FailUnless[FreeQ[fullAmplitude, Xi]];

(* Before continuing, we collect all global coefficients and constants.
 *
 *)

fullAmplitudeByPrefactor = fullAmplitude//ReplaceAll[Complex[re_,im_]:>re+im*ImagI]//BracketAssociation[Ca|Cf|Na|Tf|d33|d44|Nc|Nf|Nft|gH|gs|Xi|ImagI|_flvsum|_flvsumt];

(* Now we have reduced the amplitude to master integrals.
 *
 * ## Numerical evaluation
 *
 * The final step is to insert the values of the masters. Of
 * course the masters here are known analytically, but as an
 * example let us evaluate them numerically with [pySecDec],
 * each up to $\varepsilon^0$.
 *
 * [pySecDec]: https://github.com/gudrunhe/secdec
 *)

masters = fullAmplitude  // CaseUnion[_B];
Print["Master integrals: ", masters // Length];

basesWithoutMt2 = bases//Map[(Append[#,<|"invariants"->DeleteCases[#["invariants"],mt2]|>]&)/*Map[ReplaceAll[mt2->1]]];

SecDecPrepareSum["PathToOutputDirectory", basesWithoutMt2, fullAmplitudeByPrefactor//KeyMap[InputForm/*ToString]//Map[ReplaceAll[d->4-2*eps]],Order->0];
Print["--- Done ---"];
