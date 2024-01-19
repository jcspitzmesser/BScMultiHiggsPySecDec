(* ::Package:: *)

(* # Double Higgs production using gluon fusion
 *
 *)

(* To start off, for interactive development it is convenient
 * to reduce the width of Mathematica formatted output to 65
 * characters, and force it to clear the history (the `Out[]`
 * variables) so that old expressions will not linger on in
 * the memory.
 *)

SetOptions[$Output, PageWidth -> 65];
$HistoryLength = 2;

(* Load *alibrary* and the [[QCD Feynman rules]] from `amodel-qcd-HHH.m`.
 *)

Get["PathToAlibray"];
Get[$Apath <> "/amodel-qcd-HHH.m"];

(* Define the number of loops and kinematics
 *)

NLOOPS = 1; 
ggHHsprules = {sp[q1, q1] -> 0, sp[q2, q2] -> 0, sp[q1, q2] -> s12/2, sp[p1, p1] -> s12/2 + s13/2 + s23/2, sp[q1, p1] -> s12/4 - s13/4 + s23/4, sp[q2, p1] -> s12/4 + s13/4 - s23/4};

(* ## Diagrams
 *
 * Generate the diagrams with two incoming gluons
 * and two outgoing Higgs bosons using QGraf.
 * We immediately require the conservation of 4-momentum.
 *)

diagrams = Diagrams[{"g", "g"}, {"H", "H"}, NLOOPS]/. p2 -> -p1 + q1 + q2;
Print["Loaded ", diagrams//Length, " diagrams"];

(* ## Tensor projectors
 *
 * Because the gluon propagator has open Lorentz
 * indices corresponding to the incoming and outgoing particles,
 * we need to project it to scalar values somehow.
 *
 * The necessary projectors have been taken from
 * https://arxiv.org/pdf/1608.04798.pdf
 * and have been slightly simplified using the kinematics of this process.
 *)

T1[p1_, p2_, \[Mu]_, \[Nu]_] := delta[lor[\[Mu]], lor[\[Nu]]] - momentum[p1, lor[\[Nu]]]*momentum[p2, lor[\[Mu]]]/sp[p1, p2];
T2[p1_, p2_, p3_, \[Mu]_, \[Nu]_] := delta[lor[\[Mu]], lor[\[Nu]]] + (mH2*momentum[p1, lor[\[Nu]]]*momentum[p2, lor[\[Mu]]] - 2*sp[p1, p3]*momentum[p3, lor[\[Nu]]]*momentum[p2, lor[\[Mu]]] - 2*sp[p2, p3]*momentum[p3, lor[\[Mu]]]*momentum[p1, lor[\[Nu]]] + 2*sp[p1, p2]*momentum[p3,lor[\[Nu]]]*momentum[p3, lor[\[Mu]]])\
							/(2*sp[p1, p3]*sp[p2, p3] - mH2*sp[p1, p2]);
   
projector1 = delta[adj[-1],adj[-3]]*(T1[q1, q2, -1, -3]*(d - 2)/(d - 3) - T2[q1, q2, p1, -1, -3]*(d - 4)/(d - 3)) / 4;
projector2 = delta[adj[-1],adj[-3]]*(T2[q1, q2, p1, -1, -3]*(d - 2)/(d - 3)-T1[q1, q2, -1, -3]*(d - 4)/(d - 3)) / 4;
 
 (* Apply Feynman rules to get amplitudes out of diagrams.
 *)

amplitudes11 = diagrams // Map[Amplitude[#] * projector1&];
amplitudes12 = diagrams // Map[Amplitude[#] * projector2&];

(* ## Zero integral cleanup
 * 
 * Cleanup scaleless integrals. Some of these show up as propagators
 * with zero momentum, which means that a part of the graph is
 * disconnected from the rest, and thus scaleless. We can set
 * these to zero immediately.
 *)

amplitudes21 = amplitudes11 /. den[0] -> 0 /. momentum[0,_] -> 0;
Print["Non-zero amplitudes: ", amplitudes21//Count[Except[0]], " of ", amplitudes21//Length];
amplitudes22 = amplitudes12 /. den[0] -> 0 /. momentum[0,_] -> 0;
Print["Non-zero amplitudes: ", amplitudes22//Count[Except[0]], " of ", amplitudes22//Length];

(* In this particular example there is a set of diagrams that are
 * zero by the color factors. For example, those with subdiagrams
 * where the gluons fuse into a single gluon. In principle we
 * could try to skip these during diagram generation, but we
 * don\[CloseCurlyQuote]t need to. Let's compute color factors instead, and see
 * what turns to zero.
 *)

amplitudes31 = amplitudes21 // RunThroughForm[{ "#call colorsum\n" }];
Print["Non-zero amplitudes: ", amplitudes31//Count[Except[0]], " of ", amplitudes31//Length];
amplitudes32 = amplitudes22 // RunThroughForm[{ "#call colorsum\n" }];
Print["Non-zero amplitudes: ", amplitudes32//Count[Except[0]], " of ", amplitudes32//Length];

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
externalmomenta = diagrams // CaseUnion[q1|q2|p1|p2];
Print["External momenta: ", externalmomenta];
Print["Loop momenta: ", loopmomenta];
FailUnless[Length[loopmomenta] === NLOOPS];

constantdenominators1 = (NormalizeDens[amplitudes31] // Map[CaseUnion[_den] /* Select[FreeQ[Alternatives @@ loopmomenta]]] // 
    Flatten // Union) /. den[mom_, masssq_] :> den[mom, masssq] -> 1 / (sp[mom, mom] - masssq)//.ggHHsprules;
denominatorsets1 = amplitudes31 // NormalizeDens // Map[
  CaseUnion[_den] /* Select[NotFreeQ[Alternatives@@loopmomenta]]
];
Print["Unique denominator sets: ", denominatorsets1 // DeleteCases[{}] // Union // Length];

constantdenominators2 = (NormalizeDens[amplitudes32] // Map[CaseUnion[_den] /* Select[FreeQ[Alternatives @@ loopmomenta]]] // 
    Flatten // Union) /. den[mom_, masssq_] :> den[mom, masssq] -> 1 / (sp[mom, mom] - masssq)//.ggHHsprules;
denominatorsets2 = amplitudes32 // NormalizeDens // Map[
  CaseUnion[_den] /* Select[NotFreeQ[Alternatives@@loopmomenta]]
];
Print["Unique denominator sets: ", denominatorsets2 // DeleteCases[{}] // Union // Length];

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
 
momentamaps1 = SymmetryMaps[denominatorsets1, loopmomenta];
Print["Found ", momentamaps1 // DeleteCases[{}] // Length, " momenta mappings"];

momentamaps2 = SymmetryMaps[denominatorsets2, loopmomenta];
Print["Found ", momentamaps2 // DeleteCases[{}] // Length, " momenta mappings"];

symmetrizeddenominatorsets1 =
  MapThread[ReplaceAll, {denominatorsets1, momentamaps1}] //
  NormalizeDens;

symmetrizeddenominatorsets2 =
  MapThread[ReplaceAll, {denominatorsets2, momentamaps2}] //
  NormalizeDens;
  
(* Then, the set of unique supersets of the denominator sets is
 * the set of integral families we need.
 *)

{denominatorsupersets1, supersetindices1} =
  UniqueSupertopologyMapping[symmetrizeddenominatorsets1];
Print["Total integral families: ", denominatorsupersets1 // Length];

{denominatorsupersets2, supersetindices2} =
  UniqueSupertopologyMapping[symmetrizeddenominatorsets2];
Print["Total integral families: ", denominatorsupersets2 // Length];

(* Let us then construct the IBP basis objects for each unique
 * denominator superset. These objects are just associations
 * storing denominators, and maps from scalar products into the
 * denominators.
 *
 * Also in the case when the denominator set is not sufficient
 * to form the full linear basis of scalar products, we want to
 * complete it; [[CompleteIBPBasis]] will do this for us.
 *)

bases1 = denominatorsupersets1 //
  MapIndexed[CompleteIBPBasis[
    First[#2], #1 // NormalizeDens // Sort, loopmomenta, externalmomenta, ggHHsprules]&];

bases2 = denominatorsupersets2 //
  MapIndexed[CompleteIBPBasis[
    First[#2], #1 // NormalizeDens // Sort, loopmomenta, externalmomenta, ggHHsprules]&];
  
(* ## Amplitude conversion
 * 
 * OK, now that we have the IBP bases, we can convert the
 * amplitudes to them.
 *
 * One practical thing to start with is to identify the set of
 * sectors (integral family subsets) that correspond to scaleless
 * integrals. This is also done with [Feynson].
 *)

zerosectors1 = ZeroSectors[bases1];
zerosectors2 = ZeroSectors[bases2];

(* Next, just call FORM to do all the tensor summation and
 * conversion to IBP families. Here, we use our list of numerators
 * to simplify the ampiltude beforehand.
 *)
 
amplitudesB1 =
  MapThread[ReplaceAll, {amplitudes31, momentamaps1}] /. constantdenominators1 //
  # * BID^supersetindices1 & //
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
    FormCallToB[bases1],
    "id mt1^2 = mt2;\n",
    FormCallZeroSectors[zerosectors1]
  }] //
  MapWithProgress[FasterFactor];

FailUnless[FreeQ[amplitudesB1, l1]];
 
amplitudesB2 =
  MapThread[ReplaceAll, {amplitudes32, momentamaps2}] /. constantdenominators2 //
  # * BID^supersetindices2 & //
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
    FormCallToB[bases2],
    "id mt1^2 = mt2;\n",
    FormCallZeroSectors[zerosectors2]
  }] //
  MapWithProgress[FasterFactor];

FailUnless[FreeQ[amplitudesB2, l1]];

(* alibrary has a tendency to ignore some of the scalar products that occur here.
 * So before we process the amplitude further, we substitute out all remaining scalar
 * products repeatedly until none are left.
 * At the same time, we also substitute mH2 using the expression found with kinematics.
 * These masses occur as a result of HHH vertices.
 * Furthermore, alibrary doesn't natively resolve scalar products of sums, so we also
 * need to simplify those. These occur as a result of numerators (constant denominators),
 * and we don't want them either in our final expression for our form factors.
 *)

amplitudesB1 = amplitudesB1 //. Union[ggHHsprules, {mH2 -> (s12 + s13 + s23) / 2, sp[q1 + q2, q1 + q2] -> 2*sp[q1, q2]}];
FailUnless[FreeQ[amplitudesB1, sp]];
FailUnless[FreeQ[amplitudesB1, mH2]];

amplitudesB2 = amplitudesB2 //. Union[ggHHsprules, {mH2 -> (s12 + s13 + s23) / 2, sp[q1 + q2, q1 + q2] -> 2*sp[q1, q2]}];
FailUnless[FreeQ[amplitudesB2, sp]];
FailUnless[FreeQ[amplitudesB2, mH2]];

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

amplitudesBibp1 = amplitudesB1 // KiraIBP[bases1, ReplaceByOne->mt2];

amplitudesBibp2 = amplitudesB2 // KiraIBP[bases2, ReplaceByOne->mt2];

(* Here we use ReplaceByOne to speed up the IBP reduction. Note
 * that Kira can restore the powers of mt2 after the reduction
 * is over, but this is fairly slow at the moment. For this
 * reason we'll skip that, and update the rest of the amplitude
 * to this convention too.
 *)

amplitudesBibp1 = amplitudesBibp1 // ReplaceAll[mt2->1];
amplitudesBibp2 = amplitudesBibp2 // ReplaceAll[mt2->1];

(* The full amplitude is just the sum of the diagram amplitudes.
 *)

fullAmplitude1 = amplitudesBibp1 // Apply[Plus] // Bracket[#, _B, Factor]&;
fullAmplitude2 = amplitudesBibp2 // Apply[Plus] // Bracket[#, _B, Factor]&;

(* A good correctness check is to see if there is any Xi dependence
 * left. None should remain.
 *)

FailUnless[FreeQ[fullAmplitude1, Xi]];
FailUnless[FreeQ[fullAmplitude2, Xi]];

(* Before continuing, we collect all global coefficients and constants.
 *
 *)

fullAmplitudeByPrefactor1 = fullAmplitude1 // ReplaceAll[Complex[re_,im_]:>re+im*ImagI]//BracketAssociation[Ca|Cf|Na|Tf|d33|d44|Nc|Nf|Nft|gH|gs|Xi|ImagI|_flvsum|_flvsumt];
fullAmplitudeByPrefactor2 = fullAmplitude2 // ReplaceAll[Complex[re_,im_]:>re+im*ImagI]//BracketAssociation[Ca|Cf|Na|Tf|d33|d44|Nc|Nf|Nft|gH|gs|Xi|ImagI|_flvsum|_flvsumt];

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

masters1 = fullAmplitude1  // CaseUnion[_B];
Print["Master integrals: ", masters1 // Length];

masters2 = fullAmplitude2  // CaseUnion[_B];
Print["Master integrals: ", masters2 // Length];

bases1WithoutMt2 = bases1 // Map[(Append[#,<|"invariants"->DeleteCases[#["invariants"],mt2]|>]&)/*Map[ReplaceAll[mt2->1]]];
bases2WithoutMt2 = bases2 // Map[(Append[#,<|"invariants"->DeleteCases[#["invariants"],mt2]|>]&)/*Map[ReplaceAll[mt2->1]]];

SecDecPrepareSum["PathToOutputDirectory1", bases1WithoutMt2, fullAmplitudeByPrefactor1//KeyMap[InputForm/*ToString]//Map[ReplaceAll[d->4-2*eps]],Order->0];
Print["--- Done ---"];

SecDecPrepareSum["PathToOutputDirectory2", bases2WithoutMt2, fullAmplitudeByPrefactor2//KeyMap[InputForm/*ToString]//Map[ReplaceAll[d->4-2*eps]],Order->0];
Print["--- Done ---"];
