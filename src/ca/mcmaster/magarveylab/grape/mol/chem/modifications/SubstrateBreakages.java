package ca.mcmaster.magarveylab.grape.mol.chem.modifications;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.mcss.RMap;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment;
import ca.mcmaster.magarveylab.grape.mol.chem.MoleculeModifier;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class SubstrateBreakages {

	public static void process(MoleculeModifier mod) {
		breakPeptideBonds(mod);
		breakSugarGroups(mod);
		breakSulfateGroups(mod);
		processAAMethylations(mod);
		processAACHydroxyliations(mod);
		processHalogenations(mod);
		processSecondaryAmineExtensions(mod);
	}

	/**
	 * Break peptide bonds. This method modifies monomerFragments and should be
	 * run after thiazole, oxazole, and lactone processing.
	 * 
	 * @throws CDKException
	 */
	private static void breakPeptideBonds(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		// The resulting List is ordered C to N
		// For each nrp fragment, identify all peptide bonds/beta
		// carbons/nitrogens. Then add these to the queue.

		@SuppressWarnings("unused")
		int numStandardPeptideBonds = 0;
		@SuppressWarnings("unused")
		int numNonStandardPeptideBonds = 0;
		List<IBond> peptideBonds = new ArrayList<IBond>();

		List<Set<IAtom>> smallRings = new ArrayList<Set<IAtom>>();
		List<IAtom> atomsInRings = new ArrayList<IAtom>();
		
		// Prepare peptide template 1 (amino acids with C terminus) and
					// peptide template 2 (amino acids with a nitrogen attached to C
					// terminus)
					IAtomContainer[] peptideTemplates = new IAtomContainer[5];
					IBond[] peptideTemplateBonds = new IBond[5];
					try {
						peptideTemplates[0] = SmilesIO.readSmilesTemplates("CCNC(C)=O");
						peptideTemplates[1] = SmilesIO.readSmilesTemplates("CC(=O)NC=C");
						peptideTemplates[2] = SmilesIO.readSmilesTemplates("CC\\N=C(/C)O");
						peptideTemplates[3] = SmilesIO.readSmilesTemplates("CCNC(C)OC");
						peptideTemplates[4] = SmilesIO.readSmilesTemplates("C\\C(O)=N/C=C");
						// peptideTemplates[5] = SmilesIO.readSmiles("CC(=O)NC=C");
						// peptideTemplates.add(SmilesIO.readSmiles("CC(=O)NCC=O"));
						// peptideTemplates.add(SmilesIO.readSmiles("CC(O)NCC=O"));
						// peptideTemplates.add(SmilesIO.readSmiles("CC(=O)NCCO"));
						// peptideTemplates.add(SmilesIO.readSmiles("CC(O)NCCO"));
						// peptideTemplates.add(SmilesIO.readSmiles("C\\C(O)=N\\CCO"));
						// peptideTemplates.add(SmilesIO.readSmiles("C\\C(O)=N\\CC=O"));
					} catch (Exception e) {
						e.printStackTrace();
					}
					
					// Find the template bond
					for (int i = 0; i < peptideTemplates.length; i++) {
						IAtomContainer peptideTemplate = peptideTemplates[i];
						for (int j = 0; j < peptideTemplate.getBondCount(); j++) {
							if (peptideTemplate.getBond(j).getAtom(0).getAtomicNumber() == 6
									&& peptideTemplate.getBond(j).getAtom(1)
											.getAtomicNumber() == 7
									|| peptideTemplate.getBond(j).getAtom(0)
											.getAtomicNumber() == 7
									&& peptideTemplate.getBond(j).getAtom(1)
											.getAtomicNumber() == 6) {
								IAtom currentCarbon = null;
								if (peptideTemplate.getBond(j).getAtom(0)
										.getAtomicNumber() == 6) {
									currentCarbon = peptideTemplate.getBond(j).getAtom(
											0);
								} else {
									currentCarbon = peptideTemplate.getBond(j).getAtom(
											1);
								}
								boolean carbonNextToOxygen = false;
								for (IAtom a : peptideTemplate
										.getConnectedAtomsList(currentCarbon)) {
									if (a.getAtomicNumber() == 8) {
										carbonNextToOxygen = true;
									}
								}
								if (!carbonNextToOxygen) {
									continue;
								}
								peptideTemplateBonds[i] = (peptideTemplate.getBond(j));
							}
						}
					}
		
		for (Fragment frag : monomerFragments) {
			IAtomContainer m = frag.getAtomContainer();
			// Keep track of rings with six or fewer atoms. If the C-N or C=N bond participates in this ring, break it but do not
			// change the connectivity annotation
			smallRings.addAll(ChemicalUtilities.getSmallestRings(frag.getAtomContainer()));
			
			// Keep track of atoms in rings of any size
			atomsInRings.addAll(ChemicalUtilities.getAtomsInRings(m));
			
			for (int x = 0; x < peptideTemplates.length; x++) {
				IAtomContainer peptideTemplate = peptideTemplates[x];
				IBond peptideTemplateBond = peptideTemplateBonds[x];
				List<List<RMap>> templateMatchMap = null;
				try {
					templateMatchMap = mod.getUIT()
							.getSubgraphMaps(m, peptideTemplate);
				} catch (CDKException e) {
					e.printStackTrace();
				}
				for (int i = 0; i < templateMatchMap.size(); i++) {
					for (int j = 0; j < templateMatchMap.get(i).size(); j++) {
						IBond currentTemplateMatchBond = peptideTemplate
								.getBond(templateMatchMap.get(i).get(j)
										.getId2());
						if (currentTemplateMatchBond == peptideTemplateBond) {
							IBond candidatePeptideBond = m
									.getBond(templateMatchMap.get(i).get(j)
											.getId1());

							IAtom currentNitrogen = null;
							if (candidatePeptideBond.getAtom(0)
									.getAtomicNumber() == 6) {
								currentNitrogen = candidatePeptideBond
										.getAtom(1);
							} else {
								currentNitrogen = candidatePeptideBond
										.getAtom(0);
							}
							// Check that this nitrogen has not been assigned to
							// another peptide bond already

							boolean nitrogenAlreadyPresent = false;
							for (IBond prevBond : peptideBonds) {
								if (prevBond.contains(currentNitrogen)) {
									nitrogenAlreadyPresent = true;
									break;
								}
							}
							if (nitrogenAlreadyPresent && frag.getAtomContainer().getConnectedAtomsCount(currentNitrogen) < 3) {
								break;
							}

							if (x == 0) {
								numStandardPeptideBonds++;
							} else {
								numNonStandardPeptideBonds++;
							}
							
							peptideBonds.add(candidatePeptideBond);
						}
					}
				}
			}

			if (MoleculeModifier.showBonds) {
				if (peptideBonds.isEmpty()) {
					SmilesIO.drawMolecule(m, mod.getName()
							+ "peptide_bonds");
				} else {
					SmilesIO.drawMoleculeHighlightingBonds(
							m,
							SmilesIO.getCleanFileName(mod.getName()
									+ "peptide_bonds"), peptideBonds);
				}
			}
		}
		
		// Swap the order of peptide bonds so that bonds with the carbon in a cycle are processed first.
		// This takes precedence over bonds involved in amino acid side chains 
		
		int lastUnvisitedIndex = peptideBonds.size()-1; // All elements from lastUnvisitedIndex onward have a cyclic carbon
		for(int i = 0; i < lastUnvisitedIndex; i++) {
			IAtom peptideC = null;
			if(peptideBonds.get(i).getAtom(0).getAtomicNumber() == 6) {
				peptideC = peptideBonds.get(i).getAtom(0);
			} else {
				peptideC = peptideBonds.get(i).getAtom(1);
			}
			if(!atomsInRings.contains(peptideC)) {
				Collections.swap(peptideBonds, i, lastUnvisitedIndex);
				lastUnvisitedIndex--;
				i--;
			}
		}
		IBond lactamBond = null;
		if(lastUnvisitedIndex > 0){
			lactamBond = peptideBonds.get(lastUnvisitedIndex);
		}
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		while (!q.isEmpty()) {
			Fragment subfragment = q.poll();
			// get peptide bond index or move on if there are no peptide bonds
			int peptideBondIndex = -1;
			int bondCount = subfragment.getAtomContainer().getBondCount();
			for (int i = 0; i < bondCount; i++) {
				if (peptideBonds.contains(subfragment.getAtomContainer().getBond(i))) {
					peptideBondIndex = i;
					break;
				}
			}
			if (peptideBondIndex == -1) {
				continue;
			}
			IBond peptideBond = subfragment.getAtomContainer().getBond(
					peptideBondIndex);
			if(peptideBond.equals(lactamBond)){
				for(IAtom atom : peptideBond.atoms()){
					subfragment.addLactamAtom(atom);
				}
			}
			IAtom betaCarbon = null;
			IAtom backboneNitrogen = null;
			if (peptideBond.getAtom(0).getAtomTypeName().startsWith("N.")) {
				backboneNitrogen = peptideBond.getAtom(0);
				betaCarbon = peptideBond.getAtom(1);
			} else {
				backboneNitrogen = peptideBond.getAtom(1);
				betaCarbon = peptideBond.getAtom(0);
			}
			// For the double bond isomer, set the oxygen back to sp2 and double
			// bonded.
			if (peptideBond.getOrder() == IBond.Order.DOUBLE) {
				List<IAtom> connectedBetaAtoms = subfragment.getAtomContainer()
						.getConnectedAtomsList(betaCarbon);
				for (IAtom a : connectedBetaAtoms) {
					if (a.getAtomTypeName().startsWith("O.")) {
						// This is the carbonyl oxygen
						a.setAtomTypeName("O.sp2");
						subfragment.getAtomContainer().getBond(betaCarbon, a)
								.setOrder(IBond.Order.DOUBLE);
					}
				}
			}

			subfragment.getAtomContainer().removeBond(peptideBondIndex);

			// Add -O for the OH group on the carbon

			Atom hydroxideO = new Atom("O");
			hydroxideO.setAtomTypeName("O.sp3");
			subfragment.getAtomContainer().addAtom(hydroxideO);
			subfragment.getAtomContainer().addBond(new Bond(betaCarbon, hydroxideO));

			List<Fragment> fragments = subfragment
					.partitionIntoMonomerFragments();
			if (fragments.size() == 1) {
				boolean bondParticipatesInSmallRing = false;
				for(Set<IAtom> atomsInRing : smallRings) {
					if(atomsInRing.contains(peptideBond.getAtom(0)) &&
							atomsInRing.contains(peptideBond.getAtom(1))){
						bondParticipatesInSmallRing = true;
					}
				}
				if(!bondParticipatesInSmallRing) {
					subfragment.addAminoC(betaCarbon);
					subfragment.addAminoN(backboneNitrogen);
					subfragment.setAtomAfterCTerminus(backboneNitrogen);
					subfragment.setAtomAfterNTerminus(betaCarbon);
				}
				q.add(subfragment);
			} else {
				// there must be two fragments
				Fragment CFragment = null, NFragment = null;
				if (fragments.get(0).getAtomContainer().contains(betaCarbon)) {
					CFragment = fragments.get(0);
					NFragment = fragments.get(1);
				} else {
					CFragment = fragments.get(1);
					NFragment = fragments.get(0);
				}
				

				CFragment.addAminoC(betaCarbon);
				NFragment.addAminoN(backboneNitrogen);
				
				// If this peptide connectivity does not conflict with previously assigned connections
				if(CFragment.getAtomAfterCTerminus() == null &&
						NFragment.getAtomAfterNTerminus() == null) {	
					CFragment.setAtomAfterCTerminus(backboneNitrogen);
					NFragment.setAtomAfterNTerminus(betaCarbon);
				}
				
				for(IAtom aminoNtoAdd : subfragment.getAminoNs()){
					if(NFragment.getAtomContainer().contains(aminoNtoAdd)){
						if(!NFragment.getAminoNs().contains(aminoNtoAdd)) NFragment.addAminoN(aminoNtoAdd);
					}
				}

				int indexToReplace = monomerFragments.indexOf(subfragment);

				monomerFragments.remove(indexToReplace);
				monomerFragments.add(indexToReplace, CFragment);
				monomerFragments.add(indexToReplace + 1, NFragment);
				q.add(NFragment);
				q.add(CFragment);
			}

		}
		if (MoleculeModifier.printSteps) {
			if (numStandardPeptideBonds > 0) {
				System.out.println("Processed " + numStandardPeptideBonds
						+ " standard peptide bonds");
			}
			if (numNonStandardPeptideBonds > 0) {
				System.out.println("Processed " + numNonStandardPeptideBonds
						+ " nonstandard peptide bonds");
			}
		}
	}


	/**
	 * Remove sugar groups
	 */
	private static void breakSugarGroups(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		@SuppressWarnings("unused")
		int numSugarBonds = 0;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			ArrayList<Fragment> subfragments = new ArrayList<Fragment>();
			subfragments.add(fragment);
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();
			//prepare sugar templates, and find sugar bonds
			IAtomContainer[] sugarTemplates = new IAtomContainer[4];
			IBond[] sugarTemplateBonds = new IBond[4];
			try {
				sugarTemplates[0] = SmilesIO.readSmilesTemplates("COC1CCCCO1"); 
				sugarTemplates[1] = SmilesIO.readSmilesTemplates("COC1CCCO1");
				sugarTemplates[2] = SmilesIO.readSmilesTemplates("CSC1CCCO1");
				sugarTemplates[3] = SmilesIO.readSmilesTemplates("CSC1CCCCO1");

			} catch (Exception e) {
				e.printStackTrace();
			}

			// Set template bonds
			sugarTemplateBonds[0] = sugarTemplates[0].getBond(sugarTemplates[0].getAtom(2), sugarTemplates[0].getAtom(1));
			sugarTemplateBonds[1] = sugarTemplates[1].getBond(sugarTemplates[1].getAtom(2), sugarTemplates[1].getAtom(1));
			sugarTemplateBonds[2] = sugarTemplates[2].getBond(sugarTemplates[2].getAtom(2), sugarTemplates[2].getAtom(1));
			sugarTemplateBonds[3] = sugarTemplates[3].getBond(sugarTemplates[3].getAtom(2), sugarTemplates[3].getAtom(1));
	

			// Find all sugar bonds
			
			for(int i = 0; i < sugarTemplateBonds.length; i++){
				IAtomContainer sugarTemplate = sugarTemplates[i];
				IBond sugarTemplateBond = sugarTemplateBonds[i];
				List<IBond> sugarBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						sugarTemplate, sugarTemplateBond, currentAtomContainer);

				for (IBond sugarBond : sugarBonds) {
					Fragment subfragment = null;
					for (Fragment possibleSubfragment : subfragments) {
						if (possibleSubfragment.getAtomContainer().contains(sugarBond)) {
							subfragment = possibleSubfragment;
						}
					}
					IAtom bondCarbon = null;
					IAtom bondOxygen = null;
					if (sugarBond.getAtom(0).getAtomTypeName().startsWith("C")) {
						bondCarbon = sugarBond.getAtom(0);
						bondOxygen = sugarBond.getAtom(1);
					} else {
						bondOxygen = sugarBond.getAtom(0);
						bondCarbon = sugarBond.getAtom(1);
					}
					
					//add check where the oxygen cannot be part of a small cycle

					subfragment.getAtomContainer().removeBond(sugarBond);
					// Check if this led to two pieces
					if (ConnectivityChecker.isConnected(subfragment.getAtomContainer())) {
						subfragment.getAtomContainer().addBond(sugarBond);
						continue;
					}

					mod.addSugarCarbon(bondCarbon);
					IAtom hydroxylO = new Atom("O");
					hydroxylO.setAtomTypeName("O.sp3");
					subfragment.getAtomContainer().addAtom(hydroxylO);
					subfragment.getAtomContainer().addBond(
							new Bond(hydroxylO, bondCarbon));

					List<Fragment> partitions = subfragment
							.partitionIntoMonomerFragments();

					// There must be two partitions
					Fragment sugarPart;
					Fragment nonSugarPart;
					if (partitions.get(0).getAtomContainer().contains(bondCarbon)) {
						nonSugarPart = partitions.get(0);
						sugarPart = partitions.get(1);
					} else {
						nonSugarPart = partitions.get(1);
						sugarPart = partitions.get(0);
					}
					nonSugarPart.setAtomInAttachedSugar(bondOxygen);
					sugarPart.setAtomOppositeSugarBond(bondCarbon);
					numSugarBonds++;
					int indexToReplace = subfragments.indexOf(subfragment);
					subfragments.remove(indexToReplace);
					subfragments.add(indexToReplace, nonSugarPart);
					subfragments.add(indexToReplace + 1, sugarPart);
				}
			}

			if(subfragments.size() > 1) {
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				for(int i = 0; i < subfragments.size(); i++) {
					monomerFragments.add(indexToReplace + i,
							subfragments.get(i));
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (numSugarBonds > 0) {
				System.out.println("Processed " + numSugarBonds
						+ " sugar bonds");
			}
		}
	}

	
	/**
	 * Break off sulfate groups
	 */
	private static void breakSulfateGroups(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer template = null;
		IBond templateBond = null;
		try {
			template = SmilesIO.readSmilesTemplates("COS");
		} catch (Exception e) {
			e.printStackTrace();
		}
		for (int i = 0; i < template.getBondCount(); i++) {
			IBond bond = template.getBond(i);
			if (bond.getAtom(0).getAtomicNumber() == 6
					|| bond.getAtom(1).getAtomicNumber() == 6) {
				templateBond = bond;
				break;
			}
		}
		boolean hasSulfate = false;
		for (Fragment frag : monomerFragments) {
			List<IBond> sulfateBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, frag.getAtomContainer());
			for (IBond sulfateBond : sulfateBonds) {
				IAtom sulfateO = null;
				IAtom connectingCarbon = null;
				if (sulfateBond.getAtom(0).getAtomicNumber() == 8) {
					sulfateO = sulfateBond.getAtom(0);
					connectingCarbon = sulfateBond.getAtom(1);
				} else {
					connectingCarbon = sulfateBond.getAtom(0);
					sulfateO = sulfateBond.getAtom(1);
				}
				// Remove bond
				frag.getAtomContainer().removeBond(sulfateBond);
				// Check that the bond removal has led to two pieces, of of
				// which has exactly one sulfur and four oxygens
				IAtomContainerSet partitions = ConnectivityChecker
						.partitionIntoMolecules(frag.getAtomContainer());
				if (partitions.getAtomContainerCount() != 2) {
					frag.getAtomContainer().addBond(sulfateBond);
					continue;
				}
				IAtomContainer sulfatePiece = null;
				IAtomContainer nonSulfatePiece = null;
				if (partitions.getAtomContainer(0).contains(sulfateO)) {
					sulfatePiece = partitions.getAtomContainer(0);
					nonSulfatePiece = partitions.getAtomContainer(1);
				} else {
					sulfatePiece = partitions.getAtomContainer(1);
					nonSulfatePiece = partitions.getAtomContainer(0);
				}
				// Make sure there are no carbons and at least four oxygens
				int numOxygens = 0;
				boolean hasCarbon = false;
				for (int i = 0; i < sulfatePiece.getAtomCount(); i++) {
					if (sulfatePiece.getAtom(i).getAtomicNumber() == 8) {
						numOxygens++;
					}
					if (sulfatePiece.getAtom(i).getAtomicNumber() == 6) {
						hasCarbon = true;
					}
				}
				if (hasCarbon || numOxygens < 4) {
					frag.getAtomContainer().addBond(sulfateBond);
					continue;
				}
				hasSulfate = true;
				// At this point, we conclude that this is indeed a sulfate
				// piece. Set this fragment AtomContainer the non-sulfur piece, set
				// enum.
				IAtom hydroxylO = new Atom("O");
				hydroxylO.setAtomTypeName("O.sp3");
				nonSulfatePiece.addAtom(hydroxylO);
				nonSulfatePiece.addBond(new Bond(connectingCarbon, hydroxylO));
				frag.setAtomContainer(nonSulfatePiece);
				frag.addTailoringDomain(TailoringDomainEnums.SULFOTRANSFERASE);
			}
			
			for (Fragment m : monomerFragments) { // check for *S(O)(O)O
				List<IAtom> sulfurAtoms = new ArrayList<IAtom>();
				for (int i = 0; i < m.getAtomContainer().getAtomCount(); i++) {
					if (m.getAtomContainer().getAtom(i).getAtomicNumber() == 16) {
						sulfurAtoms.add(m.getAtomContainer().getAtom(i));
					}
				}
				if (sulfurAtoms.size() == 0) {
					continue;
				}
				List<IAtom> oxygens = new ArrayList<IAtom>(); 
				for(IAtom sulfur : sulfurAtoms){
					for(IAtom connectedAtom : m.getAtomContainer().getConnectedAtomsList(sulfur)){
						if(connectedAtom.getAtomicNumber() == 8 && m.getAtomContainer().getConnectedAtomsList(connectedAtom).size() == 1){
							oxygens.add(connectedAtom);
						}
					}
					if(oxygens.size() == 3){
						m.getAtomContainer().removeAtomAndConnectedElectronContainers(sulfur);
						for(IAtom oxygen : oxygens){
							m.getAtomContainer().removeAtomAndConnectedElectronContainers(oxygen);
						}
						hasSulfate = true;
						m.getTailoringDomains().add(TailoringDomainEnums.SULFOTRANSFERASE);
					}
				}
			}
			
		}
		if (MoleculeModifier.printSteps) {
			if (hasSulfate) {
				System.out.println("Processed sulfate");
			}
		}
	}

	/**
	 * Find, annotate, and remove N, C, and O methylations to amino acids
	 * TODO fix so will loop through multiple amino atoms
	 */
	private static void processAAMethylations(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		@SuppressWarnings("unused")
		int numNMethylations = 0;
		@SuppressWarnings("unused")
		int numOMethylations = 0;
		@SuppressWarnings("unused")
		int numCMethylations = 0;
		for (Fragment m : monomerFragments) {
			
			//guess amino N if it wasn't cleaved (ie N terminus)
			IAtomContainer mol = m.getAtomContainer();
			if(m.getAminoNs().size() == 0 && m.getAminoCs().size() == 1){
				IAtom aminoC =  m.getAminoCs().get(0);
				for(IAtom atom : mol.getConnectedAtomsList(aminoC)){
					if(atom.getAtomicNumber() == 6){
						for(IAtom atom2 : mol.getConnectedAtomsList(atom)){
							if(atom2.getAtomicNumber() == 7){
								m.addAminoN(atom2);
							}
						}
					}
				}
			}
			
			for (int i = 0; i < m.getAtomContainer().getBondCount(); i++) {
				IBond currentBond = m.getAtomContainer().getBond(i);
				
				// N methylations
				if (m.getAminoNs().size() > 0) {
					IAtom candidateMethylC = null;
					if (currentBond.getAtom(0) == m.getAminoNs().get(0)
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(1);
					} else if (currentBond.getAtom(1) == m.getAminoNs().get(0)
							&& currentBond.getAtom(0).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(0);
					}
					if (candidateMethylC != null
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateMethylC) == 1) {
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateMethylC);
						m.getTailoringDomains().add(
								TailoringDomainEnums.N_METHYLTRANSFERASE);
						numNMethylations++;
						i --;
						continue;
					}
				}
				// C methylations
				if (m.getAminoNs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						if (m.getAtomContainer()
										.getConnectedAtomsList(
												currentBond.getAtom(0))
										.contains(m.getAminoNs().get(0))) {
							candidateAlphaCarbon = currentBond.getAtom(0);
							candidateMethylC = currentBond.getAtom(1);
						}
						if (m.getAtomContainer()
										.getConnectedAtomsList(
												currentBond.getAtom(1))
										.contains(m.getAminoNs().get(0))) {
							candidateAlphaCarbon = currentBond.getAtom(1);
							candidateMethylC = currentBond.getAtom(0);
						}
						if (m.getAtomContainer().getConnectedAtomsCount(
								candidateAlphaCarbon) == 4
								&& m.getAtomContainer().getConnectedAtomsCount(
										candidateMethylC) == 1) {
							m.getAtomContainer().removeBond(currentBond);
							m.getAtomContainer().removeAtom(candidateMethylC);
							m.getTailoringDomains().add(
									TailoringDomainEnums.C_METHYLTRANSFERASE);
							numCMethylations++;
							i --;
							continue;
						}
					}
				}
				if (m.getAminoCs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						if (m.getAtomContainer()
								.getConnectedAtomsList(currentBond.getAtom(0))
								.contains(m.getAminoCs().get(0))
								) {
							candidateAlphaCarbon = currentBond.getAtom(0);
							candidateMethylC = currentBond.getAtom(1);
						}
						if (m.getAtomContainer()
								.getConnectedAtomsList(currentBond.getAtom(1))
								.contains(m.getAminoCs().get(0))
								) {
							candidateAlphaCarbon = currentBond.getAtom(1);
							candidateMethylC = currentBond.getAtom(0);
						}
						if (m.getAtomContainer().getConnectedAtomsCount(
								candidateAlphaCarbon) == 4
								&& m.getAtomContainer().getConnectedAtomsCount(
										candidateMethylC) == 1) {
							m.getAtomContainer().removeBond(currentBond);
							m.getAtomContainer().removeAtom(candidateMethylC);
							m.getTailoringDomains().add(
									TailoringDomainEnums.C_METHYLTRANSFERASE);
							numCMethylations++;
							i --;
							continue;
						}
					}
				}
				// O methylations
				if (m.getAminoCs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateO = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 8
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(1);
						candidateO = currentBond.getAtom(0);
					} else if (currentBond.getAtom(1).getAtomicNumber() == 8
							&& currentBond.getAtom(0).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(0);
						candidateO = currentBond.getAtom(1);
					}
					if (m.getAtomContainer().getConnectedAtomsList(candidateO)
							.contains(m.getAminoCs().get(0))
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateMethylC) == 1) {
						candidateO.setAtomTypeName("O.sp3");
						m.getAtomContainer().getBond(candidateO, m.getAminoCs().get(0))
								.setOrder(IBond.Order.DOUBLE);
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateMethylC);
						m.getTailoringDomains().add(
								TailoringDomainEnums.O_METHYLTRANSFERASE);
						numOMethylations++;
						i --;
						continue;
					}
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (numNMethylations > 0) {
				System.out.println("Processed " + numNMethylations
						+ " N-methylated amino acids");
			}
			if (numOMethylations > 0) {
				System.out.println("Processed " + numOMethylations
						+ " O-methylated amino acids");
			}
			if (numCMethylations > 0) {
				System.out.println("Processed " + numCMethylations
						+ " C-methylated amino acids");
			}
		}
	}
	

	private static void processAACHydroxyliations(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		for (Fragment m : monomerFragments) {
			for (int i = 0; i < m.getAtomContainer().getBondCount(); i++) {
				IBond currentBond = m.getAtomContainer().getBond(i);
				if (m.getAminoNs().size() > 0) {
					IAtom candidateHydroxylO = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 8) {
						candidateAlphaCarbon = currentBond.getAtom(0);
						candidateHydroxylO = currentBond.getAtom(1);
					}else if(currentBond.getAtom(0).getAtomicNumber() == 8
							&& currentBond.getAtom(1).getAtomicNumber() == 6){
						candidateAlphaCarbon = currentBond.getAtom(1);
						candidateHydroxylO = currentBond.getAtom(0);
					}
					if (m.getAtomContainer().getConnectedAtomsCount(
							candidateAlphaCarbon) == 4
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateHydroxylO) == 1) {
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateHydroxylO);
						m.getTailoringDomains().add(
								TailoringDomainEnums.C_HYDROXYLATION);
						continue;
					}
				}
				if (m.getAminoCs().size() > 0) {
					IAtom candidateHydroxylO = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 8) {
						candidateAlphaCarbon = currentBond.getAtom(0);
						candidateHydroxylO = currentBond.getAtom(1);
					}else if(currentBond.getAtom(0).getAtomicNumber() == 8
							&& currentBond.getAtom(1).getAtomicNumber() == 6){
						candidateAlphaCarbon = currentBond.getAtom(1);
						candidateHydroxylO = currentBond.getAtom(0);
					}
					if (m.getAtomContainer().getConnectedAtomsCount(
							candidateAlphaCarbon) == 4
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateHydroxylO) == 1) {
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateHydroxylO);
						m.getTailoringDomains().add(
								TailoringDomainEnums.C_HYDROXYLATION);
						continue;
					}
				}
			}
		}
	}


	/**
	 * Find, annotate, and remove chlorine groups
	 */
	private static void processHalogenations(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		boolean hasHalogen = false;
		for (Fragment m : monomerFragments) {
			List<IAtom> halogenAtoms = new ArrayList<IAtom>();
			for (int i = 0; i < m.getAtomContainer().getAtomCount(); i++) {
				if (m.getAtomContainer().getAtom(i).getAtomicNumber() == 17
						|| m.getAtomContainer().getAtom(i).getAtomicNumber() == 35) {
					halogenAtoms.add(m.getAtomContainer().getAtom(i));
				}
			}
			if (halogenAtoms.size() == 0) {
				continue;
			}
			m.getTailoringDomains().add(TailoringDomainEnums.HALOGENATION);
			
			if(halogenAtoms.size() >= 3) {
				for(IAtom chlorine : halogenAtoms){
					IAtom connectedToChlorine = null;
					int numConnectedToSameAtom = 0;
					if(m.getAtomContainer().getConnectedAtomsCount(chlorine) != 1) continue;
					connectedToChlorine = m.getAtomContainer().getConnectedAtomsList(chlorine).get(0);
					for(IAtom atom : m.getAtomContainer().getConnectedAtomsList(connectedToChlorine)){
						if(atom.getAtomicNumber() == 17
								|| atom.getAtomicNumber() == 35){
							numConnectedToSameAtom ++;
						}
					}
					if(numConnectedToSameAtom == 3){
						Atom firstO = new Atom("O");
						firstO.setAtomTypeName("O.sp3");
						Atom secondO = new Atom("O");
						secondO.setAtomTypeName("O.sp2");
						Atom carbonylCarbon = new Atom("C");
						carbonylCarbon.setAtomTypeName("C.sp2");
						m.getAtomContainer().addAtom(firstO);
						m.getAtomContainer().addAtom(secondO);
						m.getAtomContainer().addAtom(carbonylCarbon);
						m.getAtomContainer().addBond(new Bond(carbonylCarbon, firstO, IBond.Order.SINGLE));
						m.getAtomContainer().addBond(new Bond(carbonylCarbon, secondO, IBond.Order.DOUBLE));
						m.getAtomContainer().addBond(new Bond(carbonylCarbon, connectedToChlorine, IBond.Order.SINGLE));
						break;
					}
				}
			}
			for (IAtom halogen : halogenAtoms) {
				// Check that there is exactly one connection
				if (m.getAtomContainer().getConnectedAtomsCount(halogen) != 1) {
					continue;
				}
				hasHalogen = true;
				m.getAtomContainer().removeBond(
						m.getAtomContainer().getBond(
								halogen,
								m.getAtomContainer().getConnectedAtomsList(halogen)
										.get(0)));
				m.getAtomContainer().removeAtom(halogen);
			}
		}
		if (MoleculeModifier.printSteps) {
			if (hasHalogen) {
				System.out.println("Processed halogenation");
			}
		}
	}
	

	private static void processSecondaryAmineExtensions(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			Set<IBond> bondsToBreak = new HashSet<IBond>();
			List<Set<IAtom>> smallestRings = null;
			//find extended secondary amines
			for(IAtom atom : mol.atoms()){
				
				if(smallestRings == null){
					smallestRings = ChemicalUtilities.getSmallestRings(mol); //make sure if in ring it's 6+ size
				}

				boolean inSmallRing = false; // cannot be in a ring < 7 (then is potentially a thiazol)
				for(Set<IAtom> ring : smallestRings){
					if(ring.contains(atom) && ring.size() < 7){
						inSmallRing = true;
						break;
					}
				}
				if(inSmallRing == true){
					continue;
				}
				
				//must be a nitrogen
				boolean possible = true;
				if(atom.getAtomicNumber() != 7){
					continue;
				}
				List<IAtom> connectedCarbons = new ArrayList<IAtom>();
				for(IAtom connectedAtom : mol.getConnectedAtomsList(atom)){
					//can only be connected to carbons
					if(connectedAtom.getAtomicNumber() == 6){
						//carbons cannot be a carbonyl
						for(IAtom connectedToCarbon : mol.getConnectedAtomsList(connectedAtom)){
							if(connectedToCarbon.getAtomicNumber() == 8 && mol.getBond(connectedAtom, connectedToCarbon).getOrder().equals(Order.DOUBLE)){
								possible = false;
								break;
							}
						}
						if(!possible){
							break;
						}else{
							connectedCarbons.add(connectedAtom);
						}
					}
				}
				//make sure it's still possible and is connected to 2 carbons
				if(!possible || connectedCarbons.size() != 2){
					continue;
				}
				IAtom aminoCarbon = null;
				for(IAtom carbon : connectedCarbons){
					//determine which is connected to a COOH
					int numConnections = 0;
					for(IAtom connectedCarbon : mol.getConnectedAtomsList(carbon)){
						//ensure it's a carbon
						if(connectedCarbon.getAtomicNumber() != 1){
							numConnections ++;
						}
						
						if(connectedCarbon.getAtomicNumber() != 6){
							continue;
						}
						
						boolean connectedToOH = false; // OH
						boolean connectedToDO = false; // =O
						for(IAtom connectedAtom : mol.getConnectedAtomsList(connectedCarbon)){
							if(connectedAtom.getAtomicNumber() == 8){
								if(mol.getBond(connectedCarbon, connectedAtom).getOrder().equals(Order.DOUBLE)){
									connectedToDO = true;
								}else if(mol.getBond(connectedCarbon, connectedAtom).getOrder().equals(Order.SINGLE)){
									connectedToOH = true;
								}
							}
						}
						
						if(connectedToDO && connectedToOH){
							aminoCarbon = carbon;
							break;
						}						
					}
					if(numConnections < 2){
						possible = false;
						break;
					}
					if(aminoCarbon != null){
						break;
					}
					
				}
				if(aminoCarbon != null && possible){
					for(IAtom carbon : connectedCarbons){
						if(!carbon.equals(aminoCarbon)){
							bondsToBreak.add(mol.getBond(carbon, atom));
						}
					}
					
				}
			}
			//break the secondary amines
			for(IBond bondToBreak : bondsToBreak){
				mol.removeBond(bondToBreak);
				
				//get the carbon
				IAtom carbon = null;
				for(IAtom atom : bondToBreak.atoms()){
					
					if(atom.getAtomicNumber() == 6){
						carbon = atom;
						break;
					}					
				}
				
				for(IAtom atom : mol.getConnectedAtomsList(carbon)){
					mol.removeBond(atom, carbon);
					IBond bond = new Bond(atom, carbon, Order.SINGLE);
					mol.addBond(bond);
				}
				
				//add carboxylic acid to the carbon
				Atom firstO = new Atom("O");
				firstO.setAtomTypeName("O.sp3");
				Atom secondO = new Atom("O");
				secondO.setAtomTypeName("O.sp2");
				mol.addAtom(firstO);
				mol.addAtom(secondO);
				mol.addAtom(carbon);
				mol.addBond(new Bond(carbon, firstO, IBond.Order.SINGLE));
				mol.addBond(new Bond(carbon, secondO, IBond.Order.DOUBLE));
				
			}
			
			//set new fragments if not connected
			IAtomContainerSet fragments = ConnectivityChecker
					.partitionIntoMolecules(mol);
			if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
				monomerFragments.remove(fragment);
				for (Fragment newFragment : fragment.partitionIntoMonomerFragments()) {
					monomerFragments.add(newFragment);
				}
			}		
		}
	}


}
