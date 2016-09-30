package ca.mcmaster.magarveylab.grape.mol.chem.modifications;

import java.io.IOException;
import java.util.ArrayList;
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
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ca.mcmaster.magarveylab.grape.mol.chem.Fragment;
import ca.mcmaster.magarveylab.grape.mol.chem.MoleculeModifier;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class UniqueSubstructureBreakages{
	
	public static void process(MoleculeModifier mod) {
 
		breakAromaticEthers(mod);
		breakAdjoinedAromaticRings(mod);
		processBetaLactamLike(mod);
		processMonobactams(mod);
		processKendoLikeSubstructures(mod);
		processPieriLikeSubstructures(mod);
		processAnthramycinLikeSubstructures(mod);
		breakHybridDiCystines(mod);
		processSpectinomycinLikedoubleGlycosidicBond(mod);
		processHygromycinLikeGlycosidicBond(mod);
		processThiazs(mod);
		processOxazs(mod);
	}
	

	/**
	 * Break aromatic carbon rings connected ba an ether
	 */
	private static void breakAromaticEthers(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}

		int numBonds = 0;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();

			IAtomContainer[] templates = new IAtomContainer[2];
			IBond[] templateBonds = new IBond[2];
			try {
				// Benzene rings connected to each other
				templates[0] = SmilesIO.readSmilesTemplates("OC1=CC=CC=C1OC1=CC=CC=C1");
				templates[1] = SmilesIO
						.readSmilesTemplates("OC1=C(OC2=CC=CC=C2)C=CC=C1");
				// Benzene rings connected to another carbon chain
				// templates[2] = SmilesIO.readSmiles("CCOC1=CC=C(C)C=C1");
				// templates[3] = SmilesIO.readSmiles("CCOC1=CC(O)=CC(C)=C1");
				// templates[4] = SmilesIO.readSmiles("CCOC1=CC(C)=CC(O)=C1");

			} catch (Exception e) {
				e.printStackTrace();
			}
			// Find the template bonds for the first two templates
			for (int i = 0; i <= 1; i++) {
				// Walk from the hydroxyl oxygen to the ester bond
				IAtom hydroxylO = null;
				for (int j = 0; j < templates[i].getAtomCount(); j++) {
					if (templates[i].getAtom(j).getAtomicNumber() == 8
							|| templates[i].getConnectedAtomsCount(templates[i]
									.getAtom(j)) == 1) {
						hydroxylO = templates[i].getAtom(j);
						break;
					}
				}
				IAtom carbon1 = templates[i].getConnectedAtomsList(hydroxylO)
						.get(0);
				IAtom carbon2 = null;
				for (IAtom a : templates[i].getConnectedAtomsList(carbon1)) {
					if (templates[i].getConnectedAtomsCount(a) == 3) {
						carbon2 = a;
						break;
					}
				}
				for (IAtom a : templates[i].getConnectedAtomsList(carbon2)) {
					if (a.getAtomicNumber() == 8) {
						templateBonds[i] = templates[i].getBond(carbon2, a);
					}
				}
			}
			/*
			 * // Find the template bonds for templates 3-5 for(int i = 2; i <=
			 * 4; i++) { for(int j = 0; j < templates[i].getBondCount(); j++) {
			 * IBond bond = templates[i].getBond(j);
			 * if(bond.getAtom(0).getAtomicNumber() == 8 &&
			 * bond.getAtom(1).getAtomTypeName().equals("C.sp3")) {
			 * templateBonds[i] = bond; continue; }
			 * if(bond.getAtom(1).getAtomicNumber() == 8 &&
			 * bond.getAtom(0).getAtomTypeName().equals("C.sp3")) {
			 * templateBonds[i] = bond; continue; } } }
			 */
			for (int x = 0; x < templates.length; x++) {
				List<IBond> bondsToBreak = ChemicalUtilities.findMatchingBondsFromTemplate(
						templates[x], templateBonds[x], currentAtomContainer);
				for (IBond b : bondsToBreak) {
					if (currentAtomContainer.contains(b)) {
						currentAtomContainer.removeBond(b);
						numBonds++;
					}
				}
			}

			IAtomContainerSet fragments = ConnectivityChecker
					.partitionIntoMolecules(currentAtomContainer);
			if (fragments.getAtomContainerCount() > 1) {
				monomerFragments.remove(fragment);
				for (int i = 0; i < fragments.getAtomContainerCount(); i++) {
					monomerFragments.add(new Fragment(fragments
							.getAtomContainer(i)));
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (numBonds > 0) {
				System.out.println("Processed " + numBonds
						+ " aromatic ether(s)");
			}
		}
	}
	

	/**
	 * Break bonds which connect two phenyl rings
	 */
	private static void breakAdjoinedAromaticRings(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		boolean hasAdjoinedAromaticRings = false;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();
			IAtomContainer template = null;
			try {
				template = SmilesIO.readSmilesTemplates("C1=CC=C(C=C1)C1=CC=CC=C1");
			} catch (IOException | CDKException e) {
				e.printStackTrace();
			}
			IBond templateBond = null;
			for (int i = 0; i < template.getBondCount(); i++) {
				// The bond of interest is the bond connective two Carbons with
				// three connections
				if (template.getConnectedAtomsCount(template.getBond(i)
						.getAtom(0)) == 3
						&& template.getConnectedAtomsCount(template.getBond(i)
								.getAtom(1)) == 3) {
					templateBond = template.getBond(i);
					break;
				}
			}
			List<IBond> matchedBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, currentAtomContainer);

			// At this point, conditions have been met. Break the bond if it
			// leads to one piece.
			for (IBond currentBond : matchedBonds) {
				currentAtomContainer.removeBond(currentBond);
				if (ConnectivityChecker.partitionIntoMolecules(currentAtomContainer)
						.getAtomContainerCount() > 1) {
					currentAtomContainer.addBond(currentBond);
				} else {
					hasAdjoinedAromaticRings = true;
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (hasAdjoinedAromaticRings) {
				System.out.println("Processed adjoined aromatic rings");
			}
		}
	}


	/**
	 * Open up the rings of sulfur-containing beta-lactams: penams, penems, and
	 * cephems
	 */
	private static void processBetaLactamLike(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		// If this is a cephem, convert to a penam
		IAtomContainer cephemTemplate = null;
		try {
			cephemTemplate = SmilesIO.readSmilesTemplates("C1CN2C=CCSC12");
		} catch (Exception e) {
			e.printStackTrace();
		}
		// Store the double bond in the cephem template
		IBond cephemTemplateDoubleBond = null;
		// Store the cephem-carbon single bond that is to be broken
		IBond cephemTemplateSulfurBondToBreak = null;
		for (int i = 0; i < cephemTemplate.getBondCount(); i++) {
			if (cephemTemplate.getBond(i).getOrder() == IBond.Order.DOUBLE) {
				cephemTemplateDoubleBond = cephemTemplate.getBond(i);
				continue;
			}
			IBond candidateBond = cephemTemplate.getBond(i);
			if (candidateBond.getAtom(0).getAtomicNumber() != 16
					&& candidateBond.getAtom(1).getAtomicNumber() != 16) {
				continue;
			}
			IAtom currentSulfur = null;
			IAtom currentCarbon = null;
			if (candidateBond.getAtom(0).getAtomicNumber() == 16) {
				currentSulfur = candidateBond.getAtom(0);
				currentCarbon = candidateBond.getAtom(1);
			} else {
				currentSulfur = candidateBond.getAtom(1);
				currentCarbon = candidateBond.getAtom(0);
			}
			if (cephemTemplate.getConnectedAtomsCount(currentCarbon) != 2) {
				continue;
			}
			cephemTemplateSulfurBondToBreak = candidateBond;
		}
		ArrayList<IBond> cephemTemplateBonds = new ArrayList<IBond>();
		cephemTemplateBonds.add(cephemTemplateDoubleBond);
		cephemTemplateBonds.add(cephemTemplateSulfurBondToBreak);

		for (Fragment frag : monomerFragments) {

			List<List<IBond>> matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
					cephemTemplate, cephemTemplateBonds, frag.getAtomContainer());
			for (int i = 0; i < matchingBonds.size(); i++) {
				matchingBonds.get(i).get(0).setOrder(IBond.Order.SINGLE);
				matchingBonds.get(i).get(0).getAtom(0).setAtomTypeName("C.sp3");
				matchingBonds.get(i).get(0).getAtom(1).setAtomTypeName("C.sp3");
				frag.getAtomContainer().removeBond(matchingBonds.get(i).get(1));
				// Create new bond between the sulfur and the carbon in the
				// double bond to form a five-membered ring
				IAtom carbonPreviouslyConnectedToSulfur = null;
				IAtom carbonNewlyConnectingToSulfur = null;
				IAtom cephemSulfur = null;

				if (matchingBonds.get(i).get(1).getAtom(0).getAtomicNumber() == 6) {
					carbonPreviouslyConnectedToSulfur = matchingBonds.get(i)
							.get(1).getAtom(0);
					cephemSulfur = matchingBonds.get(i).get(1).getAtom(1);
				} else {
					carbonPreviouslyConnectedToSulfur = matchingBonds.get(i)
							.get(1).getAtom(1);
					cephemSulfur = matchingBonds.get(i).get(1).getAtom(0);
				}
				for (IAtom candidateNewlyConnectedCarbon : frag.getAtomContainer()
						.getConnectedAtomsList(
								carbonPreviouslyConnectedToSulfur)) {
					// If this candidate carbon is connected to one of the
					// previously double bonded carbons
					if (frag.getAtomContainer()
							.getConnectedAtomsList(
									candidateNewlyConnectedCarbon)
							.contains(matchingBonds.get(i).get(0).getAtom(0))
							|| frag.getAtomContainer()
									.getConnectedAtomsList(
											candidateNewlyConnectedCarbon)
									.contains(
											matchingBonds.get(i).get(0)
													.getAtom(1))) {
						carbonNewlyConnectingToSulfur = candidateNewlyConnectedCarbon;
					}
				}
				IBond newCephemBond = new Bond(cephemSulfur,
						carbonNewlyConnectingToSulfur);
				frag.getAtomContainer().addBond(newCephemBond);
				frag.setBetaLactamS(cephemSulfur);
			}
		}

		IAtomContainer sulfurBetaLactamTemplate = null;
		try {
			sulfurBetaLactamTemplate = SmilesIO.readSmilesTemplates("CSC1CC(=O)N1");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		// Store the bond in the four-membered ring that needs breaking
		IBond sulfurBetaLactamTemplateBond1 = null;
		// Store the bond in the five-membered ring that needs breaking
		IBond sulfurBetaLactamTemplateBond2 = null;
		for (int i = 0; i < sulfurBetaLactamTemplate.getBondCount(); i++) {
			IBond bond = sulfurBetaLactamTemplate.getBond(i);
			// Check if this is templateBond1
			if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(0), 16)) {
					sulfurBetaLactamTemplateBond1 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(1), 16)) {
					sulfurBetaLactamTemplateBond1 = bond;
				}
			} // Chcek if this is templateBond2
			else if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 16) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(0), 7)) {
					sulfurBetaLactamTemplateBond2 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 16) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(1), 7)) {
					sulfurBetaLactamTemplateBond2 = bond;
				}
			}
		}
		ArrayList<IBond> sulfurBetaLactamTemplateBonds = new ArrayList<IBond>();
		sulfurBetaLactamTemplateBonds.add(sulfurBetaLactamTemplateBond1);
		sulfurBetaLactamTemplateBonds.add(sulfurBetaLactamTemplateBond2);
		
		//check for clauvalnicAcidLike, these are the same as sulphur but oxygen replaces the sulfur
		IAtomContainer clavualnicAcidLikeTemplate = null;
		try {
			clavualnicAcidLikeTemplate = SmilesIO.readSmilesTemplates("COC1CC(=O)N1");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		// Store the bond in the four-membered ring that needs breaking
		IBond clavualnicAcidLikeTemplateBond1 = null;
		// Store the bond in the five-membered ring that needs breaking
		IBond clavualnicAcidLikeTemplateBond2 = null;
		for (int i = 0; i < clavualnicAcidLikeTemplate.getBondCount(); i++) {
			IBond bond = clavualnicAcidLikeTemplate.getBond(i);
			// Check if this is templateBond1
			if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(0), 8)) {
					clavualnicAcidLikeTemplateBond1 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(1), 8)) {
					clavualnicAcidLikeTemplateBond1 = bond;
				}
			} // Chcek if this is templateBond2
			else if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 8) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(0), 7)) {
					clavualnicAcidLikeTemplateBond2 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 8) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(1), 7)) {
					clavualnicAcidLikeTemplateBond2 = bond;
				}
			}
		}
		ArrayList<IBond> clavualnicAcidLikeTemplateBonds = new ArrayList<IBond>();
		clavualnicAcidLikeTemplateBonds.add(clavualnicAcidLikeTemplateBond1);
		clavualnicAcidLikeTemplateBonds.add(clavualnicAcidLikeTemplateBond2);
		
		
		IAtomContainer gammaLactamBetaLactoneTemplate = null;
		try {
			gammaLactamBetaLactoneTemplate = SmilesIO.readSmilesTemplates("O=C1OC2CC(=O)NC12");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		// Store the bond in the four-membered ring that needs breaking
		IBond gammaLactamBetaLactoneTemplateBond1 = gammaLactamBetaLactoneTemplate.getBond(gammaLactamBetaLactoneTemplate.getAtom(3),gammaLactamBetaLactoneTemplate.getAtom(8));
		// Store the bond in the five-membered ring that needs breaking
		IBond gammaLactamBetaLactoneTemplateBond2 = gammaLactamBetaLactoneTemplate.getBond(gammaLactamBetaLactoneTemplate.getAtom(3),gammaLactamBetaLactoneTemplate.getAtom(4));
		
		ArrayList<IBond> gammaLactamBetaLactoneTemplateBonds = new ArrayList<IBond>();
		gammaLactamBetaLactoneTemplateBonds.add(gammaLactamBetaLactoneTemplateBond1);
		gammaLactamBetaLactoneTemplateBonds.add(gammaLactamBetaLactoneTemplateBond2);		
		
		boolean foundMatch = false;
		for (Fragment frag : monomerFragments) {
			List<List<IBond>> matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
					sulfurBetaLactamTemplate, sulfurBetaLactamTemplateBonds, frag.getAtomContainer());
			if(matchingBonds.size() < 1){
				matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						gammaLactamBetaLactoneTemplate, gammaLactamBetaLactoneTemplateBonds, frag.getAtomContainer());
			}
			if(matchingBonds.size() < 1){
				matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						clavualnicAcidLikeTemplate, clavualnicAcidLikeTemplateBonds, frag.getAtomContainer());
			}
			for (int i = 0; i < matchingBonds.size(); i++) {
				frag.getAtomContainer().removeBond(matchingBonds.get(i).get(0));
				frag.getAtomContainer().removeBond(matchingBonds.get(i).get(1));
				// If the AtomContainer is no longer connected, then re-add these
				// bonds.
				if (!ConnectivityChecker.isConnected(frag.getAtomContainer())) {
					frag.getAtomContainer().addBond(matchingBonds.get(i).get(0));
					frag.getAtomContainer().addBond(matchingBonds.get(i).get(1));
				} else {
					foundMatch = true;
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (foundMatch) {
				System.out.println("Processed beta lactam like");
			}
		}
	}


	private static void processMonobactams(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer monobactamTemplate = null;
		try {
			monobactamTemplate = SmilesIO.readSmilesTemplates("O=C1CCN1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = monobactamTemplate.getBond(monobactamTemplate.getAtom(3),monobactamTemplate.getAtom(4));
		IBond templateBondToBreakSulfate = monobactamTemplate.getBond(monobactamTemplate.getAtom(4),monobactamTemplate.getAtom(1));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(monobactamTemplate, templateBondToBreak, mol);
			List<IBond> templateBondMatchesSulfate = ChemicalUtilities.findMatchingBondsFromTemplate(monobactamTemplate, templateBondToBreakSulfate, mol);
			if(templateBondMatches.size() == 0 && templateBondMatchesSulfate.size() == 0){
				continue;
			}
			boolean containsSulfate = false;
			for(IBond bond : templateBondMatchesSulfate){ //check if there is a sulfate
				for(IAtom atom : bond.atoms()){
					if(atom.getAtomicNumber() == 7){
						for(IAtom connectedAtom : mol.getConnectedAtomsList(atom)){
							if(connectedAtom.getAtomicNumber() == 16){
								containsSulfate = true;
							}
						}
					}
				}
			}
			if(containsSulfate){
				for(IBond bond : templateBondMatchesSulfate){
					for(IAtom atom : bond.atoms()){
						if(atom.getAtomicNumber() == 6){
							mol.removeBond(bond);
							Atom oxygen = new Atom("O");
							oxygen.setAtomTypeName("O.sp3");
							mol.addAtom(oxygen);
							mol.addBond(
									fragment.getAtomContainer().getAtomNumber(atom),
									fragment.getAtomContainer().getAtomNumber(oxygen),
									IBond.Order.SINGLE);
							fragment.addAminoC(atom);
						}else {
							fragment.addAminoN(atom);
						}
					}
				}
				continue;
			}
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(atom.getAtomicNumber() == 6){
						Atom oxygen = new Atom("O");
						oxygen.setAtomTypeName("O.sp3");
						mol.addAtom(oxygen);
						mol.addBond(
								fragment.getAtomContainer().getAtomNumber(atom),
								fragment.getAtomContainer().getAtomNumber(oxygen),
								IBond.Order.SINGLE);
					}
				}
			}
		}
	}

	

	private static void processHygromycinLikeGlycosidicBond(MoleculeModifier mod) { //CN[C@H]1C[C@@H](N)[C@H](O)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@@H]3OC4(O[C@H]23)O[C@H]([C@@H](N)CO)[C@H](O)[C@@H](O)[C@H]4O)[C@@H]1O
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer tempalte = null;
		try {
			tempalte = SmilesIO.readSmilesTemplates("COC1(C)OCCO1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = tempalte.getBond(tempalte.getAtom(2),tempalte.getAtom(4));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(tempalte, templateBondToBreak, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				IAtomContainerSet molFrags = ConnectivityChecker.partitionIntoMolecules(mol);
				if (molFrags.getAtomContainerCount() > 1){
					mol.addBond(bond);
					break;
				}
			}
		}
	}

	private static void processSpectinomycinLikedoubleGlycosidicBond(MoleculeModifier mod) { //CN[C@@H]1[C@H](O)[C@H](NC)[C@H]2O[C@]3(O)[C@@H](O[C@H](C)CC3=O)O[C@@H]2[C@H]1O
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer tempalte = null;
		try {
			tempalte = SmilesIO.readSmilesTemplates("OC1COCCO1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = tempalte.getBond(tempalte.getAtom(1),tempalte.getAtom(6));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(tempalte, templateBondToBreak, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				IAtomContainerSet molFrags = ConnectivityChecker.partitionIntoMolecules(mol);
				if (molFrags.getAtomContainerCount() > 1){
					mol.addBond(bond);
					break;
				}
			}
		}
	}
	

	private static void breakHybridDiCystines(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("CC1=NC(=CS1)C1=NC(C)=CS1");} catch (Exception e) {}
		IBond templateToBreakBond = template.getBond(template.getAtom(8),template.getAtom(9));
		IAtomContainer template2 = null;
		try {template2 = SmilesIO.readSmilesTemplates("CC1CSC(N1)C1CSC(C)=N1");} catch (Exception e) {}
		IBond templateToBreakBond2 = template2.getBond(template2.getAtom(1),template2.getAtom(0));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateToBreakBond, fragment.getAtomContainer());
			templateDoubleBondMatches.addAll(ChemicalUtilities.findMatchingBondsFromTemplate(template2, templateToBreakBond2, fragment.getAtomContainer()));
			if(templateDoubleBondMatches.size() < 1) continue;
			IAtomContainer mol = fragment.getAtomContainer();
			for(IBond bond : templateDoubleBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					IAtomContainer carboxylicAcid = null;
					try {carboxylicAcid = SmilesIO.readSmilesTemplates("C(O)=O");} catch (Exception e) {}
					IAtom carboxilicC = carboxylicAcid.getAtom(0);
					mol.add(carboxylicAcid);
					mol.addBond(
							mol.getAtomNumber(atom),
							mol.getAtomNumber(carboxilicC),
							IBond.Order.SINGLE);
					try {
						AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(mol);
					} catch (CDKException e) {
						e.printStackTrace();
					}
				}
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
	
	private static void processAnthramycinLikeSubstructures(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer[] templates = new IAtomContainer[2];
		try {
			templates[0] = SmilesIO.readSmilesTemplates("O=C1NCC=NC=C1");
			templates[1] = SmilesIO.readSmilesTemplates("O=C1CCN=CCN1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond[] templatesBond = new IBond[]{
				templates[0].getBond(templates[0].getAtom(4),templates[0].getAtom(5)),
				templates[1].getBond(templates[1].getAtom(4),templates[1].getAtom(5))
		};

		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = new ArrayList<IBond>();
			int templateIndex = 0; //template index to check
			boolean match = false; //change to true if a match comes up in the next statement
			while(!match && templates.length > templateIndex){
				templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(templates[templateIndex], templatesBond[templateIndex], mol);
				if(templateBondMatches.size() == 0){
					templateIndex ++;
				}else{
					match = true;
				}
			}
			if(!match){ //if no matches breaks the loop
				break;
			}

			for(IBond bond : templateBondMatches){
				bond.setOrder(Order.SINGLE);
				IAtom carbon = null;
				if(bond.getAtom(0).getAtomicNumber() == 6){
					carbon = bond.getAtom(0);
				}else{
					carbon = bond.getAtom(1);
				}
				Atom oxygen = new Atom("O");
				oxygen.setAtomTypeName("O.sp2");
				mol.addAtom(oxygen);
				mol.addBond(new Bond(carbon, oxygen, IBond.Order.DOUBLE));
			}
		}
	}

	private static void processPieriLikeSubstructures(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("CC1=C(C)C(=O)C(O)=C(O)N1");} catch (Exception e) {}
		IBond templateToBreakBond = template.getBond(template.getAtom(0),template.getAtom(1));
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			IAtomContainer molClone = null;
			try {molClone = fragment.getAtomContainer().clone();}catch (CloneNotSupportedException e) {}
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateToBreakBond, molClone);
			
			if(templateDoubleBondMatches.size() < 1) continue;
			
			for(IBond bond : templateDoubleBondMatches){
				molClone.removeBond(bond);
				
				if(!ConnectivityChecker.isConnected(molClone)){
					IAtomContainerSet splitFragments = ConnectivityChecker.partitionIntoMolecules(molClone);
					if(splitFragments.getAtomContainerCount() == 2){
						if(splitFragments.getAtomContainer(0).getAtomCount() > splitFragments.getAtomContainer(1).getAtomCount()){
							molClone.addBond(bond);
							continue;
						}
						if(splitFragments.getAtomContainer(0).getAtomCount() > splitFragments.getAtomContainer(1).getAtomCount()){
							molClone = splitFragments.getAtomContainer(0);
						}else{
							molClone = splitFragments.getAtomContainer(1);
						}
						IAtom terminalCarbon = null;
						if(molClone.contains(bond.getAtom(0))){
							terminalCarbon = bond.getAtom(0);
						}else if(molClone.contains(bond.getAtom(1))){
							terminalCarbon = bond.getAtom(1);
						}else{
							System.err.println("Something went wrong with method processPieriLikeSubstructures -- This should never happen, this method was ignored for the final breakdown of this compound");
							break;
						}
						IAtomContainer pieceToAdd = null;
						try {pieceToAdd = SmilesIO.readSmilesTemplates("CC(C=O)C(=O)CC(O)=O");} catch (Exception e) {}
						IAtom atomToConnect = pieceToAdd.getAtom(2);
						molClone.add(pieceToAdd);
						molClone.addBond(molClone.getAtomNumber(terminalCarbon), molClone.getAtomNumber(atomToConnect), IBond.Order.SINGLE);
						fragment.setAtomContainer(molClone);
					}else{
						molClone.addBond(bond);
						continue;
					}					
				}else{ //Shouldn't be connected
					molClone.addBond(bond);
					continue;
				}	
			}
			
		}
		
	}

	private static void processKendoLikeSubstructures(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("CC1=C2OC(C)(O)C=C2C(C)=C(O)C1=O");} catch (Exception e) {}
		IBond templateDoubleBond = template.getBond(template.getAtom(7),template.getAtom(8));
		IBond templateEtherBond = template.getBond(template.getAtom(3),template.getAtom(4));
		IBond templateHydroxylBond = template.getBond(template.getAtom(4),template.getAtom(6));
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		// go through each nrp fragment. In current implementation there should
		// be only one.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateDoubleBond, fragment.getAtomContainer());
			
			if(templateDoubleBondMatches.size() < 1) continue;
			
			List<IBond> templateEtherBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateEtherBond, fragment.getAtomContainer());
			List<IBond> templateHydroxylBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateHydroxylBond, fragment.getAtomContainer());
			
			for(IBond bond : templateDoubleBondMatches){
				fragment.getAtomContainer().removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(fragment.getAtomContainer().getConnectedAtomsCount(atom) == 2){ //Add ketone to carbon in the starter unit
						Atom ketone = new Atom("O");
						ketone.setAtomTypeName("O.sp2");
						fragment.getAtomContainer().addAtom(ketone);
						fragment.getAtomContainer().addBond(
								fragment.getAtomContainer().getAtomNumber(ketone),
								fragment.getAtomContainer().getAtomNumber(atom),
								IBond.Order.DOUBLE);
					}else if(fragment.getAtomContainer().getConnectedAtomsCount(atom) == 1){ //Add carboxcylic acid to carbon not in starter unit
						IAtomContainer carboxylicAcid = null;
						try {carboxylicAcid = SmilesIO.readSmilesTemplates("C(O)=O");} catch (Exception e) {}
						IAtom carboxilicC = carboxylicAcid.getAtom(0);
						fragment.getAtomContainer().add(carboxylicAcid);
						fragment.getAtomContainer().addBond(
								fragment.getAtomContainer().getAtomNumber(atom),
								fragment.getAtomContainer().getAtomNumber(carboxilicC),
								IBond.Order.SINGLE);
						
					}else{
						System.err.println("Something went wrong with method processKendoLikeSubstructures (Double bond match doesn't have appropriate carbons): reattached broken bond"); //TODO make exception
						fragment.getAtomContainer().addBond(bond);
					}				
				}		
			}
			for(IBond bond : templateEtherBondMatches){
				fragment.getAtomContainer().removeBond(bond);
				if(!ConnectivityChecker.isConnected(fragment.getAtomContainer())){
					fragment.getAtomContainer().addBond(bond);
					System.err.println("Something went wrong with method processKendoLikeSubstructures (Breaking ether makes AtomContainer not connected): reattached broken bond");
				}
			}
			for(IBond bond : templateHydroxylBondMatches){
				bond.setOrder(IBond.Order.DOUBLE);
			}
		}
	}
	

	/**
	 * Process thiazoles and thiazolines. Searches for thiazole/thizoline rings and 'undoes' the reaction.
	 * This method modifies monomerFragments and updates the field
	 * thiazoleSList.
	 */
	private static void processThiazs(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		// First, handle the situation corresponding to the five-membered S and
		// N containing ring
		// as in: CNC(=O)C1=CSC(CN)=N1
		List<Set<IAtom>> smallestRings = null;
		boolean hasThiaz = false;
		for (Fragment m : monomerFragments) {
			IAtomContainer mol = m.getAtomContainer();
			List<IAtom> carbonylCAtoms = new ArrayList<IAtom>();
			List<IAtom> thiazoneNAtoms = new ArrayList<IAtom>();
			List<IAtom> chainC1Atoms = new ArrayList<IAtom>();
			List<IAtom> chainC2Atoms = new ArrayList<IAtom>();
			List<IAtom> thiazSAtoms = new ArrayList<IAtom>();
			List<Boolean> thiazoles = new ArrayList<Boolean>();
			for (int i = 0; i < m.getAtomContainer().getAtomCount(); i++) {
				if (mol.getAtom(i).getAtomicNumber() != 16) {
					continue;
				}
				// check if this sulfur is connected to a carbon (C) double
				// bonded with a nitrogen
				IAtom thiazS = null;
				IAtom carbonylC = null;
				IAtom thiazN = null;
				IAtom chainC_1 = null;
				IAtom chainC_2 = null;
				boolean thiazole = false;

				for (IAtom c : mol.getConnectedAtomsList(mol
						.getAtom(i))) {
					if (mol.getBond(mol.getAtom(i), c).getOrder() != IBond.Order.SINGLE) {
						continue;
					}
					if (!c.getAtomTypeName().equals("C.sp2")) {
						continue;
					}
					for (IAtom n : mol.getConnectedAtomsList(c)) {
						if (n.getAtomTypeName().startsWith("C.")
								&& mol.getBond(c, n).getOrder() == IBond.Order.DOUBLE) {
							chainC_1 = c;
							chainC_2 = n;
							thiazole = true;
							continue;
						}
						if (n.getAtomTypeName().startsWith("N.")
								&& mol.getBond(c, n).getOrder() == IBond.Order.DOUBLE) {
							carbonylC = c;
							thiazN = n;
						}
					}
				}
				
				if(smallestRings == null){
					smallestRings = ChemicalUtilities.getSmallestRings(mol); //make sure if in ring it's 6+ size
				}
				
				boolean inRingSizeFive = false; // must be in a ring size 5
				
				for(Set<IAtom> ring : smallestRings){
					if(ring.contains(thiazN) && ring.size() == 5){
						inRingSizeFive = true;
						break;
					}
				}				
				if (carbonylC == null || !inRingSizeFive) {
					continue;
				}
				
				if (chainC_1 != null && chainC_2 != null) {
					IBond bond = mol.getBond(chainC_1, chainC_2);
					bond.setOrder(IBond.Order.SINGLE);
					chainC_1.setAtomTypeName("C.sp3");
					chainC_2.setAtomTypeName("C.sp3");
				}

				thiazS = mol.getAtom(i);

				carbonylCAtoms.add(carbonylC);
				thiazoneNAtoms.add(thiazN);
				chainC1Atoms.add(chainC_1);
				chainC2Atoms.add(chainC_2);
				thiazSAtoms.add(thiazS);
				thiazoles.add(thiazole);
			}
			for (int i = 0; i < thiazSAtoms.size(); i++) {
				IAtom thiazS = thiazSAtoms.get(i);
				IAtom carbonylC = carbonylCAtoms.get(i);
				IAtom thiazN = thiazoneNAtoms.get(i);
				boolean thiazole = thiazoles.get(i);

				// Cleave the C-S bond. If this alters connectivity, then it
				// wasn't part of a cycle so re-add it.
				IBond removedBond = mol.removeBond(carbonylC, thiazS);
				if (!ConnectivityChecker.isConnected(mol)) {
					mol.addBond(removedBond);
					continue;
				}

				// Turn the C=N bond to C-N
				mol.getBond(carbonylC, thiazN).setOrder(
						IBond.Order.SINGLE);

				// Add a =0 to that C
				Atom atom = new Atom("O");
				atom.setAtomTypeName("O.sp2");
				mol.addAtom(atom);
				mol.addBond(new Bond(carbonylC, atom, IBond.Order.DOUBLE));

				// Set carbonylC to sp2

				carbonylC.setAtomTypeName("C.sp2");
				thiazN.setAtomTypeName("N.amide");
				
				if(thiazole){
					mod.addThiazoleS(thiazS);
				}else{
					mod.addThiazolineS(thiazS);
				}

				hasThiaz = true;
			}
		}
		if (MoleculeModifier.printSteps) {
			if (hasThiaz) {
				System.out.println("Processed thiaz");
			}
		}

		// Next, generally look for any cysteine whose sulfur is connected to a
		// carbon, with that S-C bond a part of a cycle.
		boolean hasCyclizedNonThiazoleCysteine = false;
		IAtomContainer template = null;
		try {
			template = SmilesIO.readSmilesTemplates("CSCC(N)C=O");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		
		IBond templateBond = template.getBond(template.getAtom(0), template.getAtom(1));

		for (Fragment frag : monomerFragments) {
			List<IBond> matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, frag.getAtomContainer());
			for (IBond matchingBond : matchingBonds) {
				IAtom thiazS = null;
				IAtom carbonylC = null;
				if (matchingBond.getAtom(0).getAtomicNumber() == 16) {
					thiazS = matchingBond.getAtom(0);
					carbonylC = matchingBond.getAtom(1);
				} else {
					carbonylC = matchingBond.getAtom(0);
					thiazS = matchingBond.getAtom(1);
				}
				
				boolean inRingSizeFive = false; // must be in a ring size 5
				
				for(Set<IAtom> ring : smallestRings){
					if(ring.contains(thiazS) && ring.size() == 5){
						inRingSizeFive = true;
						break;
					}
				}
				if(!inRingSizeFive){
					continue;
				}
				
				frag.getAtomContainer().removeBond(matchingBond);
				if (ConnectivityChecker.isConnected(frag.getAtomContainer())) {
					// Set all bonds connecting to carbonylC as single
					for (int i = 0; i < frag.getAtomContainer().getBondCount(); i++) {
						IBond b = frag.getAtomContainer().getBond(i);
						if (b.contains(carbonylC)) {
							b.setOrder(IBond.Order.SINGLE);
						}
					}
					Atom atom = new Atom("O");
					atom.setAtomTypeName("O.sp2");
					frag.getAtomContainer().addAtom(atom);
					frag.getAtomContainer().addBond(
							new Bond(carbonylC, atom, IBond.Order.DOUBLE));
					mod.addThiazoleS(thiazS);
					hasCyclizedNonThiazoleCysteine = true;
				} else {
					frag.getAtomContainer().addBond(matchingBond);
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (hasCyclizedNonThiazoleCysteine) {
				System.out.println("Processed non-thiazole cyclized cysteine");
			}
		}
	}

	/**
	 * Process oxazoles and oxazolines. Searches for oxazole/oxazoline rings and 'undoes' the reaction.
	 * Looks for CMT on the hydroxyl carbon
	 * This method modifies monomerFragments and updates the field oxazoleOList.
	 */
	private static void processOxazs(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		boolean hasOxazole = false;
		for (Fragment m : monomerFragments) {
			IAtomContainer mol = m.getAtomContainer();
			IAtomContainer[] oxazoleTemplates = new IAtomContainer[2];
			
			try {
				oxazoleTemplates[0] = SmilesIO.readSmilesTemplates("O1C=CN=C1");
				oxazoleTemplates[1] = SmilesIO.readSmilesTemplates("C1CN=CO1");
			} catch(Exception e) {
				e.printStackTrace();
			}
			
			// ArrayList of template bonds corresponding to each oxazole template.
			//The first bond is the double bond with nitrogen, the second bond is the carbon-carbon bond, and the third bond is the oxygen-carbon bond which is broken.
			List<List<IBond>> templateBonds = new ArrayList<List<IBond>>();
			for(int i = 0; i < oxazoleTemplates.length; i++) {
				// Initialize the size of the template bonds
				templateBonds.add(new ArrayList<IBond>());
				for(int j = 0; j < 3; j++) {
					templateBonds.get(i).add(null);
				}
				for(int j = 0; j < oxazoleTemplates[i].getBondCount(); j++) {
					IBond currentBond = oxazoleTemplates[i].getBond(j);
					if(currentBond.getOrder() == IBond.Order.DOUBLE &&
						(currentBond.getAtom(0).getAtomicNumber() == 7 || currentBond.getAtom(1).getAtomicNumber() == 7)) {
						templateBonds.get(i).set(0, currentBond);
						continue;
					}
					if(currentBond.getAtom(0).getAtomicNumber() == 6 && currentBond.getAtom(1).getAtomicNumber() == 6) {
						templateBonds.get(i).set(1, currentBond);
						continue;
					}
					if(currentBond.getAtom(0).getAtomicNumber() == 8 || currentBond.getAtom(1).getAtomicNumber() == 8) {
						@SuppressWarnings("unused")
						IAtom currentO = null;
						IAtom currentC = null;
						if(currentBond.getAtom(0).getAtomicNumber() == 8) {
							currentO = currentBond.getAtom(0);
							currentC = currentBond.getAtom(1);
						}
						else {
							currentO = currentBond.getAtom(1);
							currentC = currentBond.getAtom(0);
						}
						// Check that the C is connected to a nitrogen
						boolean connectedToNitrogen = false;
						for(IAtom connectedAtom : oxazoleTemplates[i].getConnectedAtomsList(currentC)) {
							if(connectedAtom.getAtomicNumber() == 7) {
								connectedToNitrogen = true;
							}
						}
						if(connectedToNitrogen) {
							templateBonds.get(i).set(2, currentBond);
						}
					}
				}
			}
			
			for(int i = 0; i < oxazoleTemplates.length; i++) {
				List<List<IBond>> matchingBondsList = ChemicalUtilities.findMatchingBondsFromTemplate(oxazoleTemplates[i], templateBonds.get(i), m.getAtomContainer());
				for(List<IBond> matchingBonds : matchingBondsList) {
					matchingBonds.get(0).setOrder(IBond.Order.SINGLE);
					if(matchingBonds.get(0).getAtom(0).getAtomicNumber() == 6) {
						matchingBonds.get(0).getAtom(0).setAtomTypeName("C.sp3");
						matchingBonds.get(0).getAtom(1).setAtomTypeName("N.sp3");
					}
					else {
						matchingBonds.get(0).getAtom(0).setAtomTypeName("N.sp3");
						matchingBonds.get(0).getAtom(1).setAtomTypeName("C.sp3");
					}
					matchingBonds.get(1).setOrder(IBond.Order.SINGLE);
					matchingBonds.get(1).getAtom(0).setAtomTypeName("C.sp3");
					matchingBonds.get(1).getAtom(1).setAtomTypeName("C.sp3");
					m.getAtomContainer().removeBond(matchingBonds.get(2));
					IAtom oxazoleO = null;
					IAtom newKetoneC = null;
					if(matchingBonds.get(2).getAtom(0).getAtomicNumber() == 8) {
						oxazoleO = matchingBonds.get(2).getAtom(0);
						newKetoneC = matchingBonds.get(2).getAtom(1);
					}
					else {
						oxazoleO = matchingBonds.get(2).getAtom(1);
						newKetoneC = matchingBonds.get(2).getAtom(0);
					}
					
					//check if the oxazole carbon is methylated
					for(IAtom atom : mol.getConnectedAtomsList(oxazoleO)){
						if(mol.getConnectedBondsCount(atom) == 3){
							for(IAtom atom2 : mol.getConnectedAtomsList(atom)){
								if(atom2.getAtomicNumber() == 6 
										&& ChemicalUtilities.getConnectedAtomsCountNonHydrogen(mol, atom2) == 1){
									mol.removeAtomAndConnectedElectronContainers(atom2);
									mod.addMethylatedCarbon(atom);
								}
							}
						}
					}
					
					Atom ketoneO = new Atom("O");
					ketoneO.setAtomTypeName("O.sp2");
					newKetoneC.setAtomTypeName("C.sp2");
					mol.addAtom(ketoneO);
					mol.addBond(new Bond(newKetoneC, ketoneO, IBond.Order.DOUBLE));
					
					hasOxazole = true;
					if(i == 0) { //oxazole
						mod.addOxazoleO(oxazoleO);
					}else { //oxazoline
						mod.addOxazolineO(oxazoleO);
					}
				}
			}
			
			}
		if (MoleculeModifier.printSteps) {
			if (hasOxazole) {
				System.out.println("Processed oxazole");
			}
		}
	}


}
