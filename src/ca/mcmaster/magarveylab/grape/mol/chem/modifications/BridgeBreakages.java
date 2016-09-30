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

import ca.mcmaster.magarveylab.grape.mol.chem.Fragment;
import ca.mcmaster.magarveylab.grape.mol.chem.MoleculeModifier;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class BridgeBreakages {

	public static void process(MoleculeModifier mod) {
		breakDisulfideBridges(mod);
		breakLinearEthers(mod);
		modifyImine(mod); 
		breakThioesters(mod);
		openLactoneRings(mod);
		processUreido(mod);
		breakEsterLinkages(mod);
	}
	

	/**
	 * Break sulfur-sulfur bonds
	 */
	private static void breakDisulfideBridges(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		@SuppressWarnings("unused")
		int numBonds = 0;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();

			IAtomContainer template = null;
			try {
				template = SmilesIO.readSmilesTemplates("CSSC");

			} catch (Exception e) {
				e.printStackTrace();
			}

			IBond templateBond = null;
			for (int i = 0; i < template.getBondCount(); i++) {
				if (template.getBond(i).getAtom(0).getAtomicNumber() == 16
						&& template.getBond(i).getAtom(1).getAtomicNumber() == 16) {
					templateBond = template.getBond(i);
					break;
				}
			}
			
			List<IBond> bondsToBreak = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, currentAtomContainer);
			
			if(bondsToBreak.size() > 0){
				for (IBond b : bondsToBreak) {
					if (currentAtomContainer.contains(b)) {
						currentAtomContainer.removeBond(b);
						numBonds++;
					}
				}
				IAtomContainerSet fragments = ConnectivityChecker
						.partitionIntoMolecules(currentAtomContainer);
				if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
					monomerFragments.remove(fragment);
					for (int i = 0; i < fragments.getAtomContainerCount(); i++) {
						monomerFragments.add(new Fragment(fragments
								.getAtomContainer(i)));
					}
				}
				
			}
		}
		if (MoleculeModifier.printSteps) {
			if (numBonds > 0) {
				System.out.println("Processed " + numBonds
						+ " disulfide bridge(s)");
			}
		}
	}

	

	private static void breakLinearEthers(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IAtom> cyclicAtoms = ChemicalUtilities.getAtomsInRings(mol);
			List<IAtom> linearEtherOxygens = new ArrayList<IAtom>(); 
			for(IAtom oxiCandidate : mol.atoms()){
				if(oxiCandidate.getAtomicNumber() != 8 
						|| mol.getConnectedAtomsCount(oxiCandidate) != 2 
						|| cyclicAtoms.contains(oxiCandidate)
						|| mol.getBondOrderSum(oxiCandidate) != 2
						) continue;
				boolean properEther = true;
				for(IAtom connectedC : mol.getConnectedAtomsList(oxiCandidate)){
					int numOxygens = 0;
					for(IAtom connectedAtom : mol.getConnectedAtomsList(connectedC)){
						if(connectedAtom.getAtomicNumber() == 8) numOxygens ++;
					}
					if(numOxygens > 1){
						properEther = false; 
						break; //Stop checking the atoms
					}
					if(connectedC.getAtomicNumber() != 6 //Must be connected to two carbons
					 || mol.getConnectedAtomsCount(connectedC) < 2){
						properEther = false; //Carbons must not be terminal
						break; //Stop checking the atoms
					}
				}
				if(properEther) linearEtherOxygens.add(oxiCandidate);
			}
			for(IAtom etherOxygen : linearEtherOxygens){
				List<IAtom>connectedCs = new ArrayList<IAtom>(); 
				for(IAtom connectedC : mol.getConnectedAtomsList(etherOxygen)){
					connectedCs.add(connectedC);
				}
				mol.removeAtomAndConnectedElectronContainers(etherOxygen);
				for(IAtom connectedC : connectedCs){
					Atom oxygen = new Atom("O");
					oxygen.setAtomTypeName("O.sp3");
					mol.addAtom(oxygen);
					mol.addBond(
							fragment.getAtomContainer().getAtomNumber(connectedC),
							fragment.getAtomContainer().getAtomNumber(oxygen),
							IBond.Order.SINGLE);
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


	private static void modifyImine(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		// go through each nrp fragment. In current implementation there should
		// be only one.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();

			for (int j = 0; j < fragment.getAtomContainer().getAtomCount(); j++) {
				if (!(fragment.getAtomContainer().getAtom(j).getAtomTypeName()
						.equals("N.sp2"))) {
					continue;
				}
				IAtom imineN = fragment.getAtomContainer().getAtom(j);
				List<Set<IAtom>> rings = ChemicalUtilities.getSmallestRings(fragment.getAtomContainer());
				boolean inSmallRing = false;
				for(Set<IAtom> ringAtoms : rings){
					if(ringAtoms.contains(imineN)){
						inSmallRing = true;
						break;
					}
				}
				if(inSmallRing) continue;
				
				IAtom imineC = null;
				for (IAtom c : fragment.getAtomContainer().getConnectedAtomsList(imineN)) {
					if (c.getAtomTypeName().equals("C.sp2")) {

						if (fragment.getAtomContainer()
										.getBond(imineN, c)
										.getOrder() == IBond.Order.DOUBLE
										&& fragment.getAtomContainer().getConnectedAtomsList(c).size() == 2){
							imineC = c;
							
						}
					}
				}
				if(imineC == null) continue;
				// Then check if breaking the imineN and imineC leaves
				// one piece.
				IBond removedBond = fragment.getAtomContainer().removeBond(
						imineN, imineC);

				if (ConnectivityChecker.isConnected(fragment.getAtomContainer())) {
					//Do the removal of double bond and create single
					fragment.getAtomContainer().addBond(
							fragment.getAtomContainer().getAtomNumber(imineN),
							fragment.getAtomContainer().getAtomNumber(imineC),
							IBond.Order.SINGLE);
					Atom ketone = new Atom("O");
					ketone.setAtomTypeName("O.sp2");
					fragment.getAtomContainer().addAtom(ketone);
					fragment.getAtomContainer().addBond(
							fragment.getAtomContainer().getAtomNumber(ketone),
							fragment.getAtomContainer().getAtomNumber(imineC),
							IBond.Order.DOUBLE);
					
				}else{
					//re-add the original bond
					fragment.getAtomContainer().addBond(removedBond);
				}
			}
		}
		
	}

	private static void breakThioesters(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer template = null;
		try {
			template = SmilesIO.readSmilesTemplates("CSC(C)=O");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		IBond templateBond = null;
		for (int i = 0; i < template.getBondCount(); i++) {
			// The bond of interest involves a sulfur and a carbon with three
			// connections
			if (template.getConnectedAtomsCount(template.getBond(i).getAtom(0)) == 3
					&& template.getBond(i).getAtom(1).getAtomicNumber() == 16
					|| template.getConnectedAtomsCount(template.getBond(i)
							.getAtom(1)) == 3
					&& template.getBond(i).getAtom(0).getAtomicNumber() == 16) {
				templateBond = template.getBond(i);
				break;
			}
		}

		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}

		while (!q.isEmpty()) {
			Fragment currentFragment = q.poll();
			List<IBond> matchedBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, currentFragment.getAtomContainer());

			for (IBond currentBond : matchedBonds) {
				IAtom thioesterCarbon = null;
				if (currentBond.getAtom(0).getAtomicNumber() == 6) {
					thioesterCarbon = currentBond.getAtom(0);
				} else {
					thioesterCarbon = currentBond.getAtom(1);
				}
				currentFragment.addThioEsterAtom(currentBond.getAtom(0));
				currentFragment.addThioEsterAtom(currentBond.getAtom(1));
				currentFragment.getAtomContainer().removeBond(currentBond);
				Atom carboxylO = new Atom("O");
				carboxylO.setAtomTypeName("O.sp3");
				currentFragment.getAtomContainer().addAtom(carboxylO);
				Bond newBond = new Bond(thioesterCarbon, carboxylO);
				currentFragment.getAtomContainer().addBond(newBond);

				List<Fragment> fragments = currentFragment
						.partitionIntoMonomerFragments();

				if (fragments.size() == 2) {
					Fragment sFrag = null;
					Fragment cFrag = null;
					if (fragments.get(0).getAtomContainer()
							.contains(thioesterCarbon)) {
						cFrag = fragments.get(0);
						sFrag = fragments.get(1);
					} else {
						cFrag = fragments.get(1);
						sFrag = fragments.get(0);
					}
					int indexToReplace = monomerFragments
							.indexOf(currentFragment);
					monomerFragments.remove(indexToReplace);
					monomerFragments.add(indexToReplace, sFrag);
					monomerFragments.add(indexToReplace + 1, cFrag);
					q.add(sFrag);
					q.add(cFrag);
				}
			}
		}
	}


	/**
	 * Open lactone rings, labelling Fragment adjacent lactone fragments
	 * when possible. This method modifies monomerFragments and also updates the
	 * field lactoneCarboxylCList.
	 */
	private static void openLactoneRings(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		@SuppressWarnings("unused")
		int numLactoneRings = 0;
		// go through each nrp fragment. In current implementation there should
		// be only one.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			ArrayList<Fragment> subfragments = new ArrayList<Fragment>();
			subfragments.add(fragment);
			// Find all cyclic esters. First find the oxygen in the ring, check
			// that it has two carbon neighbours, check that one of the carbon
			// neighbours
			// is double bonded to another O, then check that cleaving the C-O
			// bond gives a connected product

			ArrayList<IAtom> hydroxylOList = new ArrayList<IAtom>();
			ArrayList<IAtom> hydroxylCList = new ArrayList<IAtom>();
			ArrayList<IAtom> carboxylCList = new ArrayList<IAtom>();

			for (int j = 0; j < fragment.getAtomContainer().getAtomCount(); j++) {

				if (!(fragment.getAtomContainer().getAtom(j).getAtomTypeName()
						.equals("O.sp3"))) {
					continue;
				}

				IAtom hydroxyl_O = fragment.getAtomContainer().getAtom(j);
				// check if it is linked to two carbons, one of which has a
				// double bond O.

				IAtom carboxylC = null;
				IAtom hydroxylC = null;

				for (IAtom c : fragment.getAtomContainer().getConnectedAtomsList(
						hydroxyl_O)) {
					// the carboxyl carbon must be sp2 hybridized and connected
					// to another oxygen via double bond
					if (c.getAtomTypeName().equals("C.sp2")) {
						for (IAtom candidateOtherO : fragment.getAtomContainer()
								.getConnectedAtomsList(c)) {
							if (candidateOtherO != hydroxyl_O
									&& fragment.getAtomContainer()
											.getBond(candidateOtherO, c)
											.getOrder() == IBond.Order.DOUBLE
									&& candidateOtherO.getAtomicNumber() == 8) {
								carboxylC = c;
							}
						}
					} else if (c.getAtomTypeName().equals("C.sp3")) {
						hydroxylC = c;
					}
				}

				if (carboxylC == null || hydroxylC == null) {
					continue;
				}

				// Then check if breaking the hydroxyl_O and carboxylC leaves
				// one piece.
				IBond removedBond = fragment.getAtomContainer().removeBond(
						hydroxyl_O, carboxylC);

				if (ConnectivityChecker.isConnected(fragment.getAtomContainer())) {
					hydroxylOList.add(hydroxyl_O);
					hydroxylCList.add(hydroxylC);
					carboxylCList.add(carboxylC);
					// Add this carboxyl C to the global list
					mod.addLactoneCarboxylicC(carboxylC);
				}
				fragment.getAtomContainer().addBond(removedBond);
			}

			if (MoleculeModifier.showBonds) {
				// For testing: draw the AtomContainer with detected bonds
				// highlighted
				ArrayList<IBond> bonds = new ArrayList<IBond>();
				for (int j = 0; j < hydroxylOList.size(); j++) {
					bonds.add(fragment.getAtomContainer().getBond(
							hydroxylOList.get(j), carboxylCList.get(j)));
				}
				SmilesIO.drawMoleculeHighlightingBonds(
						fragment.getAtomContainer(),
						SmilesIO.getCleanFileName(mod.getName()
								.replace(' ', '_') + "Lactone"), bonds);
			}

			// Break the cyclic esters

			for (int j = 0; j < hydroxylOList.size(); j++) {
				numLactoneRings++;
				int index = -1;
				Fragment currentFragment = null;
				for (int k = 0; k < subfragments.size(); k++) {
					if (subfragments.get(k).getAtomContainer()
							.contains(hydroxylOList.get(j))) {
						if (!subfragments.get(k).getAtomContainer()
								.contains(carboxylCList.get(j))) {
							System.err
									.println("error: fragment does not contain both lactone O and C");
						}
						index = k;
						break;
					}
				}

				currentFragment = subfragments.get(index);
				IAtom currentHydroxylO = hydroxylOList.get(j);
				IAtom currentHydroxylC = hydroxylCList.get(j);
				IAtom currentCarboxylC = carboxylCList.get(j);
				currentFragment.getAtomContainer().removeBond(currentHydroxylO,
						currentCarboxylC);
				Atom atom = new Atom("O");
				atom.setAtomTypeName("O.sp3");
				currentFragment.getAtomContainer().addAtom(atom);
				currentFragment.getAtomContainer().addBond(
						new Bond(carboxylCList.get(j), atom));

				List<Fragment> partitions = currentFragment
						.partitionIntoMonomerFragments();
				if (partitions.size() == 1) {
					currentFragment.setLactoneHydroxylO(currentHydroxylO);
					currentFragment.setLactoneHydroxylC(currentHydroxylC);
					currentFragment.setLactoneCarboxylC(currentCarboxylC);
					currentFragment
							.setAtomAfterLactoneHydroxyl(currentCarboxylC);
					currentFragment
							.setAtomAfterLactoneCarboxyl(currentHydroxylO);
				} else {
					// there must be two partitions
					Fragment carboxylPartiton = null;
					Fragment hydroxylPartition = null;

					if (partitions.get(0).getAtomContainer()
							.contains(carboxylCList.get(j))) {
						carboxylPartiton = partitions.get(0);
						hydroxylPartition = partitions.get(1);
					} else {
						carboxylPartiton = partitions.get(1);
						hydroxylPartition = partitions.get(0);
					}

					carboxylPartiton.setLactoneCarboxylC(currentCarboxylC);
					hydroxylPartition.setLactoneHydroxylO(currentHydroxylO);
					hydroxylPartition.setLactoneHydroxylC(currentHydroxylC);
					carboxylPartiton
							.setAtomAfterLactoneCarboxyl(currentHydroxylO);
					hydroxylPartition
							.setAtomAfterLactoneHydroxyl(currentCarboxylC);

					subfragments.remove(index);
					subfragments.add(index, carboxylPartiton);
					subfragments.add(index + 1, hydroxylPartition);
				}
			}
			if (subfragments.size() > 1) {
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				for (int i = 0; i < subfragments.size(); i++) {
					monomerFragments.add(indexToReplace + i,
							subfragments.get(i));
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (numLactoneRings > 0) {
				System.out.println("Processed " + numLactoneRings
						+ " lactone rings");
			}
		}
	}


	/**
	 * Process ureido linkages
	 */
	private static void processUreido(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}

		boolean hasUreido = false;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();

			for (int i = 0; i < fragment.getAtomContainer().getAtomCount(); i++) {

				IAtom mainC = null;
				IAtom main_O = null;
				IAtom n_1 = null;
				IAtom c_1 = null;
				IAtom n_2 = null;
				IAtom c_2 = null;

				if (!fragment.getAtomContainer().getAtom(i).getAtomTypeName()
						.equals("C.sp2")) {
					continue;
				}
				mainC = fragment.getAtomContainer().getAtom(i);
				for (int j = 0; j < fragment.getAtomContainer()
						.getConnectedAtomsCount(mainC); j++) {
					IAtom connectedAtom = fragment.getAtomContainer()
							.getConnectedAtomsList(mainC).get(j);

					if (connectedAtom.getAtomTypeName().startsWith("O.")
							&& fragment.getAtomContainer()
									.getBond(mainC, connectedAtom).getOrder() == IBond.Order.DOUBLE) {
						main_O = connectedAtom;
					} else if (connectedAtom.getAtomTypeName().startsWith("N.")
							&& n_1 == null) {
						n_1 = connectedAtom;
					} else if (connectedAtom.getAtomTypeName().startsWith("N.")) {
						n_2 = connectedAtom;
					}

				}
				if (mainC == null || main_O == null || n_1 == null
						|| n_2 == null) {
					continue;
				}

				for (int j = 0; j < fragment.getAtomContainer()
						.getConnectedAtomsCount(n_1); j++) {
					IAtom connectedAtom = fragment.getAtomContainer()
							.getConnectedAtomsList(n_1).get(j);
					if (connectedAtom == mainC) {
						continue;
					} else if (connectedAtom.getAtomTypeName().startsWith("C.")) {
						c_1 = connectedAtom;
						break;
					}
				}
				for (int j = 0; j < fragment.getAtomContainer()
						.getConnectedAtomsCount(n_2); j++) {
					IAtom connectedAtom = fragment.getAtomContainer()
							.getConnectedAtomsList(n_2).get(j);
					if (connectedAtom == mainC) {
						continue;
					} else if (connectedAtom.getAtomTypeName().startsWith("C.")) {
						c_2 = connectedAtom;
						break;
					}
				}
				if (c_1 == null || c_2 == null) {
					continue;
				}

				// At this point we know that we have a -C-N-(C=O)-N-C-
				// structure). Now we break bonds.
				List<IBond> removedBonds = new ArrayList<IBond>();

				removedBonds.add(fragment.getAtomContainer().getBond(mainC, main_O));
				removedBonds.add(fragment.getAtomContainer().getBond(mainC, n_1));
				removedBonds.add(fragment.getAtomContainer().getBond(mainC, n_2));

				fragment.getAtomContainer().removeBond(
						fragment.getAtomContainer().getBond(mainC, main_O));
				fragment.getAtomContainer().removeAtom(main_O);
				fragment.getAtomContainer().removeBond(
						fragment.getAtomContainer().getBond(mainC, n_1));
				fragment.getAtomContainer().removeBond(
						fragment.getAtomContainer().getBond(mainC, n_2));
				fragment.getAtomContainer().removeAtom(mainC);

				IAtomContainerSet fragments = ConnectivityChecker
						.partitionIntoMolecules(fragment.getAtomContainer());
				if (fragments.getAtomContainerCount() != 2) {
					// We must have had a ring, so re-add what we removed.
					fragment.getAtomContainer().addAtom(main_O);
					fragment.getAtomContainer().addAtom(mainC);
					for (IBond bond : removedBonds) {
						fragment.getAtomContainer().addBond(bond);
					}
					continue;
				}

				if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
					monomerFragments.remove(fragment);
					for (Fragment newFragment : fragment.partitionIntoMonomerFragments()) {
						monomerFragments.add(newFragment);
					}
				}
				hasUreido = true;
				break;
			}
		}
		if (MoleculeModifier.printSteps) {
			if (hasUreido) {
				System.out.println("Processed ureido linkage");
			}
		}
	}

	private static void breakEsterLinkages(MoleculeModifier mod) {
		List<Fragment> monomerFragments = mod.getFragments();
		IAtomContainer template = null;
		IBond templateBond = null;
		try {
			template = SmilesIO.readSmilesTemplates("COC=O");
		} catch (Exception e) {
			e.printStackTrace();
		}
		templateBond = template.getBond(template.getAtom(1),template.getAtom(2)); //ester bond
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		boolean foundEster = false;
		while (!q.isEmpty()) {
			Fragment frag = q.poll();
			List<IBond> esterBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, frag.getAtomContainer());
			for (IBond esterBond : esterBonds) {
				foundEster = true;
				IAtom esterC = null;
				if (esterBond.getAtom(0).getAtomicNumber() == 6) {
					esterC = esterBond.getAtom(0);
				} else {
					esterC = esterBond.getAtom(1);
				}
				// Remove bond
				frag.getAtomContainer().removeBond(esterBond);
				IAtom carboxylO = new Atom("O");
				carboxylO.setAtomTypeName("O.sp3");
				frag.getAtomContainer().addAtom(carboxylO);
				frag.getAtomContainer().addBond(new Bond(esterC, carboxylO));
			}
			if (!ConnectivityChecker.isConnected(frag.getAtomContainer())) {
				List<Fragment> subfragments = frag
						.partitionIntoMonomerFragments();
				monomerFragments.remove(frag);
				for (Fragment subfrag : subfragments) {
					monomerFragments.add(subfrag);
				}
			}
		}
		if (MoleculeModifier.printSteps) {
			if (foundEster) {
				System.out.println("Processed non-cyclic ester(s)");
			}
		}
	}


}
