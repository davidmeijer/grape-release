package ca.mcmaster.magarveylab.grape.util.io;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.XMLIsotopeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.smiles.FixBondOrdersTool;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smsd.tools.ExtAtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.renderer.generators.*;

/**
 * Utility class for SMILES input/output operations, refactored from Mike's refacoring of Lian's iSNAP implementation.
 * @author gmchen
 *
 */
public class SmilesIO {
	
	/**
	 * Generate an IAtomContainer from a SMILES string. Ignores stereochemistry
	 * @param smiles	The structure to generate an IAtomContainer instance from, as a SMILES string
	 * @return			An IAtomContainer object corresponding to the input SMILES
	 * @throws IOException 
	 * @throws CDKException 
	 */
	public static IAtomContainer readSmiles(String smiles) throws IOException, CDKException {
		IAtomContainer mol = null;
		synchronized (SmilesIO.class) {
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
			SmilesParser parser = new SmilesParser(builder);
	        parser.kekulise(true);
			mol = parser.parseSmiles(smiles);
			mol = parser.parseSmiles(generateSmiles(mol));
			mol = AtomContainerManipulator.removeNonChiralHydrogens(mol);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		}
		return mol;
	}
	
	/**
	 * This method reads the smiles, but also converts to canonical, so no stereochemistry is considered
	 * @param smiles
	 * @return
	 * @throws IOException
	 * @throws CDKException
	 */
	public static IAtomContainer readSmilesTemplates(String smiles) throws IOException, CDKException {
		IAtomContainer mol = null;
		synchronized (SmilesIO.class) {
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
			SmilesParser parser = new SmilesParser(builder);
	        parser.kekulise(true);
			mol = parser.parseSmiles(smiles);
			// Remove hydrogens
			mol = AtomContainerManipulator.removeNonChiralHydrogens(mol);
			AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(mol);
		}
		return mol;
	}

	/**
	 * Generate a SMILES string from an IAtomContainer, this is canonical with no stereochemistry
	 * @param imol	The molecule to be parsed
	 * @return		The molecule's structures as a SMILES string
	 * @throws CDKException 
	 */
	public static String generateSmiles(IAtomContainer mol) throws CDKException {
		IAtomContainer molCopy = null;
		try {
			molCopy = mol.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		new SmilesGenerator();
		SmilesGenerator sg = SmilesGenerator.unique();
		String smiles = sg.create(molCopy);
		return smiles;
	}
	
	/**
	 * Generate png images of a molecule with coloured atoms
	 * @param molecule
	 * @param name
	 */
	public static void drawMolecule(IAtomContainer molecule, String name) {
		drawMolecule(molecule, name, false);
	}
	
	public static void drawMolecule(String outputPath, IAtomContainer molecule, String name) {
		drawMolecule(outputPath, molecule, name, false, new ArrayList<IAtom>(), new ArrayList<IBond>(), true);
	}

	/**
	 * Generate png images of a molecule with coloured atoms or numbers. Note: atom numbers are one-indexed
	 * @param molecule
	 * @param name
	 */
	public static void drawMolecule(IAtomContainer molecule, String name, boolean atomNumbers) {
		drawMolecule(molecule, name, atomNumbers, new ArrayList<IAtom>(), new ArrayList<IBond>());
	}
	
	/**
	 * Generate png images of a molecule with coloured atoms and highlighted atoms
	 * @param molecule
	 * @param name
	 */
	public static void drawMoleculeHighlightingAtoms(IAtomContainer molecule, String name, List<IAtom> highlightedAtoms) {
		drawMolecule(molecule, name, highlightedAtoms, new ArrayList<IBond>());
	}
	/**
	 * Generate png images of a molecule with uncoloured atoms and highlighted atoms
	 * @param molecule
	 * @param name
	 */
	public static void drawMoleculeHighlightingBonds(IAtomContainer molecule, String name, IBond highlightedBond) {
		ArrayList<IBond> highlightedBonds = new ArrayList<IBond>();
		highlightedBonds.add(highlightedBond);
		drawMolecule("image/", molecule, name, false, new ArrayList<IAtom>(), highlightedBonds, false);
	}
	/**
	 * Generate png images of a molecule with uncoloured atoms and highlighted atoms
	 * @param molecule
	 * @param name
	 */
	public static void drawMoleculeHighlightingBonds(IAtomContainer molecule, String name, List<IBond> highlightedBonds) {
		drawMolecule("image/", molecule, name, false, new ArrayList<IAtom>(), highlightedBonds, false);
	}
	
	/**
	 * Generate png images of a molecule with optionally coloured atoms and highlighted atoms
	 * @param molecule
	 * @param name
	 */
	public static void drawMoleculeHighlightingBonds(IAtomContainer molecule, String name, List<IBond> highlightedBonds, boolean colorByAtomType) {
		drawMolecule("image/", molecule, name, false, new ArrayList<IAtom>(), highlightedBonds, colorByAtomType);
	}
	
	/**
	 * Generate png images of a molecule with coloured atoms and highlighted atoms or bonds
	 * @param molecule
	 * @param name
	 */
	public static void drawMolecule(IAtomContainer molecule, String name, List<IAtom> highlightedAtoms, List<IBond> highlightedBonds) {
		drawMolecule(molecule, name, false, highlightedAtoms, highlightedBonds);
	}
	
	/**
	 * Generate png images of molecules with coloured atoms and highlighted atoms or bonds
	 * @param molecule
	 * @param name
	 */
	public static void drawMolecules(List<IAtomContainer> molecules, List<String> names) {
		for(int i = 0; i < molecules.size(); i++) {
			System.out.println("Drawing molecule " + i);
			drawMolecule(molecules.get(i), names.get(i));
		}
	}
	
	/**
	 * Generate png images of a molecule with coloured atoms or numbers. Note: atom numbers are one-indexed
	 * @param molecule
	 * @param name
	 * @throws CloneNotSupportedException 
	 */
	private static void drawMolecule(IAtomContainer molecule, String name, boolean atomNumbers, List<IAtom> highlightedAtoms, List<IBond> highlightedBonds) {
		drawMolecule("image/", molecule, name, atomNumbers, highlightedAtoms, highlightedBonds, true);
	}
	
	/**
	 * Generate png images of a molecule with optionally coloured atoms or numbers. Note: atom numbers are one-indexed
	 * @param molecule
	 * @param name
	 * @throws CloneNotSupportedException 
	 */
	private static void drawMolecule(String outputPath, IAtomContainer inMolecule, String name, boolean atomNumbers, List<IAtom> highlightedAtoms, List<IBond> highlightedBonds, boolean colorByAtomType) {
		// Derive from snippet from http://chem-bla-ics.blogspot.ca/2010/06/cdk-jchempaint-6-rendering-atom-numbers.html
		   int WIDTH = 1000;
	       int HEIGHT = 1000;
	       IAtomContainer molecule = AtomContainerManipulator.removeHydrogens(inMolecule);
	       // the draw area and the image should be the same size
	       Rectangle drawArea = new Rectangle(WIDTH, HEIGHT);
	       Image image = new BufferedImage(
	               WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
	       
	       // set the ID of highlighted atoms and bonds to the string "HIGHLIGHT"
	       IAtomContainer newMolecule = null;
		try {
			newMolecule = molecule.clone();
		} catch (CloneNotSupportedException e2) {
			e2.printStackTrace();
		}
	       for(int i = 0; i < molecule.getAtomCount(); i++) {
	    	   if(highlightedAtoms.contains(molecule.getAtom(i))) {
	    		   newMolecule.getAtom(i).setID("HIGHLIGHT");
	    	   }
	       }
	       for(int i = 0; i < molecule.getBondCount(); i++) {
	    	   if(highlightedBonds.contains(molecule.getBond(i))) {
	    		   newMolecule.getBond(i).setID("HIGHLIGHT");
	    	   }
	       }
	       
	       
	       StructureDiagramGenerator sdg = new StructureDiagramGenerator();
	       sdg.setMolecule(newMolecule);
           try {
			sdg.generateCoordinates();
		} catch (CDKException e1) {
			e1.printStackTrace();
		}
	       newMolecule = sdg.getMolecule();
	       
	       
	       // generators make the image elements
	       List<IGenerator<IAtomContainer>> generators = new ArrayList<IGenerator<IAtomContainer>>();
	       generators.add(new BasicSceneGenerator());
	       generators.add(new BasicBondGenerator());
	       
	       if(atomNumbers == true) {
	    	   generators.add(new AtomNumberGenerator());
	       } else {
	    	   generators.add(new BasicAtomGenerator());
	       }
	       
	       // the renderer needs to have a toolkit-specific font manager 
	       AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
	       
	       RendererModel model = renderer.getRenderer2DModel(); 
	       
	       model.set(BasicBondGenerator.BondWidth.class, 1.0);
	       if(atomNumbers == true) {
		       model.set( AtomNumberGenerator.ColorByType.class, true);
	       }
	       else {
	    	 //model.set(BasicAtomGenerator.CompactAtom.class, true);
		     model.set(BasicAtomGenerator.AtomRadius.class, 50.0);
		     model.set(BasicAtomGenerator.CompactShape.class, BasicAtomGenerator.Shape.OVAL);
		     model.set(BasicAtomGenerator.KekuleStructure.class, true);
		     if(!colorByAtomType) {
		    	 model.set(BasicAtomGenerator.ColorByType.class, false);
		     }
	       }
	       
	       //Highlight atoms and bonds
	       
	       model.getParameter(RendererModel.ExternalHighlightColor.class).setValue(Color.RED);
	      // model.set(BasicAtomGenerator.ColorByType.class, false);
	       
	       
	       Map<IChemObject, Color> map = new HashMap<IChemObject,Color>();
	       
	       
	       for(int i = 0; i < newMolecule.getAtomCount(); i++) {
	    	   if(newMolecule.getAtom(i).getID() == "HIGHLIGHT") {
	    		   // Note: this currently doesn't change the atom colour! Let gmchen know if you figure this out!
		    	   map.put(newMolecule.getAtom(i), model.getParameter(RendererModel.ExternalHighlightColor.class).getValue());
	    	   
	    	   }
	       }
	       for(int i = 0; i < newMolecule.getBondCount(); i++) {
	    	   if(newMolecule.getBond(i).getID() == "HIGHLIGHT") {
		    	   map.put(newMolecule.getBond(i), model.getParameter(RendererModel.ExternalHighlightColor.class).getValue());
	    	   }
	       }
	      
	       model.set(RendererModel.ColorHash.class, map);
	       
	       
	       
	       //This is for displaying the atom with the number offset
	       //model.set(
	    	//	   AtomNumberGenerator.Offset.class,
	    	//	   new javax.vecmath.Vector2d(10,10)
	    	//	 );
	       
	       
	       // the call to 'setup' only needs to be done on the first paint
	       //System.out.println(SmilesIO.generateSmiles(molecule));
	    
	       renderer.setup(newMolecule, drawArea);
	       
	       // paint the background
	       Graphics2D g2 = (Graphics2D)image.getGraphics();
	       g2.setColor(Color.WHITE);
	       g2.fillRect(0, 0, WIDTH, HEIGHT);
	       
	       // the paint method also needs a toolkit-specific renderer
	       renderer.paint(newMolecule, new AWTDrawVisitor(g2), new Rectangle2D.Double(0, 0, WIDTH, HEIGHT), true);
	       
       try {
    	   new File(outputPath).mkdirs();
    	   ImageIO.write((RenderedImage)image, "PNG", new File(outputPath + name + ".png"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static String getCleanFileName(String fileName) {
		String name = fileName.replaceAll("[^A-Za-z0-9\\-/,_]+", "_");
		if(name.length() > 100) {
			name = name.substring(0, 100);
		}
		return name;
	}
}
