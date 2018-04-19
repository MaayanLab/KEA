package edu.mssm.pharm.maayanlab.KEA;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.logging.Logger;

import pal.statistics.FisherExact;
import edu.mssm.pharm.maayanlab.common.core.FileUtils;
import edu.mssm.pharm.maayanlab.common.math.SetOps;
import edu.mssm.pharm.maayanlab.common.core.Settings;
import edu.mssm.pharm.maayanlab.common.core.SettingsChanger;

public class KEA implements SettingsChanger {
	
	static Logger log = Logger.getLogger(KEA.class.getSimpleName());

	private LinkedList<Kinase> kinases;

	private final String PROTEIN_BACKGROUND = "res/kinase-protein_interactions.csv";
	private final String PHOSPHO_BACKGROUND = "res/phosphorylation_reactions.csv";
	private final String KEA_2018_BACKGROUND = "TODO";
	private final String ARCHS4_BACKGROUND = "TODO";
	private final String IPTMNET_BACKGROUND = "res/iPTMnet_kinome_interactions.txt";
	private final String NETWORK_IN_BACKGROUND = "TODO";
	private final String PHOSPHO_ELM_BACKGROUND = "TODO";
	private final String PHOSPHOPOINT_BACKGROUND = "TODO";
	private final String PHOSPHOPLUS_BACKGROUND = "TODO";
	private final String MINT_BACKGROUND = "TODO";

	private final String BACKGROUND_RANKS = "res/kea_ranks.txt";
	private final String KEA_2018_BACKGROUND_RANKS = "TODO";
	private final String ARCHS4_BACKGROUND_RANKS = "TODO";
	private final String IPTMNET_BACKGROUND_RANKS = "res/iptmnet_ranks.txt";
	private final String NETWORK_IN_BACKGROUND_RANKS = "TODO";
	private final String PHOSPHO_ELM_BACKGROUND_RANKS = "TODO";
	private final String PHOSPHOPOINT_BACKGROUND_RANKS = "TODO";
	private final String PHOSPHOPLUS_BACKGROUND_RANKS = "TODO";
	private final String MINT_BACKGROUND_RANKS = "TODO";

	// Output header
	protected final String HEADER = "Kinase,Substrates in Input,Substrates in Database,Input Fraction,Database Fraction,Difference,P-value,Z-score,Combined Score,Substrates";
	
	private final Settings settings = new Settings() {
		{
			// String: kinase interactions to include as background database. [kinase-protein interactions only/phosphorylation reactions only/both types]
			set(KEA.KINASE_INTERACTIONS, KEA.BOTH_TYPES);
			// String: rank the kinases by the Fisher Exact test's p-value, rank against the background of random genes, or combined score of the two. [combined_score/p-value/rank]
			set(KEA.SORT_BY, KEA.COMBINED_SCORE);
			// String: level of kinase resolution. [kinase-group/kinase-family/kinase]
			set(KEA.RESOLUTION_LEVEL, KEA.KINASE_LEVEL);
		}
	};
	
	// Setting keys
	public final static String KINASE_INTERACTIONS = "kinase interactions to include";
	public final static String SORT_BY = "sort kinases by";
	public final static String COMBINED_SCORE = "combined score";
	public final static String RESOLUTION_LEVEL = "resolve kinases down to";
	
	// Setting values
	public final static String KINASE_PROTEIN = "kinase-protein interactions only";
	public final static String PHOSPHORYLATION = "phosphorylation reactions only";
	public final static String KEA_2018 = "kea 2018";
	public final static String ARCHS4 = "ARCHS4";
	public final static String IPTMNET = "iPTMnet";
	public final static String NETWORK_IN = "NetworkIN";
	public final static String PHOSPHO_ELM = "Phospho.ELM";
	public final static String PHOSPHOPOINT = "Phosphopoint";
	public final static String PHOSPHOPLUS = "PhosphoPlus";
	public final static String MINT = "MINT";
	public final static String BOTH_TYPES = "both types";
	public final static String PVALUE = "p-value";
	public final static String RANK = "rank";	
	public final static String KINASE_LEVEL = "kinase";
	public final static String KINASE_FAMILY_LEVEL = "kinase-family";
	public final static String KINASE_GROUP_LEVEL = "kinase-group";

	public static void main(String[] args) {
		if (args.length == 2) {
			KEA kea = new KEA();
			kea.run(args[0]);
			kea.writeFile(args[1]);
		}
		else if (args.length == 3) {
			KEA kea = new KEA();
			kea.run(args[0], args[1]);
			kea.writeFile(args[2]);
		}		
		else
			log.severe("Usage: java -jar kea.jar [background] genelist output");
	}
	
	public KEA() {
		settings.loadSettings();
	}
	
	public KEA(Settings externalSettings) {
		settings.loadSettings(externalSettings);
	}
	
	public void setSetting(String key, String value) {
		settings.set(key, value);
	}
	
	public void run(String background, String geneList) {
		ArrayList<String> inputList = FileUtils.readFile(geneList);
		
		try {
			if (FileUtils.validateList(inputList))
				computeEnrichment(FileUtils.readFile(background), inputList);
		} catch (ParseException e) {
			if (e.getErrorOffset() == -1)
				log.warning("Invalid input: " + "Input list is empty.");
			else
				log.warning("Invalid input: " + e.getMessage() + " at line " + (e.getErrorOffset() + 1) +" is not a valid Entrez Gene Symbol.");
			System.exit(-1);	
		}
	}
	
	public void run(String geneList) {
		ArrayList<String> inputList = FileUtils.readFile(geneList);
		
		try {
			if (FileUtils.validateList(inputList))
				run(inputList);
		} catch (ParseException e) {
			if (e.getErrorOffset() == -1)
				log.warning("Invalid input: " + "Input list is empty.");
			else
				log.warning("Invalid input: " + e.getMessage() + " at line " + (e.getErrorOffset() + 1) +" is not a valid Entrez Gene Symbol.");
			System.exit(-1);	
		}
	}
	
	public void run(Collection<String> genelist) {
		computeEnrichment(assembleBackgroundDatabase(), genelist);
	}
	
	public void writeFile(String filename) {
		FileUtils.writeFile(filename, HEADER, kinases);
	}

	public Collection<Kinase> getTopRanked(int ranks) {
		LinkedList<Kinase> topRanked = new LinkedList<Kinase>();
		
		Iterator<Kinase> itr = kinases.iterator();
		while (itr.hasNext() && topRanked.size() < ranks)
			topRanked.add(itr.next());
		
		return topRanked;
	}
	
	public Collection<String> getTopRankedList(int ranks) {
		LinkedList<String> topRanked = new LinkedList<String>();
		
		Iterator<Kinase> itr = kinases.iterator();
		while (itr.hasNext() && topRanked.size() < ranks)
			topRanked.add(itr.next().getName());
		
		return topRanked;
	}
	
	// For internal use to get ranked list of kinase names
	public Collection<Kinase> getRankedList() {
		return kinases;
	}
	
	@Deprecated
	// For internal use to get ranked list of kinase names
	protected Collection<Kinase> getTopRankedList() {
		return kinases;
	}
	
	private ArrayList<String> assembleBackgroundDatabase() {
		ArrayList<String> background;
		
		if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.KINASE_PROTEIN)) {
			background = FileUtils.readResource(PROTEIN_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.PHOSPHORYLATION)) {
			background = FileUtils.readResource(PHOSPHO_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.ARCHS4)) {
			background = FileUtils.readResource(ARCHS4_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.IPTMNET)) {
			background = FileUtils.readResource(IPTMNET_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.NETWORK_IN)) {
			background = FileUtils.readResource(NETWORK_IN_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.PHOSPHO_ELM)) {
			background = FileUtils.readResource(PHOSPHO_ELM_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.PHOSPHOPOINT)) {
			background = FileUtils.readResource(PHOSPHOPOINT_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.PHOSPHOPLUS)) {
			background = FileUtils.readResource(PHOSPHOPLUS_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.MINT)) {
			background = FileUtils.readResource(MINT_BACKGROUND);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.KEA_2018)) {
			background = FileUtils.readResource(KEA_2018_BACKGROUND);
		}
		else {
			background = FileUtils.readResource(PROTEIN_BACKGROUND);
			background.addAll(FileUtils.readResource(PHOSPHO_BACKGROUND));
		}

		return background;
	}
	
	private ArrayList<String> assembleBackgroundRanks() {
		ArrayList<String> kea_ranks;

		if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.KEA_2018)) {
			kea_ranks = FileUtils.readResource(KEA_2018_BACKGROUND_RANKS);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.ARCHS4)) {
			kea_ranks = FileUtils.readResource(ARCHS4_BACKGROUND_RANKS);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.IPTMNET)) {
			kea_ranks = FileUtils.readResource(IPTMNET_BACKGROUND_RANKS);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.NETWORK_IN)) {
			kea_ranks = FileUtils.readResource(NETWORK_IN_BACKGROUND_RANKS);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.PHOSPHO_ELM)) {
			kea_ranks = FileUtils.readResource(PHOSPHO_ELM_BACKGROUND_RANKS);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.PHOSPHOPOINT)) {
			kea_ranks = FileUtils.readResource(PHOSPHOPOINT_BACKGROUND_RANKS);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.PHOSPHOPLUS)) {
			kea_ranks = FileUtils.readResource(PHOSPHOPLUS_BACKGROUND_RANKS);
		}
		else if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.MINT)) {
			kea_ranks = FileUtils.readResource(MINT_BACKGROUND_RANKS);
		}
		else /* if (settings.get(KEA.KINASE_INTERACTIONS).equals(KEA.KEA_2018)) */ {
			kea_ranks = FileUtils.readResource(BACKGROUND_RANKS);
			// Phospho, Phosphopoint, PhosphositePlus, and MINT combined ??
		}

		return kea_ranks;
	}
	
	private void computeEnrichment(ArrayList<String> background, Collection<String> genes) {
		
		// Determine level of kinase resolution (0-level: kinase-group, 1-level: kinase family: 2-level: kinase)
		int level;
		if (settings.get(RESOLUTION_LEVEL).equals(KINASE_GROUP_LEVEL))
			level = 0;
		else if (settings.get(RESOLUTION_LEVEL).equals(KINASE_FAMILY_LEVEL))
			level = 1;
		else
			level = 2;
		
		// parse background
		HashMap<String, Kinase> kinaseMap = new HashMap<String, Kinase>();
		for (String record: background) {
			String[] splitRecord = record.split(",");
			String name = splitRecord[level];
			String target = splitRecord[3];
			
			if (kinaseMap.containsKey(name)) {
				kinaseMap.get(name).addSubstrate(target);
			}
			else {
				kinaseMap.put(name, new Kinase(name, target));
			}
		}
		
		// read KEA ranks
		ArrayList<String> kea_ranks = assembleBackgroundRanks();
		for (String kinase_rank : kea_ranks) {
			String[] splitLine = kinase_rank.split("\\s");
			if (kinaseMap.containsKey(splitLine[0]))
				kinaseMap.get(splitLine[0]).setRankStats(Double.parseDouble(splitLine[1]), Double.parseDouble(splitLine[2]));
		}
		
		kinases = new LinkedList<Kinase>(kinaseMap.values());
		
		// filter substrates from input list that are not associated with an upstream kinase
		Set<String> substrateInputSet = new HashSet<String>();
		for (String substrate : genes) {
			substrateInputSet.add(substrate.toUpperCase());
		}
		Set<String> substrateBgSet = new HashSet<String>();
		for (Kinase kinase : kinases) {
			substrateBgSet.addAll(kinase.getSubstrates());
		}
		substrateInputSet.retainAll(substrateBgSet);
		
		ListIterator<Kinase> kinaseIterator = kinases.listIterator();
		
		while(kinaseIterator.hasNext()) {
			Kinase currentKinase = kinaseIterator.next();
			
			Set<String> targetBgSubstrates = currentKinase.getSubstrates();
			// Target input substrates is the intersection of target background substrates and input substrates
			Set<String> targetInputSubstrates = SetOps.intersection(targetBgSubstrates, substrateInputSet);
			
			double numOfTargetBgSubstrates = targetBgSubstrates.size();
			double totalBgSubstrates = substrateBgSet.size();
			double totalInputSubstrates = substrateInputSet.size();
			double numOfTargetInputSubstrates = targetInputSubstrates.size();

			if (targetInputSubstrates.size() > 0) {

				FisherExact fisher = new FisherExact(targetInputSubstrates.size()
						+ (substrateInputSet.size() - targetInputSubstrates.size())
						+ targetBgSubstrates.size()
						+ (substrateBgSet.size() - targetBgSubstrates.size()));
				
				double pvalue = fisher.getRightTailedP(targetInputSubstrates.size(),
						(substrateInputSet.size() - targetInputSubstrates.size()), targetBgSubstrates.size(), 
						(substrateBgSet.size() - targetBgSubstrates.size()));
				
				currentKinase.setEnrichedSubstrates(targetInputSubstrates);
				currentKinase.setFractionOfSubstratesInInput(numOfTargetInputSubstrates/totalInputSubstrates);
				currentKinase.setFractionOfSubstratesInBackground(numOfTargetBgSubstrates/totalBgSubstrates);
				currentKinase.setPValue(pvalue);
			}
			else {
				kinaseIterator.remove();
			}
		}
		
		// First, sort by p-value
		Collections.sort(kinases);
		
		// Count current rank and compute z-score
		int counter = 1;
		for (Kinase kinase : kinases) {
			kinase.computeScore(counter);
			counter++;
		}
		
		if (settings.get(SORT_BY).equals(COMBINED_SCORE)) {
			// Sort by combined score
			Collections.sort(kinases, new Comparator<Kinase>() {
				public int compare(Kinase o1, Kinase o2) {
					if (o1.getCombinedScore() < o2.getCombinedScore())				
						return 1;
					else if (o1.getCombinedScore() > o2.getCombinedScore())
						return -1;
					else
						return 0;
				}
			});
		}
		else if (settings.get(SORT_BY).equals(RANK)) {
			// Sort by z-score
			Collections.sort(kinases, new Comparator<Kinase>() {
				public int compare(Kinase o1, Kinase o2) {
					if (o1.getZScore() > o2.getZScore())				
						return 1;
					else if (o1.getZScore() < o2.getZScore())
						return -1;
					else
						return 0;
				}
			});
		}
	}
	
}
