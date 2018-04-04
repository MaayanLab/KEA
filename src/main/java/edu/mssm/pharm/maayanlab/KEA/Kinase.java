package edu.mssm.pharm.maayanlab.KEA;

import java.util.Set;
import java.util.HashSet;

import com.google.gson.annotations.Expose;
public class Kinase implements Comparable<Object> {

	@Expose
	private String name;
	private Set<String> substrates = new HashSet<String>();
	
	private double mean;
	private double standardDeviation;
	
	private double fractionOfSubstratesInInput;
	private double fractionOfSubstratesInBackground;
	@Expose
	private double pvalue;
	@Expose
	private double zscore;
	@Expose
	private double combinedScore;
	
	@Expose
	private Set<String> enrichedSubstrates;
	
	public Kinase(
			String name, 
			String substrate) {
		
		this.name = name;
		substrates.add(substrate);
		
	}

	public String getName() {
		return this.name;
	}
	
	public void addSubstrate(String substrate) {
		substrates.add(substrate);
	}
	
	public Set<String> getSubstrates() {
		return substrates;
	}
	
	public void setRankStats(double mean, double standardDeviation) {
		this.mean = mean;
		this.standardDeviation = standardDeviation;
	}
	
	public Set<String> getEnrichedSubstrates() {
		return enrichedSubstrates;
	}
	
	public void setEnrichedSubstrates(Set<String> enrichedSubstrates) {
		this.enrichedSubstrates = enrichedSubstrates;
	}
	
	public void setFractionOfSubstratesInInput(double fractionOfSubstratesInInput) {
		this.fractionOfSubstratesInInput = fractionOfSubstratesInInput;
	}
	
	public void setFractionOfSubstratesInBackground(double fractionOfSubstratesInBackground) {
		this.fractionOfSubstratesInBackground = fractionOfSubstratesInBackground;
	}
	
	public double getPValue() {
		return this.pvalue;
	}
	
	public void setPValue(double pvalue) {
		this.pvalue = pvalue;
	}
	
	public double getZScore() {
		return this.zscore;
	}
	
	public double getCombinedScore() {
		return this.combinedScore;
	}
	
	/*
	1. name of the kinase
	2. number of substrates in the input gene-list
	3. number of genes that are substrates of the kinase
	4. the fraction of genes that are substrates compared to total number of genes in gene-list
	5. the fraction of genes that are substrates compared to total number of genes in background
	6. difference between the background fraction and the substrate-list fraction
	7. p-value computed using the Fisher Test
	8. rank computed using z-test
	9. combined score computed from p-value and rank
	10. list of substrates separated by a semi-colon
	 */	
	@Override
	public String toString() {
		StringBuilder outputString = new StringBuilder();
		outputString.append(name).append(",");
		outputString.append(enrichedSubstrates.size()).append(",");
		outputString.append(substrates.size()).append(",");
		outputString.append(fractionOfSubstratesInInput).append(",");
		outputString.append(fractionOfSubstratesInBackground).append(",");
		outputString.append((fractionOfSubstratesInInput - fractionOfSubstratesInBackground)).append(",");
		outputString.append(pvalue).append(",").append(zscore).append(",").append(combinedScore).append(",");
		
		boolean firstSubstrate = true;
		for (String enrichedSubstrate : enrichedSubstrates) {
			if (firstSubstrate) {
				outputString.append(enrichedSubstrate);
				firstSubstrate = false;
			}
			else
				outputString.append(";").append(enrichedSubstrate);
		}
		
		return outputString.toString();
	}

	@Override
	public int compareTo(Object o) {
		if (this.pvalue > ((Kinase) o).pvalue)			
			return 1;
		else if (this.pvalue < ((Kinase) o).pvalue)
			return -1;
		else
			return 0;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Kinase other = (Kinase) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		return true;
	}

	public void computeScore(int currentRank) {
		if (mean == 0 && standardDeviation == 0)
			zscore = 0;
		else
			zscore = (currentRank - mean)/standardDeviation;
		combinedScore = Math.log(pvalue)*zscore;
	}
}
