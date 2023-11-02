package knapsack;

import java.util.Comparator;

public class Item implements Comparable<Item>,Cloneable{
public Item(int index, int weight, int utility) {
		super();
		this.index = index;
		this.weight = weight;
		this.utility = utility;
	}
public Item(int index, double weight, double utility) {
	super();
	this.index = index;
	this.weight = weight;
	this.utility = utility;
}
private int index;
private double weight;
private double utility;
private double density;
private double rescaledWeight;

public double getRescaledWeight() {
	return rescaledWeight;
}

public void setRescaledWeight(double rescaledWeight) {
	this.rescaledWeight = rescaledWeight;
}

public double getDensity() {
	return density;
}

public void setDensity(double density) {
	this.density = density;
}

public int getId() {
	return index;
}

public void setId(int id) {
	this.index = id;
}

public double getWeight() {
	return weight;
}

public void setWeight(double weight) {
	this.weight = weight;
}

public double getUtility() {
	return utility;
}

public void setUtility(double utility) {
	this.utility = utility;
}

@Override
public int compareTo(Item o) {
	
	return Double.compare(this.getUtility(), o.getUtility());
}
public static class Comparators {//used when sorting item lists
    public static final Comparator<Item> density = (Item o1, Item o2) -> Double.compare(o1.getDensity(),o2.getDensity());
    public static final Comparator<Item> utility = (Item o1, Item o2) -> Double.compare(o1.getUtility(),o2.getUtility());
    
}

}
