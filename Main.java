package knapsack;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;



public class Main {

	private static int n = 1000;
	static int m = (int) (n / Math.exp(1)); // Number of steps for sampling
	private static int C = 10; // knapsack capacity and size of sampling set

	static double β = 0.262;
	static int k1 = 112;
	static double ratioOfCapacity = 4;

	static double c = 0.42291; // used in in small large algorithm
	static double d = 0.6457; // used in small large algorithm

	public static void main(String[] args) {
		//getResultsFromKnapsackAlgorithmFromFirstArticle();
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		//List<Item> itemList = knapsack.createWeaklyCorrelatedInstance(100, n);
		getResultsFromSecondAlgorithm();
		//getResultsFromKnapsackAlgorithmFromThirdArticle();
		//List<Item> itemList= knapsack.uncorrelatedInstancesWithSimilarWeights(100, 1000);
		//try {
		//	exportArrayOfItemsToCsv(itemList, "C:\\Users\\Kotsos\\Desktop\\uncorelatedInstanceWIthSImilarWeights.csv");
		//} catch (IOException e) {
		//	// TODO Auto-generated catch block
		//	e.printStackTrace();
		//}
		//getResultsFromSecondAlgorithmAlternative();
	}
	

	    public static void exportArrayOfItemsToCsv(List<Item> itemList, String filename) throws IOException {
	        FileWriter writer = new FileWriter(filename);

	        // Write CSV header
	        writer.append("Weight,Utility\n");

	        // Write each item as a CSV row
	        for (Item item : itemList) {
	            writer.append(Double.toString(item.getWeight())+","+Double.toString(item.getUtility()));
	            writer.append("\n");
	        }

	        writer.flush();
	        writer.close();
	    }

	    
	    
	
	public static void getResultsFromSecondAlgorithm() {
		double resultArray[][] = new double[n][3];
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		for (int i = 0; i <= 49; i++) {
			List<Item> itemList = knapsack.uncorrelatedInstancesWithSimilarWeights(100, n);
			resultArray[i][0] = solutionOfKnapsackProblemWithCostInsteadOfCapacity(itemList);
			resultArray[i][1] = calculateProfitOfFractionalOptimalSolution(itemList);
		}
		exportResultsToCsvFromSecondArticle(resultArray, "ρανδομ2");
	}
	
	public static void getResultsFromSecondAlgorithmAlternative() {
		double resultArray[][] = new double[n][2];
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		for (int i = 0; i <= n-1; i++) {
			List<Item> itemList = knapsack.createRandomInstance(10, n);
			resultArray[i][0] = solutionOfKnapsackProblemWithCostInsteadOfCapacityAlternative(itemList);
			double capacityOfKnapsack = knapsack.calculateWeightRatioOfItemList(ratioOfCapacity, itemList);
			resultArray[i][1] = knapsack.calculateBestPossibleSolutionWithGivenKnapsackLimit(itemList,capacityOfKnapsack);
			
		}
		exportResultsToCsvFromSecondArticleAlternative(resultArray, "createRandomInstance");
	}
	public static void getResultsFromThirdAlgorithm() {
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		double results[][] = new double[n][3];
		
		for(int i =0;i<n-1;i++) {
			
			List<Item> itemList =knapsack.createStronglyCorrelatedInstance(100, n);
			
			results[i] = knapsackSolutionDualAlgorithmSmallLarge(itemList,ratioOfCapacity);//uses global ratio
		}
		exportResultsToCsv(results, ratioOfCapacity, "stronglyCorrelated");
		
	}
	public static double[] knapsackSolutionDualAlgorithmSmallLarge(List<Item> itemList, double ratio) {// accepts list of items returns
																						// best possible solution
		double resultArray[] = new double[3];
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		List<Item> sampleList = new ArrayList<>();// contains items of sampling phase
		List<Item> solutionList = new ArrayList<>();// output of algorithm
		List<Item> runList = new ArrayList<>();// contains items that have been seen by algorithm

		double capacityOfKnapsack = knapsack.calculateWeightRatioOfItemList(ratio, itemList);

		double maxV = itemList.get(0).getUtility(); // initialising maxValue
		for (int l = 0; l <= n; l++) {
			runList.add(itemList.get(l));
			if (l <= c * n) {// first sampling phase
				sampleList.add(itemList.get(l));
				if (itemList.get(l).getUtility() > maxV)
					maxV = itemList.get(l).getUtility(); // keeping track of max value seen so far
			}
			if ((c * n + 1) <= l && l <= (d * n)) {
				if (algorithmForLargeItems(solutionList, capacityOfKnapsack, itemList.get(l), maxV))
					solutionList.add(itemList.get(l));
			}
			if ((d * n + 1) <= l && l <= (n)) {
				if (algorithmForSmallItems(solutionList, capacityOfKnapsack, itemList.get(l), runList))
					solutionList.add(itemList.get(l));
			}
		}
		try {
			knapsack.exportToGlpkData(itemList, capacityOfKnapsack);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		resultArray[0] = knapsack.getItemListAccumulativeUtility(solutionList);
		resultArray[1] = knapsack.calculateBestPossibleSolutionWithGivenKnapsackLimit(itemList,capacityOfKnapsack);
		knapsack.runGlpkWithGivenFileLocation("C:\\Users\\Kotsos\\Desktop\\ilpknap",
		"C:\\Users\\Kotsos\\Desktop\\filename.dat");
		resultArray[2] = knapsack.readGLPKResult();

		return resultArray;

	}

	private static boolean algorithmForSmallItems(List<Item> solutionList, double capacityOfKnapsack, Item item,
			List<Item> runList) {// judges wether or not Item should be added to solutionList
		double probabilityOfAddingItem = probabilityOfAddingItem(runList, capacityOfKnapsack, item);
		if (probabilityOfAddingItem == 0) {
			return false;
		} else if (probabilityOfAddingItem == 1) {
			return true;
		} else {
			double randNum = Math.random();
			if (randNum > probabilityOfAddingItem)
				return false;
			if (randNum < probabilityOfAddingItem)
				return true;
		}
		return false;

	}

	private static double probabilityOfAddingItem(List<Item> runList, double capacityOfKnapsack, Item newItem) {
		List<Item> fractionalOptimalSolution = new ArrayList<>();
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		Collections.sort(runList);
		Collections.reverse(runList);
		for (Item item : runList) {
			if (item.getWeight() + knapsack.calculateWeightSum(fractionalOptimalSolution) < capacityOfKnapsack) {
				fractionalOptimalSolution.add(item);
				if (item == newItem)
					return 1; // the entire new item was addedto the optimal fractional solution therefore we
								// can safely add it to final solution list

			} else {// the entire item doesnt fit as is so we must return thhe percent of the item
					// that fits
				double capacityLeft = capacityOfKnapsack - knapsack.calculateWeightSum(fractionalOptimalSolution);
				double capacityLeftRatio = (capacityLeft / capacityOfKnapsack) * 100;
				return capacityLeftRatio;
			}
		}
		return 0;

	}

	private static boolean algorithmForLargeItems(List<Item> solutionList, double capacityOfKnapsack, Item item,
			double maxV) {// judges wether or not 1/3-large Item should be added to solutionList
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);// making knapsack instance to use methods

		if (item.getWeight() < capacityOfKnapsack / 3) {
			return false;
		}
		if (item.getWeight() + knapsack.calculateWeightSum(solutionList) < capacityOfKnapsack
				&& item.getUtility() > maxV) {
			return true;
		}
		return false;

	}

	public static double solutionOfKnapsackProblemWithCostInsteadOfCapacity(List<Item> itemList) {
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		castCostIntoInt(itemList); //our cost function goes fro Z+ to Z+ therefore all costs must be positive integers
		List<Item> sampleList = new ArrayList<>();
		List<Item> solutionList = new ArrayList<>();
		double capacity = knapsack.calculateWeightRatioOfItemList(ratioOfCapacity, itemList);
		knapsack.setW(capacity);
		knapsack.setRescaledWeight(itemList, capacity);
		double randomNumber = 0 ;

		if (randomNumber <= 0.5) {// Classic dynkin secretary solution . Sample first n/e elements. and take the
									// first element after that with bigger density
			for (int i = 0; i <= n; i++) {
				if (i <= m) {
					sampleList.add(itemList.get(i));
					continue;
				}
				Collections.sort(sampleList);
				Collections.reverse(sampleList);
				if (itemList.get(i).getDensity() > sampleList.get(0).getDensity() )
						 {
					solutionList.add(itemList.get(i));
					break;
				}
			}

		} else {// solution according to 2012 article
			int k = knapsack.generateRandomNumberUsingBinomialDistribution(itemList.size());// number of elements in
																							// sample

			int counter = 0;
			double rdensity = 0;

			for (Item item : itemList) {
				if (counter++ < k) {// sampling
					sampleList.add(item);
					continue;
				}
				if (counter == k + 1) {
					rdensity = calculateR(sampleList);
				}
				if (item.getDensity() > rdensity && isItemProfitable(solutionList, item)) {
					solutionList.add(item);

				}
			}
		}

		return calculateProfit(solutionList);

	}
	public static double solutionOfKnapsackProblemWithCostInsteadOfCapacityAlternative(List<Item> itemList) {
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
		castCostIntoInt(itemList); //our cost function goes fro Z+ to Z+ therefore all costs must be positive integers
		List<Item> sampleList = new ArrayList<>();
		List<Item> solutionList = new ArrayList<>();
		double capacity = knapsack.calculateWeightRatioOfItemList(ratioOfCapacity, itemList);
		knapsack.setW(capacity);
		knapsack.setRescaledWeight(itemList, capacity);
		double randomNumber = Math.random() ;

		if (randomNumber <= 0.5) {// Classic dynkin secretary solution . Sample first n/e elements. and take the
									// first element after that with bigger density
			for (int i = 0; i <= n; i++) {
				if (i <= m) {
					sampleList.add(itemList.get(i));
					continue;
				}
				Collections.sort(sampleList);
				Collections.reverse(sampleList);
				if (itemList.get(i).getDensity() > sampleList.get(0).getDensity() )
						 {
					solutionList.add(itemList.get(i));
					break;
				}
			}

		} else {// solution according to 2012 article
			int k = knapsack.generateRandomNumberUsingBinomialDistribution(itemList.size());// number of elements in
			double w ;																			// sample

			int counter = 0;
			double rdensity = 0;

			for (Item item : itemList) {
				if (counter++ < k) {// sampling
					sampleList.add(item);
					continue;
				}
				if (counter == k + 1) {
					rdensity = calculateR(sampleList);
				}
				if (item.getDensity() > rdensity && doesItemFit(solutionList,item,knapsack)) {
					solutionList.add(item);

				}
			}
		}

		return knapsack.calculateUtilitytSum(solutionList);

	}
	private static boolean doesItemFit(List<Item> solutionList, Item item, WeightedKnapsackProblem knapsack) {
		double w =knapsack.calculateWeightSum(solutionList);
		if(w+item.getWeight()<=knapsack.getW()) return true;
		return false;
	}


	private static void castCostIntoInt(List<Item> itemList) {//casts all costs of given list into integers
	for(Item item :itemList) {
		item.setWeight((int) item.getWeight());
	}
		
	}

	private static double calculateR(List<Item> sampleList) {// calculetes density r for
																// algorithm
		double maxR;
		maxR = 0;
		List<Item> cloneList = new ArrayList<>( sampleList);
		double profitOfFractionSolution = calculateProfitOfFractionalOptimalSolution(cloneList);
		for (Item item : sampleList) {
			List<Item> listWithElementsThatHaveDensityBiggerThanmaxR = calculateSubListOfElementsWithHigherThanSpecifiedDensity(
					sampleList, item.getDensity());
			double profit = calculateProfit(listWithElementsThatHaveDensityBiggerThanmaxR);
			double x = β*(1 - (1 / k1)) * profitOfFractionSolution;
			if (profit >x && item.getDensity() > maxR) {
				maxR = item.getDensity();
			}
		}
		return maxR;
	}

	private static double calculateProfitOfFractionalOptimalSolution(List<Item> sampleList) {// finds
																								// calculate fractional
																								// optimal
																								// solution of
																								// list and
																								// returns its
																								// profit as
																								// defined by
																								// 2012 article
		List<Item> cloneList = new ArrayList<>( sampleList);
		Collections.sort(cloneList);
		Collections.reverse(cloneList);
		double profit = 0;
		List<Item> listOfFractionalOptimalSolution = new ArrayList<>();
		process:
		for (Item item : cloneList) {

			if (isItemProfitable(listOfFractionalOptimalSolution, item)) {

				profit += item.getUtility() - costFunction(item.getWeight());

				listOfFractionalOptimalSolution.add(item);
				continue;
			} else {// decrease by 1 untill item is profitable
				int counter = 0;
				
						while (true) {
						counter++; //incrementally increasing counter
						item.setWeight(item.getWeight() - counter); //decreasing weight by one temporarily
						double ratio = item.getWeight() * 100 / (item.getWeight() + counter);// calculating the ratio of the difference between original weight and substracted weight
						double utilityDifference = item.getUtility() - item.getUtility() * ratio / 100; //utility difference will be used to revert any changes that were made to the utility of the item
						item.setUtility(item.getUtility() * ratio / 100);
						if (isItemProfitable(listOfFractionalOptimalSolution, item)) {
							item.setUtility(utilityDifference + item.getUtility());
							profit = +item.getUtility();
							item.setWeight(item.getWeight() + counter);
							continue;
						} else if ((int) (item.getWeight() ) == 0) {
							item.setUtility(utilityDifference + item.getUtility());
							item.setWeight(item.getWeight() + counter);
							continue process;
						} else {
							item.setUtility(utilityDifference + item.getUtility());
							item.setWeight(item.getWeight() + counter);
							continue;
						}
					}

			}

		}
		return profit;
	}

	private static double costFunction(double x) {
		return Math.pow(x, 2);
	}

	private static boolean isItemProfitable(List<Item> itemList, Item item) {// returns true if adding item in list will
																				// result in higher profit
		List<Item> copyList = new ArrayList<>(itemList);// creating shallow copy
		double preItemProfit = calculateProfit(copyList);
		copyList.add(item);
		double PostItemProfit = calculateProfit(copyList);
		copyList.remove(item);
		if (preItemProfit < PostItemProfit) {
			return true;
		} else {
			return false;
		}

	}

	private static double calculateProfit(List<Item> itemList) {// Returns profit of given list (profit as defined by
																// 2012 online knapsack article
		double utilitySum = 0;
		double weightSum = 0;
		for (Item item : itemList) {
			utilitySum += item.getUtility();
			weightSum += item.getWeight();
		}
		double profit = utilitySum - Math.pow((int) weightSum, 2);
		return profit;
	}

	private static List<Item> calculateSubListOfElementsWithHigherThanSpecifiedDensity(List<Item> sampleList,
			double maxR) {
		List<Item> listWithElementsThatHaveDensityBiggerThanmaxR = new ArrayList<>();
		for (Item item : sampleList) {
			if (maxR < item.getDensity()) {
				listWithElementsThatHaveDensityBiggerThanmaxR.add(item);
			}
		}
		return listWithElementsThatHaveDensityBiggerThanmaxR;
	}

//Is used to export results to a csv, array consists of 
	private static void exportResultsToCsv(double[][] resultArray, double d, String name) {
		int x = resultArray.length;
		String text = "n,CapacityRatio,Solution,OptimumSolution,GLPKSolution\n";
		for (int i = 0; i < x; i++) {
			text += x + "," + d + "," + resultArray[i][0] + "," + resultArray[i][1] + "," + resultArray[i][2] + "\n";
		}

		Path fileName = Path.of("C:\\\\Users\\Kotsos\\\\Desktop\\knapsackdataThirdArticle2\\" + name + "=" + d + "n=" + x + ".xlsx");
		try {
			Files.writeString(fileName, text);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void exportResultsToCsvFromSecondArticle(double[][] resultArray, String name) {// resultsArray[0] is
																									// results from
																									// algorithm,results[1]
																									// is optimal
																									// fractional result
		int x = resultArray.length;
		String text = "n,Solution,OptimumSolution\n";
		for (int i = 0; i < 49; i++) {
			text += x + "," + resultArray[i][0] + "," + resultArray[i][1] + "\n";
		}

		Path fileName = Path
				.of("C:\\\\Users\\Kotsos\\\\Desktop\\knapsackdataSecondArticle\\" + name + "n=" + x +".csv");
		try {
			Files.writeString(fileName, text);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	private static void exportResultsToCsvFromSecondArticleAlternative(double[][] resultArray, String name) {
		int x = resultArray.length;
		String text = "n,Solution,optimum solution\n";
		for (int i = 0; i < n-1; i++) {
		text += x + "," + resultArray[i][0]  + ","+ resultArray[i][1]+ "\n";
		}
		
		Path fileName = Path
		.of("C:\\\\Users\\Kotsos\\\\Desktop\\knapsackdataSecondArticleAlternative\\" + name + "n=" + x + "RoC"+ratioOfCapacity+ ".csv");
		try {
		Files.writeString(fileName, text);
		} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
}
}
	public static double[] doALot(double R, List<Item> stronglyCorrelatedItemList, WeightedKnapsackProblem knapsack) {
		double resultArray[] = new double[3];

		List<Item> finalSetOfWeightedSolution = new ArrayList<>(); // set-solution of virtual method
		List<Item> samplingItemList = new ArrayList<>(); // Stores sampling set

		int counter = 0;

		int a = knapsack.generateA();
		int k = knapsack.determineK(a);

		// solution for stronglyCorrelatedInstance

		double C = knapsack.calculateWeightRatioOfItemList(R, stronglyCorrelatedItemList);// summing the weights of
																							// all elements we take
																							// the ratio and set it
																							// as knapsack capacity
		knapsack.setW(C);
		knapsack.setRescaledWeight(stronglyCorrelatedItemList, C);
		counter = 0;
		for (Item item : stronglyCorrelatedItemList) {
			if (counter++ < m) {
				knapsack.sampleItem(item, k, samplingItemList);
				continue;
			}

			int w = knapsack.getItemListAccumulativeRescaledWeight(finalSetOfWeightedSolution);// current weight of
																								// knapsack

			if (knapsack.evaluateItemForKnapsack(item, a, w, k, samplingItemList)) {
				finalSetOfWeightedSolution.add(item);
			}
			if (counter == n) {
				break;
			}
		}

		try {
			knapsack.exportToGlpkData(stronglyCorrelatedItemList, C);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		resultArray[0] = knapsack.getItemListAccumulativeUtility(finalSetOfWeightedSolution);
		resultArray[1] = knapsack.calculateBestPossibleSolution(stronglyCorrelatedItemList);
		knapsack.runGlpkWithGivenFileLocation("C:\\Users\\Kotsos\\Desktop\\ilpknap",
				"C:\\Users\\Kotsos\\Desktop\\filename.dat");
		resultArray[2] = knapsack.readGLPKResult();

		return resultArray;

	}

	public static void getResultsFromKnapsackAlgorithmFromFirstArticle() {
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);
//
//		 repeatAndExportForWeaklyCorrelated(3, "weaklyCorrelated", knapsack);
//		 repeatAndExportForWeaklyCorrelated(2, "weaklyCorrelated", knapsack);
//		 repeatAndExportForWeaklyCorrelated(4, "weaklyCorrelated", knapsack);
//		 repeatAndExportForWeaklyCorrelated((double) (3) / 2, "weaklyCorrelated",
//		 knapsack);
//
		 repeatAndExportForStronglyCorrelated(3, "stronglyCorrelated", knapsack);
//		 repeatAndExportForStronglyCorrelated(2, "stronglyCorrelated", knapsack);
//		 repeatAndExportForStronglyCorrelated(4, "stronglyCorrelated", knapsack);
//		 repeatAndExportForStronglyCorrelated((double) 3 / 2, "stronglyCorrelated",
//		 knapsack);
//
//		 repeatAndExportForInverseStronglyCorrelated(3, "InversestronglyCorrelated",
//		 knapsack);
//		 repeatAndExportForInverseStronglyCorrelated(2, "InversestronglyCorrelated",
//		 knapsack);
//		 repeatAndExportForInverseStronglyCorrelated(4, "InversestronglyCorrelated",
//		 knapsack);
//		 repeatAndExportForInverseStronglyCorrelated((double) 3 / 2,
//		 "InversestronglyCorrelated", knapsack);
//
//		 repeatAndExportForsubSetSum(3, "SubSetSum", knapsack);
//		 repeatAndExportForsubSetSum(2, "SubSetSum", knapsack);
//		 repeatAndExportForsubSetSum(4, "SubSetSum", knapsack);
//		 repeatAndExportForsubSetSum((double) 3 / 2, "SubSetSum", knapsack);
//
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights(3, "uncorrelatedWithSimilarWeights", knapsack);
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights(2, "uncorrelatedWithSimilarWeights", knapsack);
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights(4, "uncorrelatedWithSimilarWeights", knapsack);
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights((double) 3 / 2, "uncorrelatedWithSimilarWeights",
//				knapsack);
//
//		repeatAndExportForRandom(3, "RandomInstance", knapsack);
//		repeatAndExportForRandom(2, "RandomInstance", knapsack);
//		repeatAndExportForRandom(4, "RandomInstance", knapsack);
//		repeatAndExportForRandom((double) 3 / 2, "RandomInstance", knapsack);
	}

	public static void repeatAndExportForWeaklyCorrelated(double ratio, String name, WeightedKnapsackProblem knapsack) {
		double results[][] = new double[n][3];
		for (int i = 0; i <=50 - 1; i++) {
			List<Item> itemList = knapsack.createWeaklyCorrelatedInstance(100, n);

			results[i] = doALot(ratio, itemList, knapsack);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForuncorrelatedInstancesWithSimilarWeights(double ratio, String name,
			WeightedKnapsackProblem knapsack) {
		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.uncorrelatedInstancesWithSimilarWeights(100, n);

			results[i] = doALot(ratio, itemList, knapsack);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForsubSetSum(double ratio, String name, WeightedKnapsackProblem knapsack) {

		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.subSetSumInstance(100, n);
			results[i] = doALot(ratio, itemList, knapsack);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForStronglyCorrelated(double ratio, String name,
			WeightedKnapsackProblem knapsack) {
		double results[][] = new double[n][3];
		for (int i = 0; i <= 50-1; i++) {
			List<Item> itemList = knapsack.createStronglyCorrelatedInstance(100, n);
			try {
				knapsack.exportToGlpkData(itemList, knapsack.calculateWeightRatioOfItemList(ratio, itemList));
			} catch (IOException e) {
				
				e.printStackTrace();
			}
			for (Item item : itemList) {
				System.out.println(item.getUtility());
				System.out.println(item.getWeight());
			}
			System.out.println(knapsack.calculateWeightRatioOfItemList(ratio, itemList));
			results[i] = doALot(ratio, itemList, knapsack);
		}
		
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForInverseStronglyCorrelated(double ratio, String name,
			WeightedKnapsackProblem knapsack) {

		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.createInverseStronglyCorrelatedInstance(100, n);
			results[i] = doALot(ratio, itemList, knapsack);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForRandom(double ratio, String name, WeightedKnapsackProblem knapsack) {

		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.createRandomInstance(100, n);
			results[i] = doALot(ratio, itemList, knapsack);
		}
		exportResultsToCsv(results, ratio, name);
	}
	public static void repeatAndExportForWeaklyCorrelated3(double ratio, String name, WeightedKnapsackProblem knapsack) {
		double results[][] = new double[n][3];
		for (int i = 0; i <=50 - 1; i++) {
			List<Item> itemList = knapsack.createWeaklyCorrelatedInstance(100, n);

			results[i] = knapsackSolutionDualAlgorithmSmallLarge( itemList,ratio);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForuncorrelatedInstancesWithSimilarWeights3(double ratio, String name,
			WeightedKnapsackProblem knapsack) {
		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.uncorrelatedInstancesWithSimilarWeights(100, n);

			results[i] = knapsackSolutionDualAlgorithmSmallLarge( itemList,ratio);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForsubSetSum3(double ratio, String name, WeightedKnapsackProblem knapsack) {

		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.subSetSumInstance(100, n);
			results[i] = knapsackSolutionDualAlgorithmSmallLarge( itemList,ratio);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForStronglyCorrelated3(double ratio, String name,
			WeightedKnapsackProblem knapsack) {
		double results[][] = new double[n][3];
		for (int i = 0; i <= 50-1; i++) {
			List<Item> itemList = knapsack.createStronglyCorrelatedInstance(100, n);
			try {
				knapsack.exportToGlpkData(itemList, knapsack.calculateWeightRatioOfItemList(ratio, itemList));
			} catch (IOException e) {
				
				e.printStackTrace();
			}
			for (Item item : itemList) {
				System.out.println(item.getUtility());
				System.out.println(item.getWeight());
			}
			System.out.println(knapsack.calculateWeightRatioOfItemList(ratio, itemList));
			results[i] = knapsackSolutionDualAlgorithmSmallLarge( itemList,ratio);
		}
		
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForInverseStronglyCorrelated3(double ratio, String name,
			WeightedKnapsackProblem knapsack) {

		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.createInverseStronglyCorrelatedInstance(100, n);
			results[i] = knapsackSolutionDualAlgorithmSmallLarge( itemList,ratio);
		}
		exportResultsToCsv(results, ratio, name);
	}

	public static void repeatAndExportForRandom3(double ratio, String name, WeightedKnapsackProblem knapsack) {

		double results[][] = new double[n][3];
		for (int i = 0; i <= 50 - 1; i++) {
			List<Item> itemList = knapsack.createRandomInstance(100, n);
			results[i] = knapsackSolutionDualAlgorithmSmallLarge( itemList,ratio);
		}
		exportResultsToCsv(results, ratio, name);
	}
	public static void getResultsFromKnapsackAlgorithmFromThirdArticle() {
		WeightedKnapsackProblem knapsack = new WeightedKnapsackProblem(n, C);

//		 repeatAndExportForWeaklyCorrelated3(3, "weaklyCorrelated", knapsack);
//		 repeatAndExportForWeaklyCorrelated3(2, "weaklyCorrelated", knapsack);
//		 repeatAndExportForWeaklyCorrelated3(4, "weaklyCorrelated", knapsack);
//		 repeatAndExportForWeaklyCorrelated3((double) (3) / 2, "weaklyCorrelated",
//		 knapsack);
//
	 repeatAndExportForStronglyCorrelated3(3, "stronglyCorrelated", knapsack);
//		 repeatAndExportForStronglyCorrelated3(2, "stronglyCorrelated", knapsack);
//		 repeatAndExportForStronglyCorrelated3(4, "stronglyCorrelated", knapsack);
//		 repeatAndExportForStronglyCorrelated3((double) 3 / 2, "stronglyCorrelated",
//		 knapsack);
//
//		 repeatAndExportForInverseStronglyCorrelated3(3, "InversestronglyCorrelated",
//		 knapsack);
//		 repeatAndExportForInverseStronglyCorrelated3(2, "InversestronglyCorrelated",
//		 knapsack);
//		 repeatAndExportForInverseStronglyCorrelated3(4, "InversestronglyCorrelated",
//		 knapsack);
//		 repeatAndExportForInverseStronglyCorrelated3((double) 3 / 2,
//		 "InversestronglyCorrelated", knapsack);
//
//		 repeatAndExportForsubSetSum3(3, "SubSetSum", knapsack);
//		 repeatAndExportForsubSetSum3(2, "SubSetSum", knapsack);
//		 repeatAndExportForsubSetSum3(4, "SubSetSum", knapsack);
//		 repeatAndExportForsubSetSum3((double) 3 / 2, "SubSetSum", knapsack);
//
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights3(3, "uncorrelatedWithSimilarWeights", knapsack);
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights3(2, "uncorrelatedWithSimilarWeights", knapsack);
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights3(4, "uncorrelatedWithSimilarWeights", knapsack);
//		repeatAndExportForuncorrelatedInstancesWithSimilarWeights((double) 3 / 2, "uncorrelatedWithSimilarWeights",
//				knapsack);
//
//		repeatAndExportForRandom3(3, "RandomInstance", knapsack);
//		repeatAndExportForRandom3(2, "RandomInstance", knapsack);
//		repeatAndExportForRandom3(4, "RandomInstance", knapsack);
		//repeatAndExportForRandom3((double) 3 / 2, "RandomInstance", knapsack);
	}
}
