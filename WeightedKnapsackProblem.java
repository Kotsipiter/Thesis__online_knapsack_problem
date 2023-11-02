package knapsack;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ThreadLocalRandom;

public class WeightedKnapsackProblem {
	public WeightedKnapsackProblem(int n, int C) {
		super();
		this.n = n;
		W = C;

	}

	private int n;
	private static double W;

	public int generateA() { // generates a ∈(1..4)
		int a = ThreadLocalRandom.current().nextInt(0, 5);
		return a;
	}

	public int determineK(int a) { // Returns size of sampling set depending on a, the case of a=4 is special
		int k;
		if (a < 4) {
			k = (int) Math.pow(3, a);
			return k;
		} else {
			k = generateRandomNumberUsingBinomialDistribution(n);
			return k;
		}

	}

	public int generateRandomNumberUsingBinomialDistribution(int x) { // Returns the number of heads when spinning a
																		// fair coin x times
		int count = 1;
		for (int i = 1; i <= (x - 1); i++) {
			double randomNumber = Math.random();
			if (randomNumber >= 0.5)
				count++;
		}
		return count;
	}

	public void sampleItem(Item item, int k, List<Item> samplingItemList) {
		if (samplingItemList.size() < k) { // untill k items are sampled, all arriving items are accepted to sample set
			samplingItemList.add(item);
			Collections.sort(samplingItemList);
			Collections.reverse(samplingItemList);
		} else {
			if (samplingItemList.get(samplingItemList.size() - 1).getUtility() > item.getUtility()) { // Replaces item
																										// in sample set
																										// if more
																										// valuable
				samplingItemList.remove(samplingItemList.size() - 1);
				samplingItemList.add(item);
				Collections.sort(samplingItemList);
				Collections.reverse(samplingItemList);
			}
		}
	}

	public double calculateDensityThreshhold(List<Item> samplingItemList) {
		double counter = 0;
		Collections.sort(samplingItemList, Item.Comparators.density);
		Collections.reverse(samplingItemList);
		for (Item item : samplingItemList) {
			if (counter + item.getDensity() <= 1 / 2) {
				counter += item.getDensity();
			}
		}
		return counter;

	}

	public boolean evaluateItemForKnapsack(Item contenderItem, int a, int w, int k, List<Item> samplingItemList) {
		if (a < 4) {
			if (contenderItem.getUtility() > samplingItemList.get(samplingItemList.size() - 1).getUtility()
					&& (w + contenderItem.getRescaledWeight()) <= 1
					&& contenderItem.getRescaledWeight() < ((float) (1) / k)) { // comparing new
				// item and
				// least
				// valuable item
				// of
				// samplelist

				samplingItemList.remove(samplingItemList.size() - 1);
				samplingItemList.add(contenderItem);

				Collections.sort(samplingItemList);
				Collections.reverse(samplingItemList);// Is required since sort returns increasing list

				return true;

			}
		}
		if (a == 4) {
			double p = calculateDensityThreshhold(samplingItemList);
			if (contenderItem.getUtility() > p && (w + contenderItem.getRescaledWeight()) <= 1) { // comparing new
																									// item and
																									// least
																									// valuable item
																									// of samplelist

				samplingItemList.remove(samplingItemList.size() - 1);
				samplingItemList.add(contenderItem);

				Collections.sort(samplingItemList);
				Collections.reverse(samplingItemList);// Is required since sort returns increasing list

				return true;

			}
			return false;
		}
		return false;
	}

	public double getW() {
		return W;
	}

	public void setW(double w) {
		W = w;
	}

	public int getN() {
		return n;
	}

	public void setN(int n) {
		this.n = n;
	}

	public Item generateWeightedItem(int counter) {// creates item
		int weight = determineWeight();
		int utility = determineUtility();
		Item item = new Item(counter, weight, utility);
		item.setDensity((float) (utility) / weight);
		item.setRescaledWeight((float) (weight) / W);
		return item;
	}

	private static int determineUtility() {
		int min = 1;
		int max = 100;
		int utility = ThreadLocalRandom.current().nextInt(min, max + 1); // nextInt is exclusive of top value so add +1
		return utility;
	}

	private static int determineWeight() {
		int min = 1;
		int max = 10;
		int weight = ThreadLocalRandom.current().nextInt(min, max + 1); // nextInt is
		// exclusive of top value so add +1
		return weight;

	}

	int getItemListAccumulativeRescaledWeight(List<Item> items) {
		int counter = 0;
		if (items.size() == 0)
			return 0;
		for (Item item : items) {
			counter += item.getRescaledWeight();
		}
		return counter;
	}

	public double calculateBestPossibleSolution(List<Item> afterSamplingItemList) {
		double utilityCounter = 0;
		double weightCounter = 0;
		Collections.sort(afterSamplingItemList, Item.Comparators.density);
		Collections.reverse(afterSamplingItemList);
		for (Item item : afterSamplingItemList) {
			if (weightCounter + item.getRescaledWeight() < 1) {
				utilityCounter += item.getUtility();
				weightCounter += item.getRescaledWeight();
			}
		}
		return utilityCounter;
	}
	public double calculateBestPossibleSolutionWithGivenKnapsackLimit(List<Item> afterSamplingItemList,double knapsackLimit) {
		double utilityCounter = 0;
		double weightCounter = 0;
		Collections.sort(afterSamplingItemList, Item.Comparators.density);
		Collections.reverse(afterSamplingItemList);
		for (Item item : afterSamplingItemList) {
			if (weightCounter + item.getRescaledWeight() < knapsackLimit) {
				utilityCounter += item.getUtility();
				weightCounter += item.getRescaledWeight();
			}
		}
		return utilityCounter;
	}
	public int getItemListAccumulativeUtility(List<Item> items) {
		int sum = 0;
		for (Item item : items) {
			sum += item.getUtility();
		}
		return sum;
	}

	public List<Item> readCSVAndCreateItemSet(String csvFile) {// asks for file location of a csv and then uses csv to
																// create file
		final String delimiter = ",";
		List<Item> itemSet = new ArrayList<>();
		try {
			File file = new File(csvFile);
			FileReader fr = new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line = " ";
			String[] tempArr;
			int counter = 0;
			while ((line = br.readLine()) != null) {
				tempArr = line.split(delimiter);
				int weight = Integer.parseInt(tempArr[0]);
				int utility = Integer.parseInt(tempArr[1]);
				Item item = generateWeightedItemWithGivenValues(counter++, utility, weight);
				itemSet.add(item);
			}
			br.close();

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		return itemSet;

	}

	public static Item generateWeightedItemWithGivenValues(int counter, double u, double w) {// creates item with given
																								// values of weight and
																								// utility

		Item item = new Item(counter, w, u);
		item.setDensity((float) (u) / w);
		item.setRescaledWeight((float) (w) / W);
		return item;
	}

	public List<Item> createWeaklyCorrelatedInstance(int R, int S) {// Weights wj are chosen randomly in [1; R] and the
																	// values pj in [wj − R=10; wj + R=10] such that pj
																	// ¿ 1
		List<Item> weaklyCorrelatedItemList = new ArrayList<>();
		for (int i = 0; i <= S; i++) {
			double weight = ThreadLocalRandom.current().nextDouble(1, R + 1);
			double utility = ThreadLocalRandom.current().nextDouble(weight - ((double) (R) / 10),
					weight + ((double) (R) / 10) + 1);
			if (utility < 1)
				utility = 1;
			Item item = generateWeightedItemWithGivenValues(i, utility, weight);
			weaklyCorrelatedItemList.add(item);
		}

		return weaklyCorrelatedItemList;

	}

	public List<Item> createStronglyCorrelatedInstance(int R, int S) {// Weights wj are chosen randomly in [1; R] and
																		// the pro6ts pj in [wj − R=10; wj + R=10] such
																		// that pj ¿ 1
		List<Item> StronglyCorrelatedItemList = new ArrayList<>();
		for (int i = 0; i <= S; i++) {
			double weight = ThreadLocalRandom.current().nextDouble(1, R + 1);
			double utility = weight + (double) (R) / 10;

			Item item = generateWeightedItemWithGivenValues(i, utility, weight);
			StronglyCorrelatedItemList.add(item);
		}

		return StronglyCorrelatedItemList;

	}

	public List<Item> createInverseStronglyCorrelatedInstance(int R, int S) {// Weights wj are chosen randomly in [1; R]
																				// and the pro6ts pj in [wj − R=10; wj +
																				// R=10] such that pj ¿ 1
		List<Item> InverseStronglyCorrelatedItemList = new ArrayList<>();
		for (int i = 0; i <= S; i++) {
			double utility = ThreadLocalRandom.current().nextDouble(1, R + 1);
			double weight = utility + (double) (R) / 10;

			Item item = generateWeightedItemWithGivenValues(i, utility, weight);
			InverseStronglyCorrelatedItemList.add(item);
		}

		return InverseStronglyCorrelatedItemList;

	}

	public List<Item> subSetSumInstance(int R, int S) {// Weights wj are chosen randomly in [1; R] and the pro6ts pj in
														// [wj − R=10; wj + R=10] such that pj ¿ 1
		List<Item> subSetSumItemList = new ArrayList<>();
		for (int i = 0; i <= S; i++) {
			double weight = ThreadLocalRandom.current().nextDouble(1, R + 1);
			double utility = weight;

			Item item = generateWeightedItemWithGivenValues(i, utility, weight);
			subSetSumItemList.add(item);
		}

		return subSetSumItemList;

	}

	public List<Item> uncorrelatedInstancesWithSimilarWeights(int R, int S) {// Weights wj are chosen randomly in [1; R]
																				// and the pro6ts pj in
		// [wj − R=10; wj + R=10] such that pj ¿ 1
		List<Item> uncorrelatedInstancesWithSimilarWeightsItemList = new ArrayList<>();
		for (int i = 0; i <= S; i++) {
			// double weight = ThreadLocalRandom.current().nextDouble(100000, 100100 + 1);
			// double utility = ThreadLocalRandom.current().nextDouble(1, 1000 + 1);

			double weight = ThreadLocalRandom.current().nextDouble(2, 10 + 1);
			double utility = ThreadLocalRandom.current().nextDouble(100, 1000 + 1);
			Item item = generateWeightedItemWithGivenValues(i, utility, weight);
			uncorrelatedInstancesWithSimilarWeightsItemList.add(item);
		}

		return uncorrelatedInstancesWithSimilarWeightsItemList;

	}

	public List<Item> createRandomInstance(int R, int S) {// Weights wj are chosen randomly in [1; R]
		// and the pro6ts pj in [wj − R=10; wj +
		// R=10] such that pj ¿ 1
		List<Item> InverseStronglyCorrelatedItemList = new ArrayList<>();
		for (int i = 0; i <= S; i++) {
			double utility = ThreadLocalRandom.current().nextDouble(1, R + 1);
			double weight = ThreadLocalRandom.current().nextDouble(1, R + 1);

			Item item = generateWeightedItemWithGivenValues(i, utility, weight);
			InverseStronglyCorrelatedItemList.add(item);
		}

		return InverseStronglyCorrelatedItemList;

	}

	public void exportToGlpkData(List<Item> itemList, double C) throws IOException {// export item list given to
																					// specified
		// location , format readable to glpk solver, C is Capacity of the knapsack
		int itemListLength = itemList.size();
		// Assigning the content of the file
		String text = "data;\nparam n := " + (itemListLength) + ";\nparam C := " + C + ";\nparam v :=";
		int counter = 0;
		for (Item item : itemList) {
			int u = (int) item.getUtility();
			text = text + " " + counter++ + " " + u;
		}
		text += ";";

		text += "\nparam w :=";
		counter = 0;
		for (Item item : itemList) {
			int v = (int) item.getWeight();
			text = text + " " + counter++ + " " + v;
		}
		text += ";";
		text += "\nend;";
		// Defining the file name of the file
		Path fileName = Path.of("C:\\\\Users\\Kotsos\\\\Desktop\\filename.dat");

		// Writing into the file
		Files.writeString(fileName, text);
	}

	public void runGlpkWithGivenFileLocation(String modelLocation, String fileLocation) {// runs the default glpk solver
																							// with the model and file
																							// of the location given
		ProcessBuilder processBuilder = new ProcessBuilder();

		// Run this on Windows, cmd, /c = terminate after this run
		processBuilder.command("cmd.exe", "/c",
				"  C:\\glpk-4.65\\w64\\glpsol.exe -m C:\\Users\\Kotsos\\Desktop\\ilpknap.mod -d C:\\Users\\Kotsos\\Desktop\\filename.dat -y C:\\\\Users\\\\Kotsos\\\\Desktop\\\\results.txt --dfs --cuts --tmlim 5 ");

		try {

			Process process = processBuilder.start();

			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));

			String line;
			while ((line = reader.readLine()) != null) {
				System.out.println(line);
			}

			int exitCode = process.waitFor();
			System.out.println("\nExited with error code : " + exitCode);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	public double readGLPKResult() {
		try {
			File myObj = new File("C:\\Users\\Kotsos\\Desktop\\results.txt");
			Scanner myReader = new Scanner(myObj);
			while (myReader.hasNextLine()) {
				String data = myReader.nextLine();
				System.out.println(data);
				if (data.contains("Profit.val")) {
					double numberOnlyData = Double.parseDouble(data.replaceAll("[^0-9]", ""));
					myReader.close();
					return numberOnlyData;
				}
			}

		} catch (FileNotFoundException e) {
			System.out.println("An error occurred.");
			e.printStackTrace();
		}
		return 0;

	}

	public double calculateWeightRatioOfItemList(double x, List<Item> itemList) {// returns sumOfWeigthInList/x
		double sum = 0;
		for (Item item : itemList) {
			double w = item.getWeight();
			sum += w;
		}

		return sum / x;
	}

	public void setRescaledWeight(List<Item> itemList, double C) {
		for (Item item : itemList) {
			item.setRescaledWeight(item.getRescaledWeight() / C);
		}
	}

	public double calculateWeightSum(List<Item> itemList) {
		double sum = 0;
		for (Item item : itemList) {
			sum += item.getWeight();
		}
		return sum;
	}

	public double calculateUtilitytSum(List<Item> itemList) {
		double sum = 0;
		for (Item item : itemList) {
			sum += item.getUtility();
		}
		return sum;
	}
}
